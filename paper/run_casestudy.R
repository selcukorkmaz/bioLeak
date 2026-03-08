#!/usr/bin/env Rscript
## =================================================================
## Case study: Multi-site gene-expression classification (Section 3.2)
##
## Uses curatedOvarianData to demonstrate guarded vs leaky pipelines
## =================================================================

cat("=== bioLeak Case Study ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

library(bioLeak)
suppressPackageStartupMessages(library(curatedOvarianData))

B <- 500

## ---------------------------------------------------------------
## 1. Load and prepare data
## ---------------------------------------------------------------
cat("Loading curatedOvarianData...\n")
data(package = "curatedOvarianData")

## Get all available datasets
all_ds <- data(package = "curatedOvarianData")$results[, "Item"]
## Filter to ExpressionSets
all_ds <- gsub("\\s.*", "", all_ds)

## Load each dataset and check eligibility
study_list <- list()
for (ds_name in all_ds) {
  tryCatch({
    data(list = ds_name, package = "curatedOvarianData", envir = environment())
    eset <- get(ds_name, envir = environment())
    if (!inherits(eset, "ExpressionSet")) next

    pd <- Biobase::pData(eset)
    ## Check for overall survival and vital status
    has_os <- "days_to_death" %in% names(pd) && "vital_status" %in% names(pd)
    if (!has_os) next

    ## Create binary endpoint: 3-year OS
    os_days <- as.numeric(pd$days_to_death)
    vital   <- pd$vital_status

    ## Need both time and status
    valid <- !is.na(os_days) & !is.na(vital)
    if (sum(valid) < 50) next

    ## Binary endpoint: survived >= 3 years (1095 days) OR censored after 3 years
    os_binary <- ifelse(os_days >= 1095, 1,
                        ifelse(vital == "deceased" & os_days < 1095, 0, NA))
    valid2 <- valid & !is.na(os_binary)
    if (sum(valid2) < 50) next

    ## Both classes needed
    if (length(unique(os_binary[valid2])) < 2) next

    study_list[[ds_name]] <- list(
      eset = eset,
      valid_idx = which(valid2),
      os_binary = os_binary[valid2],
      n = sum(valid2)
    )
    cat(sprintf("  %s: n=%d (class 0: %d, class 1: %d)\n",
                ds_name, sum(valid2),
                sum(os_binary[valid2] == 0), sum(os_binary[valid2] == 1)))
  }, error = function(e) NULL)
}

cat(sprintf("\nEligible studies: %d\n", length(study_list)))

## ---------------------------------------------------------------
## 2. Combine datasets
## ---------------------------------------------------------------
cat("Combining datasets...\n")

## Find common genes
gene_lists <- lapply(study_list, function(x) Biobase::featureNames(x$eset))
common_genes <- Reduce(intersect, gene_lists)
cat(sprintf("Common genes across all studies: %d\n", length(common_genes)))

## If too few common genes, use top studies with most overlap
if (length(common_genes) < 100) {
  cat("Too few common genes, selecting studies with maximum overlap...\n")
  ## Iteratively find best subset
  study_names <- names(study_list)
  ## Start with the two studies with most common genes
  best_pair <- NULL; best_common <- 0
  for (i in seq_along(study_names)) {
    for (j in seq_along(study_names)) {
      if (j <= i) next
      cg <- length(intersect(gene_lists[[i]], gene_lists[[j]]))
      if (cg > best_common) {
        best_common <- cg
        best_pair <- c(i, j)
      }
    }
  }
  selected <- study_names[best_pair]
  common_genes <- intersect(gene_lists[[best_pair[1]]], gene_lists[[best_pair[2]]])

  ## Add studies that maintain at least 500 common genes
  for (sn in setdiff(study_names, selected)) {
    cg <- intersect(common_genes, gene_lists[[sn]])
    if (length(cg) >= 500) {
      common_genes <- cg
      selected <- c(selected, sn)
    }
  }
  study_list <- study_list[selected]
  cat(sprintf("Selected %d studies with %d common genes\n",
              length(study_list), length(common_genes)))
}

## Limit to top 2000 most variable genes for speed
if (length(common_genes) > 2000) {
  ## Compute variance across first study to rank genes
  first_eset <- study_list[[1]]$eset
  gene_var <- apply(Biobase::exprs(first_eset)[common_genes, ], 1, var, na.rm = TRUE)
  common_genes <- names(sort(gene_var, decreasing = TRUE))[1:2000]
}

## Build combined matrix
combined_X  <- list()
combined_y  <- c()
combined_study <- c()
sample_ids  <- c()

for (sn in names(study_list)) {
  sl <- study_list[[sn]]
  expr_mat <- t(Biobase::exprs(sl$eset)[common_genes, sl$valid_idx, drop = FALSE])
  combined_X[[sn]]  <- expr_mat
  combined_y         <- c(combined_y, sl$os_binary)
  combined_study     <- c(combined_study, rep(sn, sl$n))
  sample_ids         <- c(sample_ids, paste0(sn, "_", seq_len(sl$n)))
}

X_combined <- do.call(rbind, combined_X)
rownames(X_combined) <- sample_ids

N_total <- nrow(X_combined)
S_total <- length(unique(combined_study))
G_total <- ncol(X_combined)

cat(sprintf("\nCombined dataset: N=%d samples, S=%d studies, G=%d genes\n",
            N_total, S_total, G_total))
cat(sprintf("Class distribution: 0=%d (%.1f%%), 1=%d (%.1f%%)\n",
            sum(combined_y == 0), 100 * mean(combined_y == 0),
            sum(combined_y == 1), 100 * mean(combined_y == 1)))

## Create data.frame
df_combined <- data.frame(X_combined, check.names = TRUE)
df_combined$y <- factor(combined_y)
df_combined$study <- combined_study

## ---------------------------------------------------------------
## 3. Guarded pipeline (study LOOCV)
## ---------------------------------------------------------------
cat("\n=== Guarded Pipeline (Study LOOCV) ===\n")

splits_guarded <- make_split_plan(
  df_combined, outcome = "y", mode = "study_loocv",
  study = "study", stratify = FALSE, seed = 42
)

preprocess_guarded <- list(
  impute    = list(method = "median"),
  normalize = list(method = "zscore"),
  filter    = list(var_thresh = 0, iqr_thresh = 0),
  fs        = list(method = "ttest", top_k = 100)
)

fit_guarded <- fit_resample(
  df_combined, outcome = "y", splits = splits_guarded,
  preprocess = preprocess_guarded,
  learner = "glmnet", metrics = c("auc"), seed = 42
)

cat("Guarded pipeline fitted. Auditing...\n")

X_ref_guarded <- as.matrix(df_combined[, grep("^[Xx]", names(df_combined))])

audit_guarded <- audit_leakage(
  fit_guarded, metric = "auc", B = B,
  perm_refit = FALSE, perm_stratify = TRUE, seed = 42,
  batch_cols = "study",
  coldata = df_combined[, "study", drop = FALSE],
  X_ref = X_ref_guarded,
  target_scan = TRUE,
  target_scan_multivariate = TRUE,
  target_threshold = 0.9,
  sim_threshold = 0.995,
  max_pairs = 50000
)

cat(sprintf("Guarded AUC: %.3f (SD: %.3f)\n",
            audit_guarded@permutation_gap$metric_obs,
            sd(fit_guarded@metrics$auc, na.rm = TRUE)))
cat(sprintf("Guarded perm gap: %.3f, p=%.4f\n",
            audit_guarded@permutation_gap$gap,
            audit_guarded@permutation_gap$p_value))

## ---------------------------------------------------------------
## 4. Leaky pipeline (standard 5-fold, global preprocessing, injected features)
## ---------------------------------------------------------------
cat("\n=== Leaky Comparator Pipeline ===\n")

## Inject 3 leaky features
df_leaky <- df_combined
y_num <- as.numeric(as.character(df_leaky$y))

## (i) Per-study mean outcome
df_leaky$leak_study_mean <- ave(y_num, df_leaky$study, FUN = mean)

## (ii) Globally z-scored outcome
y_sd <- sd(y_num); if (y_sd == 0) y_sd <- 1
df_leaky$leak_global_y <- (y_num - mean(y_num)) / y_sd

## (iii) Globally z-scored first PC of expression matrix
X_for_pca <- as.matrix(df_combined[, grep("^[Xx]", names(df_combined))])
## Handle missing values for PCA
X_for_pca[is.na(X_for_pca)] <- 0
pc1 <- prcomp(X_for_pca, center = TRUE, scale. = TRUE)$x[, 1]
pc1_sd <- sd(pc1); if (pc1_sd == 0) pc1_sd <- 1
df_leaky$leak_global_pc1 <- (pc1 - mean(pc1)) / pc1_sd

## Use standard 5-fold CV (no study blocking)
## Create a dummy subject variable for subject_grouped mode
df_leaky$dummy_subject <- seq_len(nrow(df_leaky))

splits_leaky <- make_split_plan(
  df_leaky, outcome = "y", mode = "subject_grouped",
  group = "dummy_subject", v = 5, stratify = TRUE, seed = 42
)

## Leaky preprocessing: no feature selection — all features (including leak
## features) go directly to glmnet, which will heavily weight the leak features.
## This mirrors a sloppy pipeline that skips proper feature screening.
preprocess_leaky <- list(
  impute    = list(method = "median"),
  normalize = list(method = "zscore"),
  filter    = list(var_thresh = 0, iqr_thresh = 0),
  fs        = list(method = "none")
)

fit_leaky <- fit_resample(
  df_leaky, outcome = "y", splits = splits_leaky,
  preprocess = preprocess_leaky,
  learner = "glmnet", metrics = c("auc"), seed = 42
)

cat("Leaky pipeline fitted. Auditing...\n")

X_ref_leaky <- as.matrix(df_leaky[, grep("^[Xx]|^leak", names(df_leaky))])

audit_leaky <- audit_leakage(
  fit_leaky, metric = "auc", B = B,
  perm_refit = FALSE, perm_stratify = TRUE, seed = 42,
  batch_cols = "study",
  coldata = df_leaky[, "study", drop = FALSE],
  X_ref = X_ref_leaky,
  target_scan = TRUE,
  target_scan_multivariate = TRUE,
  target_threshold = 0.9,
  sim_threshold = 0.995,
  max_pairs = 50000
)

cat(sprintf("Leaky AUC: %.3f (SD: %.3f)\n",
            audit_leaky@permutation_gap$metric_obs,
            sd(fit_leaky@metrics$auc, na.rm = TRUE)))
cat(sprintf("Leaky perm gap: %.3f, p=%.4f\n",
            audit_leaky@permutation_gap$gap,
            audit_leaky@permutation_gap$p_value))

## ---------------------------------------------------------------
## 5. Collect and save all results
## ---------------------------------------------------------------
cat("\n=== Collecting Results ===\n")

## Target scan for leaky pipeline
ta_leaky <- audit_leaky@target_assoc
ta_guarded <- audit_guarded@target_assoc

## Batch association
ba_leaky <- audit_leaky@batch_assoc
ba_guarded <- audit_guarded@batch_assoc

## Duplicates
dup_leaky <- audit_leaky@duplicates
dup_guarded <- audit_guarded@duplicates

results <- list(
  ## Dataset info
  N = N_total,
  S = S_total,
  G = G_total,
  study_names = names(study_list),
  study_sizes = sapply(study_list, function(x) x$n),

  ## Guarded pipeline
  guarded = list(
    auc = audit_guarded@permutation_gap$metric_obs,
    auc_sd = sd(fit_guarded@metrics$auc, na.rm = TRUE),
    gap = audit_guarded@permutation_gap$gap,
    p_value = audit_guarded@permutation_gap$p_value,
    batch_assoc = ba_guarded,
    target_assoc = ta_guarded,
    duplicates = dup_guarded,
    multi_p = audit_guarded@info$target_multivariate$p_value
  ),

  ## Leaky pipeline
  leaky = list(
    auc = audit_leaky@permutation_gap$metric_obs,
    auc_sd = sd(fit_leaky@metrics$auc, na.rm = TRUE),
    gap = audit_leaky@permutation_gap$gap,
    p_value = audit_leaky@permutation_gap$p_value,
    batch_assoc = ba_leaky,
    target_assoc = ta_leaky,
    duplicates = dup_leaky,
    multi_p = audit_leaky@info$target_multivariate$p_value
  )
)

saveRDS(results, "paper/casestudy_results.rds")

## Print summary
cat("\n=== SUMMARY ===\n")
cat(sprintf("Dataset: N=%d, S=%d studies, G=%d genes\n", N_total, S_total, G_total))
cat(sprintf("\nGuarded pipeline: AUC=%.3f (SD=%.3f), gap=%.3f, p=%.4f\n",
            results$guarded$auc, results$guarded$auc_sd,
            results$guarded$gap, results$guarded$p_value))
cat(sprintf("Leaky pipeline:   AUC=%.3f (SD=%.3f), gap=%.3f, p=%.4f\n",
            results$leaky$auc, results$leaky$auc_sd,
            results$leaky$gap, results$leaky$p_value))

if (nrow(ba_guarded) > 0)
  cat(sprintf("\nBatch association (guarded): V=%.3f, p=%.4f\n",
              ba_guarded$cramer_v[1], ba_guarded$pval[1]))
if (nrow(ba_leaky) > 0)
  cat(sprintf("Batch association (leaky):   V=%.3f, p=%.4f\n",
              ba_leaky$cramer_v[1], ba_leaky$pval[1]))

## Target scan
if (nrow(ta_leaky) > 0) {
  flagged <- ta_leaky[ta_leaky$flag == TRUE, ]
  cat(sprintf("\nTarget scan (leaky): %d features flagged\n", nrow(flagged)))
  if (nrow(flagged) > 0) {
    for (r in seq_len(min(10, nrow(flagged)))) {
      cat(sprintf("  %s: score=%.3f\n", flagged$feature[r], flagged$score[r]))
    }
  }
}

cat(sprintf("\nDuplicates (guarded): %d pairs\n", nrow(dup_guarded)))
cat(sprintf("Duplicates (leaky):   %d pairs\n", nrow(dup_leaky)))

cat(sprintf("\nMultivariate scan p (guarded): %s\n",
            ifelse(is.null(results$guarded$multi_p), "NULL",
                   sprintf("%.4f", results$guarded$multi_p))))
cat(sprintf("Multivariate scan p (leaky):   %s\n",
            ifelse(is.null(results$leaky$multi_p), "NULL",
                   sprintf("%.4f", results$leaky$multi_p))))

cat("\n=== Case study complete ===\n")
cat("End time:", format(Sys.time()), "\n")
