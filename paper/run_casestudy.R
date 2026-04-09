#!/usr/bin/env Rscript
## =================================================================
## Case study: Multi-site gene-expression classification (Section 3.2)
##
## Uses curatedOvarianData to demonstrate guarded vs leaky pipelines
## =================================================================

cat("=== bioLeak Case Study ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

library(bioLeak)

B <- 500

## ---------------------------------------------------------------
## 1. Load and prepare data
## ---------------------------------------------------------------
source("paper/_helpers.R")
cs_data <- load_ovarian_casestudy()

df_combined    <- cs_data$df_combined
X_combined     <- cs_data$X_combined
combined_y     <- cs_data$combined_y
combined_study <- cs_data$combined_study
study_list     <- cs_data$study_list

N_total <- nrow(X_combined)
S_total <- length(unique(combined_study))
G_total <- ncol(X_combined)

## ---------------------------------------------------------------
## 3. Guarded pipeline (study LOOCV)
## ---------------------------------------------------------------
cat("\n=== Guarded Pipeline (Study LOOCV) ===\n")

splits_guarded <- make_split_plan(
  df_combined, outcome = "y", mode = "study_loocv",
  study = "study", stratify = FALSE, seed = 42
)

preprocess_guarded <- list(
  impute    = list(method = "median", winsor = FALSE),
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
## 3b. Guarded pipeline WITHOUT feature selection (sensitivity analysis)
## ---------------------------------------------------------------
cat("\n=== Guarded Pipeline — No Feature Selection (Sensitivity) ===\n")

preprocess_guarded_nofs <- list(
  impute    = list(method = "median", winsor = FALSE),
  normalize = list(method = "zscore"),
  filter    = list(var_thresh = 0, iqr_thresh = 0),
  fs        = list(method = "none")
)

fit_guarded_nofs <- fit_resample(
  df_combined, outcome = "y", splits = splits_guarded,
  preprocess = preprocess_guarded_nofs,
  learner = "glmnet", metrics = c("auc"), seed = 42
)

audit_guarded_nofs <- audit_leakage(
  fit_guarded_nofs, metric = "auc", B = B,
  perm_refit = FALSE, perm_stratify = TRUE, seed = 42,
  batch_cols = "study",
  coldata = df_combined[, "study", drop = FALSE],
  X_ref = X_ref_guarded,
  target_scan = FALSE,
  target_scan_multivariate = FALSE
)

cat(sprintf("Guarded (no FS) AUC: %.3f (SD: %.3f)\n",
            audit_guarded_nofs@permutation_gap$metric_obs,
            sd(fit_guarded_nofs@metrics$auc, na.rm = TRUE)))
cat(sprintf("Guarded (no FS) perm gap: %.3f, p=%.4f\n",
            audit_guarded_nofs@permutation_gap$gap,
            audit_guarded_nofs@permutation_gap$p_value))

## ---------------------------------------------------------------
## 3c. Naive-split clean pipeline (decomposition arm: split effect only)
## ---------------------------------------------------------------
cat("\n=== Naive-Split Clean Pipeline (Split-Design Sensitivity) ===\n")

## Same row-wise 5-fold as leaky, but on CLEAN data (no leak features, no FS)
df_naive_clean <- df_combined
df_naive_clean$dummy_subject <- seq_len(nrow(df_naive_clean))

splits_naive_clean <- make_split_plan(
  df_naive_clean, outcome = "y", mode = "subject_grouped",
  group = "dummy_subject", v = 5, stratify = TRUE, seed = 42
)

fit_naive_clean <- fit_resample(
  df_naive_clean, outcome = "y", splits = splits_naive_clean,
  preprocess = preprocess_guarded_nofs,
  learner = "glmnet", metrics = c("auc"), seed = 42
)

audit_naive_clean <- audit_leakage(
  fit_naive_clean, metric = "auc", B = B,
  perm_refit = FALSE, perm_stratify = TRUE, seed = 42,
  batch_cols = "study",
  coldata = df_naive_clean[, "study", drop = FALSE],
  X_ref = X_ref_guarded,
  target_scan = FALSE,
  target_scan_multivariate = FALSE
)

cat(sprintf("Naive clean AUC: %.3f (SD: %.3f)\n",
            audit_naive_clean@permutation_gap$metric_obs,
            sd(fit_naive_clean@metrics$auc, na.rm = TRUE)))
cat(sprintf("Naive clean perm gap: %.3f, p=%.4f\n",
            audit_naive_clean@permutation_gap$gap,
            audit_naive_clean@permutation_gap$p_value))

## ---------------------------------------------------------------
## 4. Leaky pipeline (standard 5-fold, global preprocessing, injected features)
## ---------------------------------------------------------------
cat("\n=== Leaky Comparator Pipeline ===\n")

## Inject 3 leaky features
df_leaky <- df_combined
y_num <- as.numeric(as.character(df_leaky$y))

## (i) Per-study mean outcome
df_leaky$leak_study_mean <- ave(y_num, df_leaky$study, FUN = mean)

## (ii) Noisy outcome encoding (y + Gaussian noise, sigma = 1.0)
set.seed(42)
df_leaky$leak_global_y <- y_num + rnorm(nrow(df_leaky), 0, 1.0)

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
  impute    = list(method = "median", winsor = FALSE),
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

  ## Guarded pipeline without feature selection (sensitivity)
  guarded_nofs = list(
    auc = audit_guarded_nofs@permutation_gap$metric_obs,
    auc_sd = sd(fit_guarded_nofs@metrics$auc, na.rm = TRUE),
    gap = audit_guarded_nofs@permutation_gap$gap,
    p_value = audit_guarded_nofs@permutation_gap$p_value
  ),

  ## Naive-split clean pipeline (split-design decomposition arm)
  naive_clean = list(
    auc = audit_naive_clean@permutation_gap$metric_obs,
    auc_sd = sd(fit_naive_clean@metrics$auc, na.rm = TRUE),
    gap = audit_naive_clean@permutation_gap$gap,
    p_value = audit_naive_clean@permutation_gap$p_value
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
cat(sprintf("Guarded (no FS):  AUC=%.3f (SD=%.3f), gap=%.3f, p=%.4f\n",
            results$guarded_nofs$auc, results$guarded_nofs$auc_sd,
            results$guarded_nofs$gap, results$guarded_nofs$p_value))
cat(sprintf("Naive clean:      AUC=%.3f (SD=%.3f), gap=%.3f, p=%.4f\n",
            results$naive_clean$auc, results$naive_clean$auc_sd,
            results$naive_clean$gap, results$naive_clean$p_value))
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
