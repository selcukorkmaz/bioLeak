#!/usr/bin/env Rscript
## =================================================================
## Delta LSI analysis for bioLeak manuscript
##
## Part A: Case study — compute delta_lsi() between guarded and leaky
##         pipelines from the curatedOvarianData case study.
##
## Part B: Small simulation — Type I error and power of delta_lsi()
##         under no-leakage and peek_norm leakage conditions.
##
## Results are appended to paper/casestudy_results.rds (Part A)
## and saved to paper/sim_results/delta_lsi_sim.rds (Part B).
## =================================================================

cat("=== Delta LSI Analysis ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

library(bioLeak)

base_dir <- "paper"

## ---------------------------------------------------------------
## Part A: Case Study Delta LSI
## ---------------------------------------------------------------
cat("=== Part A: Case Study Delta LSI ===\n\n")

## We need the fitted objects. Re-load the case study data and re-fit
## both pipelines to obtain LeakFit objects for delta_lsi().
## (The stored casestudy_results.rds has summaries but not the fit objects.)

suppressPackageStartupMessages(library(curatedOvarianData))

## --- Reconstruct case study data (same logic as run_casestudy.R) ---
data(package = "curatedOvarianData")
all_ds <- data(package = "curatedOvarianData")$results[, "Item"]
all_ds <- gsub("\\s.*", "", all_ds)

study_list <- list()
for (ds_name in all_ds) {
  tryCatch({
    data(list = ds_name, package = "curatedOvarianData", envir = environment())
    eset <- get(ds_name, envir = environment())
    if (!inherits(eset, "ExpressionSet")) next
    pd <- Biobase::pData(eset)
    has_os <- "days_to_death" %in% names(pd) && "vital_status" %in% names(pd)
    if (!has_os) next
    os_days <- as.numeric(pd$days_to_death)
    vital   <- pd$vital_status
    valid <- !is.na(os_days) & !is.na(vital)
    if (sum(valid) < 50) next
    os_binary <- ifelse(os_days >= 1095, 1,
                        ifelse(vital == "deceased" & os_days < 1095, 0, NA))
    valid2 <- valid & !is.na(os_binary)
    if (sum(valid2) < 50) next
    if (length(unique(os_binary[valid2])) < 2) next
    study_list[[ds_name]] <- list(
      eset = eset, valid_idx = which(valid2),
      os_binary = os_binary[valid2], n = sum(valid2)
    )
  }, error = function(e) NULL)
}

gene_lists <- lapply(study_list, function(x) Biobase::featureNames(x$eset))
common_genes <- Reduce(intersect, gene_lists)

if (length(common_genes) < 100) {
  study_names <- names(study_list)
  best_pair <- NULL; best_common <- 0
  for (i in seq_along(study_names)) {
    for (j in seq_along(study_names)) {
      if (j <= i) next
      cg <- length(intersect(gene_lists[[i]], gene_lists[[j]]))
      if (cg > best_common) { best_common <- cg; best_pair <- c(i, j) }
    }
  }
  selected <- study_names[best_pair]
  common_genes <- intersect(gene_lists[[best_pair[1]]], gene_lists[[best_pair[2]]])
  for (sn in setdiff(study_names, selected)) {
    cg <- intersect(common_genes, gene_lists[[sn]])
    if (length(cg) >= 500) { common_genes <- cg; selected <- c(selected, sn) }
  }
  study_list <- study_list[selected]
}

if (length(common_genes) > 2000) {
  first_eset <- study_list[[1]]$eset
  gene_var <- apply(Biobase::exprs(first_eset)[common_genes, ], 1, var, na.rm = TRUE)
  common_genes <- names(sort(gene_var, decreasing = TRUE))[1:2000]
}

combined_X <- list(); combined_y <- c(); combined_study <- c(); sample_ids <- c()
for (sn in names(study_list)) {
  sl <- study_list[[sn]]
  expr_mat <- t(Biobase::exprs(sl$eset)[common_genes, sl$valid_idx, drop = FALSE])
  combined_X[[sn]] <- expr_mat
  combined_y <- c(combined_y, sl$os_binary)
  combined_study <- c(combined_study, rep(sn, sl$n))
  sample_ids <- c(sample_ids, paste0(sn, "_", seq_len(sl$n)))
}
X_combined <- do.call(rbind, combined_X)
rownames(X_combined) <- sample_ids

df_combined <- data.frame(X_combined, check.names = TRUE)
df_combined$y <- factor(combined_y)
df_combined$study <- combined_study

cat(sprintf("Dataset: N=%d, S=%d, G=%d\n",
            nrow(df_combined), length(unique(combined_study)), ncol(X_combined)))

## --- Guarded pipeline (study-blocked CV) with 20 repeats ---
## study_loocv ignores the repeats parameter (splits are deterministic),
## so we use batch_blocked with batch="study" instead.  This keeps studies
## as atomic units while randomizing their fold assignment each repeat,
## giving R_eff = 20 for delta_lsi() inference.
cat("\nFitting guarded pipeline (20 repeats, study-blocked)...\n")
splits_guarded <- make_split_plan(
  df_combined, outcome = "y", mode = "batch_blocked",
  batch = "study", v = 5, stratify = TRUE, seed = 42,
  repeats = 20
)

preprocess_guarded <- list(
  impute = list(method = "median"), normalize = list(method = "zscore"),
  filter = list(var_thresh = 0, iqr_thresh = 0),
  fs = list(method = "ttest", top_k = 100)
)

fit_guarded <- fit_resample(
  df_combined, outcome = "y", splits = splits_guarded,
  preprocess = preprocess_guarded,
  learner = "glmnet", metrics = "auc", seed = 42
)
cat(sprintf("Guarded AUC: %.3f\n", mean(fit_guarded@metrics$auc, na.rm = TRUE)))

## --- Leaky pipeline with 20 repeats (same splits for pairing) ---
cat("Fitting leaky pipeline (20 repeats, study-blocked)...\n")
df_leaky <- df_combined
y_num <- as.numeric(as.character(df_leaky$y))
df_leaky$leak_study_mean <- ave(y_num, df_leaky$study, FUN = mean)
y_sd <- sd(y_num); if (y_sd == 0) y_sd <- 1
df_leaky$leak_global_y <- (y_num - mean(y_num)) / y_sd
X_for_pca <- as.matrix(df_combined[, grep("^[Xx]", names(df_combined))])
X_for_pca[is.na(X_for_pca)] <- 0
pc1 <- prcomp(X_for_pca, center = TRUE, scale. = TRUE)$x[, 1]
pc1_sd <- sd(pc1); if (pc1_sd == 0) pc1_sd <- 1
df_leaky$leak_global_pc1 <- (pc1 - mean(pc1)) / pc1_sd
df_leaky$dummy_subject <- seq_len(nrow(df_leaky))

## Use same fold structure for pairing — convert study_loocv to naive by
## re-running fit_resample with same splits but including leak features.
## For proper delta_lsi pairing, both fits need the same split structure.
## So we fit the leaky pipeline on the SAME guarded splits to enable pairing.
fit_leaky <- fit_resample(
  df_leaky, outcome = "y", splits = splits_guarded,
  preprocess = list(
    impute = list(method = "median"), normalize = list(method = "zscore"),
    filter = list(var_thresh = 0, iqr_thresh = 0),
    fs = list(method = "none")
  ),
  learner = "glmnet", metrics = "auc", seed = 42
)
cat(sprintf("Leaky AUC: %.3f\n", mean(fit_leaky@metrics$auc, na.rm = TRUE)))

## --- Compute delta_lsi ---
cat("\nComputing delta_lsi()...\n")
dlsi <- delta_lsi(
  fit_naive = fit_leaky,
  fit_guarded = fit_guarded,
  metric = "auc",
  M_boot = 2000,
  M_flip = 10000,
  seed = 42
)

cat("\n--- Delta LSI Summary (Case Study) ---\n")
summary(dlsi)

## Save delta_lsi result
dlsi_cs <- list(
  delta_lsi_obj = dlsi,
  delta_metric = dlsi@delta_metric,
  delta_lsi_robust = dlsi@delta_lsi,
  ci_lsi = dlsi@delta_lsi_ci,
  ci_metric = dlsi@delta_metric_ci,
  p_value = dlsi@p_value,
  tier = dlsi@tier,
  n_repeats = dlsi@R_eff
)

## Append to casestudy_results.rds
cs <- readRDS(file.path(base_dir, "casestudy_results.rds"))
cs$delta_lsi <- dlsi_cs
saveRDS(cs, file.path(base_dir, "casestudy_results.rds"))
cat("Delta LSI appended to casestudy_results.rds\n")

## ---------------------------------------------------------------
## Part B: Delta LSI Simulation
## ---------------------------------------------------------------
cat("\n=== Part B: Delta LSI Simulation ===\n\n")

## Small simulation to evaluate delta_lsi() properties:
##   - Type I error: both arms get different pure-noise features (no leakage)
##     → δ should fluctuate around 0, testing non-trivial calibration
##   - Power: naive has peek_norm leakage → δ should be significantly > 0
## Design: n=200, p=20, 20 repeats × 5-fold CV, 20 seeds
##
## The null arm uses the SAME splits for both pipelines (enabling pairing
## and full inference tier). Both arms get 3 distinct pure-noise features
## (uncorrelated with y). Since each arm has the same number of features
## and the same regularization burden, the expected delta is 0, but
## individual realizations fluctuate due to different noise draws.
## This symmetric design avoids the systematic bias that arises when only
## one arm has extra features (regularization dilution), while still
## producing non-trivially-zero deltas that exercise the inference machinery.

library(future)
library(future.apply)

physical_cores <- parallel::detectCores(logical = FALSE)
n_workers <- max(1L, min(4L, physical_cores - 1L))
plan(multisession, workers = n_workers)

n_seeds <- 20
n_obs <- 200
p_dim <- 20
n_repeats <- 20
v_folds <- 5

run_dlsi_sim <- function(seed, inject_leakage) {
  library(bioLeak)
  set.seed(seed)

  X <- matrix(rnorm(n_obs * p_dim), n_obs, p_dim)
  colnames(X) <- sprintf("x%02d", seq_len(p_dim))
  subject <- sample(seq_len(max(5, n_obs %/% 5)), n_obs, replace = TRUE)

  linpred <- rowSums(X[, 1:5, drop = FALSE])
  linpred <- as.numeric(scale(linpred))
  y_prob <- pnorm(linpred)
  y <- rbinom(n_obs, 1, y_prob)

  batch <- sample(c("a", "b", "c"), n_obs, replace = TRUE)

  ## Base data (shared)
  df_base <- data.frame(X, y = factor(y))
  df_base$subject <- subject
  df_base$batch <- batch

  if (inject_leakage) {
    ## Power arm: naive has peek_norm leakage, guarded is clean
    df_guarded <- df_base
    df_naive <- df_base
    df_naive$leak_global <- as.numeric(y) + rnorm(n_obs, 0, 0.3)
  } else {
    ## Null arm: BOTH arms get different noise features (symmetric design).
    ## Each arm gets 3 pure-noise columns uncorrelated with y.
    ## Same feature count, same regularization burden, but different noise
    ## realizations → delta fluctuates around 0 (not systematically biased).
    ## Using same splits → pairing works → full inference tier.
    df_guarded <- df_base
    df_guarded$noise_g1 <- rnorm(n_obs)
    df_guarded$noise_g2 <- rnorm(n_obs)
    df_guarded$noise_g3 <- rnorm(n_obs)
    df_naive <- df_base
    df_naive$noise_n1 <- rnorm(n_obs)
    df_naive$noise_n2 <- rnorm(n_obs)
    df_naive$noise_n3 <- rnorm(n_obs)
  }

  preprocess <- list(
    impute = list(method = "median"), normalize = list(method = "zscore"),
    filter = list(var_thresh = 0, iqr_thresh = 0), fs = list(method = "none")
  )

  ## Same splits for both arms → pairing enabled → full inference tier
  splits <- make_split_plan(
    df_guarded, outcome = "y", mode = "subject_grouped",
    group = "subject", v = v_folds, repeats = n_repeats,
    stratify = TRUE, seed = seed, progress = FALSE
  )

  fit_g <- fit_resample(
    df_guarded, outcome = "y", splits = splits,
    preprocess = preprocess, learner = "glmnet", metrics = "auc", seed = seed
  )
  fit_n <- fit_resample(
    df_naive, outcome = "y", splits = splits,
    preprocess = preprocess, learner = "glmnet", metrics = "auc", seed = seed
  )

  dlsi <- tryCatch(
    delta_lsi(fit_naive = fit_n, fit_guarded = fit_g,
              metric = "auc", M_boot = 1000, M_flip = 5000, seed = seed),
    error = function(e) NULL
  )

  if (is.null(dlsi)) {
    return(data.frame(
      seed = seed, leakage = ifelse(inject_leakage, "peek_norm", "none"),
      delta_metric = NA, delta_lsi = NA, p_value = NA, tier = NA,
      ci_lo = NA, ci_hi = NA, stringsAsFactors = FALSE
    ))
  }

  ci <- if (length(dlsi@delta_lsi_ci) == 2) dlsi@delta_lsi_ci else c(NA, NA)

  data.frame(
    seed = seed,
    leakage = ifelse(inject_leakage, "peek_norm", "none"),
    delta_metric = dlsi@delta_metric,
    delta_lsi = dlsi@delta_lsi,
    p_value = dlsi@p_value,
    tier = dlsi@tier,
    ci_lo = ci[1],
    ci_hi = ci[2],
    stringsAsFactors = FALSE
  )
}

## Run null (no leakage) and alternative (peek_norm)
cat("Running null condition (no leakage)...\n")
null_res <- future_lapply(seq_len(n_seeds), function(s) {
  tryCatch(run_dlsi_sim(s, inject_leakage = FALSE),
           error = function(e) {
             data.frame(seed = s, leakage = "none", delta_metric = NA,
                        delta_lsi = NA, p_value = NA, tier = NA,
                        ci_lo = NA, ci_hi = NA, stringsAsFactors = FALSE)
           })
}, future.seed = TRUE)
null_df <- do.call(rbind, null_res)

cat("Running alternative condition (peek_norm leakage)...\n")
alt_res <- future_lapply(seq_len(n_seeds), function(s) {
  tryCatch(run_dlsi_sim(s, inject_leakage = TRUE),
           error = function(e) {
             data.frame(seed = s, leakage = "peek_norm", delta_metric = NA,
                        delta_lsi = NA, p_value = NA, tier = NA,
                        ci_lo = NA, ci_hi = NA, stringsAsFactors = FALSE)
           })
}, future.seed = TRUE)
alt_df <- do.call(rbind, alt_res)

dlsi_sim <- rbind(null_df, alt_df)

## Summary
cat("\n--- Delta LSI Simulation Summary ---\n")
for (cond in c("none", "peek_norm")) {
  sub <- dlsi_sim[dlsi_sim$leakage == cond & !is.na(dlsi_sim$p_value), ]
  if (nrow(sub) == 0) next
  cat(sprintf("\n  %s (n=%d):\n", cond, nrow(sub)))
  cat(sprintf("    Mean delta_metric: %.4f\n", mean(sub$delta_metric, na.rm = TRUE)))
  cat(sprintf("    Mean delta_lsi:    %.4f\n", mean(sub$delta_lsi, na.rm = TRUE)))
  cat(sprintf("    Rejection rate (p<0.05): %.1f%%\n",
              mean(sub$p_value < 0.05, na.rm = TRUE) * 100))
  valid_ci <- sub[!is.na(sub$ci_lo) & !is.na(sub$ci_hi), ]
  if (nrow(valid_ci) > 0) {
    coverage <- mean(valid_ci$ci_lo <= 0 & valid_ci$ci_hi >= 0) * 100
    cat(sprintf("    CI coverage of 0: %.1f%%\n", coverage))
  }
}

## Save
out_dir <- file.path(base_dir, "sim_results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
saveRDS(dlsi_sim, file.path(out_dir, "delta_lsi_sim.rds"))
cat("\nDelta LSI simulation saved to paper/sim_results/delta_lsi_sim.rds\n")

plan(sequential)
cat("\n=== Delta LSI analysis complete ===\n")
cat("End time:", format(Sys.time()), "\n")
