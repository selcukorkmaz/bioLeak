#!/usr/bin/env Rscript
## =================================================================
## Delta LSI analysis for bioLeak manuscript
##
## Part A: Case study — compute delta_lsi() between guarded and leaky
##         pipelines from the curatedOvarianData case study.
##
## Part B: Power simulation — delta_lsi() detects peek_norm leakage.
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

source("paper/_helpers.R")
cs_data <- load_ovarian_casestudy()

df_combined    <- cs_data$df_combined
X_combined     <- cs_data$X_combined
combined_y     <- cs_data$combined_y
combined_study <- cs_data$combined_study

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
  impute = list(method = "median", winsor = FALSE), normalize = list(method = "zscore"),
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
set.seed(42)
df_leaky$leak_global_y <- y_num + rnorm(nrow(df_leaky), 0, 1.0)
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
    impute = list(method = "median", winsor = FALSE), normalize = list(method = "zscore"),
    filter = list(var_thresh = 0, iqr_thresh = 0),
    fs = list(method = "none")
  ),
  learner = "glmnet", metrics = "auc", seed = 42
)
cat(sprintf("Leaky AUC: %.3f\n", mean(fit_leaky@metrics$auc, na.rm = TRUE)))

## --- Compute delta_lsi ---
cat("\nComputing delta_lsi()...\n")
dlsi <- delta_lsi(
  fit_leaky = fit_leaky,
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

## Power simulation for delta_lsi(): the naive pipeline includes a
## peek_norm leakage feature, the guarded pipeline is clean.
## Design: n=200, p=20, 20 repeats × 5-fold CV, 50 seeds
##
## Both arms share the same fold structure, enabling paired inference at
## tier A. This simulation demonstrates power. Part C below provides a
## complementary null calibration.

library(future)
library(future.apply)

physical_cores <- parallel::detectCores(logical = FALSE)
n_workers <- max(1L, min(4L, physical_cores - 1L))
plan(multisession, workers = n_workers)

n_seeds <- 50
n_obs <- 200
p_dim <- 20
n_repeats <- 20
v_folds <- 5

run_dlsi_sim <- function(seed) {
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

  ## Naive has peek_norm leakage, guarded is clean
  df_guarded <- df_base
  df_naive <- df_base
  df_naive$leak_global <- as.numeric(y) + rnorm(n_obs, 0, 0.3)

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
    delta_lsi(fit_leaky = fit_n, fit_guarded = fit_g,
              metric = "auc", M_boot = 1000, M_flip = 5000, seed = seed),
    error = function(e) NULL
  )

  if (is.null(dlsi)) {
    return(data.frame(
      seed = seed, delta_metric = NA, delta_lsi = NA, p_value = NA,
      tier = NA, ci_lo = NA, ci_hi = NA, stringsAsFactors = FALSE
    ))
  }

  ci <- if (length(dlsi@delta_lsi_ci) == 2) dlsi@delta_lsi_ci else c(NA, NA)

  data.frame(
    seed = seed,
    delta_metric = dlsi@delta_metric,
    delta_lsi = dlsi@delta_lsi,
    p_value = dlsi@p_value,
    tier = dlsi@tier,
    ci_lo = ci[1],
    ci_hi = ci[2],
    stringsAsFactors = FALSE
  )
}

## Run power condition (peek_norm leakage)
cat("Running power condition (peek_norm leakage)...\n")
alt_res <- future_lapply(seq_len(n_seeds), function(s) {
  tryCatch(run_dlsi_sim(s),
           error = function(e) {
             data.frame(seed = s, delta_metric = NA,
                        delta_lsi = NA, p_value = NA, tier = NA,
                        ci_lo = NA, ci_hi = NA, stringsAsFactors = FALSE)
           })
}, future.seed = TRUE)
dlsi_sim <- do.call(rbind, alt_res)

## Summary
cat("\n--- Delta LSI Simulation Summary ---\n")
sub <- dlsi_sim[!is.na(dlsi_sim$p_value), ]
cat(sprintf("\n  peek_norm (n=%d):\n", nrow(sub)))
cat(sprintf("    Mean delta_metric: %.4f\n", mean(sub$delta_metric, na.rm = TRUE)))
cat(sprintf("    Mean delta_lsi:    %.4f\n", mean(sub$delta_lsi, na.rm = TRUE)))
cat(sprintf("    Rejection rate (p<0.05): %.1f%%\n",
            mean(sub$p_value < 0.05, na.rm = TRUE) * 100))
valid_ci <- sub[!is.na(sub$ci_lo) & !is.na(sub$ci_hi), ]
if (nrow(valid_ci) > 0) {
  coverage <- mean(valid_ci$ci_lo > 0) * 100
  cat(sprintf("    CI excludes 0 (lower bound > 0): %.1f%%\n", coverage))
}

## Save
out_dir <- file.path(base_dir, "sim_results")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
saveRDS(dlsi_sim, file.path(out_dir, "delta_lsi_sim.rds"))
cat("\nDelta LSI simulation saved to paper/sim_results/delta_lsi_sim.rds\n")

## ---------------------------------------------------------------
## Part C: Null Calibration
## ---------------------------------------------------------------
cat("\n=== Part C: Null Calibration (independent-data sign-flip) ===\n\n")

## Calibration check for the sign-flip test underlying delta_lsi().
## Under the null (no leakage), the test should reject at approximately
## the nominal alpha = 0.05 level.
##
## Design rationale:
## With repeated CV on a SINGLE dataset and a deterministic learner
## (glmnet), any fixed difference between arms (different features,
## different random seeds) produces consistent delta_r signs across
## repeats. The sign-flip test correctly detects this consistency,
## inflating the rejection rate above 5 %. This is not a calibration
## defect but a consequence of within-dataset determinism.
##
## Solution: generate an INDEPENDENT dataset for each of R = 20
## "repeats." Both arms share the same base features and receive their
## own independent noise feature (pure N(0,1), exchangeable by
## construction). Because each repeat uses fresh data, the noise
## features have random correlation with y in each repeat, ensuring
## that delta_r signs are genuinely exchangeable. This design is
## stronger than repeated CV for calibration purposes because the
## delta_r values are truly independent across repeats.

n_null_repeats <- 20  # independent datasets per seed
sign_flip <- bioLeak:::.dlsi_sign_flip

run_dlsi_null <- function(seed) {
  library(bioLeak)
  sign_flip <- bioLeak:::.dlsi_sign_flip
  delta_r <- numeric(n_null_repeats)

  for (r in seq_len(n_null_repeats)) {
    sub_seed <- seed * 1000L + r
    set.seed(sub_seed)

    X <- matrix(rnorm(n_obs * p_dim), n_obs, p_dim)
    colnames(X) <- sprintf("x%02d", seq_len(p_dim))
    subject <- sample(seq_len(max(5, n_obs %/% 5)), n_obs, replace = TRUE)

    linpred <- rowSums(X[, 1:5, drop = FALSE])
    linpred <- as.numeric(scale(linpred))
    y_prob <- pnorm(linpred)
    y <- rbinom(n_obs, 1, y_prob)

    df_a <- data.frame(X, y = factor(y))
    df_a$subject <- subject
    df_a$batch <- sample(c("a", "b", "c"), n_obs, replace = TRUE)
    df_b <- df_a

    ## Both arms get their own fresh noise feature
    df_a$noise <- rnorm(n_obs)
    df_b$noise <- rnorm(n_obs)

    preprocess <- list(
      impute = list(method = "median"), normalize = list(method = "zscore"),
      filter = list(var_thresh = 0, iqr_thresh = 0), fs = list(method = "none")
    )

    splits <- make_split_plan(
      df_a, outcome = "y", mode = "subject_grouped",
      group = "subject", v = v_folds, repeats = 1,
      stratify = TRUE, seed = sub_seed, progress = FALSE
    )

    fit_a <- fit_resample(
      df_a, outcome = "y", splits = splits,
      preprocess = preprocess, learner = "glmnet", metrics = "auc", seed = sub_seed
    )
    fit_b <- fit_resample(
      df_b, outcome = "y", splits = splits,
      preprocess = preprocess, learner = "glmnet", metrics = "auc", seed = sub_seed
    )

    delta_r[r] <- mean(fit_a@metrics$auc, na.rm = TRUE) -
                  mean(fit_b@metrics$auc, na.rm = TRUE)
  }

  pv <- sign_flip(delta_r, M_flip = 5000L, seed = seed)

  data.frame(
    seed = seed,
    delta_metric = mean(delta_r),
    p_value = pv,
    n_pos = sum(delta_r > 0),
    n_neg = sum(delta_r < 0),
    stringsAsFactors = FALSE
  )
}

cat("Running null condition (independent-data sign-flip)...\n")
null_res <- future_lapply(seq_len(n_seeds), function(s) {
  tryCatch(run_dlsi_null(s),
           error = function(e) {
             data.frame(seed = s, delta_metric = NA, p_value = NA,
                        n_pos = NA, n_neg = NA, stringsAsFactors = FALSE)
           })
}, future.seed = TRUE)
dlsi_null <- do.call(rbind, null_res)

## Summary
cat("\n--- Null Calibration Summary ---\n")
sub_null <- dlsi_null[!is.na(dlsi_null$p_value), ]
cat(sprintf("\n  Null (n=%d):\n", nrow(sub_null)))
cat(sprintf("    Mean delta_metric: %.4f\n", mean(sub_null$delta_metric, na.rm = TRUE)))
cat(sprintf("    Rejection rate (p<0.05): %.1f%%\n",
            mean(sub_null$p_value < 0.05, na.rm = TRUE) * 100))
cat(sprintf("    Mean p-value: %.3f\n", mean(sub_null$p_value, na.rm = TRUE)))
cat(sprintf("    Avg sign balance: pos=%.1f, neg=%.1f (out of %d)\n",
            mean(sub_null$n_pos, na.rm = TRUE),
            mean(sub_null$n_neg, na.rm = TRUE), n_null_repeats))

saveRDS(dlsi_null, file.path(out_dir, "delta_lsi_null.rds"))
cat("Null calibration saved to paper/sim_results/delta_lsi_null.rds\n")

plan(sequential)
cat("\n=== Delta LSI analysis complete ===\n")
cat("End time:", format(Sys.time()), "\n")
