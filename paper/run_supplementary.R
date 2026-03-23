#!/usr/bin/env Rscript
## =================================================================
## Supplementary simulations for bioLeak manuscript (Section 3.1.3)
##
## 1. Splitting mode robustness: 4 modes x 5 leakage at (n=500, p=20, s=1.0)
## 2. Target leakage scan: same subset with target_scan=TRUE
##
## NOTE: Like run_simulation.R, this script uses a CUSTOM data-generation
## pipeline that differs from the package-level simulate_leakage_suite().
## See run_simulation.R for a full description of the differences.
## Do NOT replace this script with simulate_leakage_suite() calls.
## =================================================================

cat("=== bioLeak Supplementary Simulations ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

library(bioLeak)
library(future)
library(future.apply)

physical_cores <- parallel::detectCores(logical = FALSE)
n_workers <- max(1L, min(4L, physical_cores - 1L))
plan(multisession, workers = n_workers)
cat(sprintf("Parallel backend: %d workers\n\n", n_workers))

B       <- 50
n_seeds <- 20
seeds   <- seq_len(n_seeds)
n_fix   <- 500
p_fix   <- 20
s_fix   <- 1.0

leakage_types <- c("none", "subject_overlap", "batch_confounded",
                   "peek_norm", "lookahead")
split_modes   <- c("subject_grouped", "batch_blocked", "study_loocv",
                   "time_series")

out_dir <- "paper/sim_results"

## ---------------------------------------------------------------
## Part 1: Splitting mode robustness
## ---------------------------------------------------------------
cat("=== Part 1: Splitting Mode Robustness ===\n")

grid_modes <- expand.grid(
  mode    = split_modes,
  leakage = leakage_types,
  stringsAsFactors = FALSE
)
cat("Configs:", nrow(grid_modes), "\n\n")

run_one_mode <- function(seed, mode_str, leakage_type, B, s_val) {
  library(bioLeak)
  set.seed(seed)

  n <- n_fix; p <- p_fix; s <- s_val

  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- sprintf("x%02d", seq_len(p))
  subject  <- sample(seq_len(max(5, n %/% 5)), n, replace = TRUE)

  ## Generate outcome: s=0 → pure noise, s>0 → signal + AR noise
  if (s > 0) {
    linpred <- rowSums(X[, seq_len(min(5, p)), drop = FALSE])
    linpred <- as.numeric(scale(linpred)) * s
    ar_noise <- as.numeric(arima.sim(model = list(ar = 0.9), n = n,
                                     sd = 0.3))
    linpred <- linpred + ar_noise
  } else {
    linpred <- rep(0, n)
  }

  y_prob  <- pnorm(linpred)
  y       <- rbinom(n, 1, y_prob)

  ## Batch: independent of outcome by default
  batch <- sample(c("a","b","c"), n, replace = TRUE)
  study    <- sample(seq_len(max(3, n %/% 80)), n, replace = TRUE)
  time_var <- seq_len(n)

  ## Inject leakage features
  if (leakage_type == "subject_overlap") {
    X <- cbind(X, leak_subj = ave(y, subject, FUN = mean))
  } else if (leakage_type == "batch_confounded") {
    ## Make batch outcome-dependent, then leak the group mean
    batch <- ifelse(y == 1,
                    sample(c("a","b","c"), n, replace = TRUE, prob = c(0.6, 0.2, 0.2)),
                    sample(c("a","b","c"), n, replace = TRUE, prob = c(0.15, 0.5, 0.35)))
    X <- cbind(X, leak_batch = ave(y, batch, FUN = mean))
  } else if (leakage_type == "peek_norm") {
    X <- cbind(X, leak_global = as.numeric(y) + rnorm(n, 0, 0.3))
  } else if (leakage_type == "lookahead") {
    biomarker <- linpred + rnorm(n, 0, 0.5)
    X <- cbind(X, leak_future = c(biomarker[-1], biomarker[n]))
  }

  df <- data.frame(X, y = factor(y))
  df$subject <- subject; df$batch <- batch
  df$study <- study; df$time <- time_var

  ## Determine split args based on mode
  split_args <- list(
    x = df, outcome = "y", mode = mode_str,
    batch = "batch", study = "study", time = "time",
    v = 5, stratify = TRUE, seed = seed,
    progress = FALSE
  )
  if (mode_str == "subject_grouped") split_args$group <- "subject"
  if (mode_str == "time_series") {
    split_args$horizon <- 1
  }

  splits <- do.call(make_split_plan, split_args)

  preprocess <- list(
    impute = list(method = "median"), normalize = list(method = "zscore"),
    filter = list(var_thresh = 0, iqr_thresh = 0), fs = list(method = "none")
  )

  fit <- fit_resample(
    df, outcome = "y", splits = splits,
    preprocess = preprocess,
    learner = "glmnet", metrics = "auc", seed = seed
  )

  aud <- audit_leakage(
    fit, metric = "auc", B = B,
    perm_refit = FALSE, perm_stratify = TRUE,
    seed = seed, target_scan = FALSE, return_perm = FALSE
  )

  data.frame(
    seed = seed, s = s, mode = mode_str, leakage = leakage_type,
    metric_obs = aud@permutation_gap$metric_obs,
    gap = aud@permutation_gap$gap,
    p_value = aud@permutation_gap$p_value,
    stringsAsFactors = FALSE
  )
}

## Run at both s=0 (leakage-specific) and s=1.0 (signal present)
mode_signals <- c(0, 1.0)
mode_results <- list()
idx <- 0
for (s_mode in mode_signals) {
  cat(sprintf("\n--- s = %.1f ---\n", s_mode))
  for (i in seq_len(nrow(grid_modes))) {
    cfg <- grid_modes[i, ]
    idx <- idx + 1
    t0 <- proc.time()
    cat(sprintf("  s=%.1f [%d/%d] mode=%-18s leakage=%-20s ... ",
                s_mode, i, nrow(grid_modes), cfg$mode, cfg$leakage))

    res_list <- future_lapply(seeds, function(seed) {
      tryCatch(
        run_one_mode(seed, cfg$mode, cfg$leakage, B, s_mode),
        error = function(e) {
          data.frame(seed = seed, s = s_mode, mode = cfg$mode,
                     leakage = cfg$leakage,
                     metric_obs = NA, gap = NA, p_value = NA,
                     stringsAsFactors = FALSE)
        }
      )
    }, future.seed = TRUE)

    mode_results[[idx]] <- do.call(rbind, res_list)
    cat(sprintf("done (%.1fs)\n", (proc.time() - t0)[3]))
  }
}

mode_all <- do.call(rbind, mode_results)
saveRDS(mode_all, file.path(out_dir, "supplementary_modes.rds"))
cat("Part 1 saved.\n\n")

## ---------------------------------------------------------------
## Part 2: Target leakage scan
## ---------------------------------------------------------------
cat("=== Part 2: Target Leakage Scan ===\n")

## Run at both s=0 and s=1.0:
##   s=0  → true null for real features; leakage features are the ONLY
##           source of target association, so multivariate scan can
##           properly evaluate leakage-specific detection.
##   s=1.0 → real signal present; shows that the multivariate scan
##           cannot separate leakage from genuine signal.
scan_signals <- c(0, 1.0)
scan_results <- list()
idx <- 0

for (s_scan in scan_signals) {
  cat(sprintf("\n--- s = %.1f ---\n", s_scan))
  for (lk in leakage_types) {
    t0 <- proc.time()
    cat(sprintf("  s=%.1f  leakage=%-20s ... ", s_scan, lk))

    res_list <- future_lapply(seeds, function(seed) {
      library(bioLeak)
      set.seed(seed)
      n <- n_fix; p <- p_fix; s <- s_scan

      X <- matrix(rnorm(n * p), n, p)
      colnames(X) <- sprintf("x%02d", seq_len(p))
      subject <- sample(seq_len(max(5, n %/% 5)), n, replace = TRUE)

      ## Generate outcome: s=0 → pure noise, s>0 → signal + AR noise
      if (s > 0) {
        linpred <- rowSums(X[, seq_len(min(5, p)), drop = FALSE])
        linpred <- as.numeric(scale(linpred)) * s
        ar_noise <- as.numeric(arima.sim(model = list(ar = 0.9), n = n,
                                         sd = 0.3))
        linpred <- linpred + ar_noise
      } else {
        linpred <- rep(0, n)
      }

      y_prob  <- pnorm(linpred)
      y       <- rbinom(n, 1, y_prob)

      ## Batch: independent of outcome by default
      batch <- sample(c("a","b","c"), n, replace = TRUE)
      study   <- sample(seq_len(max(3, n %/% 80)), n, replace = TRUE)
      time_var <- seq_len(n)

      ## Inject leakage features
      if (lk == "subject_overlap") {
        X <- cbind(X, leak_subj = ave(y, subject, FUN = mean))
      } else if (lk == "batch_confounded") {
        batch <- ifelse(y == 1,
                        sample(c("a","b","c"), n, replace = TRUE, prob = c(0.6, 0.2, 0.2)),
                        sample(c("a","b","c"), n, replace = TRUE, prob = c(0.15, 0.5, 0.35)))
        X <- cbind(X, leak_batch = ave(y, batch, FUN = mean))
      } else if (lk == "peek_norm") {
        X <- cbind(X, leak_global = as.numeric(y) + rnorm(n, 0, 0.3))
      } else if (lk == "lookahead") {
        biomarker <- linpred + rnorm(n, 0, 0.5)
        X <- cbind(X, leak_future = c(biomarker[-1], biomarker[n]))
      }

      df <- data.frame(X, y = factor(y))
      df$subject <- subject; df$batch <- batch
      df$study <- study; df$time <- time_var

      splits <- make_split_plan(
        df, outcome = "y", mode = "subject_grouped",
        group = "subject", batch = "batch", study = "study", time = "time",
        v = 5, stratify = TRUE, seed = seed,
        progress = FALSE
      )

      preprocess <- list(
        impute = list(method = "median"), normalize = list(method = "zscore"),
        filter = list(var_thresh = 0, iqr_thresh = 0), fs = list(method = "none")
      )

      fit <- fit_resample(
        df, outcome = "y", splits = splits,
        preprocess = preprocess,
        learner = "glmnet", metrics = "auc", seed = seed
      )

      ## Get X_ref (numeric features only)
      X_ref <- as.matrix(df[, grep("^x|^leak", names(df))])

      aud <- audit_leakage(
        fit, metric = "auc", B = B,
        perm_refit = FALSE, perm_stratify = TRUE,
        seed = seed,
        X_ref = X_ref,
        target_scan = TRUE,
        target_scan_multivariate = TRUE,
        target_threshold = 0.9,
        return_perm = FALSE
      )

      ## Extract target scan results
      ta <- aud@target_assoc
      leak_flagged <- any(grepl("^leak", ta$feature[ta$flag == TRUE]))

      ## Multivariate scan p-value
      mv <- aud@info$target_multivariate
      multi_p <- if (!is.null(mv) && nrow(mv) > 0 && "p_value" %in% names(mv))
        mv$p_value[1] else NA

      data.frame(
        seed = seed, s = s, leakage = lk,
        n_flagged = sum(ta$flag, na.rm = TRUE),
        leak_feature_flagged = leak_flagged,
        multivariate_p = multi_p,
        stringsAsFactors = FALSE
      )
    }, future.seed = TRUE)

    idx <- idx + 1
    scan_results[[idx]] <- do.call(rbind, res_list)
    cat(sprintf("done (%.1fs)\n", (proc.time() - t0)[3]))
  }
}

scan_all <- do.call(rbind, scan_results)
saveRDS(scan_all, file.path(out_dir, "supplementary_target_scan.rds"))
cat("Part 2 saved.\n\n")

cat("=== Supplementary simulations complete ===\n")
cat("End time:", format(Sys.time()), "\n")
