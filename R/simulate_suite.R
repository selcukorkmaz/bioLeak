#' Simulate leakage scenarios and audit results
#'
#' @description
#' Runs a Monte Carlo simulation suite for evaluating leakage detection
#' across multiple leakage types, CV modes, and learners.
#'
#' For each seed, this function:
#' 1. Simulates a dataset with a defined leakage mechanism.
#' 2. Generates leakage-resistant CV splits using `make_splits()`.
#' 3. Fits a machine learning model using `fit_resample()`.
#' 4. Audits leakage with the permutation-gap diagnostic (`audit_leakage()`).
#'
#' @details
#' The function can optionally run in parallel using `future.apply`.
#' Set a future plan (e.g., `future::plan("multisession")`) before calling.
#' \code{signal_strength} scales the linear predictor before converting to
#' outcome probabilities; larger values yield stronger separation.
#'
#' @param n Integer. Number of samples.
#' @param p Integer. Number of predictors.
#' @param prevalence Numeric. Prevalence (probability of class 1).
#' @param mode Character. Cross-validation mode:
#'   `"subject_grouped"`, `"batch_blocked"`, `"study_loocv"`, `"time_series"`.
#' @param learner Character. Base learner: `"glmnet"` or `"ranger"`.
#' @param leakage Character. Leakage type:
#'   `"none"`, `"subject_overlap"`, `"batch_confounded"`, `"peek_norm"`, `"lookahead"`.
#' @param rho Numeric. Autocorrelation coefficient (0 = independent).
#' @param K Integer. Number of folds.
#' @param repeats Integer. Number of repeats.
#' @param horizon Numeric. Time gap for time-series CV.
#' @param B Integer. Number of permutations for `audit_leakage`.
#' @param seeds Integer vector. Random seeds for Monte Carlo replicates.
#' @param parallel Logical. If TRUE, uses `future.apply` for multi-seed execution.
#' @param signal_strength Numeric. Scaling factor for signal-to-noise ratio.
#' @param verbose Logical. If TRUE, prints progress messages for each seed and mode.
#'
#' @return
#' A `LeakSimResults` data frame with columns:
#' - `seed`: simulation seed
#' - `metric_obs`: observed model performance
#' - `gap`: permutation-gap statistic
#' - `p_value`: permutation-based significance
#' - `leakage`: leakage scenario
#' - `mode`: CV mode used
#'
#' @examples
#' \dontrun{
#' res <- simulate_leakage_suite(
#'   n = 300, p = 10, mode = "subject_grouped",
#'   learner = "ranger", leakage = "subject_overlap",
#'   seeds = 1:5, parallel = FALSE, verbose = TRUE
#' )
#' head(res)
#' }
#'
#' @export
simulate_leakage_suite <- function(
    n = 500, p = 20, prevalence = 0.5,
    mode = c("subject_grouped", "batch_blocked", "study_loocv", "time_series"),
    learner = c("glmnet", "ranger"),
    leakage = c("none", "subject_overlap", "batch_confounded", "peek_norm", "lookahead"),
    rho = 0, K = 5, repeats = 1, horizon = 0,
    B = 1000, seeds = 1:10, parallel = FALSE, signal_strength = 1,
    verbose = FALSE
) {
  mode <- match.arg(mode)
  learner <- match.arg(learner)
  leakage <- match.arg(leakage)

  FUN <- function(s) {
    if (verbose)
      message(sprintf("[Seed %d] Mode=%s, Leakage=%s", s, mode, leakage))

    set.seed(s)
    sim <- .simulate_dataset(n, p, prevalence, mode, leakage, rho, signal_strength)

    splits <- make_splits(
      sim$data, outcome = "y", mode = mode,
      group = sim$group_col, batch = sim$batch_col,
      study = sim$study_col, time = sim$time_col,
      v = K, repeats = repeats, stratify = TRUE,
      horizon = horizon, seed = s
    )

    # Automatically select metric based on outcome type
    metrics <- if (is.factor(sim$data$y)) "auc" else "rmse"

    fit <- fit_resample(
      sim$data, outcome = "y", splits = splits,
      learner = learner, metrics = metrics, seed = s
    )

    aud <- audit_leakage(fit, metric = metrics, B = B, seed = s)

    data.frame(
      seed = s,
      metric_obs = aud@permutation_gap$metric_obs,
      gap = aud@permutation_gap$gap,
      p_value = aud@permutation_gap$p_value,
      leakage = leakage,
      mode = mode
    )
  }

  results <- if (parallel) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      warning("parallel=TRUE requires the 'future.apply' package; running sequentially.")
      lapply(seeds, FUN)
    } else {
      future.apply::future_lapply(seeds, FUN, future.seed = TRUE)
    }
  } else {
    lapply(seeds, FUN)
  }

  out <- do.call(rbind, results)
  class(out) <- c("LeakSimResults", class(out))
  out
}


.simulate_dataset <- function(n, p, prevalence, mode, leakage, rho, signal_strength = 1) {
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- sprintf("x%02d", seq_len(p))
  subject <- sample(seq_len(max(5, n %/% 5)), n, replace = TRUE)
  batch <- sample(letters[1:max(3, pmin(6, n %/% 50))], n, replace = TRUE)
  study <- sample(seq_len(max(3, n %/% 80)), n, replace = TRUE)
  time <- seq_len(n)
  if (rho != 0) {
    for (j in seq_len(p)) {
      for (i in 2:n) {
        X[i, j] <- rho * X[i - 1, j] + sqrt(1 - rho^2) * X[i, j]
      }
    }
  }
  linpred <- rowSums(X[, seq_len(min(5, p)), drop = FALSE])
  linpred <- scale(linpred) * signal_strength
  thr <- stats::qnorm(prevalence)
  y_prob <- stats::pnorm(linpred - thr)
  y <- rbinom(n, 1, y_prob)
  if (leakage == "subject_overlap") {
    X <- cbind(X, leak_subj = ave(y, subject, FUN = mean))
  } else if (leakage == "batch_confounded") {
    batch_mean <- ave(y, batch, FUN = mean)
    X <- cbind(X, leak_batch = batch_mean)
  } else if (leakage == "peek_norm") {
    global_mean <- mean(y)
    X <- cbind(X, leak_global = global_mean)
  } else if (leakage == "lookahead") {
    lead <- c(y[-1], y[length(y)])
    X <- cbind(X, leak_future = lead)
  }
  df <- data.frame(X, y = factor(y))
  df$subject <- subject
  df$batch <- batch
  df$study <- study
  df$time <- time
  list(data = df, group_col = "subject", batch_col = "batch",
       study_col = "study", time_col = "time")
}
