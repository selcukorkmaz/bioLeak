#' Simulate leakage scenarios and audit results
#'
#' @description
#' Simulates synthetic binary classification datasets with optional leakage
#' mechanisms, fits a model using a leakage-aware cross-validation scheme, and
#' summarizes the permutation-gap audit for each Monte Carlo seed. The suite is
#' designed to surface validation failures such as subject overlap across folds,
#' batch-confounded outcomes, global normalization/summary leakage, and
#' time-series look-ahead. The output is a per-seed summary of observed CV
#' performance and its gap versus a label-permutation null; it does not return
#' fitted models or the full audit object. Results are limited to the built-in
#' data generator and leakage types implemented here, and should be interpreted
#' as a simulation-based sanity check rather than a comprehensive leakage
#' detector for real data.
#'
#' @details
#' The generator draws \code{p} standard normal predictors, builds a linear
#' predictor from the first \code{min(5, p)} features, scales it by
#' \code{signal_strength}, and samples a binary outcome to achieve the requested
#' \code{prevalence}. Outcomes are returned as a two-level factor, so the audited
#' metric is AUC. Simulated metadata include subject, batch, study, and time
#' fields used by \code{mode} to create leakage-aware splits. Leakage mechanisms
#' are injected by adding a single extra predictor as described in
#' \code{leakage}. Parallel execution uses \code{future.apply} when installed and
#' does not change results.
#'
#' @param n Integer scalar. Number of samples to simulate (default 500). Larger
#'   values stabilize the Monte Carlo summary but increase runtime.
#' @param p Integer scalar. Number of baseline predictors before any leakage
#'   feature is added (default 20). Increasing \code{p} changes the signal-to-noise
#'   ratio and increases fitting time.
#' @param prevalence Numeric scalar in (0, 1). Target prevalence of class 1 in
#'   the simulated outcome (default 0.5). Changing this alters class imbalance
#'   and can affect AUC and the permutation gap.
#' @param mode Character scalar. Cross-validation scheme passed to
#'   \code{make_split_plan()}; one of \code{"subject_grouped"},
#'   \code{"batch_blocked"}, \code{"study_loocv"}, \code{"time_series"}.
#'   Defaults to \code{"subject_grouped"}. This controls how samples are grouped
#'   into folds (by subject, batch, study, or time) and therefore which leakage
#'   mechanisms are realistically challenged.
#' @param learner Character scalar. Base learner, \code{"glmnet"} (default) or
#'   \code{"ranger"}. Requires the corresponding package in \code{Suggests}.
#'   Switching learners changes the fitted model, runtime, and performance.
#' @param leakage Character scalar. Leakage mechanism to inject; one of
#'   \code{"none"}, \code{"subject_overlap"}, \code{"batch_confounded"},
#'   \code{"peek_norm"}, \code{"lookahead"}. Leakage is added as an extra
#'   predictor: \code{"subject_overlap"} adds per-subject mean outcome,
#'   \code{"batch_confounded"} adds per-batch mean outcome, \code{"peek_norm"}
#'   adds the globally normalized (z-scored) outcome, and \code{"lookahead"} adds the next-time
#'   outcome. Changing this controls whether and how leakage is present.
#' @param preprocess Optional preprocessing list or recipe passed to
#'   [fit_resample()]. When NULL (default), the simulator uses the
#'   fit_resample defaults; for \code{"peek_norm"} leakage, normalization is
#'   set to \code{"none"} to avoid attenuating the constant leakage feature.
#' @param rho Numeric scalar in [-1, 1]. AR(1)-style autocorrelation applied to
#'   each predictor across row order (default 0). Higher absolute values increase
#'   serial correlation and make time-ordered leakage more pronounced.
#' @param K Integer scalar. Number of folds/partitions (default 5). Used as the
#'   fold count for \code{"subject_grouped"} and \code{"batch_blocked"}, and as
#'   the number of rolling partitions for \code{"time_series"}. Ignored for
#'   \code{"study_loocv"} (folds equal the number of studies).
#' @param repeats Integer scalar >= 1. Number of repeated CV runs for
#'   \code{"subject_grouped"} and \code{"batch_blocked"} (default 1). Increasing
#'   \code{repeats} increases the number of folds and runtime. Ignored for
#'   \code{"study_loocv"} and \code{"time_series"}.
#' @param horizon Numeric scalar >= 0. Minimum time gap enforced between train
#'   and test for \code{"time_series"} splits (default 0). Larger values make the
#'   split more conservative and can reduce leakage from temporal proximity.
#' @param B Integer scalar >= 1. Number of permutations used by
#'   \code{audit_leakage()} to compute the permutation gap and p-value (default
#'   1000). Larger values yield more stable p-values but increase runtime.
#' @param seeds Integer vector. Monte Carlo seeds (default \code{1:10}). One row
#'   of output is produced per seed; changing \code{seeds} changes the simulated
#'   datasets and splits.
#' @param parallel Logical scalar. If \code{TRUE}, evaluates seeds in parallel
#'   using \code{future.apply} (if installed). Results are identical to sequential
#'   execution; only runtime changes.
#' @param signal_strength Numeric scalar. Scales the linear predictor before
#'   sampling outcomes (default 1). Larger values increase class separation and
#'   tend to increase AUC; smaller values make the task harder.
#' @param verbose Logical scalar. If \code{TRUE}, prints progress messages for
#'   each seed. Does not affect results.
#'
#' @return
#' A \code{LeakSimResults} data frame with one row per seed and columns:
#' \itemize{
#'   \item \code{seed}: seed used for data generation, splitting, and auditing.
#'   \item \code{metric_obs}: observed CV performance (AUC for this simulation).
#'   \item \code{gap}: permutation-gap statistic (observed minus permutation mean).
#'   \item \code{p_value}: permutation p-value for the gap.
#'   \item \code{leakage}: leakage scenario used.
#'   \item \code{mode}: CV mode used.
#' }
#' Only the permutation-gap summary is returned; fitted models, predictions, and
#' other audit components are not included.
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("glmnet", quietly = TRUE)) {
#'   set.seed(1)
#'   res <- simulate_leakage_suite(
#'     n = 120, p = 6, prevalence = 0.4,
#'     mode = "subject_grouped",
#'     learner = "glmnet",
#'     leakage = "subject_overlap",
#'     K = 3, repeats = 1,
#'     B = 50, seeds = 1,
#'     parallel = FALSE
#'   )
#'   # One row per seed with observed AUC, permutation gap, and p-value
#'   res
#' }
#'
#' }
#' @export
simulate_leakage_suite <- function(
    n = 500, p = 20, prevalence = 0.5,
    mode = c("subject_grouped", "batch_blocked", "study_loocv", "time_series"),
    learner = c("glmnet", "ranger"),
    leakage = c("none", "subject_overlap", "batch_confounded", "peek_norm", "lookahead"),
    preprocess = NULL,
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

    splits <- make_split_plan(
      sim$data, outcome = "y", mode = mode,
      group = sim$group_col, batch = sim$batch_col,
      study = sim$study_col, time = sim$time_col,
      v = K, repeats = repeats, stratify = TRUE,
      horizon = horizon, seed = s
    )

    # Automatically select metric based on outcome type
    metrics <- if (is.factor(sim$data$y)) {
      if (nlevels(sim$data$y) > 2) "accuracy" else "auc"
    } else {
      "rmse"
    }

    preprocess_use <- preprocess
    if (is.null(preprocess_use) && leakage == "peek_norm") {
      preprocess_use <- list(
        impute = list(method = "median"),
        normalize = list(method = "none"),
        filter = list(var_thresh = 0, iqr_thresh = 0),
        fs = list(method = "none")
      )
    }

    fit <- fit_resample(
      sim$data, outcome = "y", splits = splits,
      preprocess = preprocess_use,
      learner = learner, metrics = metrics, seed = s
    )

    aud <- audit_leakage(
      fit,
      metric = metrics,
      B = B,
      perm_refit = "auto",
      perm_refit_auto_max = B,
      seed = s
    )

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
    global_sd <- stats::sd(y)
    if (!is.finite(global_sd) || global_sd == 0) global_sd <- 1
    leak_global <- (y - global_mean) / global_sd
    X <- cbind(X, leak_global = leak_global)
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
