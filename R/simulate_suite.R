# Simulation framework for leakage diagnostics ---------------------------------

#' Simulate leakage scenarios and audit results
#' @export
simulate_leakage_suite <- function(
  n = 500,
  p = 20,
  prevalence = 0.5,
  mode = c("subject_grouped", "batch_blocked", "study_loocv", "time_series"),
  learner = c("glmnet", "ranger", "xgboost"),
  leakage = c("none", "subject_overlap", "batch_confounded", "peek_norm", "lookahead"),
  rho = 0,
  K = 5, repeats = 1,
  horizon = 0,
  B = 1000,
  seeds = 1:10,
  parallel = FALSE
) {
  mode <- match.arg(mode)
  learner <- match.arg(learner)
  leakage <- match.arg(leakage)
  seeds <- as.integer(seeds)
  results <- lapply(seeds, function(s) {
    set.seed(s)
    sim <- .simulate_dataset(n, p, prevalence, mode, leakage, rho)
    splits <- make_splits(sim$data, outcome = "y", mode = mode,
                          group = sim$group_col, batch = sim$batch_col,
                          study = sim$study_col, time = sim$time_col,
                          v = K, repeats = repeats, stratify = TRUE,
                          horizon = horizon, seed = s)
    fit <- fit_resample(sim$data, outcome = "y", splits = splits,
                        learner = learner, metrics = c("auc"), seed = s)
    aud <- audit_leakage(fit, metric = "auc", B = B, seed = s)
    data.frame(seed = s,
               metric_obs = aud@permutation_gap$metric_obs,
               gap = aud@permutation_gap$gap,
               p_value = aud@permutation_gap$p_value,
               leakage = leakage,
               mode = mode)
  })
  do.call(rbind, results)
}

.simulate_dataset <- function(n, p, prevalence, mode, leakage, rho) {
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
  linpred <- scale(linpred)
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
