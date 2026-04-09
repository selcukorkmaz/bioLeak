## Smoke test: run 1 seed for a representative sample of configs
## Run in RStudio before committing to the full simulation

library(bioLeak)
suppressMessages(library(glmnet))

leakage_types <- c("none", "subject_overlap", "batch_confounded",
                   "peek_norm", "lookahead")
## Test edge cases: smallest and largest n/p/s
test_configs <- expand.grid(
  leakage = leakage_types,
  n = c(100, 1000),
  p = c(10, 100),
  s = c(0, 0.5, 2.0),
  stringsAsFactors = FALSE
)

cat(sprintf("Testing %d configs (1 seed each, B=10)...\n\n", nrow(test_configs)))

results <- data.frame()
for (i in seq_len(nrow(test_configs))) {
  cfg <- test_configs[i, ]
  t0 <- proc.time()
  cat(sprintf("[%2d/%d] leakage=%-18s n=%-5d p=%-4d s=%.1f ... ",
              i, nrow(test_configs), cfg$leakage, cfg$n, cfg$p, cfg$s))

  res <- tryCatch({
    set.seed(42)
    n <- cfg$n; p <- cfg$p; s <- cfg$s; leakage_type <- cfg$leakage

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

    if (leakage_type == "subject_overlap") {
      X <- cbind(X, leak_subj = ave(y, subject, FUN = mean))
    } else if (leakage_type == "batch_confounded") {
      ## Make batch outcome-dependent, then leak the group mean
      batch <- ifelse(y == 1,
                      sample(c("a","b","c"), n, replace = TRUE, prob = c(0.6, 0.2, 0.2)),
                      sample(c("a","b","c"), n, replace = TRUE, prob = c(0.15, 0.5, 0.35)))
      X <- cbind(X, leak_batch = ave(y, batch, FUN = mean))
    } else if (leakage_type == "peek_norm") {
      ## Continuous feature encoding y via global statistics (simulates
      ## normalisation using full dataset incl. test fold)
      X <- cbind(X, leak_global = as.numeric(y) + rnorm(n, 0, 0.3))
    } else if (leakage_type == "lookahead") {
      ## Continuous biomarker (noisy proxy of latent process);
      ## lookahead = next time point's measurement (future information)
      biomarker <- linpred + rnorm(n, 0, 0.5)
      X <- cbind(X, leak_future = c(biomarker[-1], biomarker[n]))
    }

    df <- data.frame(X, y = factor(y))
    df$subject <- subject; df$batch <- batch
    df$study <- study; df$time <- time_var

    splits <- make_split_plan(
      df, outcome = "y", mode = "subject_grouped",
      group = "subject", batch = "batch", study = "study", time = "time",
      v = 5, stratify = TRUE, seed = 42,
      progress = FALSE
    )

    preprocess <- list(
      impute = list(method = "median"), normalize = list(method = "zscore"),
      filter = list(var_thresh = 0, iqr_thresh = 0), fs = list(method = "none")
    )

    fit <- fit_resample(
      df, outcome = "y", splits = splits,
      preprocess = preprocess,
      learner = "glmnet", metrics = "auc", seed = 42
    )

    aud <- audit_leakage(
      fit, metric = "auc", B = 10,
      perm_refit = FALSE, perm_stratify = TRUE,
      seed = 42, target_scan = FALSE, return_perm = FALSE
    )

    data.frame(
      leakage = leakage_type, n = n, p = p, s = s,
      metric_obs = aud@permutation_gap$metric_obs,
      gap = aud@permutation_gap$gap,
      p_value = aud@permutation_gap$p_value,
      status = "OK",
      error = "",
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(
      leakage = cfg$leakage, n = cfg$n, p = cfg$p, s = cfg$s,
      metric_obs = NA, gap = NA, p_value = NA,
      status = "FAIL",
      error = conditionMessage(e),
      stringsAsFactors = FALSE
    )
  })

  elapsed <- (proc.time() - t0)[3]
  cat(sprintf("%s (%.1fs)", res$status, elapsed))
  if (res$status == "FAIL") cat(sprintf(" -- %s", res$error))
  cat("\n")
  results <- rbind(results, res)
}

cat("\n=== SUMMARY ===\n")
cat(sprintf("Passed: %d / %d\n", sum(results$status == "OK"), nrow(results)))
cat(sprintf("Failed: %d / %d\n", sum(results$status == "FAIL"), nrow(results)))

if (any(results$status == "FAIL")) {
  cat("\nFailed configs:\n")
  print(results[results$status == "FAIL", c("leakage", "n", "p", "s", "error")])
} else {
  cat("\nAll configs passed! Safe to run the full simulation.\n")
}
