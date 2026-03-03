## Diagnose why future workers fail
## Run this in RStudio

library(bioLeak)
library(future)
library(future.apply)
plan(multisession, workers = 2)

## Test 1: Can workers load the package?
cat("=== Test 1: pkgload::load_all in worker ===\n")
res1 <- future::future({
  tryCatch({
    library(bioLeak)
    "load_all OK"
  }, error = function(e) paste("FAIL:", conditionMessage(e)))
}, seed = TRUE)
cat(future::value(res1), "\n")

## Test 2: Can workers run the full pipeline?
cat("\n=== Test 2: Full pipeline in worker ===\n")
res2 <- future::future({
  tryCatch({
    library(bioLeak)
    suppressMessages(library(glmnet))
    set.seed(1)
    n <- 100; p <- 10; s <- 1.0
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- sprintf("x%02d", seq_len(p))
    subject <- sample(seq_len(20), n, replace = TRUE)
    batch <- sample(letters[1:3], n, replace = TRUE)
    study <- sample(1:3, n, replace = TRUE)
    linpred <- rowSums(X[, 1:5, drop = FALSE])
    linpred <- scale(linpred) * s
    y <- rbinom(n, 1, pnorm(linpred - qnorm(0.5)))
    X <- cbind(X, leak_batch = ave(y, batch, FUN = mean))
    df <- data.frame(X, y = factor(y))
    df$subject <- subject; df$batch <- batch
    df$study <- study; df$time <- seq_len(n)

    splits <- make_split_plan(df, outcome = "y", mode = "subject_grouped",
                              group = "subject", v = 5, stratify = TRUE, seed = 1,
                              progress = FALSE)
    fit <- fit_resample(df, outcome = "y", splits = splits,
                        preprocess = list(impute = list(method = "median"),
                                          normalize = list(method = "zscore"),
                                          filter = list(var_thresh = 0, iqr_thresh = 0),
                                          fs = list(method = "none")),
                        learner = "glmnet", metrics = "auc", seed = 1)
    aud <- audit_leakage(fit, metric = "auc", B = 10,
                         perm_refit = FALSE, perm_stratify = TRUE,
                         seed = 1, target_scan = FALSE, return_perm = FALSE)
    paste("OK! metric_obs =", aud@permutation_gap$metric_obs)
  }, error = function(e) paste("FAIL:", conditionMessage(e)))
}, seed = TRUE)
cat(future::value(res2), "\n")

## Test 3: Same but with future_lapply (closer to actual sim)
cat("\n=== Test 3: future_lapply with 2 seeds ===\n")
results <- future_lapply(1:2, function(seed) {
  tryCatch({
    library(bioLeak)
    suppressMessages(library(glmnet))
    set.seed(seed)
    n <- 100; p <- 10; s <- 1.0
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- sprintf("x%02d", seq_len(p))
    subject <- sample(seq_len(20), n, replace = TRUE)
    batch <- sample(letters[1:3], n, replace = TRUE)
    linpred <- rowSums(X[, 1:5, drop = FALSE])
    linpred <- scale(linpred) * s
    y <- rbinom(n, 1, pnorm(linpred - qnorm(0.5)))
    X <- cbind(X, leak_batch = ave(y, batch, FUN = mean))
    df <- data.frame(X, y = factor(y))
    df$subject <- subject; df$batch <- batch
    df$study <- sample(1:3, n, replace = TRUE); df$time <- seq_len(n)

    splits <- make_split_plan(df, outcome = "y", mode = "subject_grouped",
                              group = "subject", v = 5, stratify = TRUE, seed = seed,
                              progress = FALSE)
    fit <- fit_resample(df, outcome = "y", splits = splits,
                        preprocess = list(impute = list(method = "median"),
                                          normalize = list(method = "zscore"),
                                          filter = list(var_thresh = 0, iqr_thresh = 0),
                                          fs = list(method = "none")),
                        learner = "glmnet", metrics = "auc", seed = seed)
    aud <- audit_leakage(fit, metric = "auc", B = 10,
                         perm_refit = FALSE, perm_stratify = TRUE,
                         seed = seed, target_scan = FALSE, return_perm = FALSE)
    data.frame(seed = seed, metric_obs = aud@permutation_gap$metric_obs)
  }, error = function(e) {
    data.frame(seed = seed, metric_obs = NA, error = conditionMessage(e))
  })
}, future.seed = TRUE)

print(do.call(rbind, results))

plan(sequential)
cat("\n=== Done ===\n")
