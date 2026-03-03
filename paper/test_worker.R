## Quick test: can workers see bioLeak functions?
library(bioLeak)
library(future)
library(future.apply)
plan(multisession, workers = 2)

## Test 1: direct call
cat("Direct test:\n")
set.seed(1)
n <- 100; p <- 10; s <- 1.0
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- sprintf("x%02d", seq_len(p))
subject <- sample(1:20, n, replace = TRUE)
linpred <- rowSums(X[, 1:5])
linpred <- scale(linpred) * s
y <- rbinom(n, 1, pnorm(linpred - qnorm(0.5)))
df <- data.frame(X, y = factor(y), subject = subject)

splits <- make_split_plan(df, outcome = "y", mode = "subject_grouped",
                          group = "subject", v = 5, stratify = TRUE, seed = 1)
fit <- fit_resample(df, outcome = "y", splits = splits,
                    learner = "glmnet", metrics = "auc", seed = 1)
aud <- audit_leakage(fit, metric = "auc", B = 50, perm_refit = FALSE,
                     perm_stratify = TRUE, seed = 1, target_scan = FALSE,
                     return_perm = FALSE)
cat(sprintf("AUC=%.3f, gap=%.3f, p=%.4f\n",
            aud@permutation_gap$metric_obs,
            aud@permutation_gap$gap,
            aud@permutation_gap$p_value))

## Test 2: worker call
cat("\nWorker test:\n")
res <- future_lapply(1:2, function(seed) {
  tryCatch({
    library(bioLeak)
    set.seed(seed)
    n <- 100; p <- 10; s <- 1.0
    X <- matrix(rnorm(n * p), n, p)
    colnames(X) <- sprintf("x%02d", seq_len(p))
    subject <- sample(1:20, n, replace = TRUE)
    linpred <- rowSums(X[, 1:5])
    linpred <- scale(linpred) * s
    y <- rbinom(n, 1, pnorm(linpred - qnorm(0.5)))
    df <- data.frame(X, y = factor(y), subject = subject)
    splits <- make_split_plan(df, outcome = "y", mode = "subject_grouped",
                              group = "subject", v = 5, stratify = TRUE, seed = seed)
    fit <- fit_resample(df, outcome = "y", splits = splits,
                        learner = "glmnet", metrics = "auc", seed = seed)
    aud <- audit_leakage(fit, metric = "auc", B = 50, perm_refit = FALSE,
                         perm_stratify = TRUE, seed = seed, target_scan = FALSE,
                         return_perm = FALSE)
    data.frame(seed = seed, auc = aud@permutation_gap$metric_obs,
               gap = aud@permutation_gap$gap, p = aud@permutation_gap$p_value)
  }, error = function(e) {
    data.frame(seed = seed, auc = NA, gap = NA, p = NA, error = conditionMessage(e))
  })
}, future.seed = TRUE)

print(do.call(rbind, res))
cat("Done.\n")
