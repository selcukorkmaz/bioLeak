## Quick diagnostic: check which configs have NAs in sim_results_all.rds
sim <- readRDS("paper/sim_results_all.rds")

cat("Total rows:", nrow(sim), "\n")
cat("Total NAs:", sum(is.na(sim$metric_obs)), "\n\n")

## NA count by leakage x n x p x s
na_tab <- aggregate(metric_obs ~ leakage + n + p + s, data = sim,
                    FUN = function(x) sum(is.na(x)))
names(na_tab)[5] <- "n_na"
na_tab <- na_tab[na_tab$n_na > 0, ]
cat("Configs with NAs (", nrow(na_tab), " configs):\n", sep = "")
print(na_tab)

## Also try reproducing one failure with actual error message
cat("\n=== Reproducing one failure to get error message ===\n")
library(bioLeak)
set.seed(1)
n <- 100; p <- 10; s <- 1.0
X <- matrix(rnorm(n * p), n, p)
colnames(X) <- sprintf("x%02d", seq_len(p))
subject <- sample(seq_len(max(5, n %/% 5)), n, replace = TRUE)
batch <- sample(letters[1:max(3, min(6, n %/% 50))], n, replace = TRUE)
study <- sample(seq_len(max(3, n %/% 80)), n, replace = TRUE)
time_var <- seq_len(n)

linpred <- rowSums(X[, seq_len(min(5, p)), drop = FALSE])
linpred <- scale(linpred) * s
y_prob <- pnorm(linpred - qnorm(0.5))
y <- rbinom(n, 1, y_prob)

X <- cbind(X, leak_batch = ave(y, batch, FUN = mean))

df <- data.frame(X, y = factor(y))
df$subject <- subject; df$batch <- batch
df$study <- study; df$time <- time_var

tryCatch({
  splits <- make_split_plan(
    df, outcome = "y", mode = "subject_grouped",
    group = "subject", v = 5, stratify = TRUE, seed = 1
  )
  cat("make_split_plan: OK\n")

  fit <- fit_resample(
    df, outcome = "y", splits = splits,
    preprocess = list(
      impute = list(method = "median"),
      normalize = list(method = "zscore"),
      filter = list(var_thresh = 0, iqr_thresh = 0),
      fs = list(method = "none")
    ),
    learner = "glmnet", metrics = "auc", seed = 1
  )
  cat("fit_resample: OK\n")

  aud <- audit_leakage(
    fit, metric = "auc", B = 200,
    perm_refit = FALSE, perm_stratify = TRUE,
    seed = 1, target_scan = FALSE, return_perm = FALSE
  )
  cat("audit_leakage: OK\n")
  cat("metric_obs:", aud@permutation_gap$metric_obs, "\n")
}, error = function(e) {
  cat("ERROR:", conditionMessage(e), "\n")
  cat("Call:", deparse(conditionCall(e)), "\n")
})
