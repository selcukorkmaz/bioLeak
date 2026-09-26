# Summarize a LeakFit object

Prints a compact console report for a \[LeakFit\] object created by
\[fit_resample()\]. The report lists task/outcome metadata, learners,
total folds, and cross-validated metrics summarized as mean and standard
deviation across completed folds, plus a small audit table with per-fold
train/test sizes and retained feature counts.

## Usage

``` r
# S3 method for class 'LeakFit'
summary(object, digits = 3, ...)
```

## Arguments

- object:

  A \[LeakFit\] object returned by \[fit_resample()\]. It should contain
  \`metric_summary\` and \`audit\` slots; missing entries result in
  empty sections in the printed report.

- digits:

  Integer scalar. Number of decimal places to print in numeric summary
  tables. Defaults to 3; affects printed output only, not the returned
  data.

- ...:

  Unused. Included for S3 method compatibility; changing these values
  has no effect.

## Value

Invisibly returns \`object@metric_summary\`, a data frame of per-learner
metric means and standard deviations computed across folds. This
function does not recompute metrics.

## Details

This summary is meant for quick sanity checks of the resampling setup
and performance. It does not run leakage diagnostics and will not detect
target leakage, duplicate samples, or batch/study confounding; use
\[audit_leakage()\] or \`summary()\` on a \[LeakAudit\] object for those
checks.

## Examples

``` r
set.seed(1)
df <- data.frame(
  subject = rep(1:6, each = 2),
  outcome = factor(rep(c(0, 1), each = 6)),
  x1 = rnorm(12),
  x2 = rnorm(12)
)
splits <- make_split_plan(
  df,
  outcome = "outcome",
  mode = "subject_grouped",
  group = "subject",
  v = 3,
  stratify = TRUE,
  progress = FALSE
)
custom <- list(
  glm = list(
    fit = function(x, y, task, weights, ...) {
      stats::glm(y ~ ., data = data.frame(y = y, x),
                 family = stats::binomial(), weights = weights)
    },
    predict = function(object, newdata, task, ...) {
      as.numeric(stats::predict(object,
                                newdata = as.data.frame(newdata),
                                type = "response"))
    }
  )
)
fit <- fit_resample(df, outcome = "outcome", splits = splits,
                    learner = "glm", custom_learners = custom,
                    metrics = "auc", seed = 1)
#>   |                                                                              |                                                                      |   0%
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#>   |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
summary_df <- summary(fit)
#> 
#> ===========================
#>  bioLeak Model Fit Summary
#> ===========================
#> 
#> Task: binomial
#> Outcome: outcome
#> Positive class: 1
#> Learners: glm
#> Total folds: 3
#> Fold status: 3 success, 0 skipped, 0 failed
#> Refit performed: Yes
#> Hash: h91c169ae
#> 
#> Cross-validated metrics (mean ± SD):
#>   learner auc_mean auc_sd auc_ci_lo auc_ci_hi
#> 1     glm    0.375  0.375    -1.098     1.848
#> 
#> Audit overview:
#>  fold n_train n_test learner features_final
#>     1       8      4     glm              2
#>     2       8      4     glm              2
#>     3       8      4     glm              2
#> 
summary_df
#>   learner auc_mean auc_sd auc_ci_lo auc_ci_hi
#> 1     glm    0.375  0.375 -1.097912  1.847912
```
