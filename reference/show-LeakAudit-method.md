# Display summary for LeakAudit objects

Prints a brief one-screen summary of a LeakAudit, including task and
outcome, the permutation-gap statistic, and counts of batch-association
rows, target-leakage features, and duplicate pairs. Use
[`summary()`](https://rdrr.io/r/base/summary.html) for the full
diagnostic report.

## Usage

``` r
# S4 method for class 'LeakAudit'
show(object)
```

## Arguments

- object:

  A `LeakAudit` object.

## Value

No return value, called for side effects (prints a brief summary to the
console). Returns `object` invisibly.

## Examples

``` r
set.seed(1)
df <- data.frame(
  subject = rep(1:6, each = 2),
  outcome = rbinom(12, 1, 0.5),
  x1 = rnorm(12),
  x2 = rnorm(12)
)
splits <- make_split_plan(df, outcome = "outcome",
                      mode = "subject_grouped", group = "subject", v = 3,
                      progress = FALSE)
custom <- list(
  glm = list(
    fit = function(x, y, task, weights, ...) {
      stats::glm(y ~ ., data = as.data.frame(x),
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
                    metrics = "auc", refit = FALSE, seed = 1)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#>   |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
aud <- audit_leakage(fit, metric = "auc", B = 10,
                     X_ref = df[, c("x1", "x2")])
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#>   |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#>   |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#>   |                                                                              |=======================                                               |  33%
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#>   |                                                                              |===============================================                       |  67%
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#>   |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#>   |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#>   |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#>   |                                                                              |======================================================================| 100%
show(aud)
#> A LeakAudit object
#>   Task:               binomial
#>   Outcome:            outcome
#>   Permutation-gap:    metric=0.611, gap=0.169, p=0.273
#>   Batch association:  0 row(s)
#>   Target leakage:     2 feature(s)
#>   Duplicate pairs:    1
#> Use summary(<obj>) for the full diagnostic report.
```
