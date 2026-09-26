# Display summary for LeakFit objects

Prints a brief one-screen summary of a LeakFit, including task and
outcome, fold count and status (successful, skipped, failed), and the
headline cross-validated metric. Use
[`summary()`](https://rdrr.io/r/base/summary.html) for the full per-fold
diagnostic report.

## Usage

``` r
# S4 method for class 'LeakFit'
show(object)
```

## Arguments

- object:

  A `LeakFit` object.

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
show(fit)
#> A LeakFit object
#>   Task:           binomial
#>   Outcome:        outcome
#>   Learners:       glm
#>   Folds:          3  (mode = subject_grouped)
#>   Fold status:    3 success, 0 skipped, 0 failed
#> Use summary(<obj>) for the full diagnostic report.
```
