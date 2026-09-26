# Calibration diagnostics for binomial predictions

Computes reliability curve summaries and calibration metrics for a
binomial \[LeakFit\] using out-of-fold predictions.

## Usage

``` r
calibration_summary(fit, bins = 10, min_bin_n = 5, learner = NULL)
```

## Arguments

- fit:

  A \[LeakFit\] object from \[fit_resample()\].

- bins:

  Integer number of probability bins for the calibration curve.

- min_bin_n:

  Minimum samples per bin used in plotting; bins smaller than this are
  retained in the output but can be filtered by the caller.

- learner:

  Optional character scalar. When predictions include multiple learners,
  selects the learner to summarize.

## Value

A list with a \`curve\` data.frame and a one-row \`metrics\` data.frame
containing ECE, MCE, and Brier score.

## Examples

``` r
set.seed(42)
df <- data.frame(
  subject = rep(1:15, each = 2),
  outcome = factor(rep(c(0, 1), 15)),
  x1 = rnorm(30),
  x2 = rnorm(30)
)
splits <- make_split_plan(df, outcome = "outcome",
                          mode = "subject_grouped", group = "subject",
                          v = 3, progress = FALSE)
custom <- list(
  glm = list(
    fit = function(x, y, task, weights, ...) {
      stats::glm(y ~ ., data = as.data.frame(x),
                 family = stats::binomial(), weights = weights)
    },
    predict = function(object, newdata, task, ...) {
      as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
                                type = "response"))
    }
  )
)
fit <- fit_resample(df, outcome = "outcome", splits = splits,
                    learner = "glm", custom_learners = custom,
                    metrics = "auc", refit = FALSE, seed = 1)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
cal <- calibration_summary(fit, bins = 5)
cal$metrics
#>    n bins min_bin_n       ece       mce     brier
#> 1 30    5         5 0.3569617 0.9200306 0.3951334
```
