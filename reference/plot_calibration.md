# Plot calibration curve for binomial predictions

Visualizes observed outcome rates versus predicted probabilities across
bins to diagnose calibration (binomial tasks only). Requires ggplot2.

## Usage

``` r
plot_calibration(fit, bins = 10, min_bin_n = 5, learner = NULL)
```

## Arguments

- fit:

  LeakFit.

- bins:

  Number of probability bins to use.

- min_bin_n:

  Minimum samples per bin shown in the plot.

- learner:

  Optional character scalar. When predictions include multiple learners,
  selects the learner to summarize.

## Value

A list containing the calibration curve, metrics, and a ggplot object.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
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
  plot_calibration(fit, bins = 5)
}
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
```
