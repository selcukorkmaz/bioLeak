# Confounder sensitivity summaries

Computes performance metrics within confounder strata to surface
potential confounding. Requires aligned metadata in \`coldata\`.

## Usage

``` r
confounder_sensitivity(
  fit,
  confounders = NULL,
  metric = NULL,
  min_n = 10,
  coldata = NULL,
  numeric_bins = 4,
  learner = NULL,
  strict_align = FALSE
)
```

## Arguments

- fit:

  A \[LeakFit\] object from \[fit_resample()\].

- confounders:

  Character vector of columns in \`coldata\` to evaluate. Defaults to
  common batch/study identifiers when available.

- metric:

  Metric name to compute within each stratum. Defaults to the first
  metric used in the fit (or task defaults if unavailable).

- min_n:

  Minimum samples per stratum; smaller strata return NA metrics.

- coldata:

  Optional data.frame of sample metadata. Defaults to
  \`fit@splits@info\$coldata\` when available.

- numeric_bins:

  Integer number of quantile bins for numeric confounders with many
  unique values.

- learner:

  Optional character scalar. When predictions include multiple learners,
  selects the learner to summarize.

- strict_align:

  Logical scalar. If TRUE, errors when coldata cannot be aligned by row
  names or IDs and would fall back to row-order matching. Default is
  FALSE.

## Value

A data.frame with per-confounder, per-level metrics and counts.

## Examples

``` r
set.seed(42)
df <- data.frame(
  subject = rep(1:15, each = 2),
  outcome = factor(rep(c(0, 1), 15)),
  batch = factor(rep(c("A", "B", "C"), 10)),
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
confounder_sensitivity(fit, confounders = "batch", coldata = df)
#>   confounder level metric direction value  n positive_rate
#> 1      batch     A    auc    higher  0.08 10           0.5
#> 2      batch     B    auc    higher  0.24 10           0.5
#> 3      batch     C    auc    higher  0.28 10           0.5
```
