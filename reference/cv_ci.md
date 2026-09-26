# Confidence intervals for cross-validated metrics

Computes per-learner confidence intervals for each metric column in a
per-fold metrics data.frame. Supports the standard normal/t approach and
the Nadeau-Bengio (2003) corrected variance for repeated K-fold CV.

## Usage

``` r
cv_ci(
  metrics_df,
  level = 0.95,
  method = c("normal", "nadeau_bengio"),
  n_train = NULL,
  n_test = NULL
)
```

## Arguments

- metrics_df:

  Data.frame with columns `fold`, `learner`, and one or more numeric
  metric columns.

- level:

  Confidence level (default 0.95).

- method:

  One of `"normal"` or `"nadeau_bengio"`.

- n_train:

  Average number of training samples per fold. Used only when
  `method = "nadeau_bengio"`. NULL to use fallback variance.

- n_test:

  Average number of test samples per fold. Used only when
  `method = "nadeau_bengio"`. NULL to use fallback variance.

## Value

A data.frame with `learner` and, for each metric, columns
`<metric>_mean`, `<metric>_sd`, `<metric>_ci_lo`, and `<metric>_ci_hi`.
