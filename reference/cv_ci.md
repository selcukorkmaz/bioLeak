# Confidence intervals for cross-validated metrics

Computes per-learner confidence intervals for each metric column in a
per-fold metrics data.frame. Supports the standard normal/t approach and
the Nadeau-Bengio (2003) corrected variance for repeated K-fold CV.

## Usage

``` r
cv_ci(
  metrics_df,
  level = 0.95,
  method = c("nadeau_bengio", "normal"),
  n_train = NULL,
  n_test = NULL,
  folds_per_repeat = NULL
)
```

## Arguments

- metrics_df:

  A
  [`LeakFit`](https://selcukorkmaz.github.io/bioLeak/reference/LeakClasses.md)
  or a `LeakTune` (recommended: fold sizes, folds per repeat and
  repeated test sets are then taken from the object's splits), or a
  data.frame with columns `fold`, `learner`, and one or more numeric
  metric columns, such as `fit@metrics`. An optional character column
  `test_set` identifies each fold's test set: folds sharing a test set
  (repeats that re-use the same partition) are averaged and counted
  once.

- level:

  Confidence level (default 0.95).

- method:

  `"nadeau_bengio"` (default) or `"normal"`. `"normal"` treats folds as
  independent and is far too narrow for cross-validation (about 37%
  coverage of a nominal 95% interval for repeated 5-fold CV in
  simulation); it is kept for comparison.

- n_train:

  Average number of training samples per fold. Used only when
  `method = "nadeau_bengio"`. NULL to use fallback variance.

- n_test:

  Average number of test samples per fold. Used only when
  `method = "nadeau_bengio"`. Without `n_train` and `n_test` the
  correction cannot be applied: a warning is issued and the uncorrected
  variance `var / K` is used.

- folds_per_repeat:

  Number of folds in one CV repeat (e.g. 5 for 5-fold CV repeated 10
  times). Used only when `method = "nadeau_bengio"` to set the degrees
  of freedom. The resulting intervals are slightly conservative in
  simulation (about 98% coverage at the 95% level for repeated 5-fold
  CV). Intervals are not clipped to a metric's range (e.g. an AUC bound
  can exceed 1).

## Value

A data.frame with `learner` and, for each metric, columns
`<metric>_mean`, `<metric>_sd`, `<metric>_ci_lo`, and `<metric>_ci_hi`.

## Details

Fold metrics from cross-validation are correlated: training sets
overlap, and repeats re-partition the same data. `"normal"` ignores this
(`sd / sqrt(K)` over all folds) and is too narrow for repeated CV.
`"nadeau_bengio"` inflates the variance to
`(1/K + n_test/n_train) * var` and uses a t reference on
`folds_per_repeat - 1` degrees of freedom (calibrated degrees of
freedom, Bouckaert 2003): repeats sharpen the estimate but add no
independent information. Without `folds_per_repeat` all folds are
treated as one repeat.
