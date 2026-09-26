# Apply a fitted GuardFit transformer to new data

Applies the preprocessing steps stored in a `GuardFit` object to new
data without refitting any statistics. This is designed to prevent
validation leakage that would occur if imputation, scaling, filtering,
or feature selection were recomputed on evaluation data. It enforces the
training schema by aligning columns and factor levels, and it errors
when a numeric-only training fit receives non-numeric predictors. It
does not detect label leakage, duplicate samples, or train/test
contamination.

`predict.GuardFit()` is the canonical S3 method — callers can use
`predict(fit, newdata)` on a `GuardFit` object and the right method is
dispatched. `predict_guard()` is retained as a backward-compatible thin
alias that simply forwards to the S3 method, so existing code that calls
`predict_guard(fit, x)` continues to work.

## Usage

``` r
# S3 method for class 'GuardFit'
predict(object, newdata, ...)

predict_guard(fit, newdata)
```

## Arguments

- object, fit:

  A `GuardFit` object created by \[guard_fit()\]. Contains the
  training-time preprocessing settings and statistics. Changing the
  object (for example, a different imputation method or feature
  selection step) changes the output columns and values. `object` is the
  canonical name (matching the S3
  [`predict()`](https://rdrr.io/r/stats/predict.html) generic); `fit` is
  the legacy name accepted only by `predict_guard()`.

- newdata:

  A matrix or data.frame of predictors with one row per sample. This
  required argument (no default) is transformed using the training-time
  parameters in the fit only. Missing columns are added and filled,
  extra columns are dropped, and factor levels are aligned to the
  training levels; if the training fit was numeric-only, non-numeric
  columns in `newdata` trigger an error.

- ...:

  Ignored. Present so that the S3 method signature matches the
  \[stats::predict()\] generic; additional arguments are silently
  dropped.

## Value

A data.frame of transformed predictors with the same number of rows as
`newdata`. Column order and content match the training pipeline and may
include derived features (one-hot encodings, missingness indicators, or
PCA components). This output is not a prediction; it is intended as
input to a downstream model and assumes the training-time preprocessing
is valid for the new data.

## Examples

``` r
x_train <- data.frame(a = c(1, 2, NA, 4), b = c(10, 11, 12, 13))
fit <- guard_fit(
  x_train,
  y = c(0.1, 0.2, 0.3, 0.4),
  steps = list(impute = list(method = "median")),
  task = "gaussian"
)
x_new <- data.frame(a = c(NA, 5), b = c(9, 14))
## Canonical: dispatch through the predict() generic.
out <- predict(fit, x_new)
out
#>            a         b
#> 1 -0.1986799 -1.936492
#> 2  2.1854784  1.936492
## Equivalent legacy form (kept for backward compatibility).
identical(out, predict_guard(fit, x_new))
#> [1] TRUE
```
