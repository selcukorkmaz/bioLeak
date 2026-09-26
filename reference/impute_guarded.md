# Leakage-safe data imputation via guarded preprocessing

Fits imputation parameters on the training data only, then applies the
same guarded transformation to the test data. This function is a thin
wrapper around the guarded preprocessing used by
[`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md).
Output is the transformed feature matrix used by the guarded pipeline
(categorical variables are one-hot encoded).

## Usage

``` r
impute_guarded(
  train,
  test,
  method = c("median", "knn", "missForest", "none"),
  constant_value = 0,
  k = 5,
  seed = 123,
  winsor = TRUE,
  winsor_thresh = 3,
  parallel = FALSE,
  return_outliers = FALSE,
  vars = NULL
)
```

## Arguments

- train:

  data frame (training set)

- test:

  data frame (test set)

- method:

  one of "median", "knn", "missForest", or "none"

- constant_value:

  unused; retained for backward compatibility

- k:

  number of neighbors for kNN imputation (if method = "knn")

- seed:

  unused; retained for backward compatibility. Set seed before calling
  this function if reproducibility is needed.

- winsor:

  logical; apply MAD-based winsorization before imputation

- winsor_thresh:

  numeric; MAD cutoff (default = 3)

- parallel:

  logical; unused (kept for compatibility)

- return_outliers:

  logical; unused (outlier flags not returned)

- vars:

  optional character vector; impute only selected variables

## Value

A list (S3 class `"LeakImpute"`) with elements `train`, `test`, `model`,
`method`, `summary`, and `outliers`.

## See also

\[fit_resample()\], \[predict_guard()\]

## Examples

``` r
train <- data.frame(x = c(1, 2, NA, 4), y = c(NA, 1, 1, 0))
test <- data.frame(x = c(NA, 5), y = c(1, NA))
imp <- impute_guarded(train, test, method = "median", winsor = FALSE)
imp$train
#>   x y
#> 1 1 1
#> 2 2 1
#> 3 2 1
#> 4 4 0
imp$test
#>   x y
#> 1 2 1
#> 2 5 1
```
