# Fit leakage-safe preprocessing pipeline

Builds and fits a guarded preprocessing pipeline on training data, then
constructs a transformer for consistent application to new data.

## Usage

``` r
guard_fit(
  X,
  y = NULL,
  steps = list(),
  task = c("binomial", "multiclass", "gaussian", "survival")
)
```

## Arguments

- X:

  matrix/data.frame of predictors (training).

- y:

  Optional outcome for supervised feature selection.

- steps:

  List of configuration options (see Details).

- task:

  "binomial", "multiclass", "gaussian", or "survival".

## Value

An object of class "GuardFit" with elements \`transform\`, \`state\`,
\`p_out\`, and \`steps\`.

## Details

The pipeline applies, in order:

- Winsorization (optional) to limit outliers.

- Imputation learned on training data only.

- Normalization (z-score or robust).

- Variance/IQR filtering. Thresholds are compared with the variance and
  IQR of the imputed training data on its original scale (before
  normalization), so `var_thresh` and `iqr_thresh` keep their meaning
  under z-score or robust scaling.

- Feature selection (optional; t-test, lasso, PCA).

All statistics are estimated on the training data and re-used for new
data.

## See also

\[predict_guard()\]

## Examples

``` r
x <- data.frame(a = c(1, 2, NA), b = c(3, 4, 5))
fit <- guard_fit(x, y = c(1, 2, 3),
                 steps = list(impute = list(method = "median")),
                 task = "gaussian")
fit$transform(x)
#>    a  b
#> 1 -1 -1
#> 2  1  0
#> 3  0  1
```
