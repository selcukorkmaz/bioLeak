# Convert guard preprocessing steps to a recipes recipe

Maps bioLeak guard preprocessing steps (impute, normalize, filter, fs)
to their closest recipes equivalents. Requires the recipes package.
Steps that have no direct recipe equivalent are skipped with a warning.

## Usage

``` r
guard_to_recipe(steps, formula, training_data)
```

## Arguments

- steps:

  A named list of guard preprocessing steps, e.g.,
  `list(impute = list(method = "median"), normalize = list(method = "zscore"))`.

- formula:

  A model formula (e.g., `outcome ~ .`).

- training_data:

  A data.frame used to initialize the recipe.

## Value

A
[`recipes::recipe`](https://recipes.tidymodels.org/reference/recipe.html)
object with the mapped steps added.

## Details

Mapping:

- `impute$method = "median"`:
  `step_impute_median(all_numeric_predictors())`

- `impute$method = "knn"`:
  `step_impute_knn(all_predictors(), neighbors = k)`

- `impute$method = "missForest"` or `"mice"`: Warning +
  `step_impute_median()` fallback

- `normalize$method = "zscore"`:
  `step_normalize(all_numeric_predictors())`

- `normalize$method = "robust"`: Warning + `step_normalize()` fallback

- `normalize$method = "none"`: No step added

- `filter$var_thresh > 0`: `step_nzv(all_numeric_predictors())`

- `fs$method = "pca"`:
  `step_pca(all_numeric_predictors(), num_comp = ncomp)`

- `fs$method = "ttest"` or `"lasso"`: Warning, skipped (no recipe
  equivalent)
