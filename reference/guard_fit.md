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

- Imputation learned on training data only. The requested method is
  fitted even when the training rows have no missing values, so missing
  values in new rows are still imputed (for `"missForest"` with the
  training medians, as no per-column model exists then).

- Normalization (`"zscore"`, `"robust"`, frozen reference `"quantile"`,
  `"log_cpm"`, `"clr"`, or `"none"`).

- Optional frozen train-only batch adjustment (`steps$batch_adjust` with
  `method` in `"combat_frozen"`, `"location_scale"`, or
  `"mean_center"`). Batch labels come from `batch` (a vector) or
  `batch_col` (a column of `X`, removed from the predictors unless
  `keep_batch_col = TRUE`; a kept batch column is passed through
  unadjusted). Two further options:

  - `protect`:

    `FALSE` (default) or `TRUE`. When batch and outcome are associated
    in the training data (p \< 0.01), a warning is issued either way,
    because both choices are biased: without protection, part of the
    outcome signal is removed as batch effect, so performance is
    underestimated; with `protect = TRUE` (and a non-survival `y`),
    batch effects are estimated net of the outcome, per feature, from
    `z ~ batch + y`, the role of ComBat's `mod`. That keeps the signal,
    but the estimated batch effects carry fitted outcome noise, which
    can bias cross-validated performance upward (Nygaard et al., 2016);
    in simulation with 80/20 confounding and no true signal,
    cross-validated AUC rose by about 0.03-0.10, more with more
    features. The default is the conservative choice. New rows need only
    their batch label. With `protect = TRUE` only batch locations are
    adjusted (no per-batch scaling), because rescaling would also
    rescale the protected outcome effect. Outcome effects are estimated
    from variation within batches only; outcome contrasts that do not
    vary within batches are not protected (with a warning). If every
    batch holds a single outcome class, only differences between batches
    of the same class are removed; if no class spans two batches, no
    adjustment is applied. Both cases warn.

  - `unseen`:

    What happens at
    [`transform()`](https://rdrr.io/r/base/transform.html) for rows
    whose batch was not seen in training, has a missing label, or when
    no batch labels are supplied: such rows are left unadjusted, because
    batch effects are known only for training batches. `"warn"`
    (default) warns, `"error"` stops, `"grand_mean"` accepts it
    silently. Under batch-blocked or study-held-out CV every test batch
    is unseen;
    [`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md)
    then reports this once, not per fold.

  The `"combat_frozen"` method is a simplified empirical-Bayes
  location/scale adjustment in the spirit of ComBat, not a
  reimplementation of ComBat's priors.

- Variance/IQR filtering. Thresholds are compared with the variance and
  IQR of the imputed training data on its original scale (before
  normalization), so `var_thresh` and `iqr_thresh` keep their meaning
  under z-score or robust scaling.

- Feature selection (optional; t-test, lasso, PCA).

All statistics are estimated on the training data and re-used for new
data.

## References

Nygaard, V., Rødland, E. A. and Hovig, E. (2016). Methods that remove
batch effects while retaining group differences may lead to exaggerated
confidence in downstream analyses. *Biostatistics*, 17(1), 29-39.

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
