# Delta Leakage Sensitivity Index (Delta LSI)

Compares a naive (potentially leaky) cross-validation pipeline against a
guarded (leakage-protected) pipeline and quantifies leakage-induced
performance inflation using the Leakage Sensitivity Index (LSI).

## Usage

``` r
delta_lsi(
  fit_leaky,
  fit_guarded,
  metric = "auc",
  exchangeability = c("iid", "by_group", "within_batch", "blocked_time"),
  learner = NULL,
  higher_is_better = NULL,
  block_size = NULL,
  M_boot = 2000L,
  M_flip = 10000L,
  strict = FALSE,
  return_details = FALSE,
  seed = 42L,
  ...
)
```

## Arguments

- fit_leaky:

  A
  [`LeakFit`](https://selcukorkmaz.github.io/bioLeak/reference/LeakClasses.md)
  object from the leaky (unprotected) evaluation pipeline.

- fit_guarded:

  A
  [`LeakFit`](https://selcukorkmaz.github.io/bioLeak/reference/LeakClasses.md)
  object from the guarded (leakage-protected) evaluation pipeline.

- metric:

  Character. Performance metric to compare. Must appear in `fit@metrics`
  of both fits (e.g., `"auc"`, `"rmse"`).

- exchangeability:

  Character. One of `"iid"` (default), `"by_group"`, `"within_batch"`,
  `"blocked_time"`. Stored and reported; all settings use the same
  corrected t-test, because the correction already accounts for the
  dependence these settings used to address. For `"by_group"` /
  `"within_batch"` the number of group / batch levels is reported in
  `info$n_groups_used` / `info$n_batches_used`. Earlier versions ran
  sign-flip tests here that were anti-conservative with real model fits.

- learner:

  Optional character. Learner name to select from multi-learner fits. If
  `NULL`, the first learner found in `fit@metrics` is used.

- higher_is_better:

  Logical or `NULL`. Whether a higher value of `metric` indicates better
  performance. When `NULL` (default), auto-detected from the metric
  name: names containing `"rmse"`, `"mse"`, `"mae"`, `"mape"`, `"mase"`,
  `"log_loss"`, `"brier"`, `"error"`, `"loss"`, `"deviance"` or `"cost"`
  are lower-is-better; names such as `"auc"`, `"accuracy"`, `"f1"`,
  `"kappa"`, `"rsq"` or `"cindex"` are higher-is-better; any other name
  warns and is treated as higher-is-better. Setting this correctly
  ensures that a positive `delta_lsi` always indicates leakage inflation
  (the naive pipeline is artificially more optimistic than the guarded
  one).

- block_size:

  Ignored; retained for compatibility (a warning is issued when
  supplied).

- M_boot:

  Ignored; retained for compatibility (no bootstrap is used).

- M_flip:

  Ignored; retained for compatibility (no sign-flip test is used).

- strict:

  Logical. If `TRUE`, error on insufficient R_eff instead of a warning.

- return_details:

  Logical. If `TRUE`, include the per-repeat \\\Delta_r\\ vector, the
  paired fold deltas and the original fit objects in the `info` slot.

- seed:

  Ignored; retained for compatibility (the method is deterministic).

- ...:

  Unused. Reserved for deprecated aliases such as `fit_naive`.

## Value

A
[`LeakDeltaLSI`](https://selcukorkmaz.github.io/bioLeak/reference/LeakClasses.md)
object. Use [`summary()`](https://rdrr.io/r/base/summary.html) to print
a formatted report. For programmatic access to inflation estimates,
confidence intervals, p-values, tier labels, and the per-repeat
fold-level deltas, use the S4 accessor methods registered for
`LeakDeltaLSI`:

- [`dlsi_metric`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_metric.md)`(dlsi)`
  – raw mean of per-repeat metric differences (the `delta_metric` slot).

- [`dlsi_robust`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_robust.md)`(dlsi)`
  – Huber-robust point estimate (the `delta_lsi` slot).

- [`dlsi_ci`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_ci.md)`(dlsi, which = "robust" | "metric")`
  – corrected-t 95% confidence interval for either the robust or the raw
  estimate.

- [`dlsi_p_value`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_p_value.md)`(dlsi)`
  – corrected repeated-CV t-test p-value.

- [`dlsi_tier`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_tier.md)`(dlsi)`
  – inference-tier label (`"A_full_inference"`, `"B_signflip_ci"`,
  `"C_signflip"`, or `"D_insufficient"`).

- [`dlsi_R_eff`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_R_eff.md)`(dlsi)`
  – effective number of independent units (folds per repeat).

- [`dlsi_repeats`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_repeats.md)`(dlsi, which = "naive" | "guarded")`
  – per-repeat metric data frame for the requested pipeline.

The full list of accessor methods for the class is available through
`methods(class = "LeakDeltaLSI")`.

## Details

### Method

For each fit, per-fold metric values are extracted from `fit@metrics`
(or recomputed from `fit@predictions` if necessary). Fold test-set sizes
are used as weights to aggregate fold metrics into per-repeat estimates
\\\mu_r\\. The repeat-level delta \\\Delta_r = s \cdot
(\mu_r^{\text{naive}} - \mu_r^{\text{guarded}})\\ captures
leakage-induced performance inflation for each CV repeat, where \\s =
+1\\ for higher-is-better metrics (e.g., AUC) and \\s = -1\\ for
lower-is-better metrics (e.g., RMSE), so that \\\Delta_r \> 0\\ always
indicates the naive pipeline is more optimistic than the guarded one.

The **delta_lsi** point estimate is the Huber M-estimator (k = 1.345)
applied to \\\\\Delta_r\\\\, which is robust to occasional outlier
repeats. **delta_metric** is the arithmetic mean of \\\\\Delta_r\\\\ for
easy interpretation in the original metric's units.

Pairing requires that `fit_leaky` and `fit_guarded` share *identical
fold structures* (same test-set membership per fold). When fold
structures differ, a warning is issued and the fits are treated as
unpaired. Folds whose metric is missing in either fit are dropped from
both before pairing, so the paired delta always compares the same test
sets.

**Inference** uses a corrected repeated cross-validation t-test on the
paired fold deltas \\d_j\\, with folds that re-use an identical test set
(e.g. repeated leave-one-batch-out) averaged into one, so \\J\\ counts
distinct test sets. Repeated CV re-partitions the same data, so fold
deltas are correlated through overlapping training sets and across
repeats; treating them (or the repeats) as independent gives badly
anti-conservative tests. The variance of the mean is corrected as
\\(1/J + n\_{\text{test}}/n\_{\text{train}})\\ s_d^2\\ (Nadeau & Bengio,
2003), and the t reference uses \\k - 1\\ degrees of freedom, where
\\k\\ is the number of folds per repeat (calibrated degrees of freedom,
Bouckaert, 2003): repeats sharpen the estimate but add no independent
information. The test statistic is `delta_metric / se`; the 95% CIs for
`delta_metric` and `delta_lsi` are the estimate \\\pm t\_{0.975, k-1}\\
se\\. In simulations with real model fits under \\H_0\\
(subject-grouped, batch-blocked and leave-one-study-out designs, 1-20
repeats) the type-I error stayed at or below the nominal 5%; with 3-fold
designs (2 df) the test is conservative.

`R_eff` is the number of folds per repeat \\k\\ (for `study_loocv`, the
number of studies), the number of effectively independent units.

### Inference tiers

Tier labels are kept from earlier versions for compatibility; every tier
except D reports a p-value and both CIs.

- `"A_full_inference"`:

  R_eff \>= 20

- `"B_signflip_ci"`:

  10 \<= R_eff \< 20

- `"C_signflip"`:

  3 \<= R_eff \< 10

- `"D_insufficient"`:

  R_eff \< 3 or unpaired: point estimate only

`inference_ok` is `TRUE` whenever a finite p-value and CI were computed.

## References

Nadeau, C. and Bengio, Y. (2003). Inference for the generalization
error. *Machine Learning*, 52, 239-281.

Bouckaert, R. R. (2003). Choosing between two learning algorithms based
on calibrated tests. *Proceedings of ICML 2003*, 51-58.

## See also

[`audit_leakage`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md),
[`fit_resample`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md),
[`LeakDeltaLSI`](https://selcukorkmaz.github.io/bioLeak/reference/LeakClasses.md)
