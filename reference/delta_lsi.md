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

  Character. Exchangeability assumption for the sign-flip test. One of
  `"iid"` (default), `"by_group"`, `"within_batch"`, `"blocked_time"`.
  `"blocked_time"` activates a block sign-flip procedure that flips
  contiguous groups of repeats together, preserving serial
  autocorrelation under the null; see `block_size`. `"by_group"` and
  `"within_batch"` are stored and reported but inference still uses the
  iid sign-flip procedure (a warning is issued; contributions welcome).
  `"iid"` (default) applies the standard independent sign-flip test.

- learner:

  Optional character. Learner name to select from multi-learner fits. If
  `NULL`, the first learner found in `fit@metrics` is used.

- higher_is_better:

  Logical or `NULL`. Whether a higher value of `metric` indicates better
  performance. When `NULL` (default), auto-detected from the metric
  name: `"rmse"`, `"mse"`, `"mae"`, `"log_loss"`, `"brier"`, `"error"`,
  `"loss"`, and `"deviance"` are treated as lower-is-better; all others
  default to higher-is-better. Setting this correctly ensures that a
  positive `delta_lsi` always indicates leakage inflation (the naive
  pipeline is artificially more optimistic than the guarded one).

- block_size:

  Integer or `NULL`. Block length for the block sign-flip test, used
  only when `exchangeability = "blocked_time"`. When `NULL` (default),
  the block size is auto-estimated from the first-order autocorrelation
  of \\\\\Delta_r\\\\ and capped at `floor(R/3)` to ensure at least
  three independent blocks. A warning is issued when the estimate is
  used with `R_eff < 20` because the AR(1) estimate is noisy at small
  sample sizes. Provide an explicit integer to override auto-estimation.

- M_boot:

  Integer. Number of bootstrap samples for BCa CI (default 2000).

- M_flip:

  Integer. Maximum Monte Carlo samples for sign-flip test when R_eff \>
  15 (default 10000).

- strict:

  Logical. If `TRUE`, error on insufficient R_eff instead of a warning.

- return_details:

  Logical. If `TRUE`, include the per-repeat \\\Delta_r\\ vector and the
  original fit objects in the `info` slot.

- seed:

  Integer. Random seed for bootstrap and sign-flip test.

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
  – BCa bootstrap confidence interval for either the robust or the raw
  estimate.

- [`dlsi_p_value`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_p_value.md)`(dlsi)`
  – sign-flip randomization-test p-value.

- [`dlsi_tier`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_tier.md)`(dlsi)`
  – inference-tier label (`"A_full_inference"`, `"B_signflip_ci"`,
  `"C_signflip"`, or `"D_insufficient"`).

- [`dlsi_R_eff`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_R_eff.md)`(dlsi)`
  – effective number of paired repeats contributing to the inference.

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
fold structures* (same test-set membership per fold) in addition to the
same number of repeats. When repeat counts match but fold structures
differ, a warning is issued and the fits are treated as unpaired.

When \\R\_{\text{eff}} \geq 5\\ (equal, paired repeats), a sign-flip
randomization test (Phipson & Smyth, 2010) is performed: under \\H_0\\
(no leakage) the sign of each \\\Delta_r\\ is exchangeable. All \\2^R\\
sign combinations are enumerated exactly for \\R \leq 15\\ (no
continuity correction); Monte Carlo sampling is used for larger \\R\\
with the Phipson & Smyth (2010) correction.

BCa bootstrap confidence intervals (Efron, 1987) require
\\R\_{\text{eff}} \geq 10\\.

### Inference tiers

- `"A_full_inference"`:

  R_eff \>= 20: point + BCa CI + sign-flip p-value;
  `inference_ok = TRUE`

- `"B_signflip_ci"`:

  10 \<= R_eff \< 20: point + sign-flip p-value + BCa CI

- `"C_signflip"`:

  5 \<= R_eff \< 10: point + sign-flip p-value (no CI). The exact
  two-sided sign-flip test has minimum achievable p-value \\2 /
  2^{R\_{\text{eff}}}\\, so at \\R\_{\text{eff}} = 5\\ the smallest
  possible p is 0.0625 and the test can never reach \\p \< 0.05\\; use
  \\R\_{\text{eff}} \geq 6\\ (minimum p = 0.03125) if a 5% test is
  needed. The attainable floor is stored in `info$min_p_achievable`.

- `"D_insufficient"`:

  R_eff \< 5 or unpaired: point estimate only

## See also

[`audit_leakage`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md),
[`fit_resample`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md),
[`LeakDeltaLSI`](https://selcukorkmaz.github.io/bioLeak/reference/LeakClasses.md)
