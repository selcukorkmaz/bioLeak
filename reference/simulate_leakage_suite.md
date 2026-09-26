# Simulate leakage scenarios and audit results

Simulates synthetic binary classification datasets with optional leakage
mechanisms, fits a model using a leakage-aware cross-validation scheme,
and summarizes the permutation-gap audit for each Monte Carlo seed. The
suite is designed to surface validation failures such as subject overlap
across folds, batch-confounded outcomes, global normalization/summary
leakage, and time-series look-ahead. The output is a per-seed summary of
observed CV performance and its gap versus a label-permutation null; it
does not return fitted models or the full audit object. Results are
limited to the built-in data generator and leakage types implemented
here, and should be interpreted as a simulation-based sanity check
rather than a comprehensive leakage detector for real data.

## Usage

``` r
simulate_leakage_suite(
  n = 500,
  p = 20,
  prevalence = 0.5,
  mode = c("subject_grouped", "batch_blocked", "study_loocv", "time_series"),
  learner = c("glmnet", "ranger"),
  leakage = c("none", "subject_overlap", "batch_confounded", "peek_norm", "lookahead"),
  preprocess = NULL,
  rho = 0,
  K = 5,
  repeats = 1,
  horizon = 0,
  B = 200,
  seeds = 1:10,
  parallel = FALSE,
  signal_strength = 1,
  verbose = FALSE
)
```

## Arguments

- n:

  Integer scalar. Number of samples to simulate (default 500). Larger
  values stabilize the Monte Carlo summary but increase runtime.

- p:

  Integer scalar. Number of baseline predictors before any leakage
  feature is added (default 20). Increasing `p` changes the
  signal-to-noise ratio and increases fitting time.

- prevalence:

  Numeric scalar in (0, 1). Target prevalence of class 1 in the
  simulated outcome (default 0.5). Changing this alters class imbalance
  and can affect AUC and the permutation gap.

- mode:

  Character scalar. Cross-validation scheme passed to
  [`make_split_plan()`](https://selcukorkmaz.github.io/bioLeak/reference/make_split_plan.md);
  one of `"subject_grouped"`, `"batch_blocked"`, `"study_loocv"`,
  `"time_series"`. Defaults to `"subject_grouped"`. This controls how
  samples are grouped into folds (by subject, batch, study, or time) and
  therefore which leakage mechanisms are realistically challenged.

- learner:

  Character scalar. Base learner, `"glmnet"` (default) or `"ranger"`.
  Requires the corresponding package in `Suggests`. Switching learners
  changes the fitted model, runtime, and performance.

- leakage:

  Character scalar. Leakage mechanism to inject; one of `"none"`,
  `"subject_overlap"`, `"batch_confounded"`, `"peek_norm"`,
  `"lookahead"`. Leakage is added as an extra predictor:
  `"subject_overlap"` adds per-subject mean outcome,
  `"batch_confounded"` adds per-batch mean outcome, `"peek_norm"` adds
  the globally normalized (z-scored) outcome, and `"lookahead"` adds the
  next-time outcome. Changing this controls whether and how leakage is
  present.

- preprocess:

  Optional preprocessing list or recipe passed to \[fit_resample()\].
  When NULL (default), the simulator uses the fit_resample defaults; for
  `"peek_norm"` leakage, normalization is set to `"none"` to avoid
  attenuating the constant leakage feature.

- rho:

  Numeric scalar in \[-1, 1\]. AR(1)-style autocorrelation applied to
  each predictor across row order (default 0). Higher absolute values
  increase serial correlation and make time-ordered leakage more
  pronounced.

- K:

  Integer scalar. Number of folds/partitions (default 5). Used as the
  fold count for `"subject_grouped"` and `"batch_blocked"`, and as the
  number of rolling partitions for `"time_series"`. Ignored for
  `"study_loocv"` (folds equal the number of studies).

- repeats:

  Integer scalar \>= 1. Number of repeated CV runs for
  `"subject_grouped"` and `"batch_blocked"` (default 1). Increasing
  `repeats` increases the number of folds and runtime. Ignored for
  `"study_loocv"` and `"time_series"`.

- horizon:

  Numeric scalar \>= 0. Minimum time gap enforced between train and test
  for `"time_series"` splits (default 0). Larger values make the split
  more conservative and can reduce leakage from temporal proximity.

- B:

  Integer scalar \>= 1. Number of permutations used by
  [`audit_leakage()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md)
  to compute the permutation gap and p-value (default 200). Larger
  values yield more stable p-values but increase runtime.

- seeds:

  Integer vector. Monte Carlo seeds (default `1:10`). One row of output
  is produced per seed; changing `seeds` changes the simulated datasets
  and splits.

- parallel:

  Logical scalar. If `TRUE`, evaluates seeds in parallel using
  `future.apply` (if installed). Results are identical to sequential
  execution; only runtime changes.

- signal_strength:

  Numeric scalar. Scales the linear predictor before sampling outcomes
  (default 1). Larger values increase class separation and tend to
  increase AUC; smaller values make the task harder.

- verbose:

  Logical scalar. If `TRUE`, prints progress messages for each seed.
  Does not affect results.

## Value

A `LeakSimResults` data frame with one row per seed and columns:

- `seed`: seed used for data generation, splitting, and auditing.

- `metric_obs`: observed CV performance (AUC for this simulation).

- `gap`: permutation-gap statistic (observed minus permutation mean).

- `p_value`: permutation p-value for the gap.

- `leakage`: leakage scenario used.

- `mode`: CV mode used.

Only the permutation-gap summary is returned; fitted models,
predictions, and other audit components are not included.

## Details

The generator draws `p` standard normal predictors, builds a linear
predictor from the first `min(5, p)` features, scales it by
`signal_strength`, and samples a binary outcome to achieve the requested
`prevalence`. Outcomes are returned as a two-level factor, so the
audited metric is AUC. Simulated metadata include subject, batch, study,
and time fields used by `mode` to create leakage-aware splits. Leakage
mechanisms are injected by adding a single extra predictor as described
in `leakage`. Parallel execution uses `future.apply` when installed and
does not change results.

## Note

This function is a general-purpose utility and its data-generation logic
intentionally differs from the custom simulation used in the bioLeak
manuscript (the manuscript's `replication.R` script, distributed with
the manuscript supplementary materials). Specific differences:

- **peek_norm leakage**: this function uses a z-scored binary outcome as
  the leak feature; the manuscript uses a noisy continuous version
  (`as.numeric(y) + rnorm(n, 0, 0.3)`).

- **lookahead leakage**: this function shifts the binary outcome
  (`c(y[-1], y[n])`); the manuscript shifts a continuous biomarker
  (`linpred + noise`).

- **signal generation**: this function applies AR correlation to
  predictors via `rho`; the manuscript adds AR(1) noise directly to the
  linear predictor.

- **audit settings**: the manuscript uses `perm_refit = FALSE` and
  `perm_stratify = TRUE`; this function uses `perm_refit = "auto"` and
  the `perm_stratify` default (`FALSE`).

Users wishing to reproduce manuscript figures should run the
manuscript-specific `replication.R` script directly rather than calling
this function.

## Examples

``` r
# \donttest{
  if (requireNamespace("glmnet", quietly = TRUE)) {
    set.seed(1)
    res <- simulate_leakage_suite(
      n = 120, p = 6, prevalence = 0.4,
      mode = "subject_grouped",
      learner = "glmnet",
      leakage = "subject_overlap",
      K = 3, repeats = 1,
      B = 50, seeds = 1,
      parallel = FALSE
    )
    # One row per seed with observed AUC, permutation gap, and p-value
    res
  }
#> subject_grouped: repeat 1/1 done.
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   seed metric_obs    gap   p_value         leakage            mode
#> 1    1   0.900857 0.4138 0.0196078 subject_overlap subject_grouped
# }
```
