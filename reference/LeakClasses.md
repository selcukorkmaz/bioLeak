# S4 Classes for bioLeak Pipeline

These classes capture splits, model fits, and audit diagnostics produced
by
[`make_split_plan()`](https://selcukorkmaz.github.io/bioLeak/reference/make_split_plan.md),
[`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md),
and
[`audit_leakage()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md).

## Usage

``` r
# S4 method for class 'LeakDeltaLSI'
show(object)
```

## Arguments

- object:

  A `LeakDeltaLSI` object.

## Value

An S4 object of the respective class.

## Slots

- `mode`:

  Splitting mode. One of "subject_grouped", "batch_blocked",
  "study_loocv", "time_series", or "combined".

- `indices`:

  List of resampling descriptors (train/test indices when available)

- `info`:

  (LeakSplits) Metadata associated with the split plan (mode, coldata,
  hash, etc.)

- `splits`:

  A \[\`LeakSplits\`\] object used for resampling

- `metrics`:

  Model performance metrics per resample

- `metric_summary`:

  Summary of metrics across resamples

- `audit`:

  Audit information per resample

- `predictions`:

  List of prediction objects

- `preprocess`:

  Preprocessing steps used during fitting

- `learners`:

  Learner definitions used in the pipeline

- `outcome`:

  Outcome variable name

- `task`:

  Modeling task name

- `feature_names`:

  Feature names included in the model

- `info`:

  (LeakFit) Metadata about the model fit (sample IDs, timings,
  provenance, etc.)

- `fit`:

  A \[\`LeakFit\`\] object used to generate the audit

- `permutation_gap`:

  Data frame summarising permutation gaps

- `perm_values`:

  Numeric vector of permutation-based scores

- `batch_assoc`:

  Data frame of batch associations

- `target_assoc`:

  Data frame of feature-wise outcome associations

- `duplicates`:

  Data frame detailing duplicate records

- `trail`:

  List capturing audit trail information

- `info`:

  (LeakAudit) Metadata about the audit (mechanism summary, settings,
  provenance, etc.)

- `metric`:

  Performance metric compared between pipelines

- `exchangeability`:

  Exchangeability assumption used for the sign-flip test

- `tier`:

  Inference tier label based on effective number of repeats

- `strict`:

  Whether strict mode was requested

- `R_eff`:

  Effective number of paired repeats available for inference

- `delta_lsi`:

  Huber-robust point estimate of repeat-level metric difference

- `delta_lsi_ci`:

  BCa 95% CI for delta_lsi (NA when R_eff \< 10)

- `delta_metric`:

  Arithmetic mean of repeat-level metric differences

- `delta_metric_ci`:

  BCa 95% CI for delta_metric (NA when R_eff \< 10)

- `p_value`:

  Sign-flip randomization test p-value (NA when R_eff \< 5 or unpaired)

- `inference_ok`:

  TRUE when tier A (R_eff \>= 20, paired, finite p and CI)

- `folds_naive`:

  Per-fold data frame for the naive pipeline

- `folds_guarded`:

  Per-fold data frame for the guarded pipeline

- `repeats_naive`:

  Per-repeat aggregate data frame for the naive pipeline

- `repeats_guarded`:

  Per-repeat aggregate data frame for the guarded pipeline

- `info`:

  (LeakDeltaLSI) Metadata including R_naive, R_guarded, paired status,
  and block details

## Accessors (LeakFit)

The S4 accessor method \[fit_metrics()\] returns the per-fold metric
data frame. It is listed under `methods(class = "LeakFit")` alongside
the `show` and `summary` methods.

## Plot method (LeakFit)

`plot(<LeakFit>)` dispatches to \[plot_fold_balance()\] by default. Use
the `which` argument to switch to one of `"overlap"`, `"calibration"`
(binary outcomes only), `"time_acf"` (time-ordered splits), or
`"confounder_sensitivity"`.

## Accessors (LeakAudit)

The S4 accessor methods \[audit_perm_gap()\], \[audit_batch_assoc()\],
\[audit_target_assoc()\], \[audit_duplicates()\], and \[audit_info()\]
return the corresponding slots. They are listed under
`methods(class = "LeakAudit")` alongside the `show` and `summary`
methods.

## Plot method (LeakAudit)

`plot(<LeakAudit>)` dispatches to \[plot_perm_distribution()\],
rendering the permutation-null distribution with the observed metric and
permuted mean marked.

## Accessors (LeakDeltaLSI)

The S4 accessor methods \[dlsi_metric()\], \[dlsi_robust()\],
\[dlsi_ci()\], \[dlsi_p_value()\], \[dlsi_tier()\], \[dlsi_R_eff()\],
and \[dlsi_repeats()\] return the corresponding components. They are
listed under `methods(class = "LeakDeltaLSI")` alongside the `show` and
`summary` methods.

## Plot method (LeakDeltaLSI)

`plot(<LeakDeltaLSI>)` dispatches to \[plot_dlsi_repeats()\], rendering
the per-repeat \\\Delta_r\\ scatter with the Huber-robust point
estimate, the arithmetic mean, and the BCa bootstrap confidence interval
band.

## See also

\[make_split_plan()\], \[fit_resample()\], \[audit_leakage()\]

\[fit_resample()\], \[fit_metrics()\], \[plot_fold_balance()\]

\[audit_leakage()\], \[audit_report()\], \[audit_perm_gap()\],
\[plot_perm_distribution()\]

\[delta_lsi()\], \[dlsi_metric()\], \[dlsi_ci()\],
\[plot_dlsi_repeats()\]
