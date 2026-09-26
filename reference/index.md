# Package index

## Splits

Leakage-resistant split policies, the splitGraph handoff, and the checks
that confirm a split actually respects the grouping axes it claims to.

- [`as_leaksplits()`](https://selcukorkmaz.github.io/bioLeak/reference/as_leaksplits.md)
  : Convert a splitGraph split_spec into bioLeak splits
- [`make_split_plan()`](https://selcukorkmaz.github.io/bioLeak/reference/make_split_plan.md)
  : Create leakage-resistant splits
- [`check_split_overlap()`](https://selcukorkmaz.github.io/bioLeak/reference/check_split_overlap.md)
  : Check split overlap invariants
- [`as_rsample()`](https://selcukorkmaz.github.io/bioLeak/reference/as_rsample.md)
  : Convert LeakSplits to an rsample resample set
- [`show(`*`<LeakDeltaLSI>`*`)`](https://selcukorkmaz.github.io/bioLeak/reference/LeakClasses.md)
  : S4 Classes for bioLeak Pipeline
- [`show(`*`<LeakSplits>`*`)`](https://selcukorkmaz.github.io/bioLeak/reference/show-LeakSplits-method.md)
  : Display summary for LeakSplits objects

## Guarded preprocessing and fitting

Imputation, normalization, filtering, and feature selection re-estimated
inside each resampling split rather than once over the whole dataset.

- [`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md)
  : Fit and evaluate with leakage guards over predefined splits
- [`tune_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/tune_resample.md)
  : Leakage-aware nested tuning with tidymodels
- [`impute_guarded()`](https://selcukorkmaz.github.io/bioLeak/reference/impute_guarded.md)
  : Leakage-safe data imputation via guarded preprocessing
- [`guard_to_recipe()`](https://selcukorkmaz.github.io/bioLeak/reference/guard_to_recipe.md)
  : Convert guard preprocessing steps to a recipes recipe
- [`guard_fit()`](https://selcukorkmaz.github.io/bioLeak/reference/guard_fit.md)
  : Fit leakage-safe preprocessing pipeline
- [`guard_ensure_levels()`](https://selcukorkmaz.github.io/bioLeak/reference/guard_ensure_levels.md)
  : Ensure consistent categorical levels for guarded preprocessing
- [`predict(`*`<GuardFit>`*`)`](https://selcukorkmaz.github.io/bioLeak/reference/predict.GuardFit.md)
  [`predict_guard()`](https://selcukorkmaz.github.io/bioLeak/reference/predict.GuardFit.md)
  : Apply a fitted GuardFit transformer to new data

## Auditing

Post-hoc evidence: permutation gaps, batch and fold association tests,
target-leakage scans, and mechanism-level risk summaries.

- [`audit_leakage()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md)
  : Audit leakage and confounding
- [`audit_leakage_by_learner()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage_by_learner.md)
  : Audit leakage per learner
- [`audit_report()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_report.md)
  : Render an HTML audit report
- [`calibration_summary()`](https://selcukorkmaz.github.io/bioLeak/reference/calibration_summary.md)
  : Calibration diagnostics for binomial predictions
- [`confounder_sensitivity()`](https://selcukorkmaz.github.io/bioLeak/reference/confounder_sensitivity.md)
  : Confounder sensitivity summaries
- [`cv_ci()`](https://selcukorkmaz.github.io/bioLeak/reference/cv_ci.md)
  : Confidence intervals for cross-validated metrics
- [`delta_lsi()`](https://selcukorkmaz.github.io/bioLeak/reference/delta_lsi.md)
  : Delta Leakage Sensitivity Index (Delta LSI)

## Audit accessors

Pull individual results out of a fitted or audited object.

- [`audit_info()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_info.md)
  : Auxiliary information from a LeakAudit
- [`audit_perm_gap()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_perm_gap.md)
  : Permutation-gap test results from a LeakAudit
- [`audit_batch_assoc()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_batch_assoc.md)
  : Batch / study association results from a LeakAudit
- [`audit_target_assoc()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_target_assoc.md)
  : Target-leakage scan results from a LeakAudit
- [`audit_duplicates()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_duplicates.md)
  : Near-duplicate pairs from a LeakAudit
- [`fit_metrics()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_metrics.md)
  : Per-fold metric data frame from a LeakFit

## Leakage severity index

Accessors for the
[`delta_lsi()`](https://selcukorkmaz.github.io/bioLeak/reference/delta_lsi.md)
result.

- [`dlsi_metric()`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_metric.md)
  : Raw delta-metric estimate from a LeakDeltaLSI
- [`dlsi_ci()`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_ci.md)
  : BCa confidence interval from a LeakDeltaLSI
- [`dlsi_p_value()`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_p_value.md)
  : Sign-flip p-value from a LeakDeltaLSI
- [`dlsi_tier()`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_tier.md)
  : Inference tier from a LeakDeltaLSI
- [`dlsi_robust()`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_robust.md)
  : Huber-robust delta_lsi point estimate from a LeakDeltaLSI
- [`dlsi_repeats()`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_repeats.md)
  : Per-repeat metric data frames from a LeakDeltaLSI
- [`dlsi_R_eff()`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_R_eff.md)
  : Effective number of paired repeats from a LeakDeltaLSI

## Plots

Visual counterparts to the audit and split diagnostics.

- [`plot(`*`<LeakAudit>`*`,`*`<missing>`*`)`](https://selcukorkmaz.github.io/bioLeak/reference/plot-LeakAudit-missing-method.md)
  : Plot method for LeakAudit
- [`plot(`*`<LeakDeltaLSI>`*`,`*`<missing>`*`)`](https://selcukorkmaz.github.io/bioLeak/reference/plot-LeakDeltaLSI-missing-method.md)
  : Plot method for LeakDeltaLSI
- [`plot(`*`<LeakFit>`*`,`*`<missing>`*`)`](https://selcukorkmaz.github.io/bioLeak/reference/plot-LeakFit-missing-method.md)
  : Plot method for LeakFit
- [`plot_calibration()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_calibration.md)
  : Plot calibration curve for binomial predictions
- [`plot_confounder_sensitivity()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_confounder_sensitivity.md)
  : Plot confounder sensitivity
- [`plot_dlsi_repeats()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_dlsi_repeats.md)
  : Plot per-repeat \\\Delta_r\\ values from a LeakDeltaLSI object
- [`plot_fold_balance()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_fold_balance.md)
  : Plot fold balance of class counts per fold
- [`plot_overlap_checks()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_overlap_checks.md)
  : Plot overlap diagnostics between train/test groups
- [`plot_perm_distribution()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_perm_distribution.md)
  : Plot permutation distribution for a LeakAudit object
- [`plot_time_acf()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_time_acf.md)
  : Plot ACF of test predictions for time-series leakage checks

## Simulation and benchmarking

Synthetic leakage scenarios for calibrating expectations and tests.

- [`benchmark_leakage_suite()`](https://selcukorkmaz.github.io/bioLeak/reference/benchmark_leakage_suite.md)
  : Simulation benchmark matrix for leakage diagnostics
- [`simulate_leakage_suite()`](https://selcukorkmaz.github.io/bioLeak/reference/simulate_leakage_suite.md)
  : Simulate leakage scenarios and audit results

## Summaries

- [`summary(`*`<LeakAudit>`*`)`](https://selcukorkmaz.github.io/bioLeak/reference/summary.LeakAudit.md)
  : Summarize a leakage audit
- [`summary(`*`<LeakDeltaLSI>`*`)`](https://selcukorkmaz.github.io/bioLeak/reference/summary.LeakDeltaLSI.md)
  : Summarize a LeakDeltaLSI object
- [`summary(`*`<LeakFit>`*`)`](https://selcukorkmaz.github.io/bioLeak/reference/summary.LeakFit.md)
  : Summarize a LeakFit object
- [`summary(`*`<LeakTune>`*`)`](https://selcukorkmaz.github.io/bioLeak/reference/summary.LeakTune.md)
  : Summarize a nested tuning result
- [`print(`*`<LeakTune>`*`)`](https://selcukorkmaz.github.io/bioLeak/reference/print.LeakTune.md)
  : Print a LeakTune object
- [`show(`*`<LeakAudit>`*`)`](https://selcukorkmaz.github.io/bioLeak/reference/show-LeakAudit-method.md)
  : Display summary for LeakAudit objects
- [`show(`*`<LeakFit>`*`)`](https://selcukorkmaz.github.io/bioLeak/reference/show-LeakFit-method.md)
  : Display summary for LeakFit objects
