# bioLeak 0.2.0

## New features

* **Leak-safe hyperparameter tuning** via `tune_resample()`: nested
  cross-validation using tidymodels `tune`/`dials` with leakage-aware outer
  splits.
* **Tidymodels interoperability**: `fit_resample()` now accepts `rsample`
  rset/rsplit objects as `splits`, `recipes::recipe` for preprocessing,
  `workflows::workflow` as `learner`, and `yardstick::metric_set` for metrics.
  `as_rsample()` converts `LeakSplits` to an `rsample` rset.
* **Parsnip model specs** accepted directly as the `learner` argument in
  `fit_resample()`.
* **Diagnostics polish**: new `calibration_summary()` and `plot_calibration()`
  for probability calibration checks; `confounder_sensitivity()` and
  `plot_confounder_sensitivity()` for sensitivity analysis.
* **Simulation utility** `simulate_leakage_suite()` for generating controlled
  leakage scenarios and benchmarking audit sensitivity.
* **HTML audit report** via `audit_report()`: renders a self-contained HTML
  summary of all audit results for sharing and review.
* **Multi-learner auditing** with `audit_leakage_by_learner()` to audit each
  learner in a multi-model fit separately.
* **Multivariate target leakage scan** enabled by default in `audit_leakage()`
  for supported tasks, complementing the existing univariate scan.
* **Refit-based permutations** (`perm_refit = TRUE` or `"auto"`) in
  `audit_leakage()` for a more powerful permutation gap test when refit data
  are available.
* **Class weights** support in `fit_resample()` for imbalanced classification
  tasks.
* New plotting functions: `plot_fold_balance()`, `plot_overlap_checks()`,
  `plot_perm_distribution()`, `plot_time_acf()`.

## Improvements

* S4 classes (`LeakSplits`, `LeakFit`, `LeakAudit`) now include `setValidity`
  checks for slot consistency.
* `summary()` methods for `LeakFit`, `LeakAudit`, and `LeakTune` improved with
  clearer console output and edge-case handling.
* `impute_guarded()` gains enhanced diagnostics and RNG safety.
* `.guard_fit()` and `.guard_ensure_levels()` made more robust with better error
  messages.
* Permutation label factory (`permute_labels`) gains verbose mode, digest-based
  caching, and improved stratification safety.
* `audit_leakage()` handles NA metrics gracefully and enriches trail metadata.
* `make_split_plan()` improved stratification logic and reproducible seeding.
* `audit_report()` now renders from a temporary copy of the Rmd template to
  avoid write failures on read-only file systems (e.g. during `R CMD check`).
* Comprehensive vignette (`bioLeak-intro`) rewritten with guided workflow and
  leaky-vs-correct comparisons.

## Bug fixes

* Fixed `fit_resample()` result aggregation when folds fail during
  preprocessing.
* Fixed `missForest` preprocessing dropping rows.
* Fixed single-level factors causing errors in guarded preprocessing.
* Fixed filter keep-column alignment by name.
* Fixed `glmnet` folds receiving non-numeric design matrices.
* Fixed constant imputation for categorical data.
* Fixed RANN self-neighbour filter in duplicate detection.
* Fixed various edge cases in outcome extraction and hashing utilities.
* Resolved multiple CRAN check issues (Rd formatting, example runtime,
  read-only file-system writes).

---

# bioLeak 0.1.0

* Initial release.
* **Core pipeline**: `make_split_plan()` for leakage-aware splitting
  (subject-grouped, batch-blocked, study leave-out, time-ordered);
  `fit_resample()` for cross-validated fitting with built-in guarded
  preprocessing (train-only imputation, normalisation, filtering, feature
  selection).
* **Leakage auditing**: `audit_leakage()` with label-permutation gap test,
  batch/study association tests, univariate target leakage scan, and
  near-duplicate detection.
* **Guarded preprocessing helpers**: `impute_guarded()`, `predict_guard()`,
  `.guard_fit()`, `.guard_ensure_levels()`.
* S4 class system: `LeakSplits`, `LeakFit`, `LeakAudit`.
* Support for binomial, multiclass, regression, and survival tasks.
* Built-in learners: `glm`, `glmnet`, `ranger`, `xgboost` (via
  `custom_learners`).
* `SummarizedExperiment` input support.
* Vignette and comprehensive documentation.
