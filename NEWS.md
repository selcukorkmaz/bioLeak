# bioLeak 0.3.1 (development)

## Breaking changes

* `delta_lsi()`: inference tier strings renamed to accurately reflect what each
  tier provides.  `"C_point_only"` → `"C_signflip"` (the sign-flip p-value is
  available at this tier, not just point estimates); `"B_ci_only"` →
  `"B_signflip_ci"` (both the sign-flip p-value and BCa CI are available).
  Code that compares `result@tier` against the old string literals must be
  updated.

## New features

* `delta_lsi()` gains a `block_size` argument and makes `exchangeability`
  actionable for `"blocked_time"` inputs.  When `exchangeability = "blocked_time"`,
  the sign-flip test now uses a block procedure that flips contiguous blocks of
  repeats together, preserving serial autocorrelation under the null.
  `block_size` is auto-estimated from the AR(1) of the repeat-level deltas when
  `NULL` (default) and capped at `floor(R/3)` to guarantee at least three
  independent blocks.  The `@info` slot gains `block_size_used` and `n_blocks`
  fields.  If the block structure yields fewer than five independent blocks,
  `@p_value` is set to `NA` and a warning is issued.
* `delta_lsi()` now emits an explicit warning when `exchangeability` is
  `"by_group"` or `"within_batch"`, informing users that those modes are stored
  but inference still uses the iid sign-flip procedure.  Previously these values
  were accepted silently without affecting computation.

## Bug fixes and improvements

* `show()` and `summary()` for `LeakDeltaLSI` now label the sign-flip p-value
  as testing `mean(Δr)` (delta_metric), not delta_lsi, making the
  estimator–inference pairing explicit.
* `summary()` prints a diagnostic note when the sign-flip p-value and BCa CI
  lead to qualitatively different conclusions (one significant, one spanning
  zero), which can occur when outlier repeats pull the arithmetic mean away from
  the Huber estimate.
* `summary()` prints the block size and number of blocks used when
  `exchangeability = "blocked_time"`.

# bioLeak 0.3.0

## New features

* Added **N-axis combined splitting** via `constraints` in `make_split_plan()`,
  generalizing beyond two-axis combined CV while preserving train/test exclusion
  across all declared axes.
* Added `compact = TRUE` split storage (fold assignments) for large datasets to
  reduce split object memory footprint.
* Added `check_split_overlap()` for explicit overlap-invariant validation across
  fold/group axes.
* Added `cv_ci()` (with Nadeau-Bengio correction) and integrated CI columns into
  `fit_resample()` and `tune_resample()` metric summaries (`*_ci_lo`, `*_ci_hi`).
* Added `guard_to_recipe()` to map guarded preprocessing configurations to
  `recipes` pipelines with explicit fallback/warning behavior.
* Added `benchmark_leakage_suite()` for reproducible modality-by-mechanism
  benchmark grids and detection-rate summaries.
* Expanded `audit_leakage()` diagnostics with mechanism taxonomy fields
  (`mechanism_class`, `taxonomy`, `mechanism_summary`) and richer risk
  attribution outputs.
* Added FDR-aware target scan outputs (`p_value_adj`, `flag_fdr`) with selectable
  multiple-testing correction (`target_p_adjust`, `target_alpha`).
* Added `feature_space` (`raw`/`rank`) and `duplicate_scope`
  (`train_test`/`all`) controls for duplicate diagnostics.
* Strengthened permutation auditing with explicit `perm_mode` handling for
  rsample-derived splits and safer `perm_refit = "auto"` behavior.
* Extended tidymodels interoperability: rsample conversion and metadata inference
  are more robust (`split_cols = "auto"`, mode/perm-mode propagation, stricter
  compatibility checks).
* Improved nested tuning safety in `tune_resample()`: final refit now aggregates
  hyperparameters across outer folds (median/majority) instead of selecting a
  single best outer fold.
* Added binomial threshold tuning support in `tune_resample()` using inner-fold
  predictions (`tune_threshold`, `threshold_grid`, `threshold_metric`).
* Added structured fold-status tracking (`fold_status`) and elapsed timing in
  both fitting and tuning paths for better failure-mode observability.
* Added strict-mode and validation-policy infrastructure (`bioLeak.strict`,
  `bioLeak.validation_mode`) with structured condition classes for safer recipe
  and workflow guardrails.
* Added provenance capture (`.bio_capture_provenance`) and attached provenance
  metadata to `LeakFit`, `LeakAudit`, and `LeakTune`.
* Improved `summary.LeakAudit()` output with explicit Mechanism Risk Assessment
  reporting.
* Hardened recipe preprocessing in `fit_resample()` to avoid fold-time failures
  when recipes reference split metadata columns (for example `subject`).
* Updated simulation defaults and audit settings for more practical runtime
  (`simulate_leakage_suite()` default `B`, auto refit cap handling).
* Updated manuscript/simulation assets under `paper/` with refreshed large-scale
  simulation outputs and case-study artifacts.

---

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
