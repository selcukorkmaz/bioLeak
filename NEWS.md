# bioLeak 0.3.7

## Documentation

* The vignette previously mixed defensive `requireNamespace()` checks
  with bare `library(<suggests_pkg>)` calls; the bare calls would
  hard-error if the suggested package was not installed, defeating
  the defensive checks elsewhere in the same vignette. The
  `tidymodels-interop` chunk now carries
  `eval = requireNamespace("recipes", quietly = TRUE) && requireNamespace("yardstick", quietly = TRUE)`
  in its chunk header, so the chunk is skipped (rather than erroring)
  during vignette build when those Suggests packages are absent. The
  `parallel-setup` chunk gains a brief comment documenting that
  `future` is a Suggests dependency. A regression test
  (`test-vignette-suggests.R`) walks the vignette and asserts that
  every chunk-level `library(<suggests_pkg>)` call is inside an
  appropriately gated chunk.

## API improvements (no behavior change)

* `predict_guard()` is now also accessible through the standard
  [stats::predict()] generic via a registered S3 method
  `predict.GuardFit()`. Calling `predict(fit, newdata)` on a `GuardFit`
  object dispatches to `predict.GuardFit()` and yields output that is
  bit-identical to the legacy `predict_guard(fit, newdata)`.
  `predict_guard()` is preserved as a thin backward-compatible alias,
  so existing code continues to work without modification.
  `methods(class = "GuardFit")` now returns `print`, `summary`, and
  `predict`, restoring the standard R idiom for transformer objects.

* Added `show()` / `print()` methods to the public result classes that
  previously only had `summary()`:
    * `LeakFit`: new `show()` (S4) — brief auto-print giving task,
      outcome, learners, fold count, and fold-status one-liner.
    * `LeakAudit`: new `show()` (S4) — brief auto-print giving task,
      outcome, permutation-gap statistics, and component row counts
      (batch association, target leakage, duplicates).
    * `LeakTune`: new `print()` (S3) — brief auto-print giving outer-fold
      success rate, tuning-grid size, selection rule, and refit status.
  Each method ends with a one-line hint pointing to `summary(<obj>)`
  for the full diagnostic report. `methods(class = ...)` now returns
  `show`/`print` alongside `summary` for all three classes.

## Renames (no behavior change)

* `.guard_fit()` is renamed to `guard_fit()` and `.guard_ensure_levels()`
  is renamed to `guard_ensure_levels()`. Leading-dot prefixes on
  exported functions are unconventional and were causing the renamed
  helpers to appear awkwardly in `help(package = "bioLeak")`. Behavior,
  arguments, and return values are unchanged; only the names move from
  the dot-prefixed form to ordinary names. Internal callers
  (`fit_resample()`, `impute_guarded()`, `predict_guard()`'s
  documentation, and the package vignette) are updated to use the new
  names.

## New features

* Added public accessor functions for the S4 result classes so that
  downstream code (replication scripts, vignettes, and end-user
  analyses) can read components of `LeakFit`, `LeakAudit`, and
  `LeakDeltaLSI` objects without reaching into S4 internals via `@`.
  The new accessors are purely additive; slot definitions are unchanged
  and existing code that uses `@` continues to work.
    * `LeakFit`: `fit_metrics()`.
    * `LeakAudit`: `audit_perm_gap()`, `audit_batch_assoc()`,
      `audit_target_assoc()`, `audit_duplicates()`, `audit_info()`.
    * `LeakDeltaLSI`: `dlsi_metric()`, `dlsi_robust()`, `dlsi_ci()`,
      `dlsi_p_value()`, `dlsi_tier()`, `dlsi_R_eff()`, `dlsi_repeats()`.
  Each accessor performs an `is(x, "<Class>")` validation and emits an
  informative error when called on the wrong object.

# bioLeak 0.3.5

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

* `fit_resample()`: compact + combined mode now correctly excludes
  constraint-axis violations from training sets.  Previously the compact
  fallback used `setdiff(all, test)`, ignoring multi-axis constraints declared
  via `make_split_plan(constraints = ...)`.  The same fix is applied in the
  `as_rsample()` conversion path for consistency.
* Guarded preprocessing: lasso and t-test feature selection now uses name-based
  column selection in the transform step, preventing index misalignment when
  constant columns are removed during fitting.
* `delta_lsi()`: `R_eff` and the inference tier are now recomputed after
  repeat-level intersection, so that dropped all-NA repeats correctly reduce the
  effective sample size and select the appropriate tier.
* `fit_resample()`: fold error messages are now correctly captured when running
  in parallel via `future.apply`.  Previously `<<-` mutations inside worker
  processes were silently lost; errors are now attached as result attributes and
  extracted after the parallel map.
* `tune_resample()`: fold-ID columns (`id`, `id2`, `.notes`) no longer leak
  into hyperparameter aggregation in the internal `select_config()` helper.
* `summary.LeakFit()` now returns `object@metric_summary` invisibly, matching
  the documented return value (previously returned the object itself).
* Fixed vignette (`bioLeak-intro`) referencing a shadowed data frame for sample
  count; now reads from `fit_safe@splits@info$coldata`.
* Fixed `audit_leakage()` roxygen documenting a `duplicates` column named
  `in_train_test`; the actual column name is `cross_fold`.
* `make_split_plan()`: time-series mode now warns and skips folds with fewer
  than 3 test samples instead of producing degenerate folds.
* `fit_resample()`: added bounds checking for `repeat_id` in compact fold
  resolution to produce a clear error instead of a cryptic index failure.
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
