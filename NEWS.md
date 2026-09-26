# bioLeak 0.3.8

## Bug fixes (results can change)

* **AUC orientation.** AUC was computed with `pROC::roc(truth, pred)` using
  pROC's default `direction = "auto"`, which silently flips the curve when
  predictions are anti-correlated with the outcome. A perfectly inverted
  model scored AUC = 1, AUC could never fall below about 0.5, and
  permutation nulls centred above 0.5 (0.51-0.57 observed). AUC is now
  always oriented, with explicit `levels = c(negative, positive)` and
  `direction = "<"` (the positive class is the second outcome level, as set
  by `positive_class`), through a single internal helper used by
  `fit_resample()` and `audit_leakage()` (observed metric and permutation
  null). The univariate and multivariate target scans already used an
  oriented rank AUC and now share the same helper. Anti-correlated
  predictions now give AUC < 0.5, and label-permutation nulls centre at 0.5.
  The target-scan flagging `score = |AUC - 0.5| * 2` is unchanged, so strong
  inverse proxies are still flagged.

* **`simulate_leakage_suite()` prevalence.** The outcome generator used
  `pnorm(linpred - qnorm(prevalence))`, which inverted the requested
  prevalence (0.2 gave about 0.72; 0.8 gave about 0.27). It now uses
  `pnorm(linpred + b0)` with `b0 = qnorm(prevalence) * sqrt(1 + signal_strength^2)`,
  which also corrects the attenuation from the latent-variable scale, so the
  marginal prevalence equals the requested value at every `signal_strength`.
  The `imaging_tabular` (0.4) and `ehr_tabular` (0.3) profiles of
  `benchmark_leakage_suite()` now have their declared prevalence.

* **`as_leaksplits()` for all splitGraph modes.** The splitGraph modes
  `site`, `region`, `platform`, `assay`, `relatedness` and `spatial` failed
  with "subscript out of bounds" (the `subject_grouped` fallback was
  unreachable), and `composite` failed with "'primary_axis' must be a list".
  These modes, including both composite strategies, now map to
  `make_split_plan(mode = "subject_grouped", group = "group_id")`, which
  keeps each splitGraph dependency group intact within a fold. Composite
  specs that require ordering, specs whose constraint collapses all samples
  into one group, and column-name clashes between `data` and the spec give
  clear errors; unrecognised future modes fall back to grouped folds with a
  warning. bioLeak now has its own `as_leaksplits()` tests.

* **Batch-confounding rule corrected for multiplicity.** The
  `confounding_alignment` rule of the audit mechanism summary took the raw
  minimum chi-square p-value over all batch columns and CV repeats, so a
  chance p = 0.038 in one of 10 repeats flagged an unrelated batch.
  `audit_leakage()` now adds a Holm-adjusted `pval_adj` column to
  `batch_assoc` (family = all batch-column x repeat rows), and the rule flags
  only when a row has `pval_adj <= 0.05` and Cramer's V >= 0.1. With a single
  repeat and one batch column the rule is unchanged.

* **Guarded variance/IQR filter uses the original scale.** `guard_fit()`
  applied `filter$var_thresh` / `filter$iqr_thresh` after z-scoring, where
  every non-constant feature has variance 1, making `var_thresh` meaningless
  under the default `normalize = "zscore"` (and distorted under `"robust"`).
  The thresholds, and the `min_keep` ranking, now use the variance and IQR of
  the imputed training data before normalization. The pipeline order and the
  default (`var_thresh = 0`, `iqr_thresh = 0`, no filtering) are unchanged.

* **Refit permutation null for stratified plans (results can change).**
  `audit_leakage(perm_refit = TRUE)` refitted on permuted outcomes but kept
  the observed folds, which `make_split_plan(stratify = TRUE)` had balanced
  on the true labels. Under permuted labels, and especially under
  group-restricted permutations, the fold class balance then varied
  widely. Each fold model's baseline tracks its training balance, which is
  anti-correlated with its test balance, so pooled AUC under the null was
  biased below 0.5. This is the stratification bias of Parker, Günter and
  Bedo (2007, *BMC Bioinformatics* 8:326). The result was an inflated gap
  and an anti-conservative p-value. In a pure-noise grouped design (40
  families of six plus 24 singletons), the fixed-fold null centred at
  0.40-0.46, with about two-thirds of draws below 0.5. On one simulated
  design, using re-drawn folds and fold-mean AUC cut the gap from 0.235 to
  0.158 and raised p from 0.0099 to 0.050. New argument
  `perm_folds = c("auto", "fixed", "redraw")`: `"redraw"` re-draws the split
  plan for each permutation with `make_split_plan()` on the permuted
  outcome. It uses the same mode, grouping columns, constraints, `v`,
  `repeats` and `stratify`, with seed `splits@info$seed + b`. The default,
  `"auto"`, re-draws whenever the plan is stratified and can be rebuilt.
  Unstratified plans keep their observed folds, so their results are
  unchanged.

* **Unrestricted refit nulls are no longer silent.** When the refit metadata
  lacks the outcome or the design column, the refit null is an
  unrestricted label shuffle that ignores the grouping. The common trigger
  is `perm_refit_spec$x` without the group column and no
  `perm_refit_spec$coldata`, which includes the refit data that
  `fit_resample()` stores by default. `audit_leakage()` now raises a
  `bioLeak_permutation_warning` in that case, naming the missing column. No
  warning is raised when the grouping has one sample per group (for example
  `group = "row_id"`).

## Behaviour changes and new arguments

* `audit_leakage()` gains `perm_summary = c("pooled", "fold_mean")`.
  `"fold_mean"` scores the observed fit and each permutation by the mean of
  per-fold metric values instead of the metric on predictions pooled across
  folds. Between-fold baseline shifts cannot bias it, and Parker et al.
  (2007) recommend it for AUC. Refit audits report both summaries in
  `audit_info(aud)$perm_gap_summaries`. The default, `"pooled"`, keeps the
  previous behaviour.

* `audit_info(aud)` now records the null used: `perm_null` is
  `"global_shuffle"`, `"refit_fixed_folds"` or `"refit_redrawn_folds"`
  (previously `"refit"` for every refit null). It also records
  `perm_folds`, `perm_folds_reason`, `perm_scheme` (`"group_restricted"`,
  `"within_batch"`, `"within_study"`, `"time_block"`, `"unrestricted"` or
  `"global_shuffle"`), `perm_summary` and `perm_gap_summaries`.
  `summary()` prints the fold handling, the permutation scheme and the
  summary statistic for refit nulls.

* `audit_leakage()` with `perm_refit = FALSE`: the permutation-gap null is,
  by design, a global label shuffle of the pooled out-of-fold predictions.
  A restricted permutation source was also built on this path but never
  used, so `perm_stratify`, `time_block` and `block_len` silently had no
  effect. The dead code is removed, a warning is now raised when any of
  these is set on the fixed-prediction path, the documentation says which
  null each argument affects, and `info$perm_null` records the null used
  (see the entry above for its values).

* `fit_resample()` gains `id_cols`, a character vector of identifier or
  metadata columns to exclude from the predictors (previously only the
  outcome and the split's group/batch/study/time columns were dropped, so a
  character `sample_id` was one-hot encoded into one predictor per row).
  With guarded preprocessing, a `bioLeak_input_warning` is raised when a
  character or factor column not listed in `id_cols` has at least 90% unique
  values. `id_cols` is stored in the fit and reused by refit-based
  permutations in `audit_leakage()`.

## Documentation

* `delta_lsi()`: the exact two-sided sign-flip test has minimum achievable
  p-value `2 / 2^R_eff`, so tier `C_signflip` at `R_eff = 5` can never
  reach p < 0.05 (minimum 0.0625); 6 or more paired repeats are needed. The
  tier documentation and vignette (which gave `1 / 2^R`) are corrected, the
  floor is stored in `info$min_p_achievable` (using `n_blocks` for
  `exchangeability = "blocked_time"`), and `summary()` notes when it exceeds
  0.05. Tier boundaries are unchanged.

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
