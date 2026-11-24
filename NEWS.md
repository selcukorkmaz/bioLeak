# bioLeak 0.4.0

## Dependency changes

* SummarizedExperiment moved from Imports to Suggests, so
  `install.packages("bioLeak")` no longer fails when the Bioconductor package
  is not installed (JSS review, J10). It is needed only for
  SummarizedExperiment input; such input without it now stops with an
  install hint.

## Bug fixes (results can change)

* **The batch confounding alarm tested the wrong table** (JSS review,
  L1.1-L1.3). `confounding_alignment` was flagged from a fold x batch
  chi-square test, which ignores the outcome: it could not detect batch
  confounding and fired for blocked designs (V = 1 by construction). It now
  flags only when both conditions for batch-driven optimism hold: batch is
  associated with the outcome (new batch-outcome test in `batch_assoc`:
  chi-square with Cramer's V for classification, Kruskal-Wallis with
  epsilon^2 for regression; run on the plan's independent units, one
  observation per subject or group when the plan records one, including a
  subject constraint of a combined plan; time-series plans without a group
  and survival outcomes are not tested; Holm across batch columns,
  p <= 0.05 and at least a small effect) and batches are shared by training
  and test rows (new
  `span_frac` > 0). Blocked designs have `span_frac = 0` and are never
  flagged. The fold x batch test is still reported and described as a
  fold-batch test. When batch varies within subjects, association within
  subjects (e.g. paired tumour and normal samples in different batches) is
  not tested; its row-level effect size is reported (`within_unit_effect`)
  and printed by `summary()`. `audit_leakage_by_learner()` prints the
  batch-outcome evidence instead of the fold-batch table. In simulations
  (200 datasets per design) the null flag rate was 0.000-0.080 (limit
  0.096) across 14 designs: one row per subject, repeated rows per subject,
  batch varying within subjects (equal, unequal, 1-2 samples per subject,
  time-ordered, pairs mixed with singletons), multiclass outcomes, rsample
  grouped folds, leave-study-out and combined plans with a subject column,
  and time series; blocked designs were never flagged; under a clear
  batch-outcome association with shared batches the flag rate was 0.98-1.00.

* **Plots were not drawn in knitted documents.** The `plot_*()` functions
  printed only in interactive sessions, so the vignette and rendered reports
  showed none of their figures; they now also print while knitting.
* `summary(<LeakAudit>)` reported "Repeats: 1" for repeated plans (it read
  the fit's metadata instead of the plan's).
* `delta_lsi()` now warns when `M_boot`, `M_flip` or `seed` are supplied:
  they belonged to the former bootstrap / sign-flip inference and have no
  effect on the corrected t-test (they were silently ignored).
* **Imputation did nothing when a training fold had no missing values.** The
  guard then recorded "no imputation" whatever method was requested, so
  missing values in the test rows stayed NA and those predictions were
  dropped silently (in leave-study-out CV, a study that never measured a
  feature disappeared from the results without a warning). The requested
  method is now always fitted; missForest has no per-column model then and
  falls back to the training medians with a warning.
* **Stratified fold assignment lost folds.** Each class was dealt to folds
  starting from fold 1, so classes with fewer groups than `v` left
  high-numbered folds empty (6 of 10 requested folds) and fold sizes
  uneven. Units are now dealt round-robin across classes, so all folds are
  kept and get an almost equal number of groups, with each class spread as
  evenly as possible (per-fold class proportions can vary slightly more
  than before, the trade-off for equal fold sizes). Stratified plans can
  differ from earlier versions for the same seed.
* **`positive_class` set to the first level failed with recipes and
  workflows** ("wrong attributes: `y`"): the outcome was releveled after the
  recipe captured its level order. Recipes and workflows now receive the
  original level order, and `tune_resample()` points tune's event level at
  the positive class. Probability columns are now looked up by the exact
  level name (`.pred_1`, `.pred_case 1`) before its `make.names()` form, and
  a missing positive-class column is an error: falling back to "the second
  column" would have scored the negative class for non-syntactic level
  names (AUC 0.09 instead of 0.91 in a test).
* **Multiclass yardstick metric sets with class probabilities failed**
  ("Failed to compute `roc_auc()`") because the probability columns were
  passed as named arguments.
* **Audit p-values were too small with repeated CV, clustered samples or
  autocorrelated series.** The fixed-prediction permutation gap
  (`perm_refit = FALSE`, or "auto" without a stored refit specification or
  with `B > 200`) shuffled pooled prediction rows independently: a sample's
  copies from repeated CV got different labels, a subject's correlated
  samples were treated as independent, and time order was ignored (type-I at
  alpha = 0.05: 0.16-0.23 with 5 repeats, 0.17-0.19 with clustered samples,
  0.27 with AR(0.9) series). Labels are now permuted per sample (copies share
  one label); in subject-grouped plans (and combined plans, via their
  primary axis, and time-series plans with a group) whole subjects swap
  label vectors with subjects of the same size and labels are shuffled
  within subjects, which keeps the class balance exact and lets paired
  designs (one case and one control sample per subject) flip within pairs; a
  warning says when few subjects can be exchanged; in time-series plans
  without a group labels move in time blocks. Type-I: 0.058, 0.042, 0.083 (limit 0.11). The multivariate
  target scan, which refits on each permuted outcome, used a free shuffle in
  grouped, batch and study designs (0.31 with clustered samples); it now uses
  the same subject-level scheme or the refit audits' restricted permutation
  (0.025). The refit null's element-wise shuffle across subjects of unique
  size (0.147 in small uneven cohorts) is replaced by the same scheme. When
  every permutation gives the same value (e.g. subject size determines the
  outcome) a warning says the p-value is uninformative.
  `perm_summary = "fold_mean"` is not calibrated for fixed predictions
  (0.11-0.14) and now falls back to "pooled" with a warning.
* **Duplicate detection flagged ordinary tabular data.** Rows were compared
  on uncentred columns, so unrelated patients looked near-identical (1803 of
  11175 pairs without duplicates). Columns are now centred and pairs must
  also lie at a similar distance from the column means (0 false pairs with 6
  tabular features; true copies are still found). Pairs must also be much
  closer than typical pairs (within 10% of the median distance), and with
  fewer than 5 features, where close pairs arise by chance, only identical
  rows are reported. Noisy technical replicates (e.g. low-depth
  re-sequencing) can fall below `sim_threshold`; lower it to look for them.
* **Leave-batch-out and leave-study-out plans were flagged as confounded.**
  The batch association test flagged the column the split blocks on, which
  aligns with the folds by design. Rows for the axis the split mode blocks on
  (batch for batch_blocked, study for study_loocv, group, constraint axes),
  and columns nested within that axis (e.g. site within study), are now
  marked `by_design` and do not trigger `confounding_alignment`; a batch
  merely recorded on a subject-grouped plan is still tested.
* **The audit summary could overstate the result.** Its interpretation used
  the gap size alone ("Strong non-random signal" at p = 0.15); it now prints
  and uses the permutation p-value, and notes that a significant gap does not
  rule out leakage. The report template said a non-significant gap "warns of
  leakage" (leakage makes the gap larger); corrected. "Batch / Study
  Association: none detected" is now "not tested" when no column was tested.
* **`tune_resample()` tuned on an inverted AUC.** Its `auc` / `roc_auc` /
  `pr_auc` metrics set `event_level = "second"` on the metric, but tune passes
  the first class's probability column, so inner AUCs were 1 minus their
  value (0.14 for a model with AUC 0.85) and the worst setting was chosen in
  every outer fold (the default metrics). The event level is now set through
  `control$event_level`, and the plain metrics are used.
* **`tune_resample()` input fixes.** A numeric 0/1 outcome (accepted by
  `fit_resample()`) now works; it failed every inner fit with a misleading
  stratification hint (with a recipe or workflow, a non-factor class outcome
  gets a clear error). Plans from `as_leaksplits()` now work without
  `nested = TRUE` (inner folds are built from the spec's design columns
  stored in the plan; before: "Column 'group_id' not found"). rsample input
  works: rsample plans record no seed, which reached `set.seed()` as NA, and
  `rsample::nested_cv()` inner resamples were never read. rsets without the
  `bioLeak_mode` attribute now use the design inferred at import (as in
  `fit_resample()`), and `nested_cv()` needs no inferable design column. When
  single-class inner folds force AUC out of tuning, the remaining metrics
  are kept (all were dropped and the call failed, e.g. with `bootstraps()`).
  Inner folds now also use the plan's `stratify_by`.
* **Split columns could reach the model through a recipe.** With
  `preprocess = <recipe>`, the recipe was prepped on data that still held the
  split columns (group, batch, study, time, `stratify_by`, sample id) and only
  exact column names were removed afterwards, so a standard
  `recipe(y ~ .) |> step_dummy(all_nominal_predictors())` passed their dummy
  columns to the learner (CV AUC 1.0 on pure noise with an outcome-derived
  stratum). Those columns now get the non-predictor role `"bioLeak_split"`
  before `prep()`, in `fit_resample()`, `tune_resample()` and for workflows
  whose recipe lists them; steps may still name them, but neither they nor
  anything derived from them is a predictor. Steps that select by type
  (`step_dummy(all_nominal(), -all_outcomes())`) bypass roles; such recipes
  are now detected after `prep()` (predictors that change when the split
  columns are shifted) and stop with an error naming the derived columns.
* **`tune_resample()` used the group column as a predictor** when given a
  parsnip model and no recipe (`outcome ~ .` on the full data), so held-out
  subjects were new factor levels and outer folds were silently skipped
  (6 of 15 in a test). Split columns are now removed, as in `fit_resample()`.
  A recipe defined on the full data (`y ~ .`), given as `preprocess` or in a
  workflow, also works in `tune_resample()` now; it failed with "Not all
  variables in the recipe are present".
* **`audit_leakage(perm_refit = TRUE)` on `as_leaksplits()` plans ignored or
  broke the design.** The refit permutation metadata came from the user's
  data, which lacks the spec's design columns (`group_id`, `batch_group`,
  `study_group`): batch and study specs failed ("Batch column not found for
  batch_blocked mode"), and subject specs fell back to an unrestricted label
  shuffle with a warning. Those columns are now joined from the plan by
  sample id (also into a user-supplied `perm_refit_spec$coldata`), giving
  the same restricted null as an equivalent `make_split_plan()` plan. When a
  design column is genuinely unavailable, batch, study and time designs now
  fall back to the documented unrestricted shuffle with a warning instead of
  erroring.
* **`perm_stratify` made batch and study refit nulls the identity.** Refit
  permutations for `batch_blocked` / `study_loocv` permute labels within each
  batch or study; stratifying that by the outcome permuted within
  batch x outcome cells, so no label moved and p was 1. `perm_stratify` is now
  ignored there with a warning.
* **`audit_leakage()` re-drawn plans kept the sample id as a predictor.** For
  plans from `as_leaksplits()`, `perm_folds = "auto"` re-drew the folds
  without the record that `sample_id` is an identifier, so every permutation
  refit encoded it as a predictor (with `glm`, all permutation values and the
  p-value were NA). The re-drawn plan now keeps it.
* **`metric_summary` intervals could belong to another learner.** With
  several learners given in non-alphabetical order (e.g.
  `c("ranger", "glmnet")`), `fit_resample()` copied `cv_ci()`'s interval
  columns by position into a summary sorted by learner name, so each learner
  could be shown another learner's interval. They are now matched by name.
* **Repeated rsample resamples had fold and repeat swapped.** rsample's
  repeated designs (`vfold_cv(repeats = )`, `group_vfold_cv(repeats = )`) put
  the repeat in `id` and the fold in `id2`; bioLeak read them the other way
  round, so `delta_lsi()` and `cv_ci()` used the wrong number of folds per
  repeat (e.g. 19 instead of 4 degrees of freedom for 5 folds x 20 repeats,
  too-small p-values and too-narrow intervals). `as_rsample()` dropped the
  repeat column for repeated plans (rsample's `manual_rset()` has none), so a
  round trip also lost the repeats. Both directions now follow rsample's
  convention; older exports with repeated `"FoldK"` ids are still read
  correctly. `tibble` is added to Suggests for the export.
* **`splits@info$stratify` now says whether stratification was applied.**
  Rolling-origin `time_series` folds are never stratified; `stratify = TRUE`
  there now warns and is recorded as `FALSE` (only the grouped panel fallback
  stratifies). A missing `stratify_by` column is also recorded as
  unstratified. A given `stratify_by` column is always recorded, applied or
  not, so `fit_resample()` / `tune_resample()` keep it out of the predictors. `as_leaksplits()` no longer requests stratification for
  ordered specs. Previously the plan claimed a stratification it did not have,
  and `audit_leakage(perm_folds = "auto")` re-drew such plans needlessly.
* **Survival C-index was inverted.** `survival::concordance()` treats higher
  predictions as longer survival, but bioLeak's survival predictions are risk
  scores, so every survival C-index (fit metrics, audits, `delta_lsi()`) was
  reported as 1 minus its value; a strong risk model scored about 0.2. The
  C-index now treats predictions as risk scores (higher = shorter survival),
  the convention custom learners are documented to follow. Survival models
  given as parsnip specs or workflows, which previously failed to predict,
  now use `linear_pred` (else `time`). Engines do not agree on the direction
  of these outputs (censored's rpart `time` rises with hazard, survreg's with
  survival), so the direction is decided on each training fold and the
  predictions are turned into a risk score; test rows are not used. When the
  direction cannot be decided (e.g. no comparable pairs in an all-censored
  training fold), that fold's predictions are NA and left out of the
  metrics, reported in one warning across folds.
* **`cv_ci()` counts repeated identical test sets once**, as `delta_lsi()`
  does, and accepts a `LeakTune`; repeated leave-one-batch-out intervals were
  too narrow.
* **`plot(<LeakDeltaLSI>)` shows deltas on the inflation scale** for
  lower-is-better metrics too (points were sign-flipped against the line and
  band), and returns the object's `R_eff`.
* **Grouped time series keeps the better of the regular folds and the panel
  fallback** (the plan that tests more rows), instead of e.g. a single fold
  training on 2 rows; tied times are named as a reason for dropped blocks.
* **`nested = TRUE`** no longer fails the whole plan when an inner split
  cannot be built for some outer fold.

* **`as_leaksplits(use_block_constraints = )` can enforce a chosen subset.**
  `spec$block_vars` lists the blocking annotations available in a splitGraph
  spec, which can be many (batch, study, site, ...). `use_block_constraints`
  now also takes a character vector naming the block_vars to enforce (`TRUE`
  still enforces all; unknown names are an error that lists the available
  ones). When excluding shared levels empties every fold's training set, the
  error now names the axes and their level counts and suggests a smaller set,
  instead of "No valid folds generated"; when only some folds survive and
  rows are never tested, a warning reports the coverage. The option is
  ignored for non-grouped specs with a warning (previously silently), and
  overriding `mode` / `constraints` through `...` warns that it no longer
  applies.

* **Batch adjustment no longer skips rows silently.** `guard_fit()`'s
  `steps$batch_adjust` now applies a policy, `unseen = c("warn", "error",
  "grand_mean")`, to rows whose batch was not seen in training, rows with a
  missing batch label, and `transform()` calls without batch labels. All of
  these were left unadjusted with no message. Under batch-blocked or
  study-held-out CV every test batch is unseen; `fit_resample()` reports this
  once ("in k of K folds ...") instead of per fold. A single training batch
  now applies no adjustment with a warning instead of failing.
* **Batch adjustment warns when batch and outcome are confounded, and can
  protect the outcome.** When batch is associated with `y` in the training
  data (p < 0.01), removing batch effects also removes part of the outcome
  signal, so performance is underestimated; this now warns. The new opt-in
  `protect = TRUE` estimates batch effects net of the outcome
  (`z ~ batch + y`, the role of ComBat's `mod`), which keeps the signal but
  can bias cross-validated performance upward (Nygaard et al., 2016): about
  +0.03 to +0.10 AUC with 80/20 confounding and no true signal in simulation
  (more with more features),
  and it warns about this too. The default stays unprotected, the
  conservative choice. Under `protect = TRUE`, outcome effects are estimated
  from within-batch variation only, only batch locations are adjusted (per-
  batch scaling would rescale the protected outcome effect), and when every
  batch holds a single class only differences between batches of the same
  class are removed. A batch that cannot be separated from the outcome at all
  is left unadjusted, with a warning. `fit_resample()` reports each kind of
  batch/outcome diagnostic once across folds.
* **`keep_batch_col = TRUE` keeps the batch column intact.** A kept batch
  column was batch-adjusted along with the features, which zeroed it out; it
  is now passed through unadjusted.
* **Batch scale is measured as in ComBat.** Each batch's scale is now taken
  relative to the pooled within-batch residual variance, not the grand
  variance that included between-batch differences. The old version
  over-inflated within-batch spread after adjustment whenever batches
  differed in location.

* **`delta_lsi()` inference replaced: corrected repeated-CV t-test.**
  Repeated cross-validation re-uses the same data, so fold and repeat
  differences are correlated. The previous sign-flip tests treated them as
  independent and were badly anti-conservative with real model fits under a
  true null: 70% type-I error at alpha = 0.05 for 5-fold x 10 repeats on 100
  samples, and 14-17% for the group-level test on leave-one-study-out designs.
  `delta_lsi()` now uses the Nadeau & Bengio (2003) corrected variance
  `(1/J + n_test/n_train) * var(fold deltas)` with a t reference on
  `folds per repeat - 1` degrees of freedom (calibrated df, Bouckaert 2003),
  for every `exchangeability` setting. Folds that re-use an identical test
  set are counted once. In real-fit simulations (subject-grouped,
  batch-blocked, leave-one-study-out and rolling-origin designs, 1-20
  repeats) type-I error stayed at or below 5%. CIs for `delta_metric` and
  `delta_lsi` use the same standard error; the BCa bootstrap is gone.
  - `R_eff` is now the number of folds per repeat (studies for
    `study_loocv`), not the number of repeats: repeats sharpen the estimate
    but add no independent units, so they no longer raise the tier.
  - Tiers: D below 3 folds per repeat or unpaired; C 3-9; B 10-19; A 20+.
    Tier labels are unchanged for compatibility; tiers A-C all report a
    p-value and CIs, and `inference_ok` is `TRUE` whenever they are finite.
  - `exchangeability` is stored and reported only. `block_size`, `M_boot`,
    `M_flip` and `seed` are ignored (a warning is issued for `block_size`).
  - `info` gains `test`, `se`, `df`, `rho`, `n_fold_deltas`, `n_test_sets`
    and `folds_per_repeat`; `return_details = TRUE` adds `fold_deltas`.
    `min_p_achievable`, `block_size_used` and `n_blocks` are now `NA`.
* **`delta_lsi()` pairs fits at the fold level.** A fold whose metric is
  missing in either fit is dropped from both, so the paired difference always
  compares the same test sets; previously one all-NA repeat in one fit made
  the comparison unpaired, and an NA fold was dropped from one side only.
  The per-pipeline summaries (`repeats_naive`, `repeats_guarded`,
  `info$metric_naive`, `info$metric_guarded`) use the same paired folds.
* **`delta_lsi()` sizes compact-split folds correctly.** Fold sizes come
  from `fold_assignments`, and training sizes are rebuilt the way
  `fit_resample()` builds them (past rows only for `time_series`, minus
  horizon / purge / embargo; secondary-constraint exclusions for
  `combined`). Compact plans previously weighted every fold equally and used
  "all non-test rows" as the training size, so the same plan could give
  different results in compact and explicit form.
* **`delta_lsi()` metric direction.** Auto-detection now treats `mape`,
  `smape`, `mase`, `*_cost`, `ece`, `fpr`, `nll`, `cross_entropy` and other
  loss/error names as lower-is-better (several were read as higher-is-better,
  flipping the sign of the result), recognises skill scores and
  explained-deviance (e.g. `brier_skill_score`) and yardstick names such as
  `kap`, `ccc` and `roc_aunp` as higher-is-better, and warns when a metric
  name is not recognised instead of silently assuming higher is better.

* **`cv_ci(method = "nadeau_bengio")` uses calibrated degrees of freedom.**
  Its t reference now uses `folds_per_repeat - 1` degrees of freedom (new
  argument) instead of all folds - 1, and no longer switches to the normal
  distribution at 30 or more folds; counting repeats as independent made the
  intervals too narrow for repeated CV (the same issue fixed in
  `delta_lsi()`); in a real-fit simulation of 5-fold x 10 CV, coverage of the
  95% interval rose from 0.92 to 0.98. `fit_resample()` and `tune_resample()`
  pass it. Invalid `n_train`, `n_test` or `folds_per_repeat` are now errors.
  Requesting
  `"nadeau_bengio"` without `n_train` / `n_test` now warns: the correction
  cannot be applied and the uncorrected variance is used.
* **`cv_ci()` defaults to `method = "nadeau_bengio"` and accepts a `LeakFit`.**
  The previous default, `"normal"`, treats folds as independent; in a
  real-fit simulation of repeated 5-fold CV its 95% interval covered the true
  expected performance only 37% of the time (vs 98% for the corrected
  method). Calls without a `method` now get the corrected interval. Passing
  the fit itself (`cv_ci(fit)`) supplies the fold sizes and folds per repeat;
  a bare metrics data frame without `n_train` / `n_test` warns. Use
  `method = "normal"` to get the old behaviour.
* **`tune_resample()` confidence intervals are corrected.** It requested the
  Nadeau-Bengio correction without train/test sizes, so the uncorrected
  variance was used silently; it now passes the outer folds' average sizes
  and folds per repeat.
* **Compact and explicit time-series plans train on the same row order.**
  Compact plans rebuilt the training rows in row order instead of time
  order, so order-dependent learners (e.g. ranger with a fixed seed) gave
  different results for the same plan.
* **Grouped time series reports lost folds.** In `make_split_plan(mode =
  "time_series", group = )`, a later time block whose training set comes out
  empty (every earlier row belongs to its own test groups, or horizon / purge
  / embargo leave nothing) was dropped silently. It now warns with the number
  of blocks lost, the reason, and how many rows are tested (e.g. 16 of 64
  where 48 were expected); rolling-origin CV's untested first block is not
  counted as a loss.
* **Grouped time series: the panel fallback says what it does.** When every
  group appears across the time range, `make_split_plan(mode = "time_series",
  group = )` uses a panel fallback that tests every group but trains only on
  rows before the earliest test time. It now warns with the training share it
  achieves (often about 20%) and suggests alternatives. (Staggered
  rolling-origin cutoffs were tried; they train on more history but leave
  rows, and on unbalanced panels whole groups, never tested.) When regular
  folds lose blocks, both plans are built and the one with the larger
  `sum(test rows x training rows)` is kept, so a fallback that tests more rows
  but trains on very few does not replace a well-trained regular plan; the
  warning names the plan kept and both plans' coverage. Fallback folds now
  follow the same minimum of 3 test rows as regular folds; the warning gives
  the actual reasons blocks were lost, and when neither plan yields a fold
  the error says why.
* `fit_resample()`'s "No successful folds were completed" error now gives the
  fold errors or the skip reason (e.g. every requested metric was NA on
  single-class test folds). The `metric_summary` docs now say that `_mean` /
  `_sd` summarize all fold rows while `_ci_lo` / `_ci_hi` count each distinct
  test set once.

## Breaking changes

* **`delta_lsi()` tier labels renamed** to describe the test they report:
  `"A_full_inference"` -> `"A_ttest_ge20_folds"`, `"B_signflip_ci"` ->
  `"B_ttest_10_19_folds"`, `"C_signflip"` -> `"C_ttest_3_9_folds"`
  (`"D_insufficient"` unchanged). The old names referred to the bootstrap and
  sign-flip inference that the corrected t-test replaced; code comparing
  `dlsi_tier()` or `@tier` with the old strings must be updated.

* **`as_leaksplits()` requires the data's sample IDs to match the spec's
  exactly.** Missing, duplicated or extra IDs (in either direction) are an
  error; previously extra rows in `data` were dropped silently by the join.

## New features

* `audit_leakage()` with the default `perm_refit = "auto"` now says what a
  refit costs before running it (`B` x folds model fits) and says when it
  falls back to fixed-prediction permutations because `B` exceeds
  `perm_refit_auto_max` (previously silent). The fallback message no longer
  calls the fixed-prediction test "slightly liberal" (JSS review, K1).
* `audit_leakage()` reports how many near-duplicate pairs are expected by
  chance and the ratio observed / expected (`info$duplicate_null`; JSS review,
  L3.4). The expected count applies the same rule to null rows drawn from a
  Gaussian copula with the features' marginals and rank correlation;
  `duplicate_null_reps` (default 3, `0` skips) sets the number of draws. On
  data without duplicates (40 datasets each) the mean expected count was
  within about 25% of the mean observed count for Gaussian, skewed and
  discrete features. `summary()` prints both counts and, in the
  interpretation, lists the flagged leakage mechanisms; the flag rules of the
  mechanism summary are documented in `?audit_leakage` (L3.5).
* `make_split_plan()` gains `stratify_method = c("groups", "rows")` (JSS
  review, J8). The default `"groups"` keeps the current assignment (groups
  labelled by majority class, dealt round-robin), which balances numbers of
  groups, not rows. `"rows"` places groups largest first in the fold that
  keeps each class's share of rows most even, as scikit-learn's
  `StratifiedGroupKFold` does. In simulations with unequal group sizes and
  mixed-outcome groups (40 groups, 5 folds, 200 seeds), the largest
  deviation of a fold's case share from the overall share averaged 0.022,
  against 0.134 for `"groups"` and 0.208 unstratified. Refit permutations,
  nested plans and `tune_resample()` inner splits re-use the recorded
  method.
* **`as_leaksplits()` covers more of the splitGraph contract.** Specs with
  `ordering_required = TRUE` map to `mode = "time_series"` with group
  separation (train strictly earlier than test, no group in both);
  `spec$stratum_var` is forwarded as `stratify_by`; additional block columns
  (`site_group`, `region_group`, `platform_group`, `assay_group`) are carried;
  and `use_block_constraints` enforces chosen `block_vars` as secondary axes
  through `make_split_plan(mode = "combined")`.
* **`make_split_plan()`** gains `stratify_by` (stratify on a column other than
  the outcome), group-disjoint `time_series` splits (`group =` with `time =`)
  with a longitudinal panel fallback, and more axis types for
  `mode = "combined"` (site, region, platform, assay, relatedness, spatial).
* **`guard_fit()`** gains train-only frozen normalisations
  (`normalize$method = "quantile"`, `"log_cpm"`, `"clr"`) and train-only batch
  adjustment (`steps$batch_adjust`; see the bug fixes above for its
  behaviour).
* **Recipe checks** flag recipes with already-trained steps and
  `step_mutate()` / `step_filter()` calls to global summaries (`mean()`,
  `sd()`, `quantile()`, ...), which leak test-set statistics at `bake()`.
* `fit_resample()` and `tune_resample()` drop split metadata columns
  (sample ID, stratification and constraint columns) from the predictors.

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
