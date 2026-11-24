# bioLeak — open issues

Defects found by verification. IDs are stable and referenced by the probes in
`.claude/skills/*/probes/`. When fixing one, keep the entry, set **Status**
to `fixed` with the commit/file, and make its probe check pass. Remove an
entry only once it is fixed and released, or disproven (say how).

Severity: **critical** = wrong results or invalid inference; **major** =
silent wrong behaviour in a supported path; **minor** = edge case, naming,
docs.

Reproduce everything below with:

```
Rscript .claude/skills/verify-change/probes/inference_calibration.R
Rscript .claude/skills/verify-change/probes/guard_invariants.R
Rscript .claude/skills/contract-check/probes/handoff_e2e.R
```

---

## BL-1 — `delta_lsi(exchangeability = "by_group")` is anti-conservative
- **Severity:** critical · **Status:** fixed (uncommitted), `R/delta_lsi.R`
  `.dlsi_cluster_deltas()` eligibility check + cluster units drive R_eff;
  tests in `test-delta-lsi.R`; probe `inference_calibration.R` passes ·
  **Found:** 2026-11-11 review
- **Where:** `R/delta_lsi.R` `.dlsi_group_deltas()` and the `R_eff`
  upgrade block in `delta_lsi()`.
- **Problem:** each group gets the delta of the fold it was tested in, so all
  groups in one fold carry the same value; the sign-flip treats them as
  independent. `R_eff` is raised to `n_groups`, which can also promote the
  inference tier to `A_full_inference`.
- **Evidence:** 50 subjects, 5 folds, 1 repeat, H0 true → type-I 0.61 at
  α = 0.05 (limit 0.088); units used 50 vs 5 independent.
- **Not affected:** `study_loocv`, where each fold is one study (probe INV-1
  passes at 0.03).
- **Fix direction:** per-group deltas from each group's own test predictions
  (metric computed on that group's rows), or restrict the upgrade to designs
  where folds and groups coincide. Add the probe's null simulation as a test.

## BL-2 — `delta_lsi(exchangeability = "within_batch")` is anti-conservative
- **Severity:** critical · **Status:** fixed (uncommitted), same change as
  BL-1. `within_batch` now clusters on batch-like columns only (batch, then
  study), pooling all repeats per batch, so (repeat, batch) cells are gone.
  It remains a cluster sign-flip, not a stratified one; the docs say so ·
  **Found:** 2026-11-11 review
- **Where:** `R/delta_lsi.R` `.dlsi_batch_deltas()`.
- **Problem:** same replication as BL-1, over (repeat × batch) cells. It is
  also not a within-batch (stratified) flip: cells are flipped independently.
- **Evidence:** type-I 0.73 at α = 0.05; 10 units used vs 5 independent.
- **Fix direction:** per-batch metrics from that batch's predictions, flips
  constrained within strata; or document as unsupported and warn.

## BL-3 — frozen batch adjustment silently skips unseen batches
- **Severity:** major · **Status:** fixed (uncommitted), `R/guards.R`
  (`unseen` policy, `.guard_unseen_batch()`), `R/fit_resample.R` (one summary
  across folds); tests in `test-guard-fit.R`; probe BL-3 passes. Three fix
  cycles with two independent reviews; the last cycle's fixes were verified
  by tests, probes and the reviewer's scripts, not by a further review ·
  **Found:** 2026-11-11 review
- **Where:** `R/guards.R` `.guard_batch_apply()` (only `known` batches are
  adjusted).
- **Problem:** under batch-blocked / study-leave-out CV every test batch is
  unseen, so train is adjusted and test is not, with no condition raised.
- **Fix direction:** `unseen = c("warn", "error", "grand_mean")` (or similar)
  with a default that warns; document in `?guard_fit`.

## BL-4 — batch adjustment removes biology when batch ~ outcome
- **Severity:** major · **Status:** fixed (uncommitted) by decision: the
  default stays unprotected (conservative, attenuates) and now warns when
  batch and outcome are associated; `protect = TRUE` is opt-in. Real-fit
  checks: default adds no null inflation (CV AUC 0.497); protect = TRUE keeps
  the signal but inflates null CV AUC by about +0.03 to +0.10 (Nygaard et al.,
  2016), documented and warned. Under protect, outcome effects come from
  within-batch variation and only locations are adjusted (matches
  `sva::ComBat(mod = ~y)` class-by-batch means in review). Same review
  history as BL-3 · **Found:** 2026-11-11 review
- **Where:** `R/guards.R` `.guard_batch_fit()`.
- **Problem:** no protected covariates (ComBat's `mod`). With batch 85/15
  confounded with the outcome and no true batch effect, `location_scale`
  cut the class mean difference from 2.04 to 1.12 (55% retained).
- **Fix direction:** fit batch parameters on residuals after the outcome (or
  user covariates) within the training fold; at minimum warn when batch and
  outcome are associated in the training fold. Also reconsider the name
  `combat_frozen`: the shrinkage is a simplified EB, not ComBat's priors.

## BL-5 — panel fallback in grouped time series trains on very little data
- **Severity:** minor (valid but weak) · **Status:** resolved by documentation (uncommitted): the fallback keeps full test coverage and now warns with its training share (`bioLeak_fold_warning`). Staggered rolling-origin cutoffs were tried and reverted after independent review: they cut test coverage (70 vs 105 of 120 rows) and silently dropped groups on unbalanced panels · **Found:** 2026-11-11 review
- **Where:** `R/make_split_plan.R`, longitudinal panel fallback in
  `mode = "time_series"` with `group`.
- **Problem:** when every subject appears at every timepoint, each fold
  trains only on other subjects' earliest timepoint (23% of available rows in
  the probe) and tests on all later visits. Leak-free, but not rolling-origin.
  Also only triggers when *all* regular folds are empty, ignores `repeats`,
  and calls `set.seed()` on the global RNG.
- **Fix direction:** rolling origin per cutoff t: train on times < t from
  groups not in the test subset, test on time t.

## BL-6 — `use_block_constraints` silently ignored outside grouped mode
- **Severity:** minor · **Status:** fixed (uncommitted) together with
  splitGraph SG-8: warns when it cannot apply; also accepts a character
  vector of block_vars and explains infeasible combinations; warns when
  rows are never tested; tests in `test-as-leaksplits.R`. One independent
  review; its findings were fixed and re-checked with the reviewer's
  scripts · **Found:** 2026-11-11 review
- **Where:** `R/splitgraph_adapter.R` (`bl_mode == "subject_grouped"` only).
- **Problem:** with `ordering_required`, batch or study specs, the argument
  has no effect and no message says so.
- **Fix direction:** warn when `use_block_constraints = TRUE` cannot apply.

## BL-8 — identical repeated partitions counted as independent repeats
- **Severity:** critical · **Status:** fixed (uncommitted),
  `.dlsi_repeat_partition_keys()` + averaging in `delta_lsi()` ·
  **Found:** 2026-11-11 independent review, cycle 1
- **Problem (pre-dated Phase 1–2, default `iid` path):** with
  `batch_blocked` and `v = n_batches`, every repeat is the same partition;
  10 such repeats were 10 "independent" units. Null rejected 95% of the time
  (4 batches x 10 repeats). Now counted once: identical test sets are
  averaged (silently; the reduced count shows in `info$n_test_sets`), in
  `delta_lsi()` and, since the 2026-11-12 verification, in `cv_ci()`.

## BL-9 — repeats that share *some* test sets are still treated as independent
- **Severity:** critical · **Status:** fixed (uncommitted). Two independent
  reviews with real fits found no calibration failure in 20+ designs (AUC and
  RMSE, all split modes, compact, nested, stratified). `delta_lsi()` now uses the Nadeau–Bengio corrected t with
  df = folds per repeat − 1 for every exchangeability setting; identical test
  sets count once. Real-fit probe (`inference_calibration.R`): all designs
  ≤ 5% type-I, including this one (0.000–0.007) ·
  **Found:** 2026-11-11 independent review, cycle 2
- **Method decision (user, 2026-11-11):** option (b) below, with calibrated
  degrees of freedom. Real-fit checks showed plain NB was still liberal for
  repeated CV (0.09–0.19), and that the earlier group-level sign-flip was
  also anti-conservative (0.14–0.17 on leave-one-study-out), so it was
  retired.
- **Where:** `delta_lsi()` repeat-level (`iid`) inference; `within_batch`
  falls back to it.
- **Evidence:** batch_blocked, 6 batches, v = 3, 20 repeats (only 15 distinct
  test sets in 60 folds): type-I 0.217, mean CI coverage 0.70.
  subject_grouped, 8 subjects, v = 4, 10 repeats: type-I 0.150. Control with
  100 subjects (no recurrence): 0.057. Script:
  `review2/t3_designs.R` in the session scratchpad.
- **Problem:** BL-8 only merges repeats whose *whole* partition matches.
  With few groups, repeats re-draw the same test sets, so their deltas are
  positively correlated. More generally, repeated-CV repeats always share
  data; the iid sign-flip assumes they don't.
- **Options:** (a) use distinct test sets as cluster units (a fold delta per
  distinct test set, each counted once); (b) a variance correction for
  repeated CV (Nadeau–Bengio style) instead of a plain sign-flip; (c) warn
  and cap the tier when distinct test sets are few relative to folds.

## BL-10 — one all-NA repeat in one fit disables pairing entirely
- **Severity:** major · **Status:** fixed (uncommitted): pairing now checks
  fold structures only, and pairs at the fold level; test in
  `test-delta-lsi.R` (pre-dates Phase 1–2) ·
  **Found:** 2026-11-11 independent review, cycle 2
- **Where:** `delta_lsi()` pairing: `same_count` requires equal repeat counts
  before the `intersect()` path, so that path and its "repeats dropped"
  warning are unreachable. 25 repeats → R_leaky 24 vs R_guarded 25 →
  `paired = FALSE`, R_eff 0, no p-value. Script: `review2/t8_na_repeat.R`.
- **Fix direction:** decide pairing from `.dlsi_splits_match()` on the
  common repeats, then intersect.

## BL-11 — compact splits lose fold-size weights in `delta_lsi()`
- **Severity:** major · **Status:** fixed (uncommitted): `.dlsi_fold_sizes()`
  and `.dlsi_fold_train_sizes()` resolve compact folds; compact and explicit
  plans now give identical results (test in `test-delta-lsi.R`) (pre-dates
  Phase 1–2) ·
  **Found:** 2026-11-11 independent review, cycle 2
- **Where:** `.dlsi_fold_sizes()` ignores `info$fold_assignments`, so every
  compact fold weighs 1. Same plan, same metrics: explicit gives
  delta_metric −0.00423 / p 0.258, compact gives −0.00523 / p 0.180, silently.
  Script: `review2/t2_misc.R`.
- **Fix direction:** size folds through `.dlsi_fold_test_idx()`.

## BL-12 — NA fold metrics are not dropped pairwise
- **Severity:** minor · **Status:** fixed (uncommitted): a fold missing in
  either fit is dropped from both, with a warning (pre-dates Phase 1–2) ·
  **Found:** 2026-11-11 independent review, cycle 2
- **Problem:** a fold NA in one fit only is dropped from that fit's repeat
  mean but kept in the other's, so the "paired" delta compares different
  test sets, with no warning. Script: `review2/t2_misc.R`.

## BL-13 — BCa CI undercovers at R_eff around 10
- **Severity:** minor · **Status:** superseded: the BCa bootstrap was removed;
  CIs now come from the corrected t, and the real-fit probe checks that the
  CI excludes 0 no more often than the test rejects (pre-dates Phase 1–2) ·
  **Found:** 2026-11-11 independent review, cycle 2
- **Evidence:** independent repeats, R_eff 10: delta_metric CI coverage 0.857
  (Huber 0.913) vs nominal 0.95; tier B shows it without a caveat. At R_eff
  20 coverage is 0.94.
- **Fix direction:** add a caveat in tier B, or raise the CI threshold.

## BL-14 — `cv_ci()` intervals too narrow for repeated CV
- **Severity:** major · **Status:** fixed (uncommitted): `cv_ci(folds_per_repeat = )`, NB uses df = folds per repeat − 1 and no normal switch; real-fit coverage of the 95% CI for 5-fold × 10: 0.980 (was 0.920) ·
  **Found:** 2026-11-11 while reusing its Nadeau–Bengio helper
- **Where:** `R/cv_uncertainty.R` `cv_ci()`.
- **Problem:** df = (all folds − 1), switching to `qnorm` at ≥ 30 folds. With
  repeated CV that counts repeats as independent units, the same flaw fixed
  in `delta_lsi()` (BL-9); in real-fit simulations plain NB with df = J − 1
  rejected 9–19% at nominal 5%.
- **Fix direction:** same calibrated df as `delta_lsi()` (folds per repeat −
  1), which needs repeat ids; `fit@metrics` does not carry them, so take the
  `LeakSplits` or a `repeat_id` column.

## BL-15 — `tune_resample()` asks for NB but gets no correction
- **Severity:** major · **Status:** fixed (uncommitted): `tune_resample()` passes outer-fold sizes and folds per repeat; `cv_ci()` warns when NB lacks sizes ·
  **Found:** 2026-11-11
- **Where:** `R/tune_resample.R` ~line 1130:
  `cv_ci(metrics_df, method = "nadeau_bengio")` without `n_train` / `n_test`,
  so `.nb_corrected_var()` silently falls back to `var / K`, the uncorrected
  variance.
- **Fix direction:** pass mean train/test sizes from the splits; make
  `cv_ci()` warn when `nadeau_bengio` is requested without them.

## BL-16 — compact and explicit time_series plans train on rows in a different order
- **Severity:** minor · **Status:** fixed (uncommitted): all compact paths rebuild time-series training rows in time order; test in `test-make-split-plan-compact.R` ·
  **Found:** 2026-11-11 independent review (round 4, NB engine)
- **Where:** `R/fit_resample.R` `resolve_fold_indices()` (and the duplicate
  inline logic in `R/tidymodels.R`): explicit plans store `train` in time
  order, compact plans rebuild it in row order. Same set, different order.
- **Effect:** learners that depend on row order (e.g. ranger with a fixed
  seed) give different results for the same plan in compact vs explicit form.
  `delta_lsi()` is unaffected (it only uses the training-set size).
- **Fix direction:** have `resolve_fold_indices()` and the tidymodels path
  call `.bio_compact_train_indices()` and order its result the same way the
  explicit builder does; that also removes the duplicated logic.

## BL-17 — `cv_ci()` default method badly undercovers for repeated CV
- **Severity:** major · **Status:** fixed (uncommitted) by decision (user, 2026-11-12): default is now `"nadeau_bengio"`; `cv_ci()` accepts a `LeakFit` and takes fold sizes and folds per repeat from it; vignette and NEWS updated; tests in `test-cv-uncertainty.R` ·
  **Found:** 2026-11-12 real-fit coverage check for BL-14
- **Evidence:** 5-fold × 10 repeats, n = 100, glm: 95% CI coverage of the
  true expected AUC was 0.373 with the default `method = "normal"` (vs 0.980
  for `"nadeau_bengio"` with calibrated df). `fit_resample()` and
  `tune_resample()` already use Nadeau–Bengio; direct `cv_ci()` calls and the
  vignette example use the default.
- **Options:** make `"nadeau_bengio"` the default (it needs `n_train` /
  `n_test`, or a `LeakFit` to derive them from); or keep the default and warn
  when the metrics come from more than one repeat.

## BL-18 — main grouped time-series path can test only a fraction of rows
- **Severity:** minor · **Status:** fixed (uncommitted): later blocks dropped for lack of training data are reported (`bioLeak_fold_warning`) with reason and tested-row count; test in `test-make-split-plan-modes.R`. Verified by tests and probes; no separate independent review for this small change ·
  **Found:** 2026-11-12 independent review (BL-5 batch)
- **Problem:** with groups that end at different times (half end at t = 6 of
  10), `make_split_plan(mode = "time_series", group = , v = 4)` keeps 1 of 4
  folds and tests 16 of 64 rows; the panel fallback is not used because one
  fold exists. Fold loss is only partly reported.
- **Fix direction:** report coverage (rows tested / rows) whenever folds are
  dropped, as `as_leaksplits(use_block_constraints = )` now does.

## BL-19 — survival C-index was inverted
- **Severity:** critical · **Status:** fixed (uncommitted) · **Found:** 2026-11-12 verification (independent review)
- **Problem:** `.cindex_survival()` used `survival::concordance()`'s default
  (higher prediction = longer survival) while survival predictions are risk
  scores; a true risk score scored 0.226 instead of 0.774, so every survival
  metric, audit and `delta_lsi()` result was flipped. The parsnip/workflow
  survival path also never worked (it requested unsupported prediction types).
- **Fix:** `reverse = TRUE`; parsnip/workflow survival predictions use
  `linear_pred` (else `time`), oriented as a risk score on the training fold
  (see BL-20, second pass); the risk-score
  convention is documented for custom learners; tests pin the orientation.

## BL-20 — verification round findings (2026-11-12), all fixed (uncommitted)
- `cv_ci()` now counts repeated identical test sets once (was 25% too narrow
  for repeated leave-one-batch-out) and accepts a `LeakTune`.
- `plot_dlsi_repeats()` shows sign-adjusted deltas (were sign-flipped for
  lower-is-better metrics), returns the object's `R_eff` and `ci`; BCa
  wording removed from plot and class docs.
- Grouped time series: when regular folds lose blocks, the panel fallback is
  also built and the plan testing more rows is kept (a real panel gave one
  fold training on 2 rows); tied times named as a drop reason; group advice
  only for grouped series; `make_split_plan()` docs updated.
- `as_leaksplits(use_block_constraints = )`: coverage warning now also for
  compact plans; advice suggests a smaller v when an axis has fewer levels.
- Batch adjustment: bias figure consistent (+0.03 to +0.10); disagreeing
  `batch` / `batch_col` warns; `keep_batch_col` on a split column warns;
  `fit_resample()`'s "No successful folds" error includes the fold errors.
- `stratify_by` missing column named in the warning. (The change to store
  `info$stratify_by` only when used was reverted in BL-21: it let the column
  into the predictors.)
- `nested = TRUE`: an inner split that cannot be built no longer fails the
  whole plan (NULL for that fold, one warning); inner warnings reported once.

### Second verification pass (2026-11-13), fixed (uncommitted)
- **major** Plan choice ignored training size: the "more rows tested" rule
  replaced a staggered-entry regular plan (training 15-79 rows) by a fallback
  training on far fewer. Plans are now scored by `sum(test_n * train_n)`;
  tests in `test-make-split-plan-modes.R`.
- **major** Survival sign: negating parsnip output assumed "higher = longer
  survival"; censored rpart `time` rises with hazard (C = 0.26 for a strong
  model). `.bio_surv_risk_pred()` decides the direction on the training fold;
  rpart/bag_tree/survreg/Cox all give C about 0.74-0.76 on the reviewer's data,
  null C stays at 0.5; tests in `test-utils-metrics.R`.
- **low** "No successful folds" gave no reason when all metrics were NA; now
  lists skip reasons (`test-fit-resample-edge.R`).
- **low** `metric_summary` `_sd` covers all folds while the CI counts distinct
  test sets; documented in `?fit_resample` and `?tune_resample`.
- Independent re-review (5 minor, fixed): plan-choice wording said "tests
  more rows" (now states the score); fallback folds skipped the 3-test-row
  minimum; "were all empty" when regular folds were only too small; survival
  direction defaulted silently when undecidable (now NA + warning); missing
  separator in the "No successful folds" message; CI centre documented.
  `censored` added to Suggests (R CMD check WARNING from the new tests).
- Third review (6 low, fixed; no correctness or leakage defects): undecided
  survival direction warned once per fold (now one summary warning; fold
  reason `survival_direction_undecided`; training C exactly 0.5 counts as
  undecided); multi-line fold errors and more than 3 errors in the "No
  successful folds" message; empty grouped time-series plan gave a generic
  error (now gives lost-block reasons and fallback skips); fallback warning
  claimed "groups recur" whatever the drop reason (now lists reasons);
  `@param group` wording. Per the three-fix-cycle rule, R4 was not re-run
  after this cycle.

## BL-21 — findings of the 2026-11-13 re-verification (older code), fixed (uncommitted)
- **major** `fit_resample()` / `tune_resample()` copied `cv_ci()` interval
  columns into `metric_summary` by position; `aggregate()` sorts learners, so
  with learners in non-alphabetical order the intervals were swapped. Now
  matched by learner (`test-cv-uncertainty.R`).
- **major** rsample repeated designs: `id` (repeat) and `id2` (fold) were read
  the wrong way round; type-I 0.127 (limit 0.088) for 5 folds x 20 repeats.
  `as_rsample()` also dropped repeats. `.bio_rsample_fold_repeat_ids()` and
  `as_rsample()` follow rsample's convention; legacy exports still read
  (`test-tidymodels-interop.R`).
- **minor** `info$stratify` was TRUE for unstratified rolling-origin
  time-series plans and for a missing `stratify_by` column; now FALSE with a
  warning, and `as_leaksplits()` no longer requests it for ordered specs
  (`test-make-split-plan-modes.R`).
- **major (regression, caught in review)** BL-20 stored `info$stratify_by`
  only when stratification was applied; fit_resample()/tune_resample() use it
  to exclude that column from the predictors, so an outcome-derived stratum
  became a predictor whenever stratification was dropped (AUC 0.52 -> 0.89 in
  the repro). `stratify_by` is again always recorded; test covers time series
  and stratify = FALSE.

## BL-22 — findings of the 2026-11-16 full re-verification (older code), fixed (uncommitted)
- **critical** Recipe preprocessing prepped on data with split columns and
  dropped only exact names afterwards; `step_dummy()` outputs of group /
  stratify_by columns reached the learner (CV AUC 1.0 on noise). Helper
  `.bio_recipe_exclude_cols()` gives them role "bioLeak_split" before prep in
  fit_resample (recipe and recipe workflows) and tune_resample
  (`test-fit-resample-edge.R`).
- **major** tune_resample with a parsnip spec and no recipe fitted
  `outcome ~ .` on the full data frame, so the group column was a predictor
  and outer folds were skipped (6/15). Split columns are removed from the
  modelling data; recipes on the full data now work (`test-tune-resample.R`).
- **medium** audit_leakage redraw lost `info$sample_id_col` from
  as_leaksplits plans; permutation refits encoded sample_id (all NA with glm).
  Carried over in `redraw_splits` (`test-audit-perm-folds.R`).
- Review of these fixes: (a) **regression** tune_resample stripped split
  columns for every workflow, breaking recipe workflows built on the full
  data (incl. `update_role(subject, "id")`); now handled like fit_resample.
  (b) **incomplete** roles do not stop type selectors
  (`step_dummy(all_nominal(), -all_outcomes())`, AUC 1.0 on noise);
  `.bio_recipe_meta_derived()` now detects predictors derived from split
  columns after prep and stops the fold. Tests in `test-fit-resample-edge.R`
  and `test-tune-resample.R`.
- **fixed (older, found in the same review)** `audit_leakage(perm_refit =
  TRUE)` errored ("Batch column not found for batch_blocked mode") for batch
  and study plans from `as_leaksplits()`, and subject plans fell back to an
  unrestricted shuffle: the refit metadata came from the user's `x`, which
  lacks the spec's design columns. `audit_leakage()` now joins the missing
  design columns from `info$coldata` by `info$sample_id_col`; p-values equal
  those of the equivalent `make_split_plan()` plans; null type-I 0.060 (batch)
  and 0.080 (subject), limit 0.103, R = 150, B = 39. The join also covers
  user-supplied coldata; a missing design column for batch/study/time now
  falls back to the documented warned unrestricted shuffle instead of
  erroring (`test-audit-perm-folds.R`).
- **fixed (older, found reviewing the above)** `perm_stratify = TRUE`/"auto"
  made batch_blocked and study_loocv refit nulls the identity (permuting
  within batch x outcome cells), so p = 1. Now ignored there with a warning
  (`R/permute_labels.R`).

## BL-23 — findings of the 2026-11-16 four-area re-verification (older code), OPEN
All earlier fixes (BL-1 to BL-22) re-verified as working. New findings, each
reproduced independently (scripts in the session scratchpad, verifyI/ verifyJ/):
- **FIXED (uncommitted)** the four `tune_resample()` items below (inverted
  AUC, 0/1 outcome, as_leaksplits plans, rsample input): event level set via
  `control$event_level`; resolved outcome written into the tuning data;
  inner splits built from `info$coldata` when it holds the design columns;
  NA seed handled; nested_cv inner resamples imported as global indices;
  plain rsets build inner folds with the design inferred at import
  (`.bio_perm_mode()`), as fit_resample() accepts them; metric sets keep
  their names (do.call() had labelled them "<prb_mtrc>", so dropping AUC
  for single-class inner folds dropped every metric and errored, e.g. with
  bootstraps()); nested_cv() no longer needs an inferable design column.
  Note: for bootstraps / mc_cv the default inner_v equals the number of
  resamples, which is fine-grained; AUC is then dropped from tuning with a
  warning (pass inner_v = 5).
  Tests in `test-tune-resample.R`.
- **fixed (uncommitted)** a recipe or workflow combined with
  `positive_class` set to the first level failed ("wrong attributes: `y`"):
  recipe/workflow data now keep the original level order; tune's event
  level follows the positive class (`test-fit-resample-edge.R`). Review
  caught that this exposed a positional fallback for the positive-class
  probability column (inverted AUC with level names like "1"/"2");
  `.bio_prob_col()` now matches exact names, else errors.
- **fixed (uncommitted)** multiclass yardstick metric sets with roc_auc
  failed: probability columns were passed as named arguments
  (`test-fit-resample-edge.R`).
- **critical** `tune_resample()` inner AUC is inverted: the tweaked
  `roc_auc`/`auc`/`pr_auc` metrics use `event_level = "second"` but tune passes
  the first class's probability column. A model with true AUC 0.85 scores 0.14
  in tuning; every outer fold picks the worst setting (largest penalty).
  `R/tune_resample.R:372-375`.
- **FIXED (uncommitted)** the audit items below, after one review cycle that
  caught three problems in the first fix (cross-size label exchange changed
  class balance: accuracy type-I 1.0; by_design covered recorded but
  unblocked batch columns; unit-variance scaling hid count-data replicates).
  Final: `.design_perm_donor()` swaps whole label vectors only between
  same-size subjects (lone sizes shuffled within), used by the
  fixed-prediction null, the scan and the refit subject null, with coverage
  and degenerate-null warnings; time-series fixed null in time blocks; scan
  uses refit-style restricted permutations for batch/study/time; fold_mean
  with fixed predictions falls back to pooled; duplicates: centred columns,
  length agreement, outlier cut-off; by_design limited to the blocked axis;
  summary/report/vignette text. Calibration (R = 120-150, real glm fits,
  limit 0.103-0.110): 5 repeats 0.058 / scan 0.042; clustered 0.042 / 0.025;
  subject size ~ outcome: accuracy 0.008, pr_auc 0.008, auc 0.033, refit
  accuracy 0.017; AR(0.9) time series 0.083; leave-batch-out 0.020 / scan
  0.073. Power: 5 repeats 1.00, leave-batch-out 0.98.
  Third cycle (second review: combined plans and time series with group
  shuffled freely, outlier cut-off for duplicates flagged clustered data,
  by_design missed nested columns; paired designs gave a degenerate null):
  unit axis via `.audit_unit_group_col()` (group, combined primary axis,
  time_series group); within-subject shuffle added; outlier cut-off removed
  (noisy count-data replicates are again only found by lowering
  sim_threshold); nested columns count as by_design. Calibration (reviewer's
  scripts, R = 150, limit 0.103): combined 0.020 / scan 0.007, time series
  with group 0.013 / 0.033. Stopped here per the three-cycle rule; no
  further independent review of this cycle.
- **major** fixed-prediction permutation p-value (`perm_refit = FALSE`, or
  "auto" without a stored refit spec / with B > 200) is anti-conservative with
  repeated CV (type-I 0.16, limit 0.115) and clustered samples (0.19); labels
  are shuffled freely over repeated / correlated predictions.
- **major** multivariate target scan (on by default) shuffles labels freely in
  grouped designs; clustered-null type-I 0.31 (limit 0.115).
- **FIXED (uncommitted)** imputation and stratified-fold items below: the
  requested imputation is always fitted (`test-guard-fit.R`); stratified
  units are dealt round-robin across classes via `.bio_stratified_deal()`
  (`test-make-split-plan-modes.R`).
- **major** guarded imputation fits `impute = "none"` when the training fold has
  no NA, leaving test NAs: a study that never measured a variable is dropped
  silently in study_loocv (no warning, summary over 3 of 4 studies).
- **major** `tune_resample()` fails on `as_leaksplits()` plans unless
  `nested = TRUE` ("Column 'group_id' not found"), and on every rsample input
  (`info$seed` NA reaches `set.seed()`; nested_cv inner folds not read).
- **major** stratified fold assignment (`sample(rep(seq_len(v), ...))` per
  class) silently loses folds when a class has fewer groups than v (10
  requested, 6 produced, `info$v` still 10) and can worsen balance.
- **medium** `tune_resample()` fails on a numeric 0/1 outcome that
  `fit_resample()` accepts; the error blames stratification.
- **medium** duplicate detection uses uncentred, unscaled cosine: ordinary
  tabular data (age, bmi, ...) gives 1803 of 11175 false "duplicate" pairs.
- **medium** batch association flags leave-batch-out / leave-study-out plans as
  `confounding_alignment` (p = 1e-110), the package's own recommended design.
- **medium** `summary(<LeakAudit>)` interpretation uses the gap size only
  ("Strong non-random signal" at p = 0.15); audit_report.Rmd says a
  non-significant gap "warns of leakage" (backwards); "Batch / Study
  Association: none detected" when nothing was tested.
- **minor** README advertises removed BCa intervals for `delta_lsi()`.
- **unreplicated** subject_grouped refit null shuffles subjects with a unique
  group size element-wise; reviewer measured type-I 0.147 (limit 0.117, R=95)
  for small cohorts of uneven sizes.

## BL-24 — fewer groups than folds gives fewer folds without a warning, OPEN
Found 2026-11-19 while checking J8 (JSS revision). `make_split_plan(mode =
"subject_grouped", v = 5)` with 3 groups returns 3 folds and no warning
(stratified or not, either `stratify_method`); `show()` prints "Total folds:
3". Predates J8. Decide: warn (as the time-series path does) or error.

## BL-25 — proxy_target_leakage flag ignores the FDR adjustment, OPEN
Found 2026-11-20 in the JSS revision review (L3.5). `.bio_mechanism_summary()`
flags `proxy_target_leakage` when any `flag_fdr` is TRUE, else when any
unadjusted `flag` is TRUE, so `target_p_adjust = "BH"` does not change the
flag (reproduced: `flag_fdr` FALSE, adjusted p 5.9e-7 at `target_alpha =
1e-9`, flag TRUE). Pre-existing; the help page now states the rule as it is.
Decide whether the flag should use `flag_fdr` when an adjustment is chosen.

## BL-7 — Phase 1–2 behaviour changes undocumented
- **Severity:** major (release blocker) · **Status:** fixed (uncommitted): NEWS has Breaking changes and New features sections for Phase 1–2 · **Found:** 2026-11-11 review
- **Problem:** `NEWS.md` has no entry for `stratify_by`, grouped time series,
  `use_block_constraints`, cluster/batch sign-flip, new normalizations,
  `batch_adjust`, or the recipe checks. `as_leaksplits()` now requires the
  data's sample IDs to equal the spec's exactly; previously extra data rows
  were allowed. That is a breaking change. All of it is uncommitted.
