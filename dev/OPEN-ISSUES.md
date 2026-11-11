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
- **Severity:** major · **Status:** open · **Found:** 2026-11-11 review
- **Where:** `R/guards.R` `.guard_batch_apply()` (only `known` batches are
  adjusted).
- **Problem:** under batch-blocked / study-leave-out CV every test batch is
  unseen, so train is adjusted and test is not, with no condition raised.
- **Fix direction:** `unseen = c("warn", "error", "grand_mean")` (or similar)
  with a default that warns; document in `?guard_fit`.

## BL-4 — batch adjustment removes biology when batch ~ outcome
- **Severity:** major · **Status:** open · **Found:** 2026-11-11 review
- **Where:** `R/guards.R` `.guard_batch_fit()`.
- **Problem:** no protected covariates (ComBat's `mod`). With batch 85/15
  confounded with the outcome and no true batch effect, `location_scale`
  cut the class mean difference from 2.04 to 1.12 (55% retained).
- **Fix direction:** fit batch parameters on residuals after the outcome (or
  user covariates) within the training fold; at minimum warn when batch and
  outcome are associated in the training fold. Also reconsider the name
  `combat_frozen`: the shrinkage is a simplified EB, not ComBat's priors.

## BL-5 — panel fallback in grouped time series trains on very little data
- **Severity:** minor (valid but weak) · **Status:** open · **Found:** 2026-11-11 review
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
- **Severity:** minor · **Status:** open · **Found:** 2026-11-11 review
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
  (4 batches x 10 repeats). Now averaged into one unit with a warning.

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
- **Severity:** major · **Status:** open (pre-dates Phase 1–2) ·
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
- **Severity:** major · **Status:** open (pre-dates Phase 1–2) ·
  **Found:** 2026-11-11
- **Where:** `R/tune_resample.R` ~line 1130:
  `cv_ci(metrics_df, method = "nadeau_bengio")` without `n_train` / `n_test`,
  so `.nb_corrected_var()` silently falls back to `var / K`, the uncorrected
  variance.
- **Fix direction:** pass mean train/test sizes from the splits; make
  `cv_ci()` warn when `nadeau_bengio` is requested without them.

## BL-16 — compact and explicit time_series plans train on rows in a different order
- **Severity:** minor · **Status:** open (pre-dates Phase 1–2) ·
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

## BL-7 — Phase 1–2 behaviour changes undocumented
- **Severity:** major (release blocker) · **Status:** open · **Found:** 2026-11-11 review
- **Problem:** `NEWS.md` has no entry for `stratify_by`, grouped time series,
  `use_block_constraints`, cluster/batch sign-flip, new normalizations,
  `batch_adjust`, or the recipe checks. `as_leaksplits()` now requires the
  data's sample IDs to equal the spec's exactly; previously extra data rows
  were allowed. That is a breaking change. All of it is uncommitted.
