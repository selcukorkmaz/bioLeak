# JSS revision tracker

Source of truth for the revision (CLAUDE.md section 7; workflow in the
`jss-revision` skill). Every change made for the revision has a row here.

Reviewer reports: `paper/submission/comments.rtf` (Reviewers J, K, L).
Changes forced by the 2026-11 package fixes: `dev/paper-updates.md` (P items).

Status: `proposed` -> `agreed` -> `in progress` -> `done` -> `in letter`
Type: `package` / `replication` / `text` / `format`
Sections and pages refer to the submitted manuscript.

## Items

### Reviewer J

| ID | Reviewer's point (short) | Type | Proposed decision + reason | Action | Location | Status |
|---|---|---|---|---|---|---|
| J1 | Mention mlr3resampling (grouped + stratified fold assignment) alongside mlr3 | text | Accept: relevant prior software | Add mlr3resampling to related work | Sec 1 / 2.2 (p. 2) | in letter |
| J2 | Can an mlr3 interface be written, or is there a barrier? Discuss | text | Accept, discussion only: no technical barrier (a `LeakSplits` plan is explicit train/test indices, as `as_rsample()` already exports) | Short paragraph: feasible via a custom resampling; not in this version | Sec 3 or 8 | in letter |
| J3 | Align decimal points and parentheses in Tables 2-5 | format | Accept | `dcolumn`/`siunitx` alignment | Tables 2-5 | in letter |
| J4 | Sec 2.2: note that StratifiedGroupKFold (scikit-learn) and mlr3resampling support grouped fold assignment | text | Accept | One sentence + references | Sec 2.2 | in letter |
| J5 | Software references for Secs 2.1, 2.3-2.5; ideally a table of software x data structure / leakage type | text | Accept: compact table | Add one table (software x leakage type) | Sec 1 or 2 | in letter |
| J6 | p. 3 text runs off the right margin | format | Accept | Fix line break | p. 3 | in letter |
| J7 | Too much code output (pp. 9, 12, 14-17); more discussion between chunks; maybe print less | text | Accept: trim printed output to the lines discussed | Shorten listings, add interpretation between them (overlaps K2) | Secs 4-5 | in letter |
| J8 | `subject_grouped` + `stratify = TRUE`: row counts per stratum and fold no better than random (unlike StratifiedGroupKFold, mlr3resampling); explain the algorithm | package + text | Accept: document the algorithm and add a row-balancing option (assign groups largest-first to the fold that keeps class proportions closest) | New option + docs + short paper text | `make_split_plan()`; Sec 3/4 | in letter |
| J9 | Discuss row-level vs group-level stratification goals; which is better when | text | Accept | Short paragraph (follows J8) | Sec 3/4 | in letter |
| J10 | `install.packages("bioLeak")` fails until SummarizedExperiment is installed from Bioconductor | package | Accept: move SummarizedExperiment from Imports to Suggests (8 `::` calls, no NAMESPACE import; only needed for SE input) | DESCRIPTION + `requireNamespace()` guards | DESCRIPTION, `R/utils_se.R`, `R/audit.R` | in letter |

### Reviewer K

| ID | Reviewer's point (short) | Type | Proposed decision + reason | Action | Location | Status |
|---|---|---|---|---|---|---|
| K1 | `perm_refit = "auto"` refits per permutation (expensive); simulations used `perm_refit = FALSE` and accepted ~10% false positives; add a clear warning | package + text | Accept: the fixed-prediction null is now calibrated (type-I 1-4%), so the trade-off text changes; add a message when "auto" selects a refit or falls back | Message in `audit_leakage()`; revise Secs 4.4, 6.2 (overlaps P2) | `R/audit.R`; Secs 4.4, 6.2 | in letter |
| K2 | Too much code; Sec 4.5 repeats 4.1-4.4 and 5 | text | Accept | Replace Sec 4.5 full repeat by a short combined listing | Sec 4.5 | in letter |
| K3 | Tables 2-3 vs figures in 6.3: standardise, prefer figures | text | Partly accept: figure for the main simulation; exact values to a compact table or appendix | Add figure; move/condense Tables 2-3 | Sec 6 | in letter |
| K4 | Features promoted but not available: survival excluded from `tune_resample()`; ΔLSI `by_group` / `within_batch` fall back to iid | text | Accept: state survival limits plainly; ΔLSI now uses one test for every `exchangeability` value, described as such (overlaps P1) | Rewrite Sec 5 and 9.2 statements | Secs 5, 9.2 | in letter |

### Reviewer L

| ID | Reviewer's point (short) | Type | Proposed decision + reason | Action | Location | Status |
|---|---|---|---|---|---|---|
| L1.1 | Batch diagnostic described inconsistently ("fold-batch" vs "batch-outcome") | text + package (docs) | Accept | One consistent description in paper and help page | Sec 4.4 (pp. 11, 13); `?audit_leakage` | in letter |
| L1.2 | Add a direct batch x outcome association test | package | Accept: confounding needs batch linked to outcome | Batch x outcome test reported alongside | `R/audit.R` | in letter |
| L1.3 | Alarm only when batch ~ outcome AND batch spans train/test; exempt blocked designs | package | Accept: current fold x batch rule cannot detect confounding and fires in the safe direction (blocked designs already exempt since 2026-11) | New `confounding_alignment` rule + calibration + fixtures | `R/audit.R` | in letter |
| L2.1 | No demonstration of real subject-overlap leakage: compare row-wise vs subject-grouped CV with identical features/preprocessing, ICC(X) > 0, ICC(Y) > 0 | replication + text | Accept: core motivation not shown | Add ablation to the worked example (also answers L2.2) | Sec 4; `replication.R` | in letter |
| L2.2 | Worked example changes two things (leak feature + preprocessing zscore -> none); clarify `leak_subject` is target-proxy leakage | text + replication | Accept | Same preprocessing in both fits; relabel the leak | Secs 4.2-4.3 | in letter |
| L2.3 | Single-seed leaky AUC 0.92 is ~2 SE above the population AUC (~0.87); report multi-seed or note variability | text | Accept | Note sampling variability (or multi-seed mean) | Sec 4.3 (p. 11) | in letter |
| L3.1 | Toy example: 19 pairs vs ~14 expected by chance in 3-D; flag rule and the "[OK]" reconciliation unspecified | package (done) + text | Accept: already changed (centred columns, distance rule, < 5 features: identical rows only; summary verdict uses the p-value) | Rewrite Sec 4.4 text and listing | Sec 4.4 (pp. 12-13) | in letter |
| L3.2 | Case study: similarities on uncentred genes dominated by the baseline profile | package (done) + replication | Accept: already changed (column centring) | Re-run case study | Sec 7.2 | in letter |
| L3.3 | Sec 7.2 interpretation of cross-fold duplicates under study_loocv is self-contradictory | text | Accept | Rewrite after re-run | Sec 7.2 (p. 30) | in letter |
| L3.4 | Report the expected number of random pairs and the enrichment ratio | package | Accept | Add expected chance count + ratio to the duplicate output | `R/audit.R` | in letter |
| L3.5 | Document the rule that sets `duplicate_overlap = TRUE` and how the interpretation text is derived | text + package (docs) | Accept | Help page + one paper sentence | `?audit_leakage`; Sec 4.4 | in letter |

### Changes forced by the 2026-11 package fixes (`dev/paper-updates.md`)

| ID | Point | Type | Proposed decision | Action | Location | Status |
|---|---|---|---|---|---|---|
| P1 | ΔLSI described as BCa CI + sign-flip test, tiers by repeats, A at 20 repeats | text | Required | Rewrite: corrected t-test, R_eff = folds per repeat, renamed tiers (overlaps K4) | Fig 1, Tab 1, Sec 5, Tab 6, Sec 9 | in letter |
| P2 | Fixed-prediction permutation test called "slightly liberal (5-11%)" | text + replication | Required | Replace with regenerated calibration (overlaps K1) | Secs 1, 6.2, 9 | in letter |
| P3 | `perm_stratify = TRUE` said to respect the grouped split | text | Required | Correct description | Sec 4.4 (l. 494) | in letter |
| P4 | Duplicate detection described as cosine on the expression matrix | text | Required | Describe current rule (overlaps L3.1, L3.5) | Sec 7.2 (ll. 1469-71) | in letter |
| P5 | Simulation, split-mode, target-scan and ΔLSI numbers no longer reproduce | replication | Required | Regenerate (Phase 3) | Tables 2-4, Figs 2-4, Secs 6-7 | in letter |
| P6 | ΔLSI null calibration (Sec 5.3) tests an inline sign-flip, not `delta_lsi()` | replication + text | Required | Rewrite Section 5.3 of `replication.R` to call `delta_lsi()` | `replication.R` Sec 5.3; paper Sec 5 | in letter |
| P7 | `replication.R` comments and arguments stale (repeats = 20 -> tier A; "bit-faithful"; BCa legend; `M_boot`/`M_flip`/`seed`; `perm_stratify` with `perm_refit = FALSE`) | replication | Required | Edit script | `replication.R` | in letter |
| P8 | Listing outputs changed (split hash, AUCs, audit summary, tune, ΔLSI) | replication + text | Required | Re-run listings, update (overlaps J7, K2) | Secs 4-5 | in letter |
| P9 | Package version in the manuscript | text | Required | Set to the frozen version | l. 1902 | in letter (0.4.0, set 2026-11-24) |
| P10 | Case-study duplicate counts (21,435 / 33,885) and "extensive near-duplicate overlap" came from uncentred cosine | replication + text | Required | Re-run; rewrite the conclusion from the new counts (overlaps L3.2, L3.3, P4) | Sec 7.2 (ll. 1467-87, 1529, 1534) | in letter |
| P11 | Text built on the old calibration: "mild liberality" paragraph, "3-13%" range, multivariate null "1 of 10", hard-coded Table 4 none column | text | Required | Rewrite from regenerated values (overlaps P2, P5) | Sec 6 (ll. 1060-83, 1312-15); Table 4 | in letter |
| P12 | Case study: reference matrix (`X_ref`) and the leak PC built with `grep("^[Xx]", ...)`, which keeps only 3 of 2,000 genes (XRCC4, XPO1, XPA); target scan, duplicate detection, multivariate scan and `leak_global_pc1` ran on 3 genes (also in the submitted version) | replication + text | Found 2026-11-24 in Phase 4, reproduced (3 columns matched); required | Use all gene columns in Sections 4 and 5.1; store duplicate chance counts; re-run Sections 4, 5.1, 6, 7; rewrite Sec 7 target-scan and duplicate text. Done 2026-11-24: Sections 1, 4, 5.1, 6, 7 re-run (54 min); AUCs and Delta-LSI unchanged (the lasso never selected the PC feature); guarded multivariate scan p 0.8713 -> 0.0198; leaky duplicate pairs 1,968 (all within study, 0 expected; none with genes only), guarded 0; target-scan top genes NUAK1, BCHE, HSD17B6 | `replication.R` Sec 4, 5.1; paper Sec 7 | in letter |
| P13 | paper-claims-check round 1 (2026-11-24): duplicate explanation false (pairs persist with genes alone; 1,965 of 1,968 in GSE51088); case-study Delta-LSI holds the split fixed (FS + leak features only); LeakTune is S3; five modes incl. combined; regression built-in metric is rmse only; 6.3 refit rationale wrong; default winsorization unmentioned; script comments/prints stale | text + replication | Found by the independent claims check, reproduced; required | Text fixed; script comments, Section 4 print and Table 5 reprint fixed (Section 6 re-run, no number changed) | Secs 3, 4.2, 6.1, 6.3, 7, 8.1; `replication.R` | in letter |

## Phase 1 evidence (package items)

| ID | Change | Evidence |
|---|---|---|
| J10 | SummarizedExperiment to Suggests; SE input without it stops with an install hint (`.bio_is_se()`) | test-utils-se.R; installed copy runs a data.frame workflow without loading SummarizedExperiment |
| J8 | `make_split_plan(stratify_method = c("groups", "rows"))`; "rows" balances class row counts (StratifiedGroupKFold approach); used by refit redraws, nested plans, tune inner splits | test-stratify-rows.R; 200 seeds, 40 groups, 5 folds: mean largest fold deviation of case share 0.022 (rows) vs 0.134 (groups) vs 0.208 (unstratified); 3000 groups in 0.6 s |
| L1.1-L1.3 | `batch_assoc` gains batch-outcome test and `span_frac`; `confounding_alignment` needs both; fold-batch test kept as description | test-batch-confounding.R, test-mechanism-summary-multiplicity.R; probe batch_confounding.R (B = 200): null flag 0.065 / 0.075 / 0.060, clustered (40 subjects x 5 rows) 0.035 (limit 0.096), batch_blocked 0.000, power 0.98 / 0.99 / 1.00 (clustered). Review round 1 found the test counted rows, not subjects (null flag 0.68 with 5 rows per subject); fixed by testing one row per unit. Review round 2 found rsample grouped plans still counted rows (null 0.735) and that collapsing hid within-subject confounding (paired design not tested); fixed: unit from the resolved design, unit-preserving permutation p-value when batch varies within units. After fixes (B = 200): null 0.035-0.075 incl. within-subject 0.070, rsample 0.035 (unequal sizes 0.065); paired power 1.00 **Review round 3 (2026-11-20) found three more L1.2 defects; stopped at the 3-cycle limit (CLAUDE.md 7.10); user approved a fourth cycle ("Fix L1.2 as proposed"):** (a) study_loocv / batch_blocked plans with a recorded `group` test rows (null flag 0.70, reproduced); (b) time_series plans without a group test autocorrelated rows (null 0.36-1.00 for AR 0.3-0.9; reviewer's result, not re-run); (c) unit permutation with unequal subject sizes loses power silently (p = 1, coverage 0, V = 0.80 confound not flagged, reproduced). Cycle 4 fix: unit = recorded group in any mode; time series without group not tested; batch varying within units: unit-level test + within-unit batch permutation, Bonferroni. Probe (B = 200): 9 null designs 0.000-0.075 (limit 0.096), 5 power designs 0.98-1.00; reviewer repros: loocv 0.02, batch_blocked 0.02, ts 0, unequal-size power 1.00. **Review round 4 (2026-11-20) found three more; stopped again, awaiting user decision:** (1) within-subject permutation ignores time order: visits processed in runs by date, AR outcome, null flag 0.175-0.245 (reviewer's result, not re-run); (2) combined plan with the subject as a constraint (not primary axis) tests rows: null flag 0.58 over 60 datasets (reproduced); (3) `print.LeakAuditList` (audit_leakage_by_learner) still shows only fold x batch p and V (reproduced). User chose the contract approach (CLAUDE.md 7.10): unit-level test only, within-unit association reported (`within_unit_effect`), subject constraint axes as units, printout fixed. Final review against the contract: 2 inside findings fixed (alphabetical tie-break biased the per-unit summary: null 0.19-0.25, pairs + singletons 0.78 -> random tie-break; survival printed a wrong reason), 1 doc fix (L3.5 proxy rule; code behaviour logged as BL-25). Probe (B = 200): 14 null designs 0.000-0.080 (limit 0.096), blocked 0.000, power 0.98-1.00, paired confound reported in 50/50, never flagged. Suite 0 failures |
| L3.4-L3.5 | `info$duplicate_null` (observed, expected, enrichment; Gaussian copula null, `duplicate_null_reps`); flag rules documented; summary lists flagged mechanisms | test-duplicate-null.R; 40 clean datasets: observed vs expected 5.1 vs 6.4 (Gaussian), 39.4 vs 44.1 (skewed), 7.1 vs 5.6 (discrete); +32 s at 2015 x 2000 with 3 draws |
| K1 | message with B x folds model fits when "auto" refits; message when it falls back for B > perm_refit_auto_max; "slightly liberal" removed | test-perm-refit-messages.R |

Found in passing, not a reviewer item: BL-24 in `dev/OPEN-ISSUES.md` (fewer groups than folds gives fewer folds without a warning).

## Freeze (Phase 2)

User will bump the version and commit at the end (2026-11-23), so numbers
are regenerated from a checksummed build instead of a commit:
`bioLeak_0.3.8.tar.gz` (DESCRIPTION still 0.3.8; the final commit must
match this build), SHA-256
`6e919f23f49d78120b9b3cee7a268f8d65d0faf7e3cf5fea7379fe7edf891753`, built
2026-11-23 from the working tree, installed into the scratch library
`freeze/lib`. Any package change after this date invalidates the freeze.

Version bump (2026-11-24, user decision 0.4.0): `bioLeak_0.4.0.tar.gz`,
SHA-256 `dfdeffe49615a716f3e7b45111e0d51e6f7f417b46ca2bade01701be082e4736`,
differs from the freeze build only in DESCRIPTION (version, packaging date),
the NEWS heading and README.md (rebuilt 2026-11-24 after a README update);
R code, man pages, tests and vignettes are identical,
so the regenerated outputs apply to 0.4.0 unchanged.

## Replication edits (Phase 3, `paper/submission/replication/`)

| ID | Edit |
|---|---|
| L2.2 | Section 1, 4.3: leaky fit uses the guarded preprocessing (median, z-score); comment names `leak_subject` as target-proxy leakage |
| L2.1, L2.3 | New Section 1b: subject-overlap ablation (row-wise vs subject-grouped CV, ICC(X) 0.7 vs 0, 20 datasets, random forest, identical features and preprocessing) and the worked example over 50 seeds; `sim_results/revision_ablation.rds`; placeholders `[ABL_*]`, `[WE_*]` |
| P6 | Section 5.3 calls `delta_lsi()` on two equally good pipelines (100 datasets, 5 x 5-fold); inline sign-flip removed |
| P7 (run) | `perm_stratify` removed from the four fixed-prediction case-study audits (ignored; re-run of Sections 4-7 gave identical values); Section 3.1 requests `stratify` only for non-time-series modes (ignored for time series; not re-run); header runtime corrected to ~3.5 h |
| P7 | `M_boot`, `M_flip`, `seed` removed from `delta_lsi()` calls; tier / BCa / sign-flip / "bit-faithful" comments corrected; Fig. 4(b) legend "95% CI (corrected t)"; case-study field `r_eff`; `ranger` added to `install_deps.R` |

## L1.2 contract (batch-outcome test and the confounding_alignment flag)

Written 2026-11-20 under CLAUDE.md 7.10; reviews check this list.

Unit of the test: the plan's own independent unit. In order: the recorded
`group`; for combined plans, a constraint axis of type `"subject"`, else the
primary axis; for rsample plans, the group of the inferred design. Each unit
counts once (its most frequent batch; its majority class, or mean outcome;
ties broken at random, separately for batch and class).
Without a unit, rows are the units.

Covered designs (each with a calibration case in
`.claude/skills/verify-change/probes/batch_confounding.R`; null flag rate
<= 0.096 at B = 200, power >= 0.8 under a clear association):

| Design | Null case | Power case |
|---|---|---|
| one row per subject (binary, regression, several batch columns, repeats) | NULL-BIN, NULL-GAUSS, NULL-MULTI | POWER-BIN, POWER-GAUSS |
| repeated rows per subject, subject-level batch | NULL-CLUSTER | POWER-CLUSTER |
| batch varying within subjects (equal and unequal sizes, 1-2 samples per subject, time-ordered runs, pairs mixed with singletons) | NULL-WITHIN, NULL-UNEQUAL, NULL-TIES, NULL-LONGITUDINAL, NULL-PAIRED-MIX | (between-subject part: POWER-CLUSTER) |
| multiclass outcome (test level) | NULL-MULTICLASS | - |
| rsample grouped folds | NULL-RSAMPLE | - |
| leave-study-out with a subject column | NULL-LOOCV-GROUP | - |
| combined plan with the subject as a constraint axis | NULL-COMBINED-SUBJECT | - |
| blocked designs (batch_blocked, study_loocv on the batch) | BLOCKED (never flagged) | - |

Outside the contract, documented, not tested:
- confounding within subjects (batch varying within a subject together with
  the outcome, e.g. paired tumour/normal): the row-level association is
  reported (`within_unit_effect`) with a note; it does not set the flag;
- time-series plans without a group column: not tested (`outcome_method =
  "not_tested"`), rows are autocorrelated;
- survival outcomes: not tested (`outcome_method = "not_tested_survival"`);
- an explicit `perm_mode` passed to `audit_leakage()` sets the design used to
  choose the unit; no calibration case covers overriding the plan's design.

## Changed values

Regenerated 2026-11-23/24 from the frozen build (SHA-256 6e919f23...),
`replication.R` of the working copy; 3 h 14 min. Old = submitted
`bioLeak_replication_v3/placeholder_values.rds`; unchanged values are not
listed. Calibration checked, not only copied: fixed-prediction type-I 6.1%
over 180 datasets (limit 9.9%; by n 5.0 / 8.3 / 5.0%, limit 13.4%);
`delta_lsi()` null 0/100 (limit 11.5%). Not as `dev/paper-updates.md`
predicted from a 10-seed probe: BB_BC and SL_SC stay 0.10. All six listing
outputs change (split hash, AUCs, audit summary, tuning, Delta-LSI); new
text in the scratchpad file `regen/listing_outputs_section1.txt`. Run
warnings: package build versions (benign); `perm_stratify` ignored in
Section 4 and `stratify` ignored for time-series folds in Section 3.1, both
removed from the script; a re-run of Sections 4-7 reproduced every value.
2026-11-24: worked example changed to one batch per subject (user
decision); Sections 1, 1b and 6 re-run: only `[WE_GUARDED_MEAN]` (0.511 ->
0.506), `[WE_GUARDED_SD]` (0.081 -> 0.065) and `[WE_LEAKY_MEAN]` (0.854 ->
0.855) moved; listings regenerated (Section 4.2 guarded AUC 0.600; 4.4
audit batch-outcome p = 0.693, not flagged). Note: `batch` is a numeric
predictor in the worked example (as submitted).

| Placeholder / listing | Old | New | Item ID |
|---|---|---|---|
| `[PLATFORM]` | macOS (4 parallel workers) | macOS (3 parallel workers) | P5 |
| `[TYPE_I_RATE]` | 10.0 | 6.1 | P2, P11 |
| `[TYPE_I_LO]` | 8.3 | 5.0 | P2, P11 |
| `[TYPE_I_HI]` | 11.7 | 8.3 | P2, P11 |
| `[TYPE_I_N100]` | 10.0 | 5.0 | P2, P11 |
| `[TYPE_I_N250]` | 11.7 | 8.3 | P2, P11 |
| `[TYPE_I_N500]` | 8.3 | 5.0 | P2, P11 |
| `[NONE_OVERALL]` | 80.8 | 80.3 | P5 |
| `[NONE_100]` | 0.61 | 0.56 | P5 |
| `[NONE_250]` | 0.82 | 0.85 | P5 |
| `[NONE_500]` | 0.99 | 1.00 | P5 |
| `[SO_100]` | 0.90 | 0.86 | P5 |
| `[BC_100]` | 0.97 | 0.96 | P5 |
| `[TL_100]` | 0.61 | 0.68 | P5 |
| `[TL_250]` | 0.96 | 0.97 | P5 |
| `[S0_SO_100]` | 0.85 | 0.68 | P5 |
| `[S0_TL_100]` | 0.13 | 0.07 | P5 |
| `[S0_TL_250]` | 0.10 | 0.08 | P5 |
| `[S0_TL_500]` | 0.03 | 0.05 | P5 |
| `[CLEAN_AUC]` | 0.74 | 0.72 | P5 |
| `[DELTA_SO]` | 0.089 | 0.102 | P5 |
| `[DELTA_PP]` | 0.253 | 0.263 | P5 |
| `[DELTA_BC]` | 0.096 | 0.107 | P5 |
| `[DELTA_TL]` | 0.024 | 0.031 | P5 |
| `[MULTI_NONE]` | 10 | 0 | P11 |
| `[GUARDED_AUC]` | 0.610 | 0.611 | P5 |
| `[G_SD]` | 0.063 | 0.081 | P5 |
| `[LEAKY_AUC]` | 0.810 | 0.808 | P5 |
| `[L_SD]` | 0.024 | 0.018 | P5 |
| `[DELTA_AUC]` | 0.200 | 0.197 | P5 |
| `[G_GAP]` | 0.101 | 0.111 | P5 |
| `[L_GAP]` | 0.302 | 0.309 | P5 |
| `[L_CRAMER]` | 0.079 | 0.078 | L1.1 (fold-batch, descriptive) |
| `[L_BATCH_P]` | 0.5302 | 0.5762 | L1.1 (fold-batch, descriptive) |
| `[MULTI_P_GUARDED]` | 0.2673 | 0.0198 | P5 |
| `[N_DUPS]` | 33885 | 1968 | P10, L3.2 |
| `[DLSI_METRIC]` | 0.180 | 0.176 | P1, P5 |
| `[DLSI_ROBUST]` | 0.181 | 0.177 | P1, P5 |
| `[DLSI_PVAL]` | < 0.001 | = 0.0016 | P1, P5 |
| `[DLSI_TIER]` | A_full_inference | C_ttest_3_9_folds | P1, P5 |
| `[DLSI_REFF]` | 20 | 5 | P1, P5 |
| `[DLSI_CI_LO]` | 0.174 | 0.112 | P1, P5 |
| `[DLSI_CI_HI]` | 0.187 | 0.241 | P1, P5 |
| `[DLSI_SIM_ALT_MEAN]` | 0.1935 | 0.1941 | P6 |
| `[DLSI_SIM_NULL_REJECT]` | 10.0 | 0.0 | P6 |
| `[DLSI_SIM_NULL_MEAN]` | 0.0006 | -0.0001 | P6 |
| `[ABL_ICC_ROW_MEAN]` | - | 0.800 | L2.1 |
| `[ABL_ICC_ROW_SD]` | - | 0.037 | L2.1 |
| `[ABL_ICC_GRP_MEAN]` | - | 0.499 | L2.1 |
| `[ABL_ICC_GRP_SD]` | - | 0.098 | L2.1 |
| `[ABL_NOICC_ROW_MEAN]` | - | 0.493 | L2.1 |
| `[ABL_NOICC_ROW_SD]` | - | 0.049 | L2.1 |
| `[ABL_NOICC_GRP_MEAN]` | - | 0.491 | L2.1 |
| `[ABL_NOICC_GRP_SD]` | - | 0.068 | L2.1 |
| `[ABL_N_DATASETS]` | - | 20 | L2.1 |
| `[WE_GUARDED_MEAN]` | - | 0.506 | L2.3 |
| `[WE_GUARDED_SD]` | - | 0.065 | L2.3 |
| `[WE_LEAKY_MEAN]` | - | 0.855 | L2.3 |
| `[WE_LEAKY_SD]` | - | 0.039 | L2.3 |
| `[WE_N_SEEDS]` | - | 50 | L2.3 |
| `[BO_V]` | - | 0.198 | L1.1-L1.3 |
| `[BO_P]` | - | <1e-10 | L1.1-L1.3 |
| `[G_SPAN]` | - | 0.00 | L1.1-L1.3 |
| `[L_SPAN]` | - | 1.00 | L1.1-L1.3 |
| `[G_CONF]` | - | FALSE | L1.1-L1.3 |
| `[L_CONF]` | - | TRUE | L1.1-L1.3 |
| `[DUPS_G_OBS]` | - | 0 | P5 |
| `[DUPS_G_EXP]` | - | 0.0 | P5 |
| `[DUPS_L_OBS]` | - | 1968 | P5 |
| `[DUPS_L_EXP]` | - | 0.0 | P5 |

## Decisions log

| Date | Item | Decision | By |
|---|---|---|---|
| 2026-11-24 | Worked example | One batch per subject (was round-robin over rows), so the batch-outcome test is meaningful in the Section 4.4 listing | User |
| 2026-11-20 | L1.2 | After review round 4: contract approach (CLAUDE.md 7.10 rule added); unit-level test only, within-unit association reported, not flagged | User |
| 2026-11-20 | L1.2 | Fourth fix cycle approved past the 3-cycle limit, as proposed (unit from recorded group; time series without group not tested; within-unit batch permutation + unit-level test) | User |
| 2026-11-19 | All items J1-J10, K1-K4, L1.1-L3.5, P1-P11 | Accepted as proposed. J8: add row-balancing option + document; J10: SummarizedExperiment to Suggests (approved); L1.2/L1.3: new batch x outcome rule; L2.1: ablation in the worked example; K3: partial (figure + compact table); J2: discussion only | User |
