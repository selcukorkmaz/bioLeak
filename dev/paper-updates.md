# JSS manuscript: changes needed after the 2026-11 fixes

Target: `paper/submission/manuscript/manuscript_jss.tex` (latest version) and
its replication script `paper/submission/bioLeak_replication_v3/replication.R`.
Line numbers refer to those files. Source: final scoped verification,
2026-11-18, against the development code (see NEWS.md, development section).

## Before regenerating anything

- Bump the package version (DESCRIPTION is still 0.3.8, the CRAN version, but
  the development code gives different results). The manuscript names the
  version at L1902.
- `replication.R` calls `library(bioLeak)`: install the development version
  first. Every bioLeak call in it still exists and Section 1 runs.
- Section 4 (case study) needs `curatedOvarianData`.

## (A) Method descriptions that are now false

| Lines | Paper says | Now |
|---|---|---|
| 767-856 | ΔLSI inference unit is the repeat; BCa CI; sign-flip test (exact 2^R or Phipson-Smyth); `blocked_time` block sign-flip; tiers by number of repeats (5/10/20) | Nadeau-Bengio corrected t-test on fold deltas over distinct test sets, df = folds per repeat - 1, CI = estimate ± t·se; `R_eff` = folds per repeat; tiers D < 3, C 3-9, B 10-19, A >= 20 folds per repeat; `exchangeability` only recorded. Huber estimate and pairing rule unchanged |
| 832-833, 1826-1829 | `by_group` / `within_batch` fall back to iid with a warning | no warning; every setting gets the same t-test |
| 205 (Fig 1), 276-277, 308 (Tab 1), 593-596, 666, 1646-1649, 1713 (Tab 6), 1849-1851 | BCa / sign-flip; "repeats = 20 gives A_full_inference" | t-test; with v = 5 tier A is unreachable at any number of repeats |
| 575-577, 1064-1072, 1180-1182 | fixed-prediction permutation test is "slightly liberal (5-11%)" because the model captures noise | the cause was the null ignoring samples/subjects; labels are now permuted per sample and per subject (time blocks for time series); type-I near nominal |
| 494 | `perm_stratify = TRUE` makes the shuffle respect the grouped split | it stratifies by outcome; grouping comes from the split mode (and it is ignored, with a warning, for fixed-prediction nulls) |
| 1469-1471 | duplicates by cosine similarity on the expression matrix | columns centred; pairs also need row lengths within 10% and distance within 10% of the median pairwise distance; with < 5 features only identical rows |

## (B) Numbers, tables, figures to regenerate

| Where | Paper | Current code (evidence) | replication.R |
|---|---|---|---|
| Tables 2-3, Fig 2(a), L1044-1046, 1066-1069, 1101 | clean baseline 10.0 / 11.7 / 8.3% | near nominal: 0.050 (n = 100), 0.043 (n = 250), 0.030 (n = 500), s = 0, p = 10 probe; subject proxy 0.82, batch 0.95 at n = 100 | Section 2 |
| Fig 3, L1312-1315 | clean multivariate false positives 1/10 | 0/10; leaky 40/40 unchanged | Section 3.2 |
| Table 4, subject-grouped row | clean 0.10 | 0/10 in probe | Section 3.1 |
| L922-927 | ΔLSI power mean 0.194, all reject (sign-flip) | 0.197, 20/20 reject, but by the t-test at tier C | Section 5.2 |
| L928-939 | "null calibration of the inference machinery" | the script codes its own sign-flip, not `delta_lsi()`; it no longer validates the package. Calibration of the shipped test with the earlier script designs: 0/120 rejections (limit 11%), conservative | Section 5.3 (rewrite to call `delta_lsi()`) |
| L1502-1513, Fig 4(b), case-study ΔLSI | tier A, R_eff = 20, BCa CI [0.174, 0.187], p < 0.001 | re-scoring the stored fold metrics: tier C, R_eff = 5, df = 4, p ~ 0.004, CI ~ [0.10, 0.26] | Section 5.1 |
| L1467-1487, 1529, 1534 | 21,435 / 33,885 duplicate pairs; "extensive near-duplicate overlap" | counts came from uncentred cosine (the removed false-positive mechanism); conclusion unsupported until re-run | Section 4 |
| Table 5 | naive / leaky AUC, SD, gap, V = 0.079; guarded multivariate p 0.2673 | will shift (stratified fold dealing; restricted scan permutation) | Section 4 |

Calibration of the replication designs with the development code: see the
section below (filled in from the calibration run).

## (C) Listings whose printed output changed

| Lines | Paper | Now |
|---|---|---|
| 355-365 | split hash 1a2865... | 7feb50... (fold sizes unchanged) |
| 419, 431 | guarded AUC 0.584 | 0.519 [0.359, 0.678] |
| 476, 488 | leaky AUC 0.92 | 0.896 |
| 509-560 | audit: "[OK]", 19 duplicate pairs flagged, quoted values | observed 0.484, gap 0.020, p = 0.412; chi-sq 12, p = 0.916; multivariate score 0.029, p = 0.802; no duplicates, `duplicate_overlap` FALSE; interpretation prints "No significant non-random signal" |
| 736-744, 748 | tune roc_auc and penalties | roc_auc 0.469 [0.361, 0.577]; penalties 0.000 / 0.063 / 0.476 (inverted tuning AUC fixed) |
| 874-912 | ΔLSI summary | tier C_ttest_3_9_folds, R_eff 5, delta_metric 0.444, CI [0.263, 0.625], t-test p = 0.0024 (df = 4); guarded arm 0.469 (not 0.561) |

Regenerate listings from the code rather than editing by hand.

## (D) replication.R

- Comments at L468, 590, 1736, 1758-1759 claim repeats = 20 reaches
  A_full_inference / R_eff = 20: false.
- L2012: the inline `sign_flip()` is described as "bit-faithful with
  `delta_lsi()`": it is not.
- L3213: Fig 4 legend hard-codes "95% BCa CI"; the band is now a t interval.
- `M_boot`, `M_flip`, `seed` passed to `delta_lsi()` now raise a warning
  (they have no effect); remove them.
- `perm_stratify = TRUE` with `perm_refit = FALSE` (Sections 2 and 4) warns
  once per task that it has no effect; remove it there.

## Calibration of the replication designs

Re-run of replication.R's designs with the development code (real glmnet
fits; R per cell as in the script unless noted). No p-value is
miscalibrated; the fixed-prediction null is now conservative.

Design differences from the older scripts: Section 2 uses n = 100/250/500,
B = 100, 20 seeds (`subject_overlap` renamed `subject_target_proxy`);
Section 3 is a new design (10 batches, 20 studies, `study_confounded` leak,
`perm_refit = TRUE`, B = 50, 10 seeds); Section 5.2 uses 20 seeds and 5.3 an
inline copy of the old sign-flip.

| Placeholder | Paper -> current code | R |
|---|---|---|
| TYPE_I_RATE / LO / HI | 10.0 / 8.3 / 11.7 -> 1.1 / 0.0 / 3.3 (limit 9.9) | 180 (60 per n) |
| TYPE_I_N100 / N250 / N500 | 10.0 / 11.7 / 8.3 -> 0.0 / 0.0 / 3.3 | 60 each |
| S0_SO_100 | 0.85 -> 0.68 | 60 |
| S0_BC_100 | 0.88 -> 0.85 | 60 |
| S0_TL_100 / S0_TL_250 | 0.13 / 0.10 -> 0.00 / 0.02 | 60 each |
| SO_100 | 0.90 -> 0.94 | 120 |
| BC_100 | 0.97 -> 0.96 | 120 |
| TL_100 / TL_250 | 0.61 / 0.96 -> 0.72 / 0.95 | 120 each |
| NONE_100 / 250 / 500 | 0.61 / 0.82 / 0.99 -> 0.63 / 0.81 / 0.98 | 120 each |
| CLEAN_AUC | 0.74 -> 0.73 | |
| DELTA_SO / BC / PP / TL | .089 / .096 / .253 / .024 -> .104 / .102 / .262 / .031 | |
| BB_BC (Table 4) | 0.10 -> 0.00 | 10 |
| SL_SC (Table 4) | 0.10 -> 0.00 | 10 |
| Table 4 none column (hard-coded) | SG 0.10 -> SG 0.0, BB 0.0, SL 0.0, TS 0.1 | 10 each |
| MULTI_NONE | 10 -> 0 | 10 |
| DLSI_SIM_ALT_MEAN | 0.1935 -> 0.1965 | 20 |
| DLSI_SIM_NULL_MEAN | 0.0006 -> -0.0002 | 10 |
| DLSI_TIER / DLSI_REFF | A_full_inference / 20 -> C_ttest_3_9_folds / 5 | from the code |
| DLSI_PVAL / DLSI_CI_LO / DLSI_CI_HI | change (t-based); stored fold metrics give p ~ 0.004, CI ~ [0.10, 0.26] | case study not run |

Unchanged: all PP_* and S0_PP_*, SO/BC/TL_500, the n = 250/500 SO/BC cells,
NONE_OVERALL, MIN_POWER_S2, 14 of 16 split-mode cells (1.00),
UNI_DETECT / SO / PP, MULTI_DETECT, DLSI_SIM_ALT_REJECT (100).
DLSI_SIM_NULL_REJECT (10.0) still reproduces, but it measures the script's
inline sign-flip, not `delta_lsi()`; the shipped test gave 0/120 under the
null (conservative).

Text consequences: the "mild liberality" paragraph (L1064-1083) and the
"3-13%" range (L1060) no longer hold; the multivariate-scan null is "0 of
10", not "1 of 10".

Not run (stopped to finish): R >= 100 null calibration of the Section 3
refit audits per split mode (10 seeds each, all <= 0.1) and of the Section
3.2 multivariate scan; the case-study (Section 4) values.

Tier labels were renamed in the package (2026-11-18): `A_full_inference` ->
`A_ttest_ge20_folds`, `B_signflip_ci` -> `B_ttest_10_19_folds`, `C_signflip`
-> `C_ttest_3_9_folds`; update every tier name in the manuscript (Fig 1,
Table 1, Section 5, Table 6, listings).
