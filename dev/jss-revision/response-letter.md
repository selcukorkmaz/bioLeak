# Response to the reviewers

*bioLeak: Leakage-Aware Modeling and Diagnostics for Machine Learning in R*

We thank the editor and the three reviewers for their careful reading. The comments led to changes in the package as well as in the paper. Page numbers refer to the revised PDF.

## Summary of main changes

- **Batch diagnostic (Reviewer L).** The confounding alarm now tests batch against the outcome, counting each subject once. It fires only when batch is associated with the outcome and batches are shared by training and test rows. Blocked designs are never flagged.
- **Subject overlap (Reviewer L).** New Section 4.4 compares row-wise and subject-grouped CV with identical features, preprocessing and learner.
- **Delta-LSI inference.** It now uses the corrected resampled *t*-test of Nadeau and Bengio (2003), replacing the BCa interval and sign-flip test, because fold- and repeat-level differences from repeated CV are correlated. A null calibration of the shipped function is added.
- **Package changes:**
  - `make_split_plan(stratify_method = "rows")` balances class row counts across folds (Reviewer J).
  - SummarizedExperiment is now optional (Reviewer J).
  - `perm_refit = "auto"` reports its cost (Reviewer K).
  - The duplicate audit reports a chance count (Reviewer L).
- **Regeneration.** All numbers, tables, figures and listings were regenerated from the revised package. Section 6 of the replication script reprints every table.
- **Length.** The paper is shorter (38 pages, from 41): the repeated workflow code is condensed and two simulation tables are merged.

---

## Reviewer J

**J1. Mention mlr3resampling alongside mlr3.**
Added to the introduction (p. 2), with the reference.

**J2. Could bioLeak interface with mlr3, or is there a barrier?**
There is no technical barrier. A `LeakSplits` plan stores explicit training and test indices, which mlr3's `ResamplingCustom` accepts. Section 8.4 (p. 32) now says so. Such a bridge is not part of this version: the audits and `delta_lsi()` would also need mlr3 predictions converted to a `LeakFit`. bioLeak already exchanges splits with rsample.

**J3. Align decimal points and parentheses in the tables.**
Done for the result tables (Tables 3–5).

**J4. Note in Section 2.2 that StratifiedGroupKFold and mlr3resampling support grouped fold assignment.**
Added to Section 2.2 (p. 4). The sentence names `GroupKFold`/`StratifiedGroupKFold`, `rsample::group_vfold_cv()`, `caret::groupKFold()`, mlr3resampling and blockCV.

**J5. Software references for Sections 2.1 and 2.3–2.5, ideally a table.**
Each subsection now names the relevant software (pp. 4–5). A new Table 1 (p. 5) shows, for eight packages, which of the five leakage mechanisms each prevents or detects with built-in functions. It replaces the earlier feature table in Section 8.4.

**J6. Text going off the right margin.**
Fixed. The revised PDF has no line running into the margin.

**J7. Too much code output; more discussion between chunks.**
- The repeated per-fold tables are omitted from the fit listings, and this is marked.
- The chained-workflow listing (Section 4.6) is reduced to the new steps.
- Each remaining output is followed by its interpretation.

The package output itself was left unchanged.

**J8. Grouped stratification balances groups, not rows; explain the algorithm.**
We agree. The default assignment labels each group by its majority class and deals the groups of each class to folds in turn. This balances numbers of groups, not rows. The revised package adds:

- `stratify_method = "rows"`, which places groups, largest first, in the fold that keeps each class's share of rows most even. This is the criterion of `StratifiedGroupKFold`.
- A description of both methods on `?make_split_plan` and in Section 4.1 (p. 9).

In simulations with unequal group sizes (40 groups, 5 folds, 200 datasets), the largest deviation of a fold's case share from the overall share averaged 0.022 with `"rows"`, against 0.134 with `"groups"` and 0.208 without stratification. These figures are in the package NEWS, not the paper.

**J9. Which goal (rows or groups) is better when?**
Discussed in Section 4.1 (p. 9):

- Row balance matters when performance is computed over rows, as for AUC pooled across samples.
- Group balance matters when each group is one unit of evaluation.
- The two coincide when groups have similar sizes and a constant outcome.

**J10. `install.packages("bioLeak")` fails without SummarizedExperiment.**
SummarizedExperiment moved from Imports to Suggests, so the revised version installs without Bioconductor. SummarizedExperiment input without the package now stops with an installation hint. Section 8.3 (p. 30) lists the dependencies.

## Reviewer K

**K1. `perm_refit = "auto"` refits for every permutation; please warn the user.**
- **New message:** `audit_leakage()` now reports the number of model fits (B × folds) before it starts refitting. It also says when "auto" falls back to fixed predictions, which was previously silent. The listing in Section 4.5 shows the message.
- **Paper text:** Section 4.5 (p. 14) explains the cost. Section 6.2 (p. 22) explains why the simulations use fixed predictions.
- **The 10% false-positive rate:** in the regenerated simulations, the type-I error of the fixed-prediction test is 6.1% over 180 null datasets (5.0–8.3% by sample size). This is within the Monte Carlo limit for a nominal 5% test.

**K2. Section 4.5 repeats Sections 4.1–4.4.**
Agreed. Section 4.6 (p. 14) now reuses the objects defined earlier and shows only the steps it adds: the repeated split plan, the audit, `delta_lsi()` and the report.

**K3. Standardise on figures.**
Partly agreed:

- **Main simulation:** shown in Figure 2. Its exact values are in one two-panel table (Table 3), which replaces the former Tables 2 and 3.
- **Split-mode results:** kept as a table (Table 4). Its 20 cells rest on 10 runs each, and the text discusses exact cell values. A figure would not show these better.

**K4. Features presented that are not available (survival tuning; Delta-LSI `by_group` / `within_batch`).**
- **Survival:** the limitations are stated where survival is mentioned (Sections 4.2 and 4.7) and in Table 6.
- **Delta-LSI:** with the corrected *t*-test, all values of `exchangeability` use the same test, and the argument is only recorded (Section 5, p. 18; Section 9.2). The paper no longer presents `by_group` or `within_batch` as distinct procedures.

## Reviewer L

**L1.1. The batch diagnostic is described inconsistently.**
The paper and the help page now describe two tests:

- **Batch–outcome test:** whether batch carries information about the outcome.
- **Fold–batch test:** only how batches are spread over folds.

The audit summary prints both, together with the share of test rows whose batch also occurs in training (Section 4.5, p. 13).

**L1.2. Add a direct batch × outcome test.**
Added. Rows of one subject are not counted as independent:

- **Units:** when the plan records a subject or group column, the test uses one observation per unit. This includes leave-study-out, combined and rsample grouped plans.
- **Test statistic:** chi-square with Cramér's V for classification; Kruskal–Wallis with ε² for regression.
- **Not tested:** time-series plans without a group column, whose rows are autocorrelated.
- **Within subjects:** when batch varies within subjects, the within-subject association is reported as an effect size but not tested.

The limits are stated in Section 9.2 (p. 33) and on the help page.

**L1.3. Alarm only when batch is associated with the outcome and spans training and test; exempt blocked designs.**
Implemented as suggested:

- `confounding_alignment` is flagged only when the batch–outcome test is significant with at least a small effect and some test row shares its batch with training.
- Blocked designs share no batch between training and test, so they are never flagged.
- **Calibration** (200 datasets per design, nominal 5%): the null flag rate was 0.000–0.080 across 14 designs, against a limit of 0.096. These designs include repeated rows per subject, batch varying within subjects, rsample folds, leave-study-out and combined plans, and time series. The flag rate was 0.98–1.00 under a clear association.
- **Case study:** study is associated with the outcome (V = 0.198). Leave-study-out is not flagged; row-wise CV is (Table 5, p. 28).

**L2.1. No demonstration of subject-overlap leakage.**
New Section 4.4 (p. 12) uses the design you proposed, with identical features, preprocessing and learner (a random forest). It has 40 subjects × 5 rows, an outcome constant within subjects and unrelated to the features, and features with ICC(X) = 0.7 or 0.

Mean AUC over 20 datasets:

| Features | Row-wise CV | Subject-grouped CV |
|---|---|---|
| ICC(X) = 0.7 | 0.800 (SD 0.037) | 0.499 (SD 0.098) |
| ICC(X) = 0 | 0.493 | 0.491 |

**L2.2. The worked example changes two things; `leak_subject` is target-proxy leakage.**
Both are corrected in Section 4.3 (p. 11):

- Both fits now use the same splits, preprocessing and learner.
- `leak_subject` is described as target-proxy leakage (Section 2.5), not subject overlap.

**L2.3. The single-seed leaky AUC of 0.92 is optimistic; report variability.**
Over 50 datasets from the same design, the leaky AUC averaged 0.855 (SD 0.039) and the guarded AUC 0.506 (SD 0.065). Sections 4.2–4.3 report these alongside the single-seed values.

**L3.1. In the 3-D toy example, 19 pairs against about 14 expected by chance; flag and "[OK]" rules unspecified.**
We agree that these were chance pairs.

- **Duplicate rule:** it now centres columns and requires similar row lengths and a small distance. With fewer than five features it reports only identical rows.
- **Worked example:** it now reports no pairs, against 0.0 expected by chance (Section 4.5).
- **Interpretation block:** it restates the permutation result and then lists any flagged mechanism. A flagged mechanism can no longer appear next to an unqualified "[OK]".

**L3.2. Case-study similarities on uncentred genes.**
- **Centring:** columns are now centred before similarities are computed.
- **A defect we found while re-running the case study:** the replication script selected the reference columns by an "X" prefix. This kept only 3 of the 2,000 genes (XRCC4, XPO1, XPA) for the target scan, the duplicate screen and the multivariate scan. The leak feature built from the first principal component was affected the same way. The counts of 21,435 and 33,885 pairs in the submitted version therefore came from 3 or 6 columns.
- **After the fix:** all audits use the 2,000 genes, and Section 7.2 (pp. 27–29) reports the new results. The case-study AUCs and Delta-LSI did not change.

**L3.3. The interpretation of cross-fold duplicates under study-LOOCV is contradictory.**
Rewritten in Section 7.2 (p. 28).

- **Guarded audit:** no near-duplicate pairs across training and test.
- **Leaky audit (row-wise CV):** 1,968 pairs, all within single studies (1,965 within GSE51088), with 0 expected by chance. Most are also found with the genes alone.
- **Why they differ:** study holdout keeps within-study pairs on one side of each fold, so they cannot inflate the guarded estimate. Under row-wise CV they cross folds.
- **Caveat:** the chance count does not model clustering by study, so the paper does not present these pairs as evidence of duplicated samples.

**L3.4. Report the expected number of chance pairs and the enrichment ratio.**
`audit_leakage()` now reports the observed count, the expected count and their ratio (`info$duplicate_null`); `summary()` prints them.

- **Expected count:** the same rule applied to null rows drawn from a Gaussian copula with the features' marginal distributions and rank correlation.
- **Calibration:** on data without duplicates (40 datasets each, Gaussian, skewed and discrete features), the mean expected count was within about 25% of the mean observed count.
- **Limitation:** clustering beyond rank correlation is not reproduced; this is stated on the help page and in Section 9.2.

**L3.5. Document the flag and interpretation rules.**
Each mechanism-summary rule is now stated on `?audit_leakage`, together with how the interpretation text is derived. Section 4.5 (p. 13) summarizes them.

---

## Other changes

- **Delta-LSI method (Section 5, pp. 17–20).** The submitted version treated repeats as independent and reported tier A at 20 repeats. The corrected *t*-test uses the number of folds per repeat as its effective size, so the case study is now tier C:

  | | Submitted | Revised |
  |---|---|---|
  | Delta-metric | 0.180 | 0.176 |
  | 95% interval | [0.174, 0.187] | [0.112, 0.241] |
  | *p* | < 0.001 | 0.0016 |

  The null calibration now calls `delta_lsi()` itself, on two equally good pipelines: 0 of 100 rejected (limit 11.5%).
- **Simulation numbers.** All simulation values changed with the package corrections listed in NEWS, including fold assignment and the permutation nulls. Text built on the earlier calibration ("mild liberality", "1 of 10") was rewritten.
- **Case-study Delta-LSI scope (Section 7.2).** The text now states that this comparison holds the split design fixed. It measures feature selection and the leak features (0.176, close to 0.074 + 0.091).
- **Smaller corrections:**
  - `LeakTune` is an S3 class.
  - A fifth split mode, `combined`, is mentioned.
  - The built-in regression metric is `rmse`.
  - The default winsorization in guarded preprocessing is stated.
- **Replication runtime.** The regenerated run took 3 h 14 min on a four-core 2017 laptop with three workers (Section 11).
- **Version.** The revised package is bioLeak 0.4.0. Every result was regenerated from its code; Computational details give the version.
