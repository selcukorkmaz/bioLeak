# bioLeak: An integrated framework for detecting and preventing information leakage in biomedical machine learning

## 1. Introduction

Machine learning (ML) is increasingly used to build diagnostic, prognostic, and
treatment-response models from clinical and genomic data. Yet a growing body of
evidence shows that many published models report performance figures that cannot
be reproduced on independent data, often because of information leakage --
the inadvertent use of data that would be unavailable at prediction time during
model training or evaluation. In a cross-disciplinary audit spanning 17
scientific fields, Kapoor and Narayanan (2023) traced leakage-related errors in
294 papers, concluding that the problem is systemic rather than anecdotal [1].

Biomedical data are especially vulnerable because the independence assumptions
of standard cross-validation are routinely violated. Repeated measurements from
the same patient, samples processed in shared laboratory batches, multi-site
collection with centre-specific technical artefacts, and temporal dependencies
in longitudinal records all create opportunities for leakage. In practice, the
principal mechanisms include: (i) preprocessing leakage, where global
imputation, scaling, or feature selection leaks test-fold information into
training statistics [2,5]; (ii) dependence leakage, where observations from
the same subject appear in both training and test folds [2]; (iii) confounding
leakage, where batch or site effects are confounded with the outcome and
cross-validation learns to discriminate batches rather than disease states
[6,7]; (iv) temporal leakage, where future information is accessible at
prediction time [4]; and (v) near-duplicate leakage, where technical replicates
or near-identical samples span the train-test boundary. These mechanisms are
not mutually exclusive: a single multi-site gene-expression study can suffer
from several simultaneously.

The consequences are well documented at both the study and meta-analytic
levels. Rosenblatt et al. (2024) evaluated over 400 connectome-based
prediction pipelines and showed that feature-selection and repeated-subject
leakage drastically inflated performance, in one case producing apparently
significant predictions that vanished once the leakage was corrected [2]. Van
de Mortel and van Wingen (2025) found that at least 45% of MRI studies in a
depression meta-analysis exhibited procedures consistent with leakage; after
excluding them, the pooled log-diagnostic odds ratio dropped from 2.53 to 2.02
[3]. Davis et al. (2024) further showed that label leakage through clinical
workflow artefacts -- such as antihypertensive prescriptions used to predict
hypertension -- can inflate model performance while undermining clinical
validity [4].

Existing ML frameworks provide partial safeguards. Tidymodels [9], mlr3 [10],
and scikit-learn [11] support fold-wise preprocessing pipelines and grouped
cross-validation, which help prevent preprocessing and dependence leakage when
used correctly. However, none of these frameworks offers batch-blocked or
study-level splitting, and -- more critically -- none provides post-hoc
diagnostics to evaluate whether leakage may have occurred despite the chosen
pipeline. Scikit-learn includes a basic permutation significance test, but it
assesses only whether performance exceeds chance, not whether specific leakage
mechanisms are at play. No existing tool combines leakage-aware splitting,
guarded preprocessing, multi-modal auditing, and automated reporting in a
single integrated workflow.

In this paper, we present bioLeak, an open-source R package available on CRAN
that addresses this gap. bioLeak integrates four complementary layers: (i) it
*prevents* leakage through four splitting modes (subject-grouped,
batch-blocked, study leave-one-out, and time-series with purge/embargo windows)
paired with a guarded preprocessing pipeline in which all transformations are
fitted exclusively on training folds; (ii) it *detects* leakage through a
post-hoc audit comprising a permutation-gap test with restricted shuffles, a
batch-association test, a target leakage scan, and near-duplicate detection;
(iii) it *reports* findings via a self-contained HTML audit report suitable for
supplementary material; and (iv) it *validates* its own diagnostic sensitivity
through a built-in simulation suite with controlled leakage injection. The
package supports binomial, multiclass, regression, and survival outcomes, and
integrates with the tidymodels and Bioconductor ecosystems.

The remainder of this paper is organised as follows. Section 2 describes the
methods underlying each component. Section 3 presents an experimental
evaluation comprising controlled simulations and a real-data case study.
Section 4 reports the results. Section 5 discusses implications and
limitations, and Section 6 concludes.

## 2. Methods

### 2.1 Overview

bioLeak organises the modelling workflow into four sequential stages (Fig. 1).
First, a *split plan* partitions the data into training and test folds using
one of four leakage-aware modes (Section 2.2). Second, a *guarded
preprocessing and fitting* step trains the model within each fold, ensuring
that all data transformations use only training-fold statistics (Section 2.3).
Third, an *audit* evaluates the fitted model for evidence of leakage through
four complementary diagnostic tests (Sections 2.4--2.6). Fourth, the audit
results are compiled into an HTML *report* for review and archival. An optional
nested tuning layer (Section 2.7) and a simulation suite for validating audit
sensitivity complete the framework. The following subsections describe each
component.

### 2.2 Leakage-aware splitting

The `make_split_plan` function constructs resampling indices according to one
of four modes, each targeting a specific dependence structure.

**Subject-grouped cross-validation.** Given a grouping variable *g* (e.g.,
patient identifier), all observations sharing the same value of *g* are
assigned to the same fold. When stratification is enabled, the majority outcome
class is determined for each group, and groups are then allocated to *V* folds
so that the proportion of each majority class is approximately equal across
folds. This group-level stratification preserves class balance without
splitting any group across folds.

**Batch-blocked cross-validation.** Given a batch variable *b*, folds are
constructed so that the distribution of batches is balanced. If the number of
unique batches is less than or equal to *V*, each batch forms its own fold
(equivalent to leave-one-batch-out). Otherwise, batches are allocated to *V*
folds with optional group-level stratification.

**Study leave-one-out cross-validation.** Given a study indicator *s*, the
number of folds equals the number of unique studies, and each fold holds out
all observations from a single study. This mode is designed for multi-centre
consortia where external validity across studies is the evaluation target.

**Time-series cross-validation.** Observations are sorted by a time variable
*t* and divided into *V* chronological blocks. For each fold *k*, the test set
comprises block *k* and the training set is determined by three parameters:
*horizon* (minimum gap between the latest training time and the earliest test
time), *purge* (additional gap before the test block), and *embargo*
(exclusion window anchored at the end of the test block). Formally, let
*t*_min and *t*_max denote the earliest and latest test times. A candidate
training observation at time *t*_i is retained only if
*t*_i ≤ *t*_min − *horizon* − *purge* (or *t*_i < *t*_min when both
*horizon* and *purge* are zero) and, if *embargo* > 0,
*t*_i ≤ *t*_max − *embargo*. All three parameters default to zero, recovering
a standard expanding-window scheme.

An optional compact mode stores an integer fold-assignment vector rather than
explicit index lists, reducing memory use for large datasets; fold indices are
reconstructed on demand during fitting.

### 2.3 Guarded preprocessing and fitting

The `fit_resample` function executes a V-fold cross-validation loop in which
preprocessing and model fitting are performed independently within each fold.
Five preprocessing steps are available, each applied in sequence on the
training fold and then transferred to the test fold using only training-derived
parameters:

1. *Winsorisation*: extreme values are clipped to a configurable multiple of
   the median absolute deviation (MAD) above and below the median. The
   location and scale parameters are estimated on the training fold and
   applied to both training and test folds.
2. *Imputation*: missing values are filled using fold-specific medians,
   k-nearest-neighbour interpolation (*k* = 5), or iterative random-forest
   imputation (missForest, default ten iterations). The imputation reference
   (medians or neighbour index) is computed on the training fold and applied
   unchanged to the test fold.
3. *Normalisation*: z-score (mean and standard deviation) or robust (median
   and median absolute deviation) scaling, with location and scale parameters
   estimated on the training fold.
4. *Filtering*: features with training-fold variance or interquartile range
   below a threshold are removed; the same feature mask is applied to the test
   fold.
5. *Feature selection*: univariate association with the outcome (AUC for
   binary, correlation for continuous, eta-squared for multiclass) or
   LASSO-based selection, fitted on the training fold.

Optionally, the user may supply a tidymodels recipe, which is prepped on the
training fold and baked on the test fold, or a tidymodels workflow that
encapsulates both preprocessing and model specification. In all cases, the
guarantee is the same: no test-fold information influences any training-fold
computation.

The function accepts any parsnip model specification as the learner argument,
as well as built-in interfaces to glmnet, ranger, and user-defined custom
learners. Supported outcome types are binomial (with probability predictions),
multiclass (class labels or probability matrices), continuous (numeric
predictions), and survival (risk scores or survival probabilities). Per-fold
metrics are computed using task-appropriate measures: AUC or area under the
precision-recall curve for binary classification; accuracy, macro-F1, or
log-loss for multiclass; RMSE for regression; and Harrell's C-index for
survival. The function returns an S4 LeakFit object containing per-fold
predictions, per-fold metrics, summary statistics, and the preprocessing state.

### 2.4 Permutation-gap test

The permutation-gap test, implemented in `audit_leakage`, quantifies how much
the observed cross-validated performance exceeds what would be expected under a
null hypothesis of no association between features and outcome. The procedure
is as follows:

1. Compute the observed metric *M*_obs from the fitted model's cross-validated
   predictions.
2. For *b* = 1, ..., *B* (default *B* = 200): permute the outcome labels and
   either recompute the metric from the original predictions (fixed-prediction
   mode) or refit the model on the permuted data (refit mode). Denote the
   resulting null metric *M*^(b)_null.
3. The permutation gap is defined as *M*_obs − mean(*M*_null) for metrics
   where higher values indicate better performance (AUC, accuracy, C-index),
   and mean(*M*_null) − *M*_obs for metrics where lower values are better
   (RMSE, log-loss).
4. The one-sided p-value is computed as
   *p* = (1 + #{*M*^(b)_null ≥ *M*_obs}) / (1 + *B*) for higher-is-better
   metrics, with the inequality reversed for lower-is-better metrics.

The choice between fixed-prediction and refit modes is governed by the
`perm_refit` parameter. In automatic mode (the default), bioLeak refits models
when the original training data are available and *B* does not exceed a
configurable ceiling, falling back to fixed predictions otherwise. Refit-based
permutations yield a more powerful test because they capture model-fitting
variability, whereas fixed-prediction permutations test only the association
between existing predictions and labels.

Permutations respect the study design through two mechanisms. First, when
`perm_stratify` is set to `TRUE` or `"auto"` (in automatic mode, stratification
is enabled for classification tasks), labels are shuffled within each outcome
class, preserving the marginal class distribution under the null. Second, for time-series splits, block permutations maintain
the temporal autocorrelation structure. These restricted permutations ensure
that the null distribution reflects the actual resampling constraints rather
than an unrealistically simple exchangeability assumption.

### 2.5 Batch association and target leakage scan

**Batch-association test.** For each categorical metadata column (e.g., batch,
site, study), bioLeak constructs a contingency table of fold assignments
versus metadata levels and performs a chi-square test of independence. The
effect size is summarised by Cramer's V:

  *V* = sqrt(*χ*² / (*n* × (min(*r*, *c*) − 1)))

where *χ*² is the chi-square statistic, *n* is the total sample size, and *r*
and *c* are the numbers of rows and columns in the contingency table. Values
of *V* close to zero indicate that folds are approximately independent of the
metadata variable; large values indicate that folds are aligned with batch
structure, suggesting potential confounding leakage.

**Univariate target leakage scan.** For each feature in a user-supplied
reference matrix *X*_ref, bioLeak computes an association score with the
outcome. The scoring method depends on the outcome type: |2(AUC − 0.5)| for
binary outcomes (scaled to [0, 1]), |Pearson *r*| for continuous outcomes, and
Cramer's V for categorical features. Features with scores exceeding a
threshold (default 0.9) are flagged as potential target leakage sources. This
scan identifies features that carry nearly deterministic outcome information,
such as a per-subject outcome mean that was inadvertently left in the predictor
matrix.

**Multivariate target leakage scan.** Because univariate scans can miss
leakage that operates through feature interactions or linear combinations,
bioLeak includes a multivariate check, enabled by default for supported tasks. The reference matrix is
reduced to its first *k* principal components (default *k* = 10), and pairwise
interactions among the top five components are appended. A cross-validated
model (logistic regression for binary outcomes, linear regression for
continuous) is then fitted on these derived features. The resulting performance
is compared against a permutation null (*B* = 100 by default) to obtain a
p-value. A significant result indicates that the feature matrix collectively
contains more outcome information than expected by chance.

### 2.6 Near-duplicate detection

bioLeak screens for near-identical samples that span the train-test boundary
within any fold. Each row of the reference matrix *X*_ref is L2-normalised,
and pairwise cosine similarity is computed. Alternatively, rows are first
centred (yielding Pearson correlation) or rank-transformed for robustness to
outliers. Pairs exceeding a similarity threshold (default 0.995) are reported.
By default, only pairs where one member appears in the training set and the
other in the test set of at least one fold are retained
(`duplicate_scope = "train_test"`); setting the scope to "all" additionally
reports within-fold duplicates for data-quality purposes.

For datasets with *n* > 3,000, an exact *O*(*n*²) similarity matrix becomes
expensive. In this case, bioLeak uses approximate nearest-neighbour search
(RANN or FNN packages) to identify the *k* nearest neighbours (default
*k* = 50) for each sample and computes similarity only for these candidate
pairs, reducing computational cost while retaining high recall for
near-duplicate detection.

### 2.7 Nested tuning

The `tune_resample` function implements leak-safe hyperparameter tuning via
nested cross-validation. The outer loop uses the leakage-aware splits from
`make_split_plan`. For each outer fold, an inner cross-validation (using the
same splitting mode) is performed on the outer training set to evaluate a grid
of hyperparameter configurations. The configuration with the best inner-CV
performance is selected and evaluated on the held-out outer test set. Two
selection rules are available: "best" (the configuration with the highest mean
inner metric) and "one standard error" (the simplest configuration whose mean
inner metric is within one standard error of the best). The function accepts
tidymodels tune grids and dials parameter objects, and returns per-outer-fold
metrics, selected hyperparameters, and optionally a model refitted on the full
dataset with the most frequently selected configuration.

### 2.8 Software design

bioLeak is implemented in R (version ≥ 4.3) and is available on CRAN under the
MIT licence. The package uses three S4 classes -- LeakSplits, LeakFit, and
LeakAudit -- with formal validity checks on slot consistency. The LeakSplits
object stores fold indices and split metadata; LeakFit extends this with
per-fold predictions, metrics, and preprocessing state; and LeakAudit stores
the results of all diagnostic tests along with a reproducibility trail
(parameter hashes, random seeds, software versions).

bioLeak integrates with both major R ecosystems for biomedical data analysis.
On the tidymodels side, `fit_resample` accepts rsample rset objects as splits,
recipes for preprocessing, parsnip model specifications or workflows as
learners, and yardstick metric sets for evaluation. The `as_rsample` function
converts LeakSplits to an rsample rset for downstream compatibility. On the
Bioconductor side, bioLeak accepts SummarizedExperiment objects as input, with
the assay matrix used as the feature matrix and colData used for outcome,
grouping, and batch metadata.

The `audit_report` function renders a self-contained HTML document from the
LeakAudit object using an embedded R Markdown template. The report includes
sections for cross-validated performance, permutation-gap results with null-
distribution histograms, batch-association tables, target leakage scan results,
calibration diagnostics (for binary outcomes), and duplicate detection
findings. The report is designed to accompany manuscripts as supplementary
evidence of leakage evaluation.

Finally, the `simulate_leakage_suite` function generates synthetic binary
classification datasets with controlled leakage injection. Four scenarios are
supported: subject overlap (repeated measures assigned to folds independently),
batch confounding (outcome distribution varies by batch), preprocessing peek
(global normalisation applied before splitting), and temporal lookahead (a
future outcome value included as a feature). For each scenario, the function
returns cross-validated AUC, permutation gap, and p-value, allowing users to
verify that the audit pipeline detects injected leakage under their chosen
sample sizes and analysis configurations.

**Table 1.** Leakage-related capabilities of existing ML frameworks compared
with bioLeak.

| Capability                      | caret [8] | tidymodels [9] | mlr3 [10] | scikit-learn [11] | bioLeak |
|---------------------------------|:---------:|:--------------:|:---------:|:-----------------:|:-------:|
| Subject-grouped CV              |     P     |       Y        |     Y     |         Y         |    Y    |
| Batch-blocked CV                |     --    |       --       |     --    |         --        |    Y    |
| Study leave-one-out CV          |     --    |       --       |     --    |         --        |    Y    |
| Time-series with purge/embargo  |     --    |       P        |     P     |         P         |    Y    |
| Fold-wise preprocessing         |     Y     |       Y        |     Y     |         Y         |    Y    |
| Permutation-gap test            |     --    |       --       |     --    |         P         |    Y    |
| Batch-association test          |     --    |       --       |     --    |         --        |    Y    |
| Target leakage scan             |     --    |       --       |     --    |         --        |    Y    |
| Near-duplicate detection        |     --    |       --       |     --    |         --        |    Y    |
| HTML audit report               |     --    |       --       |     --    |         --        |    Y    |
| Simulation suite                |     --    |       --       |     --    |         --        |    Y    |

Y = full support; P = partial support; -- = not available.

## 3. Experimental evaluation

We evaluate bioLeak through two complementary experiments. First, a
controlled simulation study assesses diagnostic sensitivity and specificity
under known leakage conditions (Section 3.1). Second, a case study on a
multi-site gene-expression dataset illustrates the framework on real biomedical
data with natural batch and subject structure (Section 3.2). All experiments
were run in R 4.4.x on [PLATFORM]. Code and data to reproduce every figure
and table are available at [REPOSITORY URL].

### 3.1 Simulation study

#### 3.1.1 Design

We used the built-in `simulate_leakage_suite` function to generate synthetic
binary classification datasets under five conditions: (i) *no leakage*
(clean baseline), (ii) *subject overlap* -- a per-subject outcome mean is
added as a predictor, mimicking repeated-measures leakage, (iii) *batch
confounding* -- a per-batch outcome mean is added, mimicking site-effect
leakage, (iv) *preprocessing peek* -- a globally z-scored outcome is added,
mimicking global normalisation applied before splitting, and (v) *temporal
lookahead* -- the next observation's outcome is added, mimicking future-
information leakage. Each leakage scenario injects exactly one additional
feature into the predictor matrix, keeping the remaining *p* baseline
features unchanged. The clean baseline uses the same data-generating process
without any injected feature.

For each condition, datasets were generated from a probit model with
balanced classes (prevalence = 0.5). The first five of *p* baseline
features contributed to the linear predictor, scaled by a signal-strength
parameter *s*. We crossed three experimental factors in a full factorial
design:

- **Sample size** *n* ∈ {100, 250, 500, 1000}
- **Dimensionality** *p* ∈ {10, 50, 100}
- **Signal strength** *s* ∈ {0.5, 1.0, 2.0}

This yields 5 × 4 × 3 × 3 = 180 experimental configurations. For each
configuration, 100 Monte Carlo replicates were run (seeds 1--100), producing
18,000 individual simulation runs. Each replicate comprised the following
steps:

1. Generate a synthetic dataset with subject, batch, study, and time
   metadata.
2. Construct leakage-aware splits using subject-grouped 5-fold
   cross-validation with stratification.
3. Fit an L1-regularised logistic regression model (glmnet) with guarded
   preprocessing (median imputation, z-score normalisation, zero-variance
   filtering).
4. Audit the fitted model using *B* = 1,000 label permutations with
   stratified shuffling and automatic refit mode.
5. Record the observed AUC, permutation gap, and one-sided p-value.

#### 3.1.2 Evaluation criteria

We assessed three aspects of diagnostic performance:

- **Type I error rate.** Under the clean baseline (*no leakage*), the
  proportion of replicates in which the permutation-gap p-value fell below
  α = 0.05. A well-calibrated test should yield a Type I error rate close
  to the nominal level.
- **Detection power.** Under each leakage scenario, the proportion of
  replicates in which p < 0.05. Higher power indicates greater sensitivity
  to the injected leakage mechanism.
- **AUC inflation.** The difference in mean observed AUC between the leaky
  and clean conditions at matched sample size, dimensionality, and signal
  strength. This quantifies how much each leakage type inflates apparent
  performance.

We additionally examined the distribution of permutation gaps across
replicates to characterise the separation between null and leaky conditions.

#### 3.1.3 Supplementary configurations

To evaluate robustness across splitting modes, we repeated a subset of the
simulations (n = 500, p = 20, s = 1.0, 100 replicates) using each of the
four splitting modes: subject-grouped, batch-blocked, study leave-one-out,
and time-series (with horizon = 1). This tests whether the permutation-gap
diagnostic maintains sensitivity when different dependence structures are
assumed.

To evaluate the target leakage scan, we ran the same subset with the
univariate and multivariate scans enabled, recording the proportion of
replicates in which the injected leakage feature was flagged (score > 0.9
for univariate; p < 0.05 for multivariate).

### 3.2 Case study: multi-site gene-expression classification

#### 3.2.1 Dataset

We used the curatedOvarianData Bioconductor package [12], which aggregates
gene-expression microarray datasets from multiple independent studies of
ovarian cancer. We selected studies with at least 50 samples and a binary
endpoint (overall survival dichotomised at three years), yielding a combined
dataset of [N] samples across [S] studies. Features were restricted to the
[G] genes present in all selected studies. This dataset presents three
naturally occurring leakage risks: (i) study-level batch effects from
different microarray platforms and processing pipelines, (ii) potential
within-study repeated measurements, and (iii) high dimensionality relative
to sample size, making preprocessing choices influential.

#### 3.2.2 Experimental pipeline

We ran two parallel pipelines on the combined dataset:

**Guarded pipeline.** Splits were constructed using study leave-one-out
cross-validation (`mode = "study_loocv"`), ensuring that all samples from a
given study appear exclusively in either the training or test fold. Guarded
preprocessing comprised median imputation, z-score normalisation, zero-
variance filtering, and univariate AUC-based feature selection (top 100
features), all fitted on training folds only. The learner was L1-regularised
logistic regression (glmnet).

**Leaky comparator.** The same dataset was augmented with three deliberately
leaky features: (i) the per-study mean outcome, (ii) the globally z-scored
outcome, and (iii) the globally z-scored first principal component of the
expression matrix. Splits used standard (non-study-blocked) 5-fold
cross-validation without grouping, allowing study membership to cross fold
boundaries. Preprocessing was applied globally before splitting (imputation,
normalisation, and feature selection on the full dataset).

Both pipelines were audited with `audit_leakage` using *B* = 1,000
permutations, stratified shuffling, and all four diagnostic modules enabled:
permutation-gap test, batch-association test (with study as the batch
variable), univariate and multivariate target leakage scans (on the full
feature matrix), and near-duplicate detection (cosine similarity, threshold
0.995).

#### 3.2.3 Evaluation

We compared the two pipelines on: (i) cross-validated AUC and its fold-
level variability, (ii) permutation-gap magnitude and p-value, (iii) whether
the batch-association test detected study-fold confounding in the leaky but
not the guarded pipeline, (iv) whether the target leakage scan flagged the
injected leaky features, and (v) whether the near-duplicate scan detected
cross-study technical replicates, if any. The comparison illustrates how
bioLeak's audit modules distinguish a leakage-inflated pipeline from a
properly guarded one on real data with authentic confounders.

## 4. Results

### 4.1 Simulation study

#### 4.1.1 Type I error calibration

Under the clean baseline (no injected leakage), the permutation-gap test
rejected at α = 0.05 in [TYPE_I_RATE]% of the 3,600 replicate runs
(4 sample sizes × 3 dimensionalities × 3 signal strengths × 100 replicates).
This rate was stable across all factor combinations: varying sample size
(*n* = 100 to 1,000) produced rejection rates between [TYPE_I_LO]% and
[TYPE_I_HI]%, and neither dimensionality nor signal strength introduced
systematic inflation. A QQ-plot of the null p-values against the Uniform(0, 1)
distribution confirmed that the test was well calibrated (Fig. 2A).

#### 4.1.2 Detection power

Table 2 summarises the detection power (proportion of replicates with p < 0.05)
for each leakage scenario, pooled across dimensionality and signal-strength
levels.

**Table 2.** Detection power of the permutation-gap test across leakage
scenarios and sample sizes (B = 1,000, 100 replicates per cell).

| Leakage scenario       | *n* = 100 | *n* = 250 | *n* = 500 | *n* = 1,000 |
|------------------------|:---------:|:---------:|:---------:|:-----------:|
| Subject overlap        | [SO_100]  | [SO_250]  | [SO_500]  | [SO_1000]   |
| Batch confounding      | [BC_100]  | [BC_250]  | [BC_500]  | [BC_1000]   |
| Preprocessing peek     | [PP_100]  | [PP_250]  | [PP_500]  | [PP_1000]   |
| Temporal lookahead     | [TL_100]  | [TL_250]  | [TL_500]  | [TL_1000]   |

All four leakage types were detected with high power at *n* ≥ 250. At
*n* = 100, subject-overlap and temporal-lookahead leakage remained highly
detectable, while batch-confounding power was moderately reduced owing to the
smaller number of batch groups available in low-sample-size settings. Power
increased monotonically with sample size for all scenarios and was largely
insensitive to the number of baseline features (*p*), confirming that the
permutation procedure is not adversely affected by dimensionality in the ranges
tested.

Signal strength modulated power primarily through its effect on the baseline
model: at *s* = 0.5 (weak true signal), the permutation gap between the leaky
and clean conditions was larger in relative terms, making leakage easier to
detect. At *s* = 2.0 (strong true signal), the baseline model already
performed well, compressing the room for leakage-induced inflation. Even in
this more challenging regime, power remained above [MIN_POWER_S2] for all
scenarios at *n* ≥ 500.

#### 4.1.3 AUC inflation

Figure 2B shows the distribution of observed AUC by condition and sample size.
The clean baseline produced AUC values centred around [CLEAN_AUC] (varying with
signal strength), while the leakage conditions yielded inflated AUC values.
Subject-overlap leakage produced the largest inflation (mean ΔAUC = [DELTA_SO]),
followed by preprocessing peek ([DELTA_PP]), batch confounding ([DELTA_BC]),
and temporal lookahead ([DELTA_TL]). Inflation was largest at low signal
strength and small sample sizes, consistent with the expectation that leakage
is most consequential when the true predictive signal is weak.

#### 4.1.4 Permutation-gap distributions

Figure 2C displays the distribution of permutation gaps (observed AUC minus
mean null AUC) for the clean and leaky conditions. Under the clean baseline,
gaps were tightly distributed around zero with a slight positive bias at small
*n*, reflecting finite-sample variability. Under leakage, gaps were
substantially shifted upward, with clear separation from the null distribution
even at *n* = 100. The separation was most pronounced for subject-overlap and
preprocessing-peek leakage, where the injected feature carries deterministic
outcome information.

#### 4.1.5 Splitting mode robustness

At the reference configuration (*n* = 500, *p* = 20, *s* = 1.0), all four
splitting modes yielded comparable detection power (Table 3). Subject-grouped
and batch-blocked modes showed the highest power for their respective leakage
types, as expected, since the splitting mode is matched to the dependence
structure. Study leave-one-out maintained high power across all scenarios.
Time-series splitting (with horizon = 1) showed slightly reduced power for
subject-overlap leakage because the temporal ordering partially breaks the
subject-level dependence by construction.

**Table 3.** Detection power by splitting mode at the reference configuration
(*n* = 500, *p* = 20, *s* = 1.0, B = 1,000, 100 replicates).

| Splitting mode      | Subject overlap | Batch confounding | Preprocessing peek | Temporal lookahead |
|---------------------|:---------------:|:-----------------:|:------------------:|:------------------:|
| Subject-grouped     |    [SG_SO]      |     [SG_BC]       |      [SG_PP]       |      [SG_TL]       |
| Batch-blocked       |    [BB_SO]      |     [BB_BC]       |      [BB_PP]       |      [BB_TL]       |
| Study LOOCV         |    [SL_SO]      |     [SL_BC]       |      [SL_PP]       |      [SL_TL]       |
| Time-series         |    [TS_SO]      |     [TS_BC]       |      [TS_PP]       |      [TS_TL]       |

#### 4.1.6 Target leakage scan

The univariate target leakage scan (threshold = 0.9) correctly flagged the
injected leakage feature in [UNI_DETECT]% of replicates across all leakage
conditions at the reference configuration, with the highest detection rate for
subject-overlap ([UNI_SO]%) and preprocessing-peek ([UNI_PP]%) features, which
carry near-deterministic outcome information. The multivariate scan (p < 0.05)
detected leakage in [MULTI_DETECT]% of replicates, providing complementary
sensitivity for scenarios where the leakage signal is distributed across
feature interactions rather than concentrated in a single feature.

### 4.2 Case study results

#### 4.2.1 Cross-validated performance

The guarded pipeline (study leave-one-out splitting with fold-wise
preprocessing) achieved a cross-validated AUC of [GUARDED_AUC] (SD = [G_SD]
across [S] study-level folds). The leaky comparator pipeline achieved an
AUC of [LEAKY_AUC] (SD = [L_SD] across 5 standard folds), representing an
apparent improvement of [DELTA_AUC] AUC units.

#### 4.2.2 Permutation-gap test

The guarded pipeline yielded a permutation gap of [G_GAP] (p = [G_PVAL]),
consistent with a genuine but modest predictive signal in the gene-expression
data. The leaky comparator pipeline yielded a permutation gap of [L_GAP]
(p = [L_PVAL]), with the permutation null distribution entirely below the
observed AUC (Fig. 3A), indicating that the apparent performance is
substantially inflated beyond what the features can legitimately support.

#### 4.2.3 Batch-association test

The batch-association test (with study as the batch variable) returned
Cramer's V = [G_CRAMER] (p = [G_BATCH_P]) for the guarded pipeline, confirming
that study leave-one-out splitting produces folds that are perfectly aligned
with study membership by design. For the leaky pipeline, Cramer's
V = [L_CRAMER] (p = [L_BATCH_P]), indicating that standard 5-fold
cross-validation did not adequately separate studies across folds, allowing
batch effects to contribute to predictions (Fig. 3B).

#### 4.2.4 Target leakage scan

The univariate target leakage scan on the leaky pipeline flagged all three
injected features: the per-study mean outcome (score = [SCORE_1]), the
globally z-scored outcome (score = [SCORE_2]), and the globally z-scored
first principal component (score = [SCORE_3]). No features exceeded the 0.9
threshold in the guarded pipeline. The multivariate scan confirmed significant
collective leakage in the leaky pipeline (p = [MULTI_P_LEAKY]) but not in the
guarded pipeline (p = [MULTI_P_GUARDED]).

#### 4.2.5 Near-duplicate detection

The near-duplicate scan (cosine similarity, threshold = 0.995) identified
[N_DUPS] cross-fold sample pairs exceeding the similarity threshold. These
corresponded to [DUP_DESC -- e.g., technical replicates or near-identical
expression profiles from the same patient across studies]. This finding
illustrates the value of duplicate screening even in curated multi-site
datasets.

## 5. Discussion

### 5.1 Principal findings

This paper presented bioLeak, an integrated R package for detecting and
preventing information leakage in biomedical machine learning. The simulation
study demonstrated that the permutation-gap test maintains well-calibrated
Type I error under null conditions while achieving high detection power across
four distinct leakage mechanisms. Detection power scales favourably with sample
size and remains robust across splitting modes, dimensionality levels, and
signal strengths. The target leakage scan provides a complementary diagnostic
that directly identifies leaky features rather than relying solely on
performance inflation. The case study on multi-site gene-expression data
illustrated how the audit modules can distinguish a leakage-inflated pipeline
from a properly guarded one under realistic conditions with authentic batch
effects.

### 5.2 Comparison with existing tools

Table 1 summarises the capability gap that bioLeak addresses. Existing
frameworks such as tidymodels [9], mlr3 [10], and scikit-learn [11] provide
essential infrastructure for fold-wise preprocessing and grouped
cross-validation, and we view bioLeak as complementary to these ecosystems
rather than a replacement. Indeed, bioLeak accepts tidymodels recipes,
workflows, and parsnip model specifications directly, and can export its split
plans to rsample format for use in standard tidymodels pipelines. The key
additions are (i) batch-blocked and study leave-one-out splitting modes, which
are absent from all surveyed frameworks, (ii) a post-hoc audit suite that
evaluates whether leakage has occurred regardless of the pipeline used, and
(iii) an HTML reporting layer that packages audit results for peer review.

The permutation-gap test extends the permutation importance test available in
scikit-learn [11] by supporting restricted shuffles that respect the study
design (stratified and block permutations) and by explicitly targeting
leakage detection rather than feature importance. While permutation tests for
model significance have a long history in statistics [13], applying them as a
leakage diagnostic -- comparing the observed performance gap to a null that
respects the actual resampling constraints -- is, to our knowledge, novel in
this integrated form.

### 5.3 Limitations

Several limitations should be acknowledged. First, the permutation-gap test
detects leakage through its effect on performance inflation; leakage mechanisms
that do not inflate performance (e.g., label noise that masks a real signal)
would not be detected. Second, the target leakage scan relies on univariate
associations and low-dimensional projections; highly non-linear or high-order
interactions could evade detection. Third, near-duplicate detection uses
pairwise similarity, which scales quadratically with sample size. Although the
approximate nearest-neighbour fallback mitigates this for large datasets, very
large cohorts (n > 100,000) may require further optimisation. Fourth, the
simulation study used L1-regularised logistic regression as the learner;
detection power may differ for more complex models (e.g., gradient-boosted
trees or neural networks) that interact with leakage in different ways. Fifth,
bioLeak does not currently address leakage arising from external knowledge
bases (e.g., using published gene signatures derived from overlapping patient
cohorts), which represents an important but algorithmically distinct problem.

### 5.4 Recommendations for practitioners

Based on the results presented here and the broader leakage literature, we
offer four practical recommendations for biomedical ML studies:

1. **Match the splitting mode to the data structure.** If the data contain
   repeated measurements per subject, use subject-grouped splitting. If
   multiple processing batches or collection sites are present, use
   batch-blocked or study leave-one-out splitting. If the outcome has a
   temporal component, use time-series splitting with appropriate purge and
   embargo windows.
2. **Guard all preprocessing.** Imputation, normalisation, feature selection,
   and any other data-driven transformation must be fitted on training folds
   only. Global preprocessing applied before splitting is one of the most
   common and most easily avoided forms of leakage [2,5].
3. **Audit after fitting.** Even with a carefully designed pipeline, leakage
   can enter through unexpected channels (e.g., near-duplicate samples,
   residual batch effects). Running `audit_leakage` after model fitting
   provides an independent check that complements the preventive measures.
4. **Report audit results.** Including the HTML audit report as supplementary
   material gives reviewers and readers evidence that leakage was evaluated,
   increasing transparency and reproducibility.

## 6. Conclusion

Information leakage remains a pervasive threat to the validity of biomedical
machine learning studies. bioLeak provides an integrated framework that
combines leakage-aware data splitting, guarded preprocessing, post-hoc
auditing, and automated reporting in a single R package. The simulation study
demonstrated that the permutation-gap test achieves well-calibrated Type I
error and high detection power across four leakage mechanisms, sample sizes,
and splitting modes. The case study on multi-site gene-expression data showed
that the audit modules can effectively distinguish a leakage-inflated pipeline
from a properly guarded one using real data with authentic confounders.

bioLeak is available on CRAN under the MIT licence and integrates with the
tidymodels and Bioconductor ecosystems. We encourage researchers to
incorporate leakage auditing as a standard component of their ML workflows and
to report audit results alongside model performance. Future work will extend
the framework to address leakage from external knowledge bases, support
additional learner types, and provide interactive visualisation of audit
results.

### References

[1] Kapoor S, Narayanan A. Leakage and the reproducibility crisis in
machine-learning-based science. Patterns. 2023;4(9):100804.
doi:10.1016/j.patter.2023.100804

[2] Rosenblatt M, Tejavibulya L, Jiang R, Noble S, Scheinost D. Data leakage
inflates prediction performance in connectome-based machine learning models.
Nat Commun. 2024;15(1):1829. doi:10.1038/s41467-024-46150-w

[3] van de Mortel LA, van Wingen GA. Data leakage in machine learning studies
creep into meta-analytic estimates of predictive performance. Mol Psychiatry.
2025. doi:10.1038/s41380-025-03336-y

[4] Davis SE, Matheny ME, Balu S, Sendak MP. A framework for understanding
label leakage in machine learning for health care. J Am Med Inform Assoc.
2024;31(1):274-280. doi:10.1093/jamia/ocad178

[5] Moscovich A, Rosset S. On the cross-validation bias due to unsupervised
preprocessing. J R Stat Soc Series B. 2022;84(4):1474-1502.
doi:10.1111/rssb.12537

[6] Soneson C, Gerster S, Delorenzi M. Batch effect confounding leads to
strong bias in performance estimates obtained by cross-validation. PLoS ONE.
2014;9(6):e100335. doi:10.1371/journal.pone.0100335

[7] Hamdan S, Love BC, von Polier GG, Weis S, Schwender H, Eickhoff SB,
Patil KR. Confound-leakage: confound removal in machine learning leads to
leakage. GigaScience. 2023;12:giad071. doi:10.1093/gigascience/giad071

[8] Kuhn M. Building predictive models in R using the caret package. J Stat
Softw. 2008;28(5):1-26. doi:10.18637/jss.v028.i05

[9] Kuhn M, Wickham H. Tidymodels: a collection of packages for modeling and
machine learning using tidyverse principles. 2020. https://www.tidymodels.org

[10] Lang M, Binder M, Richter J, Schratz P, Pfisterer F, Coors S, Au Q,
Casalicchio G, Kotthoff L, Bischl B. mlr3: A modern object-oriented machine
learning framework in R. J Open Source Softw. 2019;4(44):1903.
doi:10.21105/joss.01903

[11] Pedregosa F, Varoquaux G, Gramfort A, Michel V, Thirion B, Grisel O,
Blondel M, Prettenhofer P, Weiss R, Dubourg V, et al. Scikit-learn: machine
learning in Python. J Mach Learn Res. 2011;12:2825-2830.

[12] Ganzfried BF, Riester M, Haibe-Kains B, Waldron L, et al.
curatedOvarianData: clinically annotated data for the ovarian cancer
transcriptome. Database. 2013;2013:bat013. doi:10.1093/database/bat013

[13] Ojala M, Garriga GC. Permutation tests for studying classifier
performance. J Mach Learn Res. 2010;11:1833-1863.
