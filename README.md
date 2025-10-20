# bioLeak: Leakage-Safe Modeling and Auditing for Omics and Clinical Data

`bioLeak` is an R toolkit for building trustworthy predictive models with
biomedical data. The package focuses on preventing data leakage during
pre-processing, cross-validation, and post-hoc evaluation—issues that often lead
to overly optimistic results when modelling high-dimensional omics or
longitudinal clinical datasets.

The toolbox helps data scientists, statisticians, and translational researchers
answer a core question: *"Is my model benefitting from true signal, or from
information that leaked across folds, batches, or duplicate samples?"*

---

## Why bioLeak?

Biomedical datasets combine multiple leakage risks:

* **Grouped sampling:** multiple measurements from the same subject or patient.
* **Batch and study effects:** platforms, centers, or sequencing runs encode the
  outcome of interest.
* **Time-dependent structures:** training on future observations contaminates
  prospective validation.
* **Aggressive preprocessing:** imputation, filtering, or feature selection that
  sees the full dataset often leaks signal into validation folds.
* **Hidden duplicates and confounders:** overlapping samples across studies can
  inflate performance.

`bioLeak` provides guarded splitting, preprocessing, fitting, and auditing tools
that explicitly address these pitfalls while remaining flexible enough to plug
into existing modelling workflows.

---

## Key Capabilities

* **Leakage-aware resampling:** `make_splits()` creates grouped, batch-blocked,
  study leave-one-out, and rolling-origin time-series cross-validation splits.
* **Guarded preprocessing:** `fit_resample()` wraps common learners (e.g.
  `glmnet`, `ranger`) with fold-specific preprocessing pipelines so that
  imputation, scaling, filtering, and feature selection are learned strictly on
  training data.
* **Robust imputation:** `impute_guarded()` performs winsorised mean/median
  imputation as well as advanced methods (`knn`, `mice`, `missForest`) without
  leaking statistics.
* **Model auditing:** `audit_leakage()` quantifies permutation gaps, tests fold
  association with batch metadata, and detects duplicates or near-duplicates in
  the feature space.
* **Rich result objects:** S4 classes (`LeakSplits`, `LeakFit`, `LeakAudit`)
  capture folds, metrics, predictions, preprocessing state, and provenance for
  downstream summaries and reproducibility.

---

## Installation

`bioLeak` is distributed as an R package. Install the development version from
source using [`remotes`](https://cran.r-project.org/package=remotes):

```r
# install.packages("remotes")
remotes::install_local("/path/to/bioLeak")
```

or directly from a Git repository:

```r
remotes::install_github("owner/bioLeak")
```

Most advanced functionality relies on optional dependencies. You will be
prompted to install packages like `glmnet`, `ranger`, `pROC`, `PRROC`,
`future.apply`, `VIM`, `mice`, `missForest`, and `SummarizedExperiment` when
those features are used.

---

## Quick Start

```r
library(bioLeak)

# Assume df is a sample-by-feature data.frame with outcome, subject_id, and batch columns
splits <- make_splits(df,
                      outcome = "outcome",
                      mode    = "subject_grouped",
                      group   = "subject_id",
                      stratify = TRUE,
                      v = 5,
                      repeats = 2)

fit <- fit_resample(df,
                    outcome  = "outcome",
                    splits   = splits,
                    preprocess = list(
                      impute   = list(method = "median"),
                      normalize = list(method = "zscore")
                    ),
                    learner = "glmnet",
                    metrics = c("auc", "pr_auc", "accuracy"))

summary(fit)
```

This workflow performs repeated, subject-grouped cross-validation with guarded
imputation and z-score normalization, followed by elastic-net logistic
regression. The `summary()` method reports fold-level metrics, aggregate
performance, and preprocessing diagnostics.

---

## Core Functions and Workflows

### 1. Leakage-resistant splitting

`make_splits()` accepts data frames, matrices, or `SummarizedExperiment`
objects. Supported modes cover common biomedical scenarios:

| Mode              | Use case                                                             |
|-------------------|----------------------------------------------------------------------|
| `subject_grouped` | Multiple samples per individual; keeps patient data in a single fold |
| `batch_blocked`   | Sequencing/assay batches as blocking factors                          |
| `study_loocv`     | Leave-one-study-out validation across cohorts                         |
| `time_series`     | Rolling-origin folds with optional prediction horizon                 |

Additional options include outcome stratification, repeated resampling, nested
inner cross-validation, and progress reporting. The resulting `LeakSplits`
object stores fold indices, metadata, and provenance hashes for auditing.

### 2. Guarded modelling with `fit_resample()`

`fit_resample()` orchestrates end-to-end leakage-safe modelling:

1. Extracts features from data frames or `SummarizedExperiment` assays.
2. Applies preprocessing steps (`impute`, `normalize`, `filter`, `fs`) that are
   fit exclusively on training samples within each fold.
3. Trains one or more learners (currently `glmnet`, `ranger`) and collects
   predictions on held-out folds.
4. Computes metrics such as AUC, PR AUC, accuracy, RMSE, or user-supplied
   functions.
5. Optionally refits the final model on the full dataset after resampling.

The resulting `LeakFit` object bundles per-fold predictions, preprocessing
state, learner fits, and data summaries. You can access metrics via
`fit@metrics` or `fit@metric_summary`, and predictions via `fit@predictions`.

### 3. Leakage auditing with `audit_leakage()`

`audit_leakage()` interrogates a `LeakFit` object to quantify residual risks:

* **Permutation gap:** compares observed performance against a null distribution
  generated by permuting outcomes within folds.
* **Batch association:** evaluates whether fold assignments correlate with
  metadata columns (batch, plate, site, study, etc.) using χ² tests and
  Cramér's V.
* **Duplicate detection:** identifies near-duplicate samples using cosine or
  Pearson similarity, optionally after rank-normalisation, with scalable nearest
  neighbour search.
* **Provenance trail:** hashes fold indices and seeds to detect tampering or
  reproducibility issues.

The `LeakAudit` result contains data frames for each diagnostic and can be
summarised or plotted to document model integrity.

### 4. Standalone guarded imputation

`impute_guarded()` offers a reusable imputation engine for train/test splits:

```r
imp <- impute_guarded(train_df, test_df,
                      method = "missForest",
                      winsor = TRUE,
                      winsor_thresh = 3,
                      parallel = TRUE)

train_clean <- imp$train
test_clean  <- imp$test
```

It returns imputed training and test frames, the fitted imputation model, and
(optional) outlier flags, ensuring reproducible preprocessing outside the full
resampling workflow.

---

## End-to-End Examples

### Example 1: Minimal binary classification

```r
library(bioLeak)
library(dplyr)

# Simulated gene-expression study with subject-level replicates
set.seed(1)
df <- tibble(
  subject_id = rep(1:60, each = 2),
  batch      = sample(letters[1:5], 120, replace = TRUE),
  outcome    = rbinom(120, 1, 0.5),
  gene1      = rnorm(120),
  gene2      = rnorm(120),
  gene3      = rnorm(120)
)

splits <- make_splits(df,
                      outcome = "outcome",
                      mode    = "subject_grouped",
                      group   = "subject_id",
                      stratify = TRUE,
                      v = 5)

fit <- fit_resample(df,
                    outcome = "outcome",
                    splits  = splits,
                    learner = "glmnet")

fit@metric_summary
```

### Example 2: Multi-modal SummarizedExperiment with batch blocking

```r
library(SummarizedExperiment)

# Construct SummarizedExperiment (assay = expression matrix; colData = metadata)
expr <- matrix(rnorm(5000), nrow = 100, ncol = 50)
colnames(expr) <- paste0("gene", seq_len(ncol(expr)))
outcome <- sample(c("Responder", "NonResponder"), 100, replace = TRUE)
subject <- paste0("P", sample(1:40, 100, replace = TRUE))
batch   <- sample(paste0("Plate", 1:4), 100, replace = TRUE)
se <- SummarizedExperiment(assays = list(counts = expr),
                           colData = DataFrame(outcome = outcome,
                                               subject = subject,
                                               batch   = batch))

splits <- make_splits(se,
                      outcome = "outcome",
                      mode    = "batch_blocked",
                      batch   = "batch",
                      stratify = TRUE,
                      v = 4)

fit <- fit_resample(se,
                    outcome  = "outcome",
                    splits   = splits,
                    preprocess = list(
                      impute   = list(method = "median"),
                      normalize = list(method = "zscore"),
                      filter   = list(var_thresh = 0.1)
                    ),
                    learner = c("glmnet", "ranger"),
                    metrics = c("auc", "pr_auc"))

audit <- audit_leakage(fit, batch_cols = c("batch"))
audit@permutation_gap
```

### Example 3: Time-series regression with duplicate detection

```r
# Longitudinal lab measurements with prediction horizon
set.seed(42)
time_df <- data.frame(
  patient   = rep(1:80, each = 5),
  time      = rep(seq_len(5), times = 80),
  outcome   = rnorm(400),
  biomarker = rnorm(400),
  clinical  = rnorm(400)
)

splits <- make_splits(time_df,
                      outcome = "outcome",
                      mode    = "time_series",
                      time    = "time",
                      v       = 5,
                      horizon = 1)

fit <- fit_resample(time_df,
                    outcome = "outcome",
                    splits  = splits,
                    learner = "ranger",
                    metrics = c("rmse"))

audit <- audit_leakage(fit,
                       n_perm = 100,
                       perm_stratify = FALSE,
                       X_ref = as.matrix(time_df[, c("biomarker", "clinical")]),
                       sim_method = "cosine",
                       sim_threshold = 0.99)

audit@duplicates
```

These examples illustrate how `bioLeak` supports datasets with grouping,
batching, and temporal dependencies while providing diagnostic evidence that the
resulting models generalize beyond hidden shortcuts.

---

## Reproducibility and Extensibility

* Every `LeakSplits` object stores seeds, fold indices, and optional metadata
  for reproducible resampling.
* `LeakFit` retains preprocessing states per learner, enabling deployment on new
  cohorts with identical pipelines.
* You can extend the framework by supplying custom metric functions, integrating
  new learners inside `fit_resample()`, or post-processing predictions stored in
  `fit@predictions`.
* Because inputs support `SummarizedExperiment`, `bioLeak` integrates naturally
  with Bioconductor workflows.

---

## Getting Help

* Browse the `man/` directory for function documentation.
* Open an issue with reproducible examples if you suspect leakage risks not yet
  covered by the package.
* Contributions (new guards, learners, diagnostics) are welcome via pull
  requests.

Protect your biomedical models from hidden leakage and build analyses that stand
up to rigorous review with `bioLeak`.
