<img src="man/figures/bioLeak_hex.png" align="right" width="95"/>

# bioLeak: Leakage-Safe Modeling and Auditing for Genomic and Clinical Data

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/bioLeak)](https://CRAN.R-project.org/package=bioLeak)
[![R-CMD-check](https://github.com/selcukorkmaz/bioLeak/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/selcukorkmaz/bioLeak/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`bioLeak` is an R package for detecting, quantifying, and diagnosing data leakage in biomedical machine-learning workflows. It provides leakage-resistant resampling, guarded preprocessing, post-hoc auditing, and inference tools for cross-validation and related evaluation settings.

## Purpose and scope

**In scope:**
- **Preprocessing leakage:** global imputation, scaling, filtering, or feature selection applied before resampling.
- **Dependence leakage:** repeated measures, subject-level grouping, batch/site/study effects, and near-duplicate samples.
- **Resampling violations:** group overlap, study holdout, time-ordered evaluation, and multi-axis split constraints.
- **Diagnostic evidence:** permutation-based performance gaps, batch/fold association tests, target leakage scans, duplicate detection, and mechanism-level risk summaries.

**Out of scope:**
- Proving the absence of leakage or guaranteeing unbiased performance.
- Production deployment tooling.
- Unsupervised learning.
- Broad, non-leakage-oriented data-quality diagnostics.

## Why bioLeak is needed

Standard cross-validation assumes independent samples and exchangeable labels. Biomedical datasets often violate these assumptions due to repeated measures, site effects, batch structure, and temporal dependence. These violations can inflate performance metrics even when a model does not generalize. `bioLeak` enforces leakage-aware resampling and provides post-hoc diagnostics that estimate how much apparent performance could be driven by leakage or confounding.

## Core functionality

### Splitting and validation

| Function | Description |
|---|---|
| `make_split_plan()` | Leakage-aware splits: subject-grouped, batch-blocked, study leave-out, time-series, and N-axis combined modes with compact storage option |
| `check_split_overlap()` | Explicit overlap-invariant validation across declared grouping axes |
| `as_rsample()` | Convert `LeakSplits` to an `rsample` rset for tidymodels interoperability |

### Guarded preprocessing and modeling

| Function | Description |
|---|---|
| `fit_resample()` | Cross-validated fitting with train-only imputation, normalization, filtering, and feature selection; supports binomial, multiclass, regression, and survival tasks |
| `tune_resample()` | Nested hyperparameter tuning via tidymodels `tune`/`dials` with leakage-aware outer splits and hyperparameter aggregation across folds |
| `impute_guarded()` | Standalone train-only imputation (median, knn, missForest, mice) |
| `guard_to_recipe()` | Convert guarded preprocessing specifications to `recipes` pipelines |

### Auditing and diagnostics

| Function | Description |
|---|---|
| `audit_leakage()` | Permutation gap test, batch/study association tests, univariate and multivariate target leakage scans, near-duplicate detection, and mechanism-level risk summaries |
| `audit_leakage_by_learner()` | Multi-learner auditing for multi-model fits |
| `audit_report()` | Self-contained HTML summary of audit results |
| `calibration_summary()` | Probability calibration checks |
| `confounder_sensitivity()` | Sensitivity analysis for confounding effects |

### Inference

| Function | Description |
|---|---|
| `delta_lsi()` | Leakage sensitivity index with BCa confidence intervals and sign-flip inference, including blocked exchangeability for time-series |
| `cv_ci()` | Cross-validation confidence intervals with Nadeau-Bengio correction |

### Simulation and benchmarking

| Function | Description |
|---|---|
| `simulate_leakage_suite()` | Generate controlled leakage scenarios for benchmarking audit sensitivity |
| `benchmark_leakage_suite()` | Reproducible modality-by-mechanism benchmark grids with detection-rate summaries |

### Visualization

`plot_calibration()`, `plot_confounder_sensitivity()`, `plot_fold_balance()`, `plot_overlap_checks()`, `plot_perm_distribution()`, `plot_time_acf()`

## Installation

Requires R >= 4.3.

From CRAN:

```r
install.packages("bioLeak")
```

Development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("selcukorkmaz/bioLeak")
```

Non-obvious dependencies:
- `SummarizedExperiment` and `BiocGenerics` are Bioconductor packages (installed automatically by `remotes`, but can be installed manually with `BiocManager::install()` if needed).
- Optional packages enable specific features: `glmnet`, `ranger`, `pROC`, `PRROC`, `survival`, `future.apply`, `RANN`, `rmarkdown`, `tune`, `dials`.

## Minimal working example

```r
library(bioLeak)

set.seed(1)
n_subject <- 40
rep_per_subject <- 3
n <- n_subject * rep_per_subject

subject <- rep(seq_len(n_subject), each = rep_per_subject)
batch <- rep(seq_len(6), length.out = n)

# Subject-level latent risk creates dependence across repeated measures
subj_risk <- rnorm(n_subject, sd = 1)
x1 <- subj_risk[subject] + rnorm(n, sd = 0.5)
x2 <- rnorm(n)
x3 <- rnorm(n)
p_subj <- stats::plogis(1.5 * subj_risk)
outcome <- factor(ifelse(runif(n) < p_subj[subject], "case", "control"),
                  levels = c("control", "case"))

df <- data.frame(subject, batch, outcome, x1, x2, x3)

# Leakage-aware splits (subjects do not cross folds)
splits <- make_split_plan(
  df,
  outcome = "outcome",
  mode = "subject_grouped",
  group = "subject",
  v = 5,
  stratify = TRUE,
  seed = 1
)

# Guarded pipeline (train-only preprocessing)
spec <- parsnip::logistic_reg(mode = "classification") |>
  parsnip::set_engine("glm")

fit_guarded <- fit_resample(
  df,
  outcome = "outcome",
  splits = splits,
  learner = spec,
  metrics = "auc",
  preprocess = list(
    impute = list(method = "median"),
    normalize = list(method = "zscore"),
    filter = list(var_thresh = 0),
    fs = list(method = "none")
  ),
  refit = FALSE,
  seed = 1
)

# Leaky comparator: add a leakage feature computed on the full dataset
df_leaky <- within(df, {
  leak_subject <- ave(as.numeric(outcome == "case"), subject, FUN = mean)
})

fit_leaky <- fit_resample(
  df_leaky,
  outcome = "outcome",
  splits = splits,
  learner = spec,
  metrics = "auc",
  preprocess = list(
    impute = list(method = "none"),
    normalize = list(method = "none"),
    filter = list(var_thresh = 0),
    fs = list(method = "none")
  ),
  refit = FALSE,
  seed = 1
)

# Use unstratified permutations here to avoid warnings in small grouped data
audit_guarded <- audit_leakage(
  fit_guarded,
  metric = "auc",
  B = 30,
  perm_stratify = FALSE,
  X_ref = df[, c("x1", "x2", "x3")]
)

audit_leaky <- audit_leakage(
  fit_leaky,
  metric = "auc",
  B = 30,
  perm_stratify = FALSE,
  X_ref = df_leaky[, c("x1", "x2", "x3", "leak_subject")]
)

summary(fit_guarded)
summary(fit_leaky)
summary(audit_guarded)
summary(audit_leaky)
```

Interpretation notes:
- If the leaky comparator shows higher AUC and `leak_subject` ranks near the top of the target leakage scan, the performance gap is likely inflated by leakage.
- Similar guarded and leaky results do not prove the absence of leakage; they only reduce specific risks tested by the audit.

## Tidymodels interoperability

`bioLeak` integrates with the tidymodels ecosystem at multiple levels:

- **Splits:** `fit_resample()` and `tune_resample()` accept `rsample` rset/rsplit objects as `splits`. `as_rsample()` converts `LeakSplits` to an `rsample` rset.
- **Preprocessing:** `preprocess` can be a `recipes::recipe` (prepped on training folds, baked on test folds). `guard_to_recipe()` converts guarded preprocessing specifications into a `recipes` pipeline.
- **Learners:** `learner` accepts `parsnip::model_spec` or `workflows::workflow` objects.
- **Metrics:** `metrics` accepts `yardstick::metric_set` objects.

Note: When using recipes/workflows, the built-in guarded preprocessing list is not applied; ensure your recipe is leakage-safe.

Example (rsample + recipes + yardstick):
```r
if (requireNamespace("rsample", quietly = TRUE) &&
    requireNamespace("recipes", quietly = TRUE) &&
    requireNamespace("yardstick", quietly = TRUE)) {
  rs <- rsample::vfold_cv(df, v = 5)
  rec <- recipes::recipe(outcome ~ ., data = df) |>
    recipes::step_normalize(recipes::all_numeric_predictors())
  ys <- yardstick::metric_set(yardstick::roc_auc, yardstick::accuracy)

  fit_rs <- fit_resample(
    df,
    outcome = "outcome",
    splits = rs,
    learner = spec,
    preprocess = rec,
    metrics = ys,
    refit = FALSE
  )
}
```

## Supported tasks and learners

`bioLeak` supports four task types:
- **Binomial classification** (metrics: AUC, accuracy, Brier score, log-loss, sensitivity, specificity, F1, with optional threshold tuning)
- **Multiclass classification** (metrics: accuracy, macro-F1, log-loss)
- **Regression** (metrics: RMSE, MAE, R-squared)
- **Survival analysis** (metric: C-index)

The `learner` argument accepts parsnip model specs, workflows, or built-in learner strings (`"glm"`, `"glmnet"`, `"ranger"`, `"xgboost"`).

## Methodological assumptions

- The split mode matches the true dependence structure (subject, batch, study, or time).
- Leakage is inferred from performance gaps and diagnostic signals, not proven or ruled out.
- Permutation gaps assume the chosen resampling scheme reflects the intended evaluation setting.
- By default, permutation gaps use refit-based permutations when refit data are available (auto mode), falling back to fixed-prediction shuffles when they are not.
- Target leakage scans include univariate associations plus a multivariate/interaction check by default for supported tasks; proxies outside `X_ref` can still pass undetected.

## Interpretation guidance

- **Permutation gap:** large positive gaps indicate non-random signal; they do not by themselves indicate or refute leakage.
- **Batch/study association warnings** indicate that folds align with metadata; this can reflect leakage or study design constraints.
- **Target leakage flags** identify features overly aligned with the outcome; inspect data provenance before removing them. When enabled, FDR-adjusted columns (`p_value_adj`, `flag_fdr`) provide a more conservative screen.
- **Duplicate detection** flags near-identical samples across train/test by default; use `duplicate_scope = "all"` to include within-fold duplicates and review for data-quality issues.
- **Mechanism risk summary** (`audit@info$mechanism_summary`) provides a compact mechanism-level risk view across permutation, confounding, target-proxy, duplicate, and temporal signals.

Common misinterpretations:
- "Non-significant permutation test means no leakage": false.
- "High AUC implies good generalization": false if resampling is violated.
- "No flagged features means no leakage": false; audits are limited to available metadata and `X_ref`.

## Intended audience

Biomedical ML researchers, biostatisticians, and methodologists reviewing cross-validation and leakage risk.

## Citation

```r
citation("bioLeak")
```

When reporting results, include: split mode, grouping columns, random seeds, preprocessing steps, learner specification, metrics, and audit settings (B, target threshold, similarity method). Include both guarded and leaky comparator results when used.

## License

MIT
