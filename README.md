# bioLeak: Leakage-Safe Modeling and Auditing for Omics and Clinical Data

`bioLeak` is an R toolkit for building trustworthy predictive models with
biomedical data. It helps prevent data leakage during preprocessing,
resampling, and evaluation, and provides diagnostics to quantify residual
leakage risk.

The core question it helps answer is: "Is my model benefitting from real
signal, or from information that leaked across folds, batches, or duplicates?"

---

## Highlights

- Leakage-aware resampling for grouped, batch, study, and time series data
- Guarded preprocessing (train-only imputation, scaling, filtering, selection)
- Cross-validated fitting for common learners and custom learners
- Leakage audits: permutation gap, batch association, and duplicate detection
- One-click HTML audit report
- S4 objects with provenance for reproducibility

---

## Installation

Install the development version from GitHub:

```r
remotes::install_github("selcukorkmaz/bioLeak")
```

Optional dependencies are used when you call the relevant functionality.
You may be prompted to install:

- Modeling: `glmnet`, `ranger`
- Metrics: `pROC`, `PRROC`, `survival`
- Imputation: `VIM`, `mice`, `missForest`
- Parallel: `future.apply`
- Duplicate detection: `RANN`
- HTML report: `rmarkdown` (requires Pandoc)

---

## Quick start

```r
library(bioLeak)

# Example data.frame with outcome and subject ids
splits <- make_splits(
  df,
  outcome = "outcome",
  mode = "subject_grouped",
  group = "subject_id",
  stratify = TRUE,
  v = 5
)

fit <- fit_resample(
  df,
  outcome = "outcome",
  splits = splits,
  preprocess = list(
    impute = list(method = "median"),
    normalize = list(method = "zscore")
  ),
  learner = "glmnet",
  metrics = c("auc", "pr_auc", "accuracy")
)

summary(fit)

audit <- audit_leakage(fit, metric = "auc", B = 100)
summary(audit)

audit_report(audit, output_dir = ".")
```

---

## Core capabilities

### 1) Leakage-resistant splitting

`make_splits()` supports common biomedical scenarios:

| Mode              | Use case |
|-------------------|----------|
| `subject_grouped` | Multiple samples per individual; keep each subject in one fold |
| `batch_blocked`   | Batch or platform effects; isolate batches across folds |
| `study_loocv`     | Leave-one-study-out validation |
| `time_series`     | Rolling origin folds with optional horizon |

Splits are stored in a `LeakSplits` object with metadata, seed, and hash.

### 2) Guarded preprocessing and fitting

`fit_resample()` orchestrates guarded preprocessing and model training:

- Imputation, scaling, filtering, and feature selection are fit on training
  samples only, then applied to test samples.
- Supports `glmnet` and `ranger` out of the box.
- Supports custom learners via `custom_learners`.

Custom learner example:

```r
custom <- list(
  glm = list(
    fit = function(x, y, task, weights, ...) {
      stats::glm(y ~ ., data = as.data.frame(x),
                 family = stats::binomial(), weights = weights)
    },
    predict = function(object, newdata, task, ...) {
      as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
                                type = "response"))
    }
  )
)

fit2 <- fit_resample(
  df,
  outcome = "outcome",
  splits = splits,
  learner = "glm",
  custom_learners = custom,
  metrics = "accuracy"
)
```

Parsnip model_spec example:

```r
spec <- parsnip::boost_tree(mode = "classification", trees = 200) |>
  parsnip::set_engine("xgboost")
fit3 <- fit_resample(df, outcome = "outcome", splits = splits,
                     learner = spec, metrics = "auc")
```

### 3) Leakage auditing

`audit_leakage()` quantifies residual leakage risk:

- Permutation significance test: observed metric vs permuted null
- Batch association: chi-square tests with Cramer V
- Duplicate detection: cosine or Pearson similarity

Duplicate detection requires a reference feature matrix:

```r
audit <- audit_leakage(
  fit,
  metric = "auc",
  B = 100,
  X_ref = df[, c("x1", "x2")],
  sim_method = "pearson",
  sim_threshold = 0.99
)
```

### 4) HTML audit report

`audit_report()` renders a shareable HTML summary:

```r
# from a LeakAudit
audit_report(audit, output_dir = ".")

# or directly from a LeakFit
audit_report(fit, metric = "auc", B = 100, output_dir = ".")
```

### 5) Per-learner audits

If you fit multiple learners, audit each one separately:

```r
audits <- audit_leakage_by_learner(fit, metric = "auc", B = 100)
audits
```

### 6) Standalone guarded imputation

`impute_guarded()` performs leakage-safe imputation for train/test splits:

```r
imp <- impute_guarded(train_df, test_df,
                      method = "missForest",
                      winsor = TRUE,
                      winsor_thresh = 3)
train_clean <- imp$train
test_clean  <- imp$test
```

---

## Data types

- `data.frame` and `matrix` inputs are supported.
- `SummarizedExperiment` inputs are supported; metadata are read from `colData`.

---

## Simulation suite

`simulate_leakage_suite()` runs Monte Carlo simulations to validate leakage
signals across modes and leakage types.

```r
res <- simulate_leakage_suite(
  n = 300, p = 10, mode = "subject_grouped",
  learner = "ranger", leakage = "subject_overlap",
  seeds = 1:3, parallel = FALSE
)
head(res)
```

---

## Reproducibility

- Splits include hashes and metadata for auditing.
- Guarded preprocessing stores fold-level state.
- Audit objects include a provenance trail.

---

## Performance notes

- Duplicate detection on large cohorts can be accelerated by `RANN` and a
  higher `sim_threshold`.
- Permutation tests (`B`) can be expensive; reduce `B` for quick checks and
  increase for final reporting.
- Parallel resampling uses `future.apply` if available.

---

## Limitations

- Built-in learners are `glmnet` and `ranger`. Use `custom_learners` for others.
- Classification supports binary outcomes; multiclass is not supported yet.
- `cindex` is implemented for regression-style concordance, not survival data.

---

## Testing

Run tests locally:

```r
testthat::test_dir("tests/testthat")
```

CI uses GitHub Actions (`.github/workflows/R-CMD-check.yaml`).

---

## Getting help

- See function docs in `man/` or `?fit_resample`.
- Open issues with reproducible examples.
- Contributions are welcome.

---

## License

GPL-3
