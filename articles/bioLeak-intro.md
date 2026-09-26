# bioLeak: Leakage-Aware Biomedical Modeling

## Why bioLeak

bioLeak is a leakage-aware modeling toolkit for biomedical and
machine-learning analyses. Its purpose is to prevent and diagnose
information leakage across resampling workflows where training and
evaluation data are not truly independent because samples share
subjects, batches, studies, or time.

Standard workflows are often insufficient. Random, row-wise
cross-validation assumes samples are independent. Global preprocessing
(imputation, scaling, feature selection) done before resampling lets
test-fold information shape the training process. These choices inflate
performance and can lead to incorrect biomarker discovery, misleading
clinical signals, or models that fail in deployment.

Data leakage means any direct or indirect use of evaluation data in
training or feature engineering. In biomedical data, leakage commonly
appears as:

- repeated measurements from the same patient split across folds
- batch or study effects that align with the outcome
- duplicates or near-duplicates across train and test
- time-series lookahead or random splits that use future information
- preprocessing that uses all samples instead of training-only
  statistics

bioLeak addresses these issues with leakage-aware splitting, guarded
preprocessing that is fit only on training data, and audit diagnostics
that surface overlaps, confounding, and duplicates.

## Guided workflow

The sections below walk through a leakage-aware workflow from data setup
to audits. Each step includes a leaky failure mode and a corrected
alternative.

### Example data

``` r

library(bioLeak)

set.seed(123)
n <- 160
subject <- rep(seq_len(40), each = 4)
batch <- sample(paste0("B", 1:6), n, replace = TRUE)
study <- sample(paste0("S", 1:4), n, replace = TRUE)
time <- seq_len(n)

x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
linpred <- 0.7 * x1 - 0.4 * x2 + 0.2 * x3 + rnorm(n, sd = 0.5)
p <- stats::plogis(linpred)
outcome <- factor(ifelse(runif(n) < p, "case", "control"),
                  levels = c("control", "case"))

df <- data.frame(
  subject = subject,
  batch = batch,
  study = study,
  time = time,
  outcome = outcome,
  x1 = x1,
  x2 = x2,
  x3 = x3
)

df_leaky <- within(df, {
  leak_subject <- ave(as.numeric(outcome == "case"), subject, FUN = mean)
  leak_batch <- ave(as.numeric(outcome == "case"), batch, FUN = mean)
  leak_global <- mean(as.numeric(outcome == "case"))
})

df_time <- df
df_time$leak_future <- c(as.numeric(df_time$outcome == "case")[-1], 0)
predictors <- c("x1", "x2", "x3")


# Example data (first 6 rows)
head(df)
#>   subject batch study time outcome          x1         x2          x3
#> 1       1    B3    S4    1 control  0.21444531  1.0149432 -0.09736927
#> 2       1    B6    S3    2    case -0.32468591 -1.9927485  0.21615254
#> 3       1    B3    S1    3    case  0.09458353 -0.4272793  0.88246516
#> 4       1    B2    S3    4 control -0.89536336  0.1166373  0.20559750
#> 5       2    B2    S4    5 control -1.31080153 -0.8932076 -0.61643584
#> 6       2    B6    S2    6    case  1.99721338  0.3339029 -0.73479925

# Outcome class counts
as.data.frame(table(df$outcome))
#>      Var1 Freq
#> 1 control   91
#> 2    case   69
```

The table preview displays the metadata columns (`subject`, `batch`,
`study`, `time`), the binary outcome, and three numeric predictors
(`x1`, `x2`, `x3`).

The class count table shows the baseline prevalence. Stratified splits
try to preserve these proportions (for grouped modes, stratification is
applied at the group level).

### Create leakage-aware splits with `make_split_plan()`

[`make_split_plan()`](https://selcukorkmaz.github.io/bioLeak/reference/make_split_plan.md)
is the foundation. It returns a `LeakSplits` object with explicit
train/test indices (or compact fold assignments) and metadata. For
data.frame inputs, a unique `row_id` column is added automatically, so
`group = "row_id"` is the explicit way to request sample-wise CV. It
assumes that the grouping columns you provide are complete and that
samples sharing a group must not cross folds. Grouped stratification
uses the majority class per group and is ignored for `study_loocv`,
`time_series`, and survival outcomes. Misuse to avoid:

- using `group = "row_id"` (sample-wise CV) when subjects repeat
- using the wrong grouping column (for example batch instead of subject)
- using time-series CV with unsorted or missing time values
- relying on stratification when the outcome is missing, single-class at
  the group level, or a time/event outcome
- assuming `time_series` will always return exactly `v` folds

**Leaky example: row-wise CV when subjects repeat**

``` r

leaky_splits <- make_split_plan(
  df,
  outcome = "outcome",
  mode = "subject_grouped",
  group = "row_id",
  v = 5,
  repeats = 1,
  stratify = TRUE
)

cat("Leaky splits summary (sample-wise CV):\n")
#> Leaky splits summary (sample-wise CV):
leaky_splits
#> LeakSplits object (mode = subject_grouped, v = 5, repeats = 1)
#> Outcome: outcome | Stratified: TRUE | Nested: FALSE
#> ------------------------------------------------------
#>   fold repeat_id train_n test_n
#> 1    1         1     127     33
#> 2    2         1     128     32
#> 3    3         1     128     32
#> 4    4         1     128     32
#> 5    5         1     129     31
#> ------------------------------------------------------
#> Total folds: 5 | Hash: 54f1f1e5ad3c6f2f9cf0c7422b2f4425
```

The printed `LeakSplits` summary reports the split mode, number of
folds, and fold-level training and test set sizes.

Because `group = "row_id"`, each sample is treated as its own group. As
a result, repeated samples from the same subject may appear in both the
training and test sets, which introduces information leakage.

**Leakage-safe alternative: group by subject**

``` r

safe_splits <- make_split_plan(
  df,
  outcome = "outcome",
  mode = "subject_grouped",
  group = "subject",
  v = 5,
  repeats = 1,
  stratify = TRUE,
  seed = 10
)

cat("Leakage-safe splits summary (subject-grouped CV):\n")
#> Leakage-safe splits summary (subject-grouped CV):
safe_splits
#> LeakSplits object (mode = subject_grouped, v = 5, repeats = 1)
#> Outcome: outcome | Stratified: TRUE | Nested: FALSE
#> ------------------------------------------------------
#>   fold repeat_id train_n test_n
#> 1    1         1     124     36
#> 2    2         1     128     32
#> 3    3         1     128     32
#> 4    4         1     128     32
#> 5    5         1     132     28
#> ------------------------------------------------------
#> Total folds: 5 | Hash: e34ddb9f1415cf2a82f5b73d6e8975f6
```

Here, each subject is confined to a single fold. Test set sizes are
multiples of the number of samples per subject, confirming that subjects
were not split across folds.

The fold sizes remain comparable because `stratify = TRUE` balances
outcome proportions across folds while respecting subject-level
boundaries.

**Other leakage-aware modes**

``` r

batch_splits <- make_split_plan(
  df,
  outcome = "outcome",
  mode = "batch_blocked",
  batch = "batch",
  v = 4,
  stratify = TRUE
)

cat("Batch-blocked splits summary:\n")
#> Batch-blocked splits summary:
batch_splits
#> LeakSplits object (mode = batch_blocked, v = 4, repeats = 1)
#> Outcome: outcome | Stratified: TRUE | Nested: FALSE
#> ------------------------------------------------------
#>   fold repeat_id train_n test_n
#> 1    1         1      73     87
#> 2    2         1     138     22
#> 3    3         1     132     28
#> 4    4         1     137     23
#> ------------------------------------------------------
#> Total folds: 4 | Hash: 6815f595d5822e72c45cb3fe135574f7

study_splits <- make_split_plan(
  df,
  outcome = "outcome",
  mode = "study_loocv",
  study = "study"
)

cat("Study leave-one-out splits summary:\n")
#> Study leave-one-out splits summary:
study_splits
#> LeakSplits object (mode = study_loocv, v = 5, repeats = 1)
#> Outcome: outcome | Stratified: FALSE | Nested: FALSE
#> ------------------------------------------------------
#>   fold repeat_id train_n test_n
#> 1    1         1      96     64
#> 2    2         1     121     39
#> 3    3         1     131     29
#> 4    4         1     132     28
#> ------------------------------------------------------
#> Total folds: 4 | Hash: 51ebad5a2ade97fc502f662d7b1c9804

time_splits <- make_split_plan(
  df,
  outcome = "outcome",
  mode = "time_series",
  time = "time",
  v = 4,
  horizon = 2
)

cat("Time-series splits summary:\n")
#> Time-series splits summary:
time_splits
#> LeakSplits object (mode = time_series, v = 4, repeats = 1)
#> Outcome: outcome | Stratified: FALSE | Nested: FALSE
#> ------------------------------------------------------
#>   fold repeat_id train_n test_n
#> 1    2         1      39     40
#> 2    3         1      79     40
#> 3    4         1     119     40
#> ------------------------------------------------------
#> Total folds: 3 | Hash: dbde2d11033db6b03c4578feffd553d0

nested_splits <- make_split_plan(
  df,
  outcome = "outcome",
  mode = "subject_grouped",
  group = "subject",
  v = 3,
  nested = TRUE,
  stratify = TRUE
)

cat("Nested CV splits summary:\n")
#> Nested CV splits summary:
nested_splits
#> LeakSplits object (mode = subject_grouped, v = 3, repeats = 1)
#> Outcome: outcome | Stratified: TRUE | Nested: TRUE
#> ------------------------------------------------------
#>   fold repeat_id train_n test_n
#> 1    1         1     104     56
#> 2    2         1     104     56
#> 3    3         1     112     48
#> ------------------------------------------------------
#> Total folds: 3 | Hash: 44a1560f6d79af026e80218128c838ed
```

**RNG seed lineage**

bioLeak uses offset-based seeding so that every sub-operation gets a
deterministic, distinct seed without requiring isolated RNG streams:

| Function | Base seed | Per-unit offset |
|----|----|----|
| [`make_split_plan()`](https://selcukorkmaz.github.io/bioLeak/reference/make_split_plan.md) | `set.seed(seed)` | repeats: `seed + 1000 * r`; nested: `seed + 1` |
| [`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md) | `set.seed(seed)` | per-fold: `seed + fold_id` |
| [`audit_leakage()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md) | `set.seed(seed)` | per-permutation: `seed + b` |

These are **not** isolated RNG streams but simple offsets. Collisions do
not occur for practical parameter values (e.g., `repeats < 1000`,
`v < 1000`). Always set an explicit `seed` argument for reproducibility;
strict mode warns when `seed` is `NULL` or `NA`.

Each summary reports the number of folds and the fold sizes for the
corresponding split strategy.

- **Batch & Study:** The test set sizes vary because batches/studies are
  not evenly sized. `study_loocv` uses one study per fold, so `v` and
  `repeats` are ignored.

- **Time series:** The training set size increases across folds while
  the test set size follows the block size. Early blocks with no prior
  training data (or too-small test windows) are skipped, so you may see
  fewer than `v` folds.

**Combined (N-axis) mode**

When leakage can occur along multiple axes simultaneously (for example,
both subject and batch), `mode = "combined"` handles multi-axis
constraint splitting. The first element in `constraints` is the
**primary fold driver** (determines which groups go to test); subsequent
elements are **exclusion constraints** (training samples sharing
constraint-axis levels with the test set are removed).

``` r

# For combined mode, constraint-axis levels should not span all folds.
# Here each subject belongs to exactly one site, so site exclusion only
# removes training samples from the same site as the test subjects.
df_comb <- df
df_comb$site <- paste0("site", rep(1:8, each = 5)[df_comb$subject])

combined_splits <- make_split_plan(
  df_comb,
  outcome = "outcome",
  mode = "combined",
  constraints = list(
    list(type = "subject", col = "subject"),
    list(type = "batch",   col = "site")
  ),
  v = 4,
  stratify = TRUE,
  seed = 42
)

cat("Combined (N-axis) splits summary:\n")
#> Combined (N-axis) splits summary:
combined_splits
#> LeakSplits object (mode = combined, v = 4, repeats = 1)
#> Outcome: outcome | Stratified: TRUE | Nested: FALSE
#> ------------------------------------------------------
#>   fold repeat_id train_n test_n
#> 1    1         1     140     20
#> 2    2         1     120     40
#> 3    3         1     120     40
#> 4    4         1     100     60
#> ------------------------------------------------------
#> Total folds: 4 | Hash: a065b9e85e204119fcca3fd19c112015
```

The fold sizes may be smaller than single-axis modes because training
samples are excluded whenever their site (or any constraint-axis level)
also appears in the test set. You can add more than two axes by
appending additional elements to `constraints`.

The legacy `primary_axis`/`secondary_axis` parameters are still accepted
for backward compatibility but `constraints` supersedes them and
supports arbitrary numbers of axes.

``` r

# Verify zero overlap on all constraint axes
overlap_result <- check_split_overlap(combined_splits)
overlap_result
#>   fold repeat_id     col n_overlap pass
#> 1    1         1 subject         0 TRUE
#> 2    1         1    site         0 TRUE
#> 3    2         1 subject         0 TRUE
#> 4    2         1    site         0 TRUE
#> 5    3         1 subject         0 TRUE
#> 6    3         1    site         0 TRUE
#> 7    4         1 subject         0 TRUE
#> 8    4         1    site         0 TRUE
```

**Assumptions and intent by mode**

- `subject_grouped`: keeps all samples from a subject together
- `batch_blocked`: keeps batches/centers together (leave-one-group-out
  when `v` is at least the number of batches; otherwise multiple batches
  per fold)
- `study_loocv`: holds out each study in turn for external validation
- `time_series`: rolling-origin splits with an optional horizon to
  prevent lookahead
- `combined`: primary axis drives fold assignment; all additional
  constraint axes are enforced as exclusion constraints on the training
  set

Use `repeats` for repeated CV in grouped/batch modes (it is ignored for
`study_loocv` and `time_series`), `stratify = TRUE` to balance outcome
proportions at the group level when applicable, and `nested = TRUE` to
attach inner folds (one repeat). For large datasets, `progress = TRUE`
reports progress; storing explicit indices can be memory intensive.

### Validating splits with `check_split_overlap()`

[`check_split_overlap()`](https://selcukorkmaz.github.io/bioLeak/reference/check_split_overlap.md)
verifies that no grouping-column levels appear in both the training and
test sets of any fold. It returns a data.frame with one row per
fold-column combination showing the overlap count and a pass/fail flag.

``` r

# Run on the subject-grouped splits from earlier
overlap_safe <- check_split_overlap(safe_splits)
overlap_safe
#>   fold repeat_id     col n_overlap pass
#> 1    1         1 subject         0 TRUE
#> 2    2         1 subject         0 TRUE
#> 3    3         1 subject         0 TRUE
#> 4    4         1 subject         0 TRUE
#> 5    5         1 subject         0 TRUE
```

When `cols` is not supplied, the function auto-infers the relevant
columns from the split mode (e.g., `subject` for `subject_grouped`, all
constraint axes for `combined`). By default `stop_on_fail = TRUE`, so
the function raises a `bioLeak_overlap_error` when any overlap is
detected. Set `stop_on_fail = FALSE` to return the results without
stopping.

When strict mode is active (`options(bioLeak.strict = TRUE)`),
[`make_split_plan()`](https://selcukorkmaz.github.io/bioLeak/reference/make_split_plan.md)
automatically runs
[`check_split_overlap()`](https://selcukorkmaz.github.io/bioLeak/reference/check_split_overlap.md)
after building non-compact splits, so you do not need to call it
manually.

### Scalability

**Handling Large Datasets (Compact Mode)**

For large datasets (e.g., $`N > 50,000`$) with many repeats, storing
explicit integer indices for every fold can consume gigabytes of memory.
Use `compact = TRUE` to store a lightweight vector of fold assignments
instead.
[`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md)
will automatically reconstruct the indices on the fly. Compact mode is
not supported when `nested = TRUE` (it falls back to full indices). For
`time_series` compact splits, the time column must be present in the
stored split metadata so folds can be reconstructed.

``` r

# Efficient storage for large N
big_splits <- make_split_plan(
  df,
  outcome = "outcome",
  mode = "subject_grouped",
  group = "subject",
  v = 5,
  compact = TRUE  # <--- Saves memory
)

cat("Compact-mode splits summary:\n")
#> Compact-mode splits summary:
big_splits
#> LeakSplits object (mode = subject_grouped, v = 5, repeats = 1)
#> Outcome: outcome | Stratified: FALSE | Nested: FALSE
#> ------------------------------------------------------
#>   fold repeat_id train_n test_n
#> 1    1         1     128     32
#> 2    2         1     128     32
#> 3    3         1     128     32
#> 4    4         1     128     32
#> 5    5         1     128     32
#> ------------------------------------------------------
#> Total folds: 5 | Hash: ee2866b7282066cda850ea9995f43c10
cat(sprintf("Compact storage enabled: %s\n", big_splits@info$compact))
#> Compact storage enabled: TRUE
```

The summary is identical to a regular split, but the underlying storage
is a compact fold-assignment vector. Use the `compact` flag to confirm
the memory- saving mode is active.

### Strict leakage mode

Setting `options(bioLeak.strict = TRUE)` activates a global safety
switch that tightens validation across the entire workflow. In strict
mode:

1.  **Trained recipe/workflow detection escalates to an error.**
    Normally, passing a pre-trained recipe or workflow to
    [`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md)
    or
    [`tune_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/tune_resample.md)
    produces a warning. In strict mode, it raises an error so
    fold-global leakage cannot proceed silently.

2.  **Seed warnings.** When `seed` is `NULL` or `NA` in
    [`make_split_plan()`](https://selcukorkmaz.github.io/bioLeak/reference/make_split_plan.md),
    [`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md),
    or
    [`tune_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/tune_resample.md),
    strict mode emits a `bioLeak_validation_warning` reminding you to
    set an explicit seed for reproducibility.

3.  **Combined + nested warning.** Using `nested = TRUE` with
    `mode = "combined"` may produce empty inner folds; strict mode warns
    about this before proceeding.

4.  **Automatic overlap check.** After building non-compact splits,
    [`make_split_plan()`](https://selcukorkmaz.github.io/bioLeak/reference/make_split_plan.md)
    auto-runs
    [`check_split_overlap()`](https://selcukorkmaz.github.io/bioLeak/reference/check_split_overlap.md)
    to verify that no grouping-column levels cross fold boundaries.

``` r

# Enable strict mode temporarily
withr::with_options(list(bioLeak.strict = TRUE), {
  strict_splits <- make_split_plan(
    df,
    outcome = "outcome",
    mode = "subject_grouped",
    group = "subject",
    v = 5,
    stratify = TRUE,
    seed = 42
  )
  cat("Strict mode splits completed without error.\n")
})
#> Strict mode splits completed without error.
```

``` r

# Strict mode catches a pre-trained recipe
if (requireNamespace("recipes", quietly = TRUE)) {
  rec <- recipes::recipe(outcome ~ ., data = df[, c("outcome", predictors)]) |>
    recipes::step_normalize(recipes::all_numeric_predictors()) |>
    recipes::prep(training = df[, c("outcome", predictors)])

  withr::with_options(list(bioLeak.strict = TRUE), {
    tryCatch(
      fit_resample(
        df, outcome = "outcome",
        splits = safe_splits,
        learner = "glmnet",
        preprocess = rec
      ),
      error = function(e) cat("Strict mode error:", conditionMessage(e), "\n")
    )
  })
}
#> Strict mode error: fit_resample: received a trained recipe. Pass an untrained recipe to avoid fold-global leakage.
```

Strict mode is off by default. Use
[`withr::with_options()`](https://withr.r-lib.org/reference/with_options.html)
for temporary activation or set `options(bioLeak.strict = TRUE)` in your
`.Rprofile` for project-wide enforcement.

### Guarded preprocessing and imputation

`bioLeak` uses guarded preprocessing to prevent leakage from global
imputation, scaling, and feature selection. The low-level API is:

- [`guard_fit()`](https://selcukorkmaz.github.io/bioLeak/reference/guard_fit.md)
  to fit preprocessing on training data only
- [`predict()`](https://rdrr.io/r/stats/predict.html) (dispatching to
  [`predict.GuardFit()`](https://selcukorkmaz.github.io/bioLeak/reference/predict.GuardFit.md))
  to apply the trained guard to new data; the legacy alias
  [`predict_guard()`](https://selcukorkmaz.github.io/bioLeak/reference/predict.GuardFit.md)
  is also kept for backward compatibility
- [`guard_ensure_levels()`](https://selcukorkmaz.github.io/bioLeak/reference/guard_ensure_levels.md)
  to align factor levels across train/test
- [`impute_guarded()`](https://selcukorkmaz.github.io/bioLeak/reference/impute_guarded.md)
  as a convenience wrapper for leakage-safe imputation

[`guard_fit()`](https://selcukorkmaz.github.io/bioLeak/reference/guard_fit.md)
fits a pipeline that can winsorize, impute, normalize, filter, and
select features. It one-hot encodes non-numeric columns and carries
factor levels forward to new data. Assumptions and misuse to avoid:

- `fs = "ttest"` requires a binary outcome with enough samples per class
- `fs = "lasso"` requires `glmnet`
- `impute = "mice"` is not supported for guarded prediction
- `impute = "none"` with missing values will trigger a median fallback
  and missingness indicators to avoid errors

Supported preprocessing options include imputation (`median`, `knn`,
`missForest`, `none`), normalization (`zscore`, `robust`, `none`),
filtering by variance or IQR (optionally `min_keep`), and feature
selection (`ttest`, `lasso`, `pca`). Winsorization is controlled via
`impute$winsor` and `impute$winsor_k` in guarded steps; in
[`impute_guarded()`](https://selcukorkmaz.github.io/bioLeak/reference/impute_guarded.md),
use `winsor` and `winsor_thresh`.
[`impute_guarded()`](https://selcukorkmaz.github.io/bioLeak/reference/impute_guarded.md)
returns a `LeakImpute` object with guarded data, diagnostics, and the
fitted guard state.

**Leaky example: global scaling before CV**

``` r

df_leaky_scaled <- df
df_leaky_scaled[predictors] <- scale(df_leaky_scaled[predictors])
scaled_summary <- data.frame(
  feature = predictors,
  mean = colMeans(df_leaky_scaled[predictors]),
  sd = apply(df_leaky_scaled[predictors], 2, stats::sd)
)
scaled_summary$mean <- round(scaled_summary$mean, 3)
scaled_summary$sd <- round(scaled_summary$sd, 3)

# Leaky global scaling: means ~0 and SDs ~1 computed on all samples
scaled_summary
#>    feature mean sd
#> x1      x1    0  1
#> x2      x2    0  1
#> x3      x3    0  1
```

The summary shows that scaling used the full dataset, so test-fold
statistics influenced the training transformation. This violates the
train/test separation and can bias performance estimates.

**Leakage-safe alternative: fit preprocessing on training only**

``` r

fold1 <- safe_splits@indices[[1]]
train_x <- df[fold1$train, predictors]
test_x <- df[fold1$test, predictors]

guard <- guard_fit(
  X = train_x,
  y = df$outcome[fold1$train],
  steps = list(
    impute = list(method = "median", winsor = TRUE),
    normalize = list(method = "zscore"),
    filter = list(var_thresh = 0, iqr_thresh = 0),
    fs = list(method = "none")
  ),
  task = "binomial"
)

train_x_guarded <- predict(guard, train_x)
test_x_guarded <- predict(guard, test_x)

cat("GuardFit object:\n")
#> GuardFit object:
guard
#> Guarded preprocessing pipeline
#>  - Imputation: median
#>  - Normalization: zscore
#>  - Filter: var>=0, iqr>=0
#>  - Feature selection: none
#>  - Output features: 3
cat("\nGuardFit summary:\n")
#> 
#> GuardFit summary:
summary(guard)
#> GuardFit summary
#>  n_train=124, p_out=3, removed_by_filter=0
#>  Steps:
#>  step                 action
#>     1            winsor: k=5
#>     2   impute: none (no NA)
#>     3      normalize: zscore
#>     4 filter: var>=0, iqr>=0
#>     5               fs: none
#>  hash: 9f7f4824b80da853f5fedaaac6ac205f

# Guarded training data (first 6 rows)
head(train_x_guarded)
#>            x1           x2          x3
#> 5  -1.4287254 -0.846746034 -0.58336326
#> 6   1.9282279  0.352168157 -0.70450573
#> 7   0.5110636  0.427913734 -0.08735148
#> 8  -1.3683145 -0.006339485  0.36484181
#> 9  -0.7187394 -2.383299658 -1.01654516
#> 10 -1.3015499  2.538309224 -0.14109035

# Guarded test data (first 6 rows)
head(test_x_guarded)
#>              x1         x2          x3
#> 1   0.119085787  1.0175596 -0.05210949
#> 2  -0.428021276 -1.9210219  0.26877352
#> 3  -0.002549215 -0.3915238  0.95073054
#> 4  -1.007141255  0.1398948  0.25797066
#> 25 -1.105722181  0.5243484  2.03908390
#> 26  1.601954989 -0.5499985  0.86726497
```

The `GuardFit` object records the preprocessing steps and the number of
features retained after filtering. The summary reports how many features
were removed and the preprocessing audit trail. The guarded train/test
previews show that missing values were imputed and (if requested) scaled
using training-only statistics; the test data never influences these
values.

**Factor level guard (advanced)**

``` r

raw_levels <- data.frame(
  site = c("A", "B", "B"),
  status = c("yes", "no", "yes"),
  stringsAsFactors = FALSE
)

level_state <- guard_ensure_levels(raw_levels)

# Aligned factor data with consistent levels
level_state$data
#>   site status
#> 1    A    yes
#> 2    B     no
#> 3    B    yes

# Levels map
level_state$levels
#> $site
#> [1] "A" "B"
#> 
#> $status
#> [1] "no"  "yes"
```

The returned `data` keeps factor levels consistent across folds, while
the `levels` list records the training-time levels (including any dummy
levels added to preserve one-hot columns). Use these when transforming
new data to avoid misaligned model matrices.

**Leaky example: imputation using train and test together**

``` r

train <- data.frame(a = c(1, 2, NA, 4), b = c(NA, 1, 1, 0))
test <- data.frame(a = c(NA, 5), b = c(1, NA))

all_median <- vapply(rbind(train, test),
                     function(col) median(col, na.rm = TRUE),
                     numeric(1))
train_leaky <- as.data.frame(Map(function(col, m) { col[is.na(col)] <- m; col },
                                 train, all_median))
test_leaky <- as.data.frame(Map(function(col, m) { col[is.na(col)] <- m; col },
                                test, all_median))

# Leaky medians computed on train + test
data.frame(feature = names(all_median), median = all_median)
#>   feature median
#> a       a      3
#> b       b      1

# Leaky-imputed training data
train_leaky
#>   a b
#> 1 1 1
#> 2 2 1
#> 3 3 1
#> 4 4 0

# Leaky-imputed test data
test_leaky
#>   a b
#> 1 3 1
#> 2 5 1
```

The medians above were computed using both train and test data, so the
test set influences the imputation values. This is a classic leakage
pathway because test information directly alters the training features.

**Leakage-safe alternative:
[`impute_guarded()`](https://selcukorkmaz.github.io/bioLeak/reference/impute_guarded.md)**

``` r

imp <- impute_guarded(
  train = train,
  test = test,
  method = "median",
  winsor = FALSE
)

# Guarded-imputed training data
imp$train
#>   a b
#> 1 1 1
#> 2 2 1
#> 3 2 1
#> 4 4 0

# Guarded-imputed test data
imp$test
#>   a b
#> 1 2 1
#> 2 5 1
```

Here, the imputation statistics are computed from the training set only.
The missing value in the test set (column `a`) is replaced by 2, which
is the training median. In contrast, in the leaky example above, it was
replaced by 3, the global median.

This confirms that the test set is transformed using fixed values
learned from the training data, preserving a clean separation between
training and evaluation. The `LeakImpute` object also contains
missingness diagnostics in `imp$summary$diagnostics`, and guarded
outputs use the same encoding as
[`guard_fit()`](https://selcukorkmaz.github.io/bioLeak/reference/guard_fit.md)
(categorical variables are one-hot encoded). Use `vars` to impute only a
subset of columns when needed.

**Bridging guard steps to recipes:
[`guard_to_recipe()`](https://selcukorkmaz.github.io/bioLeak/reference/guard_to_recipe.md)**

[`guard_to_recipe()`](https://selcukorkmaz.github.io/bioLeak/reference/guard_to_recipe.md)
converts a guard preprocessing specification into a
[`recipes::recipe`](https://recipes.tidymodels.org/reference/recipe.html),
enabling a smooth transition from guard-based to tidymodels-based
workflows. It maps supported guard steps to their closest recipe
equivalents:

| Guard step | Recipe equivalent |
|----|----|
| `impute$method = "median"` | [`step_impute_median()`](https://recipes.tidymodels.org/reference/step_impute_median.html) |
| `impute$method = "knn"` | [`step_impute_knn()`](https://recipes.tidymodels.org/reference/step_impute_knn.html) |
| `normalize$method = "zscore"` | [`step_normalize()`](https://recipes.tidymodels.org/reference/step_normalize.html) |
| `filter$var_thresh > 0` | [`step_nzv()`](https://recipes.tidymodels.org/reference/step_nzv.html) |
| `fs$method = "pca"` | [`step_pca()`](https://recipes.tidymodels.org/reference/step_pca.html) |

Steps with no direct recipe equivalent (`missForest`, `mice`, `robust`
normalization, `ttest`, `lasso` feature selection) produce a
`bioLeak_fallback_warning` and either fall back to the closest available
step or are skipped.

``` r

if (requireNamespace("recipes", quietly = TRUE)) {
  guard_steps <- list(
    impute    = list(method = "median"),
    normalize = list(method = "zscore")
  )

  rec <- guard_to_recipe(
    steps = guard_steps,
    formula = outcome ~ .,
    training_data = df[, c("outcome", predictors)]
  )

  cat("Converted recipe:\n")
  rec
}
#> Converted recipe:
```

The returned recipe is untrained and can be passed directly to
[`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md)
or
[`tune_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/tune_resample.md)
as the `preprocess` argument.

### Fit and resample with `fit_resample()`

[`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md)
combines leakage-aware splits with guarded preprocessing and model
fitting. It supports:

- built-in learners (`glmnet`, `ranger`)
- `parsnip` model specifications or
  [`workflows::workflow`](https://workflows.tidymodels.org/reference/workflow.html)
  objects (recommended)
- [`recipes::recipe`](https://recipes.tidymodels.org/reference/recipe.html)
  preprocessing (prepped on training folds only)
- `rsample` resamples (`rset`/`rsplit`) via `splits`, with `split_cols`
  to map metadata
- custom learners via `custom_learners` (legacy/advanced)
- multiple metrics, class weights, positive class override, and
  yardstick metrics
- fold-level run diagnostics in `fit@info$fold_status`
  (`success`/`skipped`/`failed`)

Assumptions and misuse to avoid:

- outcome must be binary or multiclass (factor) for classification,
  numeric for regression, or a survival outcome (a `Surv` column or
  time/event columns)
- `positive_class` must be a single value that exists in the outcome
  levels
- `class_weights` only applies to binomial/multiclass tasks; workflows
  must handle weights internally
- if you supply a matrix without metadata, you must remove leakage
  columns yourself (group, batch, study, time)

Use `learner_args` to pass model-specific arguments. For custom learners
(advanced), you can supply separate `fit` and `predict` argument lists.
Set `parallel = TRUE` to use future.apply for fold-level parallelism
when available. When `learner` is a parsnip specification,
`learner_args` and `custom_learners` are ignored. When `learner` is a
workflow, `preprocess` and `learner_args` are ignored. When a recipe or
workflow is used, the built-in guarded preprocessing list is bypassed,
so ensure the recipe/workflow itself is leakage-safe.

Metrics used in this vignette:

- **AUC**: area under the ROC curve (0.5 = random, 1.0 = perfect; higher
  is better)
- **PR AUC**: area under the precision-recall curve (0 to 1; higher is
  better)
- **Accuracy**: fraction of correct predictions (0 to 1; higher is
  better)
- **RMSE**: root mean squared error (0 to infinity; lower is better)
- **Macro F1**: average F1 across classes (multiclass; higher is better)
- **Log loss**: cross-entropy for probabilities (lower is better)
- **C-index**: concordance for regression/survival (higher is better)

Always report the mean and variability across folds, not a single fold
value. When using yardstick metrics, the positive class is the second
factor level; set `positive_class` to control this.

**Parsnip model specification (recommended)**

``` r

spec <- parsnip::logistic_reg(mode = "classification") |>
  parsnip::set_engine("glm")
```

This spec uses base R `glm` under the hood and does not require extra
model packages. Use it in the examples below; custom learners are
covered in the Advanced section.

**Leaky example: leaky features and row-wise splits**

``` r

fit_leaky <- fit_resample(
  df_leaky,
  outcome = "outcome",
  splits = leaky_splits,
  learner = spec,
  metrics = c("auc", "accuracy"),
  preprocess = list(
    impute = list(method = "median"),
    normalize = list(method = "zscore"),
    filter = list(var_thresh = 0),
    fs = list(method = "none")
  )
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%

cat("Leaky fit summary:\n")
#> Leaky fit summary:
summary(fit_leaky)
#> 
#> ===========================
#>  bioLeak Model Fit Summary
#> ===========================
#> 
#> Task: binomial
#> Outcome: outcome
#> Positive class: case
#> Learners: logistic_reg/glm
#> Total folds: 5
#> Fold status: 5 success, 0 skipped, 0 failed
#> Refit performed: Yes
#> Hash: ha16bdf5a
#> 
#> Cross-validated metrics (mean ± SD):
#>            learner auc_mean auc_sd accuracy_mean accuracy_sd auc_ci_lo
#> 1 logistic_reg/glm    0.693  0.097         0.618        0.06     0.513
#>   auc_ci_hi accuracy_ci_lo accuracy_ci_hi
#> 1     0.873          0.506          0.729
#> 
#> Audit overview:
#>  fold n_train n_test          learner features_final
#>     1     127     33 logistic_reg/glm             17
#>     2     128     32 logistic_reg/glm             17
#>     3     128     32 logistic_reg/glm             17
#>     4     128     32 logistic_reg/glm             17
#>     5     129     31 logistic_reg/glm             17
metrics_leaky <- as.data.frame(fit_leaky@metric_summary)
num_cols <- vapply(metrics_leaky, is.numeric, logical(1))
metrics_leaky[num_cols] <- lapply(metrics_leaky[num_cols], round, digits = 3)

# Leaky fit: mean and SD of metrics across folds
metrics_leaky
#>            learner auc_mean auc_sd accuracy_mean accuracy_sd auc_ci_lo
#> 1 logistic_reg/glm    0.693  0.097         0.618        0.06     0.513
#>   auc_ci_hi accuracy_ci_lo accuracy_ci_hi
#> 1     0.873          0.506          0.729
```

The summary reports cross-validated performance. AUC ranges from 0.5
(random) to 1.0 (perfect), while accuracy is the fraction of correct
predictions. Because the splits are leaky, these metrics can be
artificially high and should not be used for reporting.

**Leakage-safe alternative: grouped splits and clean predictors**

``` r

fit_safe <- fit_resample(
  df,
  outcome = "outcome",
  splits = safe_splits,
  learner = spec,
  metrics = c("auc", "accuracy"),
  preprocess = list(
    impute = list(method = "median"),
    normalize = list(method = "zscore"),
    filter = list(var_thresh = 0),
    fs = list(method = "none")
  ),
  positive_class = "case",
  class_weights = c(control = 1, case = 1),
  refit = TRUE
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%

cat("Leakage-safe fit summary:\n")
#> Leakage-safe fit summary:
summary(fit_safe)
#> 
#> ===========================
#>  bioLeak Model Fit Summary
#> ===========================
#> 
#> Task: binomial
#> Outcome: outcome
#> Positive class: case
#> Learners: logistic_reg/glm
#> Total folds: 5
#> Fold status: 5 success, 0 skipped, 0 failed
#> Refit performed: Yes
#> Hash: hba29ffa5
#> 
#> Cross-validated metrics (mean ± SD):
#>            learner auc_mean auc_sd accuracy_mean accuracy_sd auc_ci_lo
#> 1 logistic_reg/glm    0.612  0.184         0.598       0.117      0.27
#>   auc_ci_hi accuracy_ci_lo accuracy_ci_hi
#> 1     0.954          0.379          0.816
#> 
#> Audit overview:
#>  fold n_train n_test          learner features_final
#>     1     124     36 logistic_reg/glm             13
#>     2     128     32 logistic_reg/glm             13
#>     3     128     32 logistic_reg/glm             13
#>     4     128     32 logistic_reg/glm             13
#>     5     132     28 logistic_reg/glm             13
metrics_safe <- as.data.frame(fit_safe@metric_summary)
num_cols <- vapply(metrics_safe, is.numeric, logical(1))
metrics_safe[num_cols] <- lapply(metrics_safe[num_cols], round, digits = 3)

# Leakage-safe fit: mean and SD of metrics across folds
metrics_safe
#>            learner auc_mean auc_sd accuracy_mean accuracy_sd auc_ci_lo
#> 1 logistic_reg/glm    0.612  0.184         0.598       0.117      0.27
#>   auc_ci_hi accuracy_ci_lo accuracy_ci_hi
#> 1     0.954          0.379          0.816

# Per-fold metrics (first 6 rows)
head(fit_metrics(fit_safe))
#>   fold          learner       auc  accuracy
#> 1    1 logistic_reg/glm 0.6284830 0.5277778
#> 2    2 logistic_reg/glm 0.7611336 0.7187500
#> 3    3 logistic_reg/glm 0.6823529 0.6562500
#> 4    4 logistic_reg/glm 0.6944444 0.6562500
#> 5    5 logistic_reg/glm 0.2944444 0.4285714
```

[`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md)
returns a `LeakFit` object. Use `summary(fit_safe)` for a compact
report, `fit_metrics(fit_safe)` for the per-fold metric data frame, and
inspect `fit_safe@predictions` for fold-level predictions. When
`refit = TRUE`, the final model and preprocessing state are stored in
`fit_safe@info$final_model` and `fit_safe@info$final_preprocess`.
Fold-level status diagnostics are available in
`fit_safe@info$fold_status`.

The mean +/- SD table is the primary performance summary to report,
while the per-fold metrics reveal variability and potential instability
across folds.

By default,
[`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md)
stores refit inputs in `fit_safe@info$perm_refit_spec` to enable
`audit_leakage(perm_refit = "auto")`. Set `store_refit_data = FALSE` to
reduce memory and provide `perm_refit_spec` manually when you want
refit-based permutations. When multiple learners are passed,
`refit = TRUE` refits only the first learner; set `refit = FALSE` if you
do not need a final model.

``` r

# Fold-level diagnostics for reproducible troubleshooting
head(fit_safe@info$fold_status, 5)
#>   fold    stage  status reason                    notes elapsed_sec
#> 1    1 fold_run success   <NA> Successful learners: 1/1       0.025
#> 2    2 fold_run success   <NA> Successful learners: 1/1       0.023
#> 3    3 fold_run success   <NA> Successful learners: 1/1       0.022
#> 4    4 fold_run success   <NA> Successful learners: 1/1       0.022
#> 5    5 fold_run success   <NA> Successful learners: 1/1       0.022
```

Interpretation of `fold_status`:

- `status = "success"` means at least one learner produced valid metrics
  on that fold.
- `status = "skipped"` typically means a data-structure issue in that
  fold (for example, single-class training data in classification).
- `status = "failed"` indicates a runtime failure (for example, model
  fitting or preprocessing error).

Use `reason` and `notes` to decide whether to:

- adjust split design (for example, stronger grouping, stratification,
  or fewer folds),
- adjust preprocessing/model complexity,
- or exclude unstable configurations before formal comparison.

If many folds are skipped/failed, do **not** trust aggregate metric
means until the instability is resolved.

The examples above use a parsnip model specification. To swap in other
models, replace `spec` with another parsnip spec (see the gradient
boosting example below).

**Multiclass classification (optional)**

``` r

if (requireNamespace("ranger", quietly = TRUE)) {
  set.seed(11)
  df_multi <- df
  df_multi$outcome3 <- factor(sample(c("A", "B", "C"),
                                     nrow(df_multi), replace = TRUE))

  multi_splits <- make_split_plan(
    df_multi,
    outcome = "outcome3",
    mode = "subject_grouped",
    group = "subject",
    v = 4,
    stratify = TRUE,
    seed = 11
  )

  fit_multi <- fit_resample(
    df_multi,
    outcome = "outcome3",
    splits = multi_splits,
    learner = "ranger",
    metrics = c("accuracy", "macro_f1", "log_loss"),
    refit = FALSE
  )

  cat("Multiclass fit summary:\n")
  summary(fit_multi)
} else {
  cat("ranger not installed; skipping multiclass example.\n")
}
#>   |                                                                              |                                                                      |   0%  |                                                                              |==================                                                    |  25%  |                                                                              |===================================                                   |  50%  |                                                                              |====================================================                  |  75%  |                                                                              |======================================================================| 100%
#> Multiclass fit summary:
#> 
#> ===========================
#>  bioLeak Model Fit Summary
#> ===========================
#> 
#> Task: multiclass
#> Outcome: outcome3
#> Learners: ranger
#> Total folds: 4
#> Fold status: 4 success, 0 skipped, 0 failed
#> Refit performed: No
#> Hash: hc88c6471
#> 
#> Cross-validated metrics (mean ± SD):
#>   learner accuracy_mean accuracy_sd macro_f1_mean macro_f1_sd log_loss_mean
#> 1  ranger         0.319       0.072         0.313       0.071         1.154
#>   log_loss_sd accuracy_ci_lo accuracy_ci_hi macro_f1_ci_lo macro_f1_ci_hi
#> 1       0.051          0.144          0.494           0.14          0.486
#>   log_loss_ci_lo log_loss_ci_hi
#> 1          1.029          1.278
#> 
#> Audit overview:
#>  fold n_train n_test learner features_final
#>     1     116     44  ranger             14
#>     2     120     40  ranger             14
#>     3     120     40  ranger             14
#>     4     124     36  ranger             14
```

Multiclass fits compute accuracy, macro-F1, and log loss when class
probabilities are available. `positive_class` is ignored for multiclass
tasks.

**Survival outcomes (optional)**

``` r

if (requireNamespace("survival", quietly = TRUE) &&
    requireNamespace("glmnet", quietly = TRUE)) {
  set.seed(12)
  df_surv <- df
  df_surv$time_to_event <- rexp(nrow(df_surv), rate = 0.1)
  df_surv$event <- rbinom(nrow(df_surv), 1, 0.7)

  surv_splits <- make_split_plan(
    df_surv,
    outcome = c("time_to_event", "event"),
    mode = "subject_grouped",
    group = "subject",
    v = 4,
    stratify = FALSE,
    seed = 12
  )

  fit_surv <- fit_resample(
    df_surv,
    outcome = c("time_to_event", "event"),
    splits = surv_splits,
    learner = "glmnet",
    metrics = "cindex",
    refit = FALSE
  )

  cat("Survival fit summary:\n")
  summary(fit_surv)
} else {
  cat("survival or glmnet not installed; skipping survival example.\n")
}
#>   |                                                                              |                                                                      |   0%
#>   |                                                                              |==================                                                    |  25%
#>   |                                                                              |===================================                                   |  50%
#>   |                                                                              |====================================================                  |  75%
#>   |                                                                              |======================================================================| 100%
#> Survival fit summary:
#> 
#> ===========================
#>  bioLeak Model Fit Summary
#> ===========================
#> 
#> Task: survival
#> Outcome: time_to_event
#>  Outcome: event
#> Learners: glmnet
#> Total folds: 4
#> Fold status: 4 success, 0 skipped, 0 failed
#> Refit performed: No
#> Hash: h041adef3
#> 
#> Cross-validated metrics (mean ± SD):
#>   learner cindex_mean cindex_sd cindex_ci_lo cindex_ci_hi
#> 1  glmnet         0.5         0          0.5          0.5
#> 
#> Audit overview:
#>  fold n_train n_test learner features_final
#>     1     120     40  glmnet             14
#>     2     120     40  glmnet             14
#>     3     120     40  glmnet             14
#>     4     120     40  glmnet             14
```

For survival tasks, supply `outcome = c(time_col, event_col)` or a
`Surv` column; `stratify` is ignored for time/event outcomes and
`class_weights` are not used.

**Passing learner-specific arguments (optional)**

``` r

if (requireNamespace("glmnet", quietly = TRUE)) {
  fit_glmnet <- fit_resample(
    df,
    outcome = "outcome",
    splits = safe_splits,
    learner = "glmnet",
    metrics = "auc",
    learner_args = list(glmnet = list(alpha = 0.5)),
    preprocess = list(
      impute = list(method = "median"),
      normalize = list(method = "zscore"),
      filter = list(var_thresh = 0),
      fs = list(method = "none")
    )
  )
  cat("GLMNET summary with learner-specific arguments:\n")
  summary(fit_glmnet)
} else {
  cat("glmnet not installed; skipping learner_args example.\n")
}
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#> GLMNET summary with learner-specific arguments:
#> 
#> ===========================
#>  bioLeak Model Fit Summary
#> ===========================
#> 
#> Task: binomial
#> Outcome: outcome
#> Positive class: case
#> Learners: glmnet
#> Total folds: 5
#> Fold status: 5 success, 0 skipped, 0 failed
#> Refit performed: Yes
#> Hash: hba29ffa5
#> 
#> Cross-validated metrics (mean ± SD):
#>   learner auc_mean auc_sd auc_ci_lo auc_ci_hi
#> 1  glmnet    0.653  0.199     0.282     1.024
#> 
#> Audit overview:
#>  fold n_train n_test learner features_final
#>     1     124     36  glmnet             13
#>     2     128     32  glmnet             13
#>     3     128     32  glmnet             13
#>     4     128     32  glmnet             13
#>     5     132     28  glmnet             13
```

This summary reflects the same guarded preprocessing but a different
model configuration (here, `alpha = 0.5`). Use the mean +/- SD metrics
to compare learner settings under identical splits.

**SummarizedExperiment inputs (optional)**

``` r

if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = t(as.matrix(df[, predictors]))),
    colData = df[, c("subject", "batch", "study", "time", "outcome"), drop = FALSE]
  )

  se_splits <- make_split_plan(
    se,
    outcome = "outcome",
    mode = "subject_grouped",
    group = "subject",
    v = 3
  )

  se_fit <- fit_resample(
    se,
    outcome = "outcome",
    splits = se_splits,
    learner = spec,
    metrics = "auc"
  )
  cat("SummarizedExperiment fit summary:\n")
  summary(se_fit)
} else {
  cat("SummarizedExperiment not installed; skipping SE example.\n")
}
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#> SummarizedExperiment fit summary:
#> 
#> ===========================
#>  bioLeak Model Fit Summary
#> ===========================
#> 
#> Task: binomial
#> Outcome: outcome
#> Positive class: case
#> Learners: logistic_reg/glm
#> Total folds: 3
#> Fold status: 3 success, 0 skipped, 0 failed
#> Refit performed: Yes
#> Hash: h3f97c442
#> 
#> Cross-validated metrics (mean ± SD):
#>            learner auc_mean auc_sd auc_ci_lo auc_ci_hi
#> 1 logistic_reg/glm    0.673  0.062     0.431     0.914
#> 
#> Audit overview:
#>  fold n_train n_test          learner features_final
#>     1     104     56 logistic_reg/glm              3
#>     2     108     52 logistic_reg/glm              3
#>     3     108     52 logistic_reg/glm              3
```

The summary is identical in structure to the `data.frame` case because
[`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md)
extracts predictors and metadata from the `SummarizedExperiment` object.

Note that `features_final` here reflects only the assay predictors
(`x1`, `x2`, `x3`), because the assay was constructed without metadata
columns.

### Tidymodels interoperability

`bioLeak` integrates with `rsample`, `recipes`, `workflows`, and
`yardstick`:

- `splits` can be an `rsample` rset/rsplit; use `split_cols` (default
  `"auto"`) or `bioLeak_mode` / `bioLeak_perm_mode` attributes to map
  group/batch/study/time metadata.
- `preprocess` can be a
  [`recipes::recipe`](https://recipes.tidymodels.org/reference/recipe.html)
  (prepped on training folds only).
- `learner` can be a
  [`workflows::workflow`](https://workflows.tidymodels.org/reference/workflow.html).
- `metrics` can be a
  [`yardstick::metric_set`](https://yardstick.tidymodels.org/reference/metric_set.html).
- [`as_rsample()`](https://selcukorkmaz.github.io/bioLeak/reference/as_rsample.md)
  converts `LeakSplits` to an `rsample` object.

``` r

## This chunk requires the recipes and yardstick packages, which are
## listed in DESCRIPTION's Suggests rather than Imports. The eval=
## option above gates the chunk on their availability, so the vignette
## still builds when those packages are absent. Users copying this
## code into their own R session must have both packages installed
## (e.g., install.packages(c("recipes", "yardstick"))).
library(bioLeak)
library(parsnip)
library(recipes)
library(yardstick)

set.seed(123)
N <- 60
df <- data.frame(
  subject = factor(rep(paste0("S", 1:20), length.out = N)), 
  outcome = factor(rep(c("ClassA", "ClassB"), length.out = N)),
  x1 = rnorm(N),
  x2 = rnorm(N),
  x3 = rnorm(N)
)

spec <- logistic_reg() |> set_engine("glm")

# Use bioLeak's native split planner to avoid conversion errors
set.seed(13)

# Use make_split_plan instead of rsample::group_vfold_cv
# This creates a subject-grouped CV directly compatible with fit_resample
splits <- make_split_plan(
  df, 
  outcome = "outcome", 
  mode = "subject_grouped", 
  group = "subject", 
  v = 3
)

rec <- recipes::recipe(outcome ~ x1 + x2 + x3, data = df) |>
  recipes::step_impute_median(recipes::all_numeric_predictors()) |>
  recipes::step_normalize(recipes::all_numeric_predictors())

metrics_set <- yardstick::metric_set(yardstick::roc_auc, yardstick::accuracy)

fit_rs <- fit_resample(
  df,
  outcome = "outcome",
  splits = splits,
  learner = spec,
  preprocess = rec,
  metrics = metrics_set,
  refit = FALSE
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%

if (exists("as_rsample", where = asNamespace("bioLeak"), mode = "function")) {
    rs_export <- as_rsample(fit_rs@splits, data = df)
    print(rs_export)
}
#> # Manual resampling 
#> # A tibble: 3 × 2
#>   splits          id   
#>   <list>          <chr>
#> 1 <split [39/21]> Fold1
#> 2 <split [39/21]> Fold2
#> 3 <split [42/18]> Fold3
```

Use `split_cols` to ensure split-defining metadata are dropped from
predictors when using rsample inputs.

``` r

if (requireNamespace("workflows", quietly = TRUE)) {
  wf <- workflows::workflow() |>
    workflows::add_model(spec) |>
    workflows::add_formula(outcome ~ x1 + x2 + x3)

  fit_wf <- fit_resample(
    df,
    outcome = "outcome",
    splits = splits,
    learner = wf,
    metrics = "auc",
    refit = FALSE
  )

  cat("Workflow fit summary:\n")
  summary(fit_wf)
} else {
  cat("workflows not installed; skipping workflow example.\n")
}
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#> Workflow fit summary:
#> 
#> ===========================
#>  bioLeak Model Fit Summary
#> ===========================
#> 
#> Task: binomial
#> Outcome: outcome
#> Positive class: ClassB
#> Learners: workflow_1
#> Total folds: 3
#> Fold status: 3 success, 0 skipped, 0 failed
#> Refit performed: No
#> Hash: h825b7ff0
#> 
#> Cross-validated metrics (mean ± SD):
#>      learner auc_mean auc_sd auc_ci_lo auc_ci_hi
#> 1 workflow_1    0.366  0.113    -0.078      0.81
#> 
#> Audit overview:
#>  fold n_train n_test    learner features_final
#>     1      39     21 workflow_1              3
#>     2      39     21 workflow_1              3
#>     3      42     18 workflow_1              3
```

When auditing rsample-based fits, pass `perm_mode = "subject_grouped"`
(or set `attr(rs, "bioLeak_perm_mode") <- "subject_grouped"`) so
restricted permutations respect the intended split design.

### Nested tuning with `tune_resample()`

[`tune_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/tune_resample.md)
runs leakage-aware nested tuning with tidymodels. It accepts a parsnip
model specification or workflow, and supports two inner-CV modes:

- precomputed inner folds (for example from
  `make_split_plan(..., nested = TRUE)`)
- on-the-fly inner split generation for `LeakSplits` using `inner_v`,
  `inner_repeats`, and `inner_seed`

For rsample inputs, inner folds must already be present. Survival tasks
are not yet supported.

``` r

if (requireNamespace("parsnip", quietly = TRUE) &&
    requireNamespace("recipes", quietly = TRUE) &&
    requireNamespace("tune", quietly = TRUE) &&
    requireNamespace("glmnet", quietly = TRUE)) {
  # --- 1. Create Data ---
  set.seed(123)
  N <- 60
  df <- data.frame(
    subject = factor(rep(paste0("S", 1:20), length.out = N)),
    outcome = factor(sample(c("ClassA", "ClassB"), N, replace = TRUE)),
    x1 = rnorm(N),
    x2 = rnorm(N),
    x3 = rnorm(N)
  )

  # --- 2. Generate Nested Splits ---
  set.seed(1)
  nested_splits <- make_split_plan(
    df,
    outcome = "outcome",
    mode = "subject_grouped",
    group = "subject",
    v = 3,
    nested = TRUE
  )

  # --- 3. Define Recipe & Model ---
  rec <- recipes::recipe(outcome ~ x1 + x2 + x3, data = df) |>
    recipes::step_impute_median(recipes::all_numeric_predictors()) |>
    recipes::step_normalize(recipes::all_numeric_predictors())

  spec_tune <- parsnip::logistic_reg(
    penalty = tune::tune(),
    mixture = 1,
    mode = "classification"
  ) |>
    parsnip::set_engine("glmnet")

  # --- 4. Run Tuning ---
  tuned <- tune_resample(
    df,
    outcome = "outcome",
    splits = nested_splits,
    learner = spec_tune,
    preprocess = rec,
    inner_v = 2,
    grid = 3,
    metrics = c("auc", "accuracy"),
    selection = "one_std_err",
    refit = TRUE,
    seed = 14
  )

  summary(tuned)
} else {
  cat("parsnip/recipes/tune/glmnet not installed; skipping nested tuning example.\n")
}
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
#> 
#> ============================
#>  bioLeak Tuning Summary
#> ============================
#> 
#> Task: binomial
#> Outcome: outcome
#> Positive class: ClassB
#> Tuning Grid: 3
#> Selection Rule: one_std_err (Metric: roc_auc)
#> Fold status: 3 success, 0 skipped, 0 failed
#> Refit performed: Yes
#> Outer Folds: 3 successful / 3 total
#> 
#> Outer Loop Metrics (mean ± SD):
#>               learner accuracy_mean accuracy_sd roc_auc_mean roc_auc_sd
#> 1 logistic_reg/glmnet         0.619       0.048        0.577      0.134
#>   accuracy_ci_lo accuracy_ci_hi roc_auc_ci_lo roc_auc_ci_hi
#> 1          0.501          0.737         0.245         0.909
#> 
#> Best Parameters (First 5 Folds):
#>  penalty fold
#>    0.763    1
#>    0.000    2
#>    0.995    3
```

``` r

if (exists("tuned")) {
  # Fold-level status from the outer loop
  head(tuned$fold_status, 5)

  # Final tuned model is available when refit = TRUE
  is.null(tuned$final_model)
} else {
  cat("Nested tuning object not available (dependencies missing).\n")
}
#> [1] FALSE
```

How to interpret this tuning output:

- `selection = "one_std_err"` chooses a model within one standard error
  of the best inner-CV score, then applies a deterministic simplicity
  tie-break. This usually favors more regularized/simpler models with
  similar expected performance.
- `tuned$metrics` and `tuned$metric_summary` remain the primary unbiased
  performance estimates (outer folds only).
- `tuned$final_model` is the deployable model fit on full data using the
  selected hyperparameters; it is for downstream prediction, not for
  estimating generalization.
- `tuned$fold_status` should be reviewed before reporting results. A
  large number of skipped/failed folds indicates unstable tuning and
  should trigger split/parameter redesign.

### Advanced: Using Gradient Boosting with Parsnip

`bioLeak` natively supports `tidymodels` specifications. You can pass a
`parsnip` model specification directly to
[`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md).
This allows you to use state-of-the-art algorithms like XGBoost,
LightGBM, or SVMs while ensuring all preprocessing remains guarded.

``` r

if (requireNamespace("parsnip", quietly = TRUE) &&
    requireNamespace("xgboost", quietly = TRUE) &&
    requireNamespace("recipes", quietly = TRUE)) {
  
  # 1. Define the model spec
  xgb_spec <- parsnip::boost_tree(
    mode = "classification",
    trees = 100,
    tree_depth = 6,
    learn_rate = 0.01
  ) |>
    parsnip::set_engine("xgboost")
  
  # 2. Define a recipe on data without 'subject' because split metadata
  # columns are excluded from predictors by fit_resample().
  df_for_rec <- df[, !names(df) %in% "subject"]
  
  rec_xgb <- recipes::recipe(outcome ~ ., data = df_for_rec) |>
    recipes::step_dummy(recipes::all_nominal_predictors()) |>
    recipes::step_impute_median(recipes::all_numeric_predictors())
  # Note: No need for step_rm(subject) because it's already gone!
  
  # 3. Fit
  fit_xgb <- fit_resample(
    df,
    outcome = "outcome",
    splits = splits,
    learner = xgb_spec,
    metrics = "auc",
    preprocess = rec_xgb 
  )
  
  cat("XGBoost parsnip fit summary:\n")
  print(summary(fit_xgb))
} else {
  cat("parsnip/xgboost/recipes not installed.\n")
}
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#> XGBoost parsnip fit summary:
#> 
#> ===========================
#>  bioLeak Model Fit Summary
#> ===========================
#> 
#> Task: binomial
#> Outcome: outcome
#> Positive class: ClassB
#> Learners: boost_tree/xgboost
#> Total folds: 3
#> Fold status: 3 success, 0 skipped, 0 failed
#> Refit performed: Yes
#> Hash: h825b7ff0
#> 
#> Cross-validated metrics (mean ± SD):
#>              learner auc_mean auc_sd auc_ci_lo auc_ci_hi
#> 1 boost_tree/xgboost    0.616  0.058     0.388     0.844
#> 
#> Audit overview:
#>  fold n_train n_test            learner features_final
#>     1      39     21 boost_tree/xgboost              3
#>     2      39     21 boost_tree/xgboost              3
#>     3      42     18 boost_tree/xgboost              3
#> 
#>              learner  auc_mean     auc_sd auc_ci_lo auc_ci_hi
#> 1 boost_tree/xgboost 0.6157407 0.05800909 0.3878946 0.8435869
```

The summary reports cross-validated AUC for a non-linear gradient
boosting model. Use the mean $`\pm`$ SD table to compare against
baseline models, and confirm that any gains do not coincide with leakage
signals in the audit diagnostics.

### Advanced: Custom learners

Custom learners are used when a model is not available as a *parsnip*
specification or when a lightweight wrapper around base R is needed.

Each custom learner must define `fit` and `predict` components. The
`fit` function must accept `x`, `y`, `task`, and `weights` as inputs,
while the `predict` function must accept `object` and `newdata`.

``` r

custom_learners <- list(
  glm = list(
    fit = function(x, y, task, weights, ...) {
      df_fit <- data.frame(y = y, x, check.names = FALSE)
      stats::glm(y ~ ., data = df_fit,
                 family = stats::binomial(), weights = weights)
    },
    predict = function(object, newdata, task, ...) {
      as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
                                type = "response"))
    }
  )
)

cat("Custom learner names:\n")
#> Custom learner names:
names(custom_learners)
#> [1] "glm"
cat("Custom learner components (fit/predict):\n")
#> Custom learner components (fit/predict):
lapply(custom_learners, names)
#> $glm
#> [1] "fit"     "predict"
```

This output confirms that each custom learner provides both a `fit` and
a `predict` function.

Custom learners can be used with
[`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md)
as follows:

``` r

fit_resample(
  ...,
  learner = "glm",
  custom_learners = custom_learners
)
```

### Visual diagnostics

**bioLeak** includes plotting helpers for quick diagnostic checks:

- [`plot_fold_balance()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_fold_balance.md)
  shows class balance per fold.  
- [`plot_overlap_checks()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_overlap_checks.md)
  highlights train/test overlap for a metadata column.  
- [`plot_time_acf()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_time_acf.md)
  checks autocorrelation in time-series predictions.
- [`plot_calibration()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_calibration.md)
  shows reliability of binomial probabilities.  
- [`plot_confounder_sensitivity()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_confounder_sensitivity.md)
  summarizes performance by batch/study strata.

Tabular helpers are also available:

- [`calibration_summary()`](https://selcukorkmaz.github.io/bioLeak/reference/calibration_summary.md)
  returns calibration curve data and ECE/MCE/Brier metrics.
- [`confounder_sensitivity()`](https://selcukorkmaz.github.io/bioLeak/reference/confounder_sensitivity.md)
  returns per-stratum performance summaries.

**Misuse to avoid**

- Plotting without predictions (fits with no successful folds).  
- Using
  [`plot_overlap_checks()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_overlap_checks.md)
  when the column is not present in the split metadata.  
- Using
  [`plot_time_acf()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_time_acf.md)
  without a time column or for non-temporal data.
- Using
  [`plot_calibration()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_calibration.md)
  outside binomial tasks.
- Using
  [`plot_confounder_sensitivity()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_confounder_sensitivity.md)
  with a metric unsupported by the task.

For classification,
[`plot_fold_balance()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_fold_balance.md)
uses the positive class recorded in `fit@info$positive_class`, or
defaults to the second factor level if this is not set. For multiclass
tasks, it shows per-class counts without a proportion line.

``` r

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_fold_balance(fit_safe)
} else {
  cat("ggplot2 not installed; skipping fold balance plot.\n")
}
```

The bar chart shows the counts of Positives (blue) and Negatives (tan)
in each fold. The dashed blue line represents the proportion of the
positive class.

In a well-stratified fit, this line remains relatively stable across
folds. Large fluctuations indicate poor stratification, which can lead
to unstable fold-level performance estimates.

``` r

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_calibration(fit_safe, bins = 10)
} else {
  cat("ggplot2 not installed; skipping calibration plot.\n")
}
```

Calibration curves compare predicted probabilities to observed event
rates. Large deviations from the diagonal indicate miscalibration. Use
`min_bin_n` to suppress tiny bins and `learner =` to select a model when
multiple learners are stored in the fit.

``` r

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_confounder_sensitivity(fit_safe, confounders = c("batch", "study"),
                              metric = "auc", min_n = 3)
} else {
  cat("ggplot2 not installed; skipping confounder sensitivity plot.\n")
}
```

Confounder sensitivity plots highlight whether performance varies
substantially across batch or study strata.

``` r

cal <- calibration_summary(fit_safe, bins = 10, min_bin_n = 5)
conf_tbl <- confounder_sensitivity(fit_safe,
                                   confounders = c("batch", "study"),
                                   metric = "auc",
                                   min_n = 3)

# Calibration metrics
cal$metrics
#>     n bins min_bin_n       ece       mce     brier
#> 1 160   10         5 0.1005285 0.9113016 0.2499484

# Confounder sensitivity table (first 6 rows)
head(conf_tbl)
#>   confounder level metric direction     value  n positive_rate
#> 1      batch    B1    auc    higher 0.7569444 26     0.3076923
#> 2      batch    B2    auc    higher 0.7196970 23     0.4782609
#> 3      batch    B3    auc    higher 0.4153846 28     0.4642857
#> 4      batch    B4    auc    higher 0.6944444 24     0.5000000
#> 5      batch    B5    auc    higher 0.6071429 22     0.3636364
#> 6      batch    B6    auc    higher 0.5058824 37     0.4594595
```

``` r

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_overlap_checks(fit_leaky, column = "subject")
  plot_overlap_checks(fit_safe, column = "subject")
} else {
  cat("ggplot2 not installed; skipping overlap plots.\n")
}
```

The overlap plots compare the number of unique subjects appearing in the
training and test sets for each fold.

**Top plot (`fit_leaky`)**  
The red line shows high overlap counts, accompanied by the explicit
*“WARNING: Overlaps detected!”* annotation. This confirms that the same
subjects appear in both the training and test sets, indicating
information leakage.

**Bottom plot (`fit_safe`)**  
The overlap line remains flat at zero across all folds. This confirms
that `subject_grouped` splitting successfully keeps subjects isolated,
ensuring that no subject appears in both sets simultaneously.

### Audit leakage with `audit_leakage()`

[`audit_leakage()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md)
combines four diagnostics in one object:

- permutation gap: tests whether model signal is non-random
- batch or study association with folds (chi-square and Cramers V)
- target leakage scan: feature-wise outcome association to flag proxy
  columns
- duplicate and near-duplicate detection across train/test (by default)
  using feature similarity

Assumptions and misuse to avoid:

- `coldata` (or stored split metadata) must align with prediction IDs;
  rownames or a `row_id` column help alignment
- `X_ref` must align to prediction IDs (rownames, `row_id`, or matching
  order)
- [`plot_perm_distribution()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_perm_distribution.md)
  requires `return_perm = TRUE`
- `target_threshold` should be set high enough to flag only strong
  proxies
- target leakage scan includes univariate associations plus an optional
  multivariate/interaction check (`target_scan_multivariate`,
  `*_components`, `*_interactions`, `*_B`); proxies outside `X_ref` can
  still pass undetected
- for rsample splits, set `perm_mode` (or `bioLeak_perm_mode`) so
  restricted permutations respect the intended design
- for time-series audits, ensure a valid time column and set
  `time_block` / `block_len` as needed

**Interpretation Note:** The label-permutation test refits models by
default when refit inputs are available (`perm_refit = "auto"` with
`store_refit_data = TRUE`) and `B <= perm_refit_auto_max`. Otherwise it
keeps predictions fixed, so the p-value quantifies prediction-label
association rather than a full refit null. A large gap indicates a
non-random association. To determine if that signal is *real* or
*leaked*, check the **Batch Association** (confounding), **Target
Leakage Scan** (proxy features), and **Duplicate Detection**
(memorization) tables. Use `perm_refit = FALSE` to force
fixed-prediction shuffles or `perm_refit = TRUE` with `perm_refit_spec`
to always refit.

Use `feature_space = "rank"` to compare samples by rank profiles, and
`sim_method` (`cosine` or `pearson`) to control similarity. For large
`n`, `nn_k` and `max_pairs` limit duplicate searches. Use
`duplicate_scope = "all"` to include within-fold duplicates
(data-quality checks). Use `ci_method = "if"` (influence-function) or
`ci_method = "bootstrap"` with `boot_B` to obtain a confidence interval
for the permutation gap; `include_z` controls whether a z-score is
reported. Set `perm_stratify = "auto"` to stratify permutations only
when outcomes support it.

**Leakage-aware audit**

``` r

# Use df_leaky (the original 160-row data) so X_ref aligns with fit_safe's predictions.
# df may have been redefined to a smaller dataset in the tidymodels section above.
X_ref <- df_leaky[, predictors]
X_ref[c(1, 5), ] <- X_ref[1, ]

audit <- audit_leakage(
  fit_safe,
  metric = "auc",
  B = 20,
  perm_stratify = TRUE,
  batch_cols = c("batch", "study"),
  X_ref = X_ref,
  sim_method = "cosine",
  sim_threshold = 0.995,
  return_perm = TRUE
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%

cat("Leakage audit summary:\n")
#> Leakage audit summary:
summary(audit)
#> 
#> ==============================
#>  bioLeak Leakage Audit Summary
#> ==============================
#> 
#> Task: binomial | Outcome: outcome | Splitting mode: subject_grouped | Positive class: case
#> Hash: hba29ffa5 | Folds: 5 | Repeats: 1
#> 
#> Label-Permutation Association Test:
#>   Method: refit per permutation (auto)
#>   Null: re-drawn per permutation | Permutation: group_restricted | Summary: pooled
#>   Observed metric: 0.614
#>   Permuted mean ± SD: 0.515 ± 0.067
#>   Gap: 0.099 (larger gap = stronger non-random signal)
#>   This test does NOT diagnose information leakage. Use the Batch Association,
#>   Target Leakage Scan, and Duplicate Detection sections to check for leakage.
#> 
#> Batch / Study Association:
#>   batch (repeat 1): χ² = 17.043 (df = 20.000), p = 0.650
#>   study (repeat 1): χ² = 9.487 (df = 12.000), p = 0.661
#> 
#> Target Leakage Scan:
#>   Features checked: 3 | Flagged (score ≥ 0.900): 0
#>   No strong proxy features detected.
#> 
#> Multivariate Target Scan:
#>   Metric: auc | Score: 0.302 | p = 0.030
#>   Features: 3 | Components: 3 | Interactions: 3 | Permutations: 100
#> 
#> Near-Duplicate Samples:
#>   Scope: train/test only
#>   20 pairs detected above cosine ≥ 0.995
#>   Example pairs:
#>    mechanism_class   i   j       sim cross_fold   cos_sim
#>  duplicate_overlap   1   5 1.0000000       TRUE 1.0000000
#>  duplicate_overlap  35  61 0.9994623       TRUE 0.9994623
#>  duplicate_overlap  80 106 0.9990700       TRUE 0.9990700
#>  duplicate_overlap 120 153 0.9978776       TRUE 0.9978776
#>  duplicate_overlap   7 107 0.9976746       TRUE 0.9976746
#> 
#> Mechanism Risk Assessment:
#>        mechanism_class flagged        evidence statistic  p_value
#>      non_random_signal    TRUE permutation_gap 0.0988852 0.047619
#>  confounding_alignment   FALSE     batch_assoc 0.1631865 1.000000
#>   proxy_target_leakage   FALSE    target_assoc 0.3556299       NA
#>      duplicate_overlap    TRUE      duplicates 1.0000000       NA
#>     temporal_lookahead   FALSE      duplicates        NA       NA
#> 
#> Interpretation:
#>   ✓ Strong non-random signal.
pg <- audit_perm_gap(audit)
if (!is.null(pg) && nrow(pg) > 0) {
  # Permutation significance results
  pg
}
#>     mechanism_class metric_obs perm_mean   perm_sd       gap       z  p_value
#> 1 non_random_signal   0.614111  0.515225 0.0674816 0.0988852 2.16674 0.047619
#>   n_perm
#> 1     20
ba <- audit_batch_assoc(audit)
if (!is.null(ba) && nrow(ba) > 0) {
  # Batch/study association with folds (Cramer's V)
  ba
} else {
  cat("No batch or study associations detected.\n")
}
#>             mechanism_class batch_col repeat_id      stat df     pval  cramer_v
#> batch confounding_alignment     batch         1 17.043088 20 0.650174 0.1631865
#> study confounding_alignment     study         1  9.487022 12 0.660865 0.1405867
#>       pval_adj
#> batch        1
#> study        1
ta <- audit_target_assoc(audit)
if (!is.null(ta) && nrow(ta) > 0) {
  # Top features by target association score
  head(ta)
} else {
  cat("No target leakage scan results available.\n")
}
#>        mechanism_class feature    type metric     value      score      p_value
#> 1 proxy_target_leakage      x1 numeric    auc 0.6778149 0.35562988 0.0001205498
#> 2 proxy_target_leakage      x2 numeric    auc 0.4207676 0.15846472 0.0868346476
#> 3 proxy_target_leakage      x3 numeric    auc 0.5446727 0.08934544 0.3346967486
#>     n n_levels  flag p_value_adj flag_fdr
#> 1 160       NA FALSE          NA    FALSE
#> 2 160       NA FALSE          NA    FALSE
#> 3 160       NA FALSE          NA    FALSE
du <- audit_duplicates(audit)
if (!is.null(du) && nrow(du) > 0) {
  # Top duplicate/near-duplicate pairs by similarity
  head(du)
} else {
  cat("No near-duplicates detected.\n")
}
#>      mechanism_class   i   j       sim cross_fold   cos_sim
#> 1  duplicate_overlap   1   5 1.0000000       TRUE 1.0000000
#> 12 duplicate_overlap  35  61 0.9994623       TRUE 0.9994623
#> 18 duplicate_overlap  80 106 0.9990700       TRUE 0.9990700
#> 25 duplicate_overlap 120 153 0.9978776       TRUE 0.9978776
#> 4  duplicate_overlap   7 107 0.9976746       TRUE 0.9976746
#> 22 duplicate_overlap 104 142 0.9968751       TRUE 0.9968751
```

The permutation table reports the observed metric, the mean under random
label permutation, the gap (difference), and a permutation *p*-value.
For metrics where higher values indicate better performance, larger gaps
reflect stronger non-random signal.

**Which permutation null was used?**

With refit-based permutations, two design choices decide whether the
null is calibrated for grouped, stratified plans such as `safe_splits`.

1.  *Restricted permutation.* Labels should be permuted in a way that
    respects the grouping. For `subject_grouped` plans, whole outcome
    vectors are swapped between groups of equal size. This needs the
    outcome and the group column in `perm_refit_spec$coldata`, or in
    `perm_refit_spec$x` when it is a data.frame. If the group column is
    missing, the null falls back to an unrestricted shuffle, and
    [`audit_leakage()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md)
    raises a `bioLeak_permutation_warning`.
2.  *Fold construction.* `make_split_plan(stratify = TRUE)` balances
    each fold on the observed labels. If those folds are reused for
    permuted labels, the permuted class balance varies between folds.
    Each fold model’s baseline follows its training balance, which moves
    opposite to its test balance, so AUC pooled across folds is biased
    below 0.5 under the null. This is the stratification bias of Parker,
    Günter and Bedo (2007, *BMC Bioinformatics* 8:326). It inflates the
    gap and makes the p-value anti-conservative. In a pure-noise grouped
    design (40 families of six plus 24 singletons), the fixed-fold
    pooled null centred between about 0.40 and 0.46, and roughly
    two-thirds of the draws fell below 0.5.

By default (`perm_folds = "auto"`), a stratified plan is re-drawn with
[`make_split_plan()`](https://selcukorkmaz.github.io/bioLeak/reference/make_split_plan.md)
on every permuted outcome. The call reuses the plan’s mode, grouping,
`v`, `stratify` and constraints, with seed `splits@info$seed + b`.
Observed and permuted statistics then come from the same fold
construction. Use `perm_folds = "fixed"` to keep the observed folds, or
`perm_folds = "redraw"` to insist on re-drawing (this errors if the plan
cannot be rebuilt). `perm_summary = "fold_mean"` scores each draw by the
mean of per-fold AUCs instead of the pooled AUC. The fold mean never
compares predictions from different fold models, so between-fold
baseline shifts cannot bias it. Refit audits report both summaries in
`audit_info(audit)$perm_gap_summaries`.

``` r

info <- audit_info(audit)
# Null used, fold handling, permutation scheme and summary statistic
info[c("perm_null", "perm_folds", "perm_scheme", "perm_summary")]
#> $perm_null
#> [1] "refit_redrawn_folds"
#> 
#> $perm_folds
#> [1] "redrawn"
#> 
#> $perm_scheme
#> [1] "group_restricted"
#> 
#> $perm_summary
#> [1] "pooled"
info$perm_folds_reason
#> [1] "auto: stratified splits, so the split plan is re-drawn on each permuted outcome."
# Pooled and per-fold-mean permutation gaps (refit nulls only)
info$perm_gap_summaries
#>     summary metric_obs perm_mean    perm_sd        gap    p_value n_perm
#> 1    pooled  0.6141105 0.5152254 0.06748161 0.09888517 0.04761905     20
#> 2 fold_mean  0.6121717 0.5232547 0.06765256 0.08891695 0.09523810     20
```

The batch association table reports chi-square statistics and Cramer’s
V. Large p-values and small V values indicate that folds are not aligned
with batch or study labels (which is the desired outcome when these
should be independent).

The target scan table lists features with the strongest associations
with the outcome. For numeric features, the score is
$`|\mathrm{AUC} - 0.5| \times 2`$ for classification or the absolute
correlation for regression. For categorical features, the score is
Cramér’s V or eta-squared. Scores closer to 1 indicate stronger outcome
association.

The duplicate table lists pairs of samples with near-identical profiles
that cross train/test folds (by default). In this setup, the
artificially duplicated pair (`X_ref[c(1, 5), ] <- X_ref[1, ]`) should
appear near the top of the list. Use `duplicate_scope = "all"` to
include within-fold duplicates.

**Mechanism risk assessment**

`summary(audit)` now includes a **Mechanism Risk Assessment** section
that classifies leakage evidence into four mechanism classes:

| Mechanism class | Evidence slot | Flagging rule |
|----|----|----|
| `non_random_signal` | permutation gap | p ≤ 0.05 and gap \> 0 |
| `confounding_alignment` | batch association | p ≤ 0.05 and Cramér’s V ≥ 0.1 |
| `proxy_target_leakage` | target scan | any feature flagged (FDR or raw) |
| `duplicate_overlap` | duplicates | any cross-fold duplicates present |
| `temporal_lookahead` | duplicates | cross-fold duplicates in time-series splits |

The summary is stored in `audit@info$mechanism_summary` and can be
accessed directly for programmatic triage:

``` r

mech <- audit@info$mechanism_summary
if (is.data.frame(mech) && nrow(mech) > 0) {
  mech
} else {
  cat("No mechanism summary available.\n")
}
#>         mechanism_class flagged        evidence statistic  p_value
#> 1     non_random_signal    TRUE permutation_gap 0.0988852 0.047619
#> 2 confounding_alignment   FALSE     batch_assoc 0.1631865 1.000000
#> 3  proxy_target_leakage   FALSE    target_assoc 0.3556299       NA
#> 4     duplicate_overlap    TRUE      duplicates 1.0000000       NA
#> 5    temporal_lookahead   FALSE      duplicates        NA       NA
```

Each row reports whether the mechanism was flagged, the evidence type,
and the associated test statistic and p-value. This provides a quick
triage overview before drilling into individual audit slots.

``` r

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_perm_distribution(audit)
} else {
  cat("ggplot2 not installed; skipping permutation plot.\n")
}
```

The histogram shows the null distribution (gray bars) of the performance
metric under random label permutation.

The blue dashed line represents the average performance of a random
model (the permuted mean). The red solid line represents the observed
performance of the fitted model.

A genuine signal is indicated when the red line lies well to the right
of the gray distribution; overlap suggests weak or unstable signal.

Use the mechanism risk assessment (`audit@info$mechanism_summary`) as a
quick triage tool: flagged mechanisms point you to the specific audit
slots that warrant deeper investigation.

**Audit per learner**

``` r

if (requireNamespace("ranger", quietly = TRUE) && 
    requireNamespace("parsnip", quietly = TRUE)) {
    
    # 1. Define specs
    # Standard GLM (no tuning)
    spec_glm <- parsnip::logistic_reg() |> 
        parsnip::set_engine("glm")
    
    # Random Forest
    spec_rf <- parsnip::rand_forest(
        mode = "classification",
        trees = 100
    ) |>
        parsnip::set_engine("ranger")
    
    # 2. Fit using the current split object
    fit_multi <- fit_resample(
        df,
        outcome = "outcome",
        splits = splits,
        learner = list(glm = spec_glm, rf = spec_rf),
        metrics = "auc"
    )
    
    # 3. Run the audit
    audits <- audit_leakage_by_learner(fit_multi, metric = "auc", B = 20)
    cat("Per-learner audit summary:\n")
    print(audits)
    
} else {
    cat("ranger/parsnip not installed.\n")
}
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#> Per-learner audit summary:
#> LeakAuditList with 2 learners
#>  learner metric_obs perm_mean   gap p_value     z n_perm batch_p_min
#>      glm      0.537     0.418 0.120    0.19 1.316     20          NA
#>       rf      0.589     0.445 0.145    0.19 1.685     20          NA
#>  batch_col_min_p batch_v_max batch_col_max_v dup_pairs dup_threshold
#>             <NA>          NA            <NA>        NA         0.995
#>             <NA>          NA            <NA>        NA         0.995
#>  dup_max_sim
#>           NA
#>           NA
```

This example uses **ranger** for the random forest specification; if it
is not installed, the code chunk is skipped.

Use `parallel_learners = TRUE` to audit learners concurrently when
**future.apply** is available.

The printed table summarizes each learner’s observed metric, permutation
gap, p-value, and key batch/duplicate summaries. Use it to compare
signal strength across models while checking for leakage risks. Pass
`learners =` to audit a subset, or `mc.cores =` to control parallel
workers.

**HTML audit report**

[`audit_report()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_report.md)
accepts either a `LeakAudit` or a `LeakFit` object. When a `LeakFit` is
provided, it first runs
[`audit_leakage()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md)
and then forwards any additional arguments to the audit step. If
multiple learners were fit, pass `learner =` via `...` to select one.
Use `open = TRUE` to open the report in a browser after rendering.

``` r

if (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available()) {
  report_path <- audit_report(audit, output_dir = ".")
  cat("HTML report written to:\n", report_path, "\n")
} else {
  cat("rmarkdown or pandoc not available; skipping audit report rendering.\n")
}
```

The report path points to a standalone HTML file containing the same
audit tables and plots, suitable for sharing with collaborators or
archiving as a quality control record.

### Time-series leakage checks

Time-series data require special handling. Random splits can leak
information from the future into the past. Use `mode = "time_series"`
with a prediction `horizon`, and audit with block permutations. Choose
`time_block = "circular"` or `"stationary"`; when `block_len` is NULL,
the audit uses a default block length (~10% of the test block size,
minimum 5).

**Leaky example: lookahead feature**

``` r

time_splits <- make_split_plan(
  df_time,
  outcome = "outcome",
  mode = "time_series",
  time = "time",
  v = 4,
  horizon = 1
)

cat("Time-series splits summary:\n")
#> Time-series splits summary:
time_splits
#> LeakSplits object (mode = time_series, v = 4, repeats = 1)
#> Outcome: outcome | Stratified: FALSE | Nested: FALSE
#> ------------------------------------------------------
#>   fold repeat_id train_n test_n
#> 1    2         1      40     40
#> 2    3         1      80     40
#> 3    4         1     120     40
#> ------------------------------------------------------
#> Total folds: 3 | Hash: a0169c15b3a0826566edf98a905fbc04

fit_time_leaky <- fit_resample(
  df_time,
  outcome = "outcome",
  splits = time_splits,
  learner = spec,
  metrics = "auc"
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%

cat("Time-series leaky fit summary:\n")
#> Time-series leaky fit summary:
summary(fit_time_leaky)
#> 
#> ===========================
#>  bioLeak Model Fit Summary
#> ===========================
#> 
#> Task: binomial
#> Outcome: outcome
#> Positive class: case
#> Learners: logistic_reg/glm
#> Total folds: 3
#> Fold status: 3 success, 0 skipped, 0 failed
#> Refit performed: Yes
#> Hash: h3b18276c
#> 
#> Cross-validated metrics (mean ± SD):
#>            learner auc_mean auc_sd auc_ci_lo auc_ci_hi
#> 1 logistic_reg/glm    0.619  0.066      0.36     0.877
#> 
#> Audit overview:
#>  fold n_train n_test          learner features_final
#>     1      40     40 logistic_reg/glm             14
#>     2      80     40 logistic_reg/glm             14
#>     3     120     40 logistic_reg/glm             14
```

Time-series splitting trains on growing windows, so performance can
differ from standard CV simply because early folds have smaller training
sets. Regardless of the score, the presence of the `leak_future` feature
makes this estimate methodologically invalid.

**Leakage-safe alternative: remove lookahead and audit with blocks**

``` r

df_time_safe <- df_time
df_time_safe$leak_future <- NULL

fit_time_safe <- fit_resample(
  df_time_safe,
  outcome = "outcome",
  splits = time_splits,
  learner = spec,
  metrics = "auc"
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%

cat("Time-series safe fit summary:\n")
#> Time-series safe fit summary:
summary(fit_time_safe)
#> 
#> ===========================
#>  bioLeak Model Fit Summary
#> ===========================
#> 
#> Task: binomial
#> Outcome: outcome
#> Positive class: case
#> Learners: logistic_reg/glm
#> Total folds: 3
#> Fold status: 3 success, 0 skipped, 0 failed
#> Refit performed: Yes
#> Hash: h3b18276c
#> 
#> Cross-validated metrics (mean ± SD):
#>            learner auc_mean auc_sd auc_ci_lo auc_ci_hi
#> 1 logistic_reg/glm    0.608  0.042     0.442     0.773
#> 
#> Audit overview:
#>  fold n_train n_test          learner features_final
#>     1      40     40 logistic_reg/glm             13
#>     2      80     40 logistic_reg/glm             13
#>     3     120     40 logistic_reg/glm             13

audit_time <- audit_leakage(
  fit_time_safe,
  metric = "auc",
  B = 20,
  time_block = "stationary",
  block_len = 5
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%

cat("Time-series leakage audit summary:\n")
#> Time-series leakage audit summary:
summary(audit_time)
#> 
#> ==============================
#>  bioLeak Leakage Audit Summary
#> ==============================
#> 
#> Task: binomial | Outcome: outcome | Splitting mode: time_series | Positive class: case
#> Hash: h3b18276c | Folds: 3 | Repeats: 1
#> 
#> Label-Permutation Association Test:
#>   Method: refit per permutation (auto)
#>   Null: observed folds reused | Permutation: time_block | Summary: pooled
#>   Observed metric: 0.601
#>   Permuted mean ± SD: 0.500 ± 0.055
#>   Gap: 0.101 (larger gap = stronger non-random signal)
#>   This test does NOT diagnose information leakage. Use the Batch Association,
#>   Target Leakage Scan, and Duplicate Detection sections to check for leakage.
#> 
#> Batch / Study Association:
#>   batch (repeat 1): χ² = 8.237 (df = 10.000), p = 0.606
#>   study (repeat 1): χ² = 2.383 (df = 6.000), p = 0.881
#> 
#> Target Leakage Scan: not available.
#> 
#> Multivariate Target Scan: not available.
#> 
#> No near-duplicates detected.
#> 
#> Mechanism Risk Assessment:
#>        mechanism_class flagged        evidence statistic   p_value
#>      non_random_signal   FALSE permutation_gap 0.1009140 0.0952381
#>  confounding_alignment   FALSE     batch_assoc 0.1852551 1.0000000
#>   proxy_target_leakage   FALSE    target_assoc        NA        NA
#>      duplicate_overlap   FALSE      duplicates        NA        NA
#>     temporal_lookahead   FALSE      duplicates        NA        NA
#> 
#> Interpretation:
#>   ✓ Strong non-random signal.
if (!is.null(audit_time@permutation_gap) && nrow(audit_time@permutation_gap) > 0) {
  # Time-series permutation significance results
  audit_time@permutation_gap
}
#>     mechanism_class metric_obs perm_mean   perm_sd      gap       z   p_value
#> 1 non_random_signal   0.600679  0.499765 0.0551572 0.100914 1.84995 0.0952381
#>   n_perm
#> 1     20

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_time_acf(fit_time_safe, lag.max = 20)
} else {
  cat("ggplot2 not installed; skipping ACF plot.\n")
}
```

The safe fit summary provides the leakage-resistant performance
estimate. Compare leaky vs. safe fits to gauge inflation risk;
`features_final` should drop when the lookahead feature is removed.

The ACF plot shows the autocorrelation of out-of-fold predictions by
lag. The dashed red lines represent the 95% confidence interval for
white noise; large bars outside the bands indicate residual temporal
dependence.
[`plot_time_acf()`](https://selcukorkmaz.github.io/bioLeak/reference/plot_time_acf.md)
requires numeric predictions and time metadata aligned to the fit.

### Cross-validation uncertainty with `cv_ci()`

Standard cross-validation yields one metric value per fold. Because
folds share training data, the per-fold values are positively
correlated: the naive standard error `sd / sqrt(K)` systematically
underestimates variability, making the resulting confidence intervals
too narrow.
[`cv_ci()`](https://selcukorkmaz.github.io/bioLeak/reference/cv_ci.md)
addresses this with two approaches:

- **`"normal"`** (default): uses the naive SEM `sd / sqrt(K)` with a
  *t*-distribution on `K − 1` degrees of freedom. Fast and appropriate
  when fold overlap is low (large datasets, many folds).

- **`"nadeau_bengio"`**: applies the Nadeau–Bengio (2003) variance
  correction `(1/K + n_test / n_train) × var(x)`, explicitly accounting
  for the fraction of training samples shared across folds. Recommended
  for small datasets or few-fold CV where fold correlation is
  non-negligible.

The function accepts any data frame with `fold` and `learner` columns —
such as `fit@metrics` from a `LeakFit` object — and returns one row per
learner with four columns per metric: `<metric>_mean`, `<metric>_sd`,
`<metric>_ci_lo`, and `<metric>_ci_hi`.

``` r

# Standard 95% CI from the leakage-safe fit
ci_std <- cv_ci(fit_safe@metrics, level = 0.95, method = "normal")
cat("Standard 95% CI:\n")
#> Standard 95% CI:
print(ci_std)
#>            learner  auc_mean    auc_sd auc_ci_lo auc_ci_hi accuracy_mean
#> 1 logistic_reg/glm 0.6121717 0.1837801 0.3839785 0.8403648     0.5975198
#>   accuracy_sd accuracy_ci_lo accuracy_ci_hi
#> 1   0.1172633      0.4519182      0.7431215

# Nadeau-Bengio corrected CI
# For K-fold CV: n_train ≈ (K-1)/K × n, n_test ≈ n/K
K_folds  <- length(unique(fit_safe@metrics$fold))
n_total  <- nrow(fit_safe@splits@info$coldata)

ci_nb <- cv_ci(
  fit_safe@metrics,
  level   = 0.95,
  method  = "nadeau_bengio",
  n_train = round(n_total * (K_folds - 1L) / K_folds),
  n_test  = round(n_total / K_folds)
)
cat("Nadeau-Bengio corrected 95% CI:\n")
#> Nadeau-Bengio corrected 95% CI:
print(ci_nb)
#>            learner  auc_mean    auc_sd auc_ci_lo auc_ci_hi accuracy_mean
#> 1 logistic_reg/glm 0.6121717 0.1837801 0.2698819 0.9544614     0.5975198
#>   accuracy_sd accuracy_ci_lo accuracy_ci_hi
#> 1   0.1172633      0.3791174      0.8159223
```

The Nadeau–Bengio interval is typically wider because it correctly
inflates the variance to account for fold overlap. Use it whenever the
training-set fraction shared between consecutive folds is large — for
example, 5-fold CV reuses 80% of training data across fold pairs, so
correcting for this correlation matters.

When multiple learners were fit,
[`cv_ci()`](https://selcukorkmaz.github.io/bioLeak/reference/cv_ci.md)
returns one row per learner. For repeated K-fold CV, average `n_train`
and `n_test` across folds before passing them in; `fit@metrics` already
stores per-fold counts if
[`audit_leakage()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md)
annotated the object.

## Delta Leakage Sensitivity Index: Quantifying Performance Inflation

### Motivation

[`audit_leakage()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md)
detects whether leakage is present via permutation tests and confounding
checks. It cannot, however, directly answer: *by how much does leakage
inflate my reported performance?* The **Delta Leakage Sensitivity
Index** (ΔLSI) addresses this by rigorously comparing a naive
(potentially leaky) cross-validation pipeline against a leakage-guarded
one run on the same data.

### Mathematical framework

#### Notation

Let $`\mathcal{D} = \{(x_i, y_i)\}_{i=1}^n`$ denote a dataset of $`n`$
observations. A repeated $`K`$-fold cross-validation design with $`R`$
independent repeats partitions the observation indices into $`RK`$
folds. Let $`\mathcal{T}_{r,k} \subset \{1, \ldots, n\}`$ denote the
test index set for fold $`k`$ in repeat $`r`$, with
$`n_{r,k} = |\mathcal{T}_{r,k}|`$. Let $`m_{r,k}^{(p)}`$ denote the
performance metric evaluated on fold $`(r,k)`$ under pipeline
$`p \in \{\text{naive},\, \text{guarded}\}`$.

#### Stage 1: Within-repeat aggregation

For each repeat $`r`$ and pipeline $`p`$, the $`K`$ fold-level metric
values are combined into a size-weighted repeat mean:

``` math
\mu_r^{(p)} = \frac{\displaystyle\sum_{k=1}^{K} n_{r,k}\, m_{r,k}^{(p)}}{\displaystyle\sum_{k=1}^{K} n_{r,k}}
```

The test-fold size $`n_{r,k}`$ is the weight, so folds with more
held-out observations contribute proportionally more to the repeat-level
estimate. Folds with missing metric values are excluded from both the
numerator and denominator. This produces the `@repeats_naive` and
`@repeats_guarded` data frames (column `metric`), one row per repeat.

#### Stage 2: Per-repeat inflation score

For each repeat $`r`$, the signed performance gap between the two
pipelines is:

``` math
\Delta_r = s \cdot \bigl(\mu_r^{\text{naive}} - \mu_r^{\text{guarded}}\bigr),
\qquad r = 1, \ldots, R
```

where

``` math
s = \begin{cases}
+1 & \text{if the metric is higher-is-better (e.g.\ AUC, accuracy)} \\
-1 & \text{if the metric is lower-is-better (e.g.\ RMSE, MAE, log-loss)}
\end{cases}
```

By construction $`\Delta_r > 0`$ if and only if the naive pipeline is
more optimistic than the guarded pipeline in repeat $`r`$, irrespective
of metric direction. Because both pipelines are evaluated on the same
test fold $`\mathcal{T}_{r,k}`$ (the paired design), any between-fold
variability common to both pipelines cancels in $`\Delta_r`$,
concentrating the signal on genuine pipeline differences rather than
partition noise.

#### Point estimates

Two summaries of the inflation vector
$`\boldsymbol{\Delta} = (\Delta_1, \ldots, \Delta_R)^\top`$ are
reported.

**Arithmetic mean — `delta_metric`:**

``` math
\hat{\delta}_{\text{metric}} = \frac{1}{R} \sum_{r=1}^{R} \Delta_r
```

This is the most directly interpretable quantity. For AUC,
$`\hat{\delta}_{\text{metric}} = 0.08`$ means the naive pipeline’s
reported AUC exceeds the guarded pipeline’s by approximately $`0.08`$
units, averaged across all repeats. This value is stored in
`@delta_metric`.

**Huber M-estimate — `delta_lsi`:**

``` math
\hat{\delta}_{\text{LSI}} = \operatorname*{arg\,min}_{\mu} \sum_{r=1}^{R}
\rho_{k}\!\left(\frac{\Delta_r - \mu}{\hat{\sigma}}\right)
```

where $`\rho_k`$ is Huber’s loss function with tuning constant
$`k = 1.345`$:

``` math
\rho_k(u) = \begin{cases}
\tfrac{1}{2}u^2 & |u| \leq k \\
k|u| - \tfrac{1}{2}k^2 & |u| > k
\end{cases}
```

and
$`\hat{\sigma} = 1.4826 \times \operatorname{MAD}(\boldsymbol{\Delta})`$
is a fixed robust scale estimate. This minimisation is solved via
iteratively reweighted least squares (IRLS):

1.  Initialise:
    $`\hat{\mu}^{(0)} = \operatorname{median}(\boldsymbol{\Delta})`$,
    $`\hat{\sigma} = 1.4826 \times \operatorname{MAD}(\boldsymbol{\Delta})`$
2.  Compute standardised residuals:
    $`u_r^{(t)} = (\Delta_r - \hat{\mu}^{(t)}) / \hat{\sigma}`$
3.  Assign Huber influence weights:
    $`w_r^{(t)} = \min\!\left(1,\; k / |u_r^{(t)}|\right)`$
4.  Update location:
    $`\hat{\mu}^{(t+1)} = \displaystyle\frac{\sum_{r=1}^{R} w_r^{(t)} \Delta_r}{\sum_{r=1}^{R} w_r^{(t)}}`$
5.  Repeat until
    $`|\hat{\mu}^{(t+1)} - \hat{\mu}^{(t)}| < 10^{-8}(1 + |\hat{\mu}^{(t)}|)`$

Note that $`\hat{\sigma}`$ is **fixed** at its initial value throughout
iteration; only the location $`\hat{\mu}`$ is updated. Outlier repeats
with $`|u_r^{(t)}| > k`$ receive weight $`w_r < 1`$, reducing their
contribution to the update. When $`\hat{\sigma} = 0`$ (all $`\Delta_r`$
identical), the estimator returns the common value without iteration.
This quantity is stored in `@delta_lsi` and is the primary reported
point estimate.

The two estimates coincide when $`\{\Delta_r\}`$ are symmetric and free
of outliers. When they differ substantially, `delta_lsi` is more
trustworthy because it is not pulled by a single anomalous repeat.

**Unpaired fallback.** When the two fits do not share the same fold
structure (`@info$paired = FALSE`), the paired vector
$`\boldsymbol{\Delta}`$ cannot be formed. The point estimates are
instead computed from the separately summarised repeat-mean vectors
$`\boldsymbol{\mu}^{\text{naive}}`$ and
$`\boldsymbol{\mu}^{\text{guarded}}`$:

``` math
\hat{\delta}_{\text{metric}}^{\text{unpaired}}
  = s \cdot \bigl(\bar{\mu}^{\text{naive}} - \bar{\mu}^{\text{guarded}}\bigr),
\qquad
\hat{\delta}_{\text{LSI}}^{\text{unpaired}}
  = s \cdot \bigl(\hat{H}(\boldsymbol{\mu}^{\text{naive}})
               - \hat{H}(\boldsymbol{\mu}^{\text{guarded}})\bigr)
```

where $`\bar{\mu}`$ is the arithmetic mean and $`\hat{H}`$ is the Huber
M-estimator applied to each pipeline’s repeat-means vector
independently. No inference (sign-flip test or BCa interval) is
available in the unpaired case ($`R_{\text{eff}} = 0`$).

### Setting up a two-pipeline comparison

A ΔLSI analysis requires two `LeakFit` objects:

1.  **Naive fit** — standard cross-validation without leakage safeguards
    (e.g., random row splits, global preprocessing, or leaky features).
2.  **Guarded fit** — leakage-aware cross-validation using
    subject-grouped splits and preprocessing fit only on training folds.

For the most powerful **paired** design, generate a single `LeakSplits`
object and pass it to both
[`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md)
calls. Pairing is validated automatically by comparing the sorted
test-set indices of every fold; if the structures differ despite equal
repeat counts, a warning is issued and analysis falls back to unpaired
mode.

``` r

set.seed(100)
n_d    <- 120L
subj_d <- rep(paste0("P", seq_len(30L)), each = 4L)   # 30 subjects, 4 obs each

df_dlsi_guarded <- data.frame(
  subject = subj_d,
  outcome = factor(sample(c("case", "control"), n_d, replace = TRUE)),
  x1      = rnorm(n_d),
  x2      = rnorm(n_d)
)

# Leaky feature: subject-level mean of outcome
# With subject-grouped splits, test subjects' leak value encodes their outcome
df_dlsi_naive <- df_dlsi_guarded
df_dlsi_naive$leak <- ave(
  as.numeric(df_dlsi_guarded$outcome == "case"),
  subj_d, FUN = mean
)

# Shared splits: 5-fold × 5 repeats → R_eff = 5 (tier C: point + p-value)
splits_dlsi <- make_split_plan(
  df_dlsi_guarded,
  outcome = "outcome",
  mode    = "subject_grouped",
  group   = "subject",
  v       = 5L,
  repeats = 5L
)
```

``` r

# Reuse the parsnip spec (base R glm) defined earlier
# Naive fit: uses leaky feature (leak encodes test-subject outcome label)
fit_dlsi_naive <- fit_resample(
  df_dlsi_naive,
  outcome = "outcome",
  splits  = splits_dlsi,
  learner = spec,
  metrics = "auc",
  seed    = 1L
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%

# Guarded fit: same splits, clean features only
fit_dlsi_guarded <- fit_resample(
  df_dlsi_guarded,
  outcome = "outcome",
  splits  = splits_dlsi,
  learner = spec,
  metrics = "auc",
  seed    = 1L
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%

cat("Naive   AUC:", round(mean(fit_dlsi_naive@metrics$auc,   na.rm = TRUE), 3), "\n")
#> Naive   AUC: 0.716
cat("Guarded AUC:", round(mean(fit_dlsi_guarded@metrics$auc, na.rm = TRUE), 3), "\n")
#> Guarded AUC: 0.502
```

### Running `delta_lsi()`

``` r

result_dlsi <- delta_lsi(
  fit_leaky   = fit_dlsi_naive,
  fit_guarded = fit_dlsi_guarded,
  metric      = "auc",
  M_boot      = 300L,
  M_flip      = 300L,
  seed        = 1L
)

summary(result_dlsi)
#> 
#> =====================================
#>  bioLeak Delta LSI Summary
#> =====================================
#> 
#> Metric:            auc
#> Exchangeability:   iid
#> Inference tier:    C_signflip
#> R_eff:             5  (R_leaky=5, R_guarded=5, paired=TRUE)
#> 
#> Leaky pipeline:    0.716
#> Guarded pipeline:  0.502
#> 
#> Point estimates (leaky - guarded; positive = leakage inflation):
#>   delta_metric:  0.213  (raw metric difference)
#>   delta_lsi:     0.216  (Huber-robust)
#> 
#> Hypothesis test  [H0: no systematic inflation; paired repeat signs exchangeable]:
#>   Sign-flip p:   0.0625 
#>   [Note: with R_eff = 5 the smallest attainable p is 0.0625 > 0.05;
#>    use R_eff >= 6 repeats for a 5% sign-flip test.]
#> 
#> Inference valid:   NO (C_signflip; need paired R_eff >= 20)
```

[`summary()`](https://rdrr.io/r/base/summary.html) reports seven
quantities; here is how to read each one:

- **Naive / Guarded pipeline values**: the Huber-robust mean AUC across
  repeats for each pipeline in isolation. Compare these directly to see
  the raw performance level of each pipeline before computing the gap.

- **`delta_metric`**: the arithmetic mean of the per-repeat inflation
  scores $`\Delta_r`$. Interpret this in the original metric’s units: a
  value of 0.10 means the naive pipeline reports AUC ~0.10 higher than
  the guarded pipeline, averaged across repeats. This is the most
  intuitive summary for reporting.

- **`delta_lsi`**: the Huber M-estimate of $`\{\Delta_r\}`$. Because it
  down-weights outlier repeats, it is the more reliable point estimate
  when repeat-to-repeat variation is high. Values very close to
  `delta_metric` indicate a consistent, non-skewed inflation pattern.

- **95% BCa CI** (shown when `R_eff ≥ 10`): a bootstrap confidence
  interval for `delta_lsi`. If the lower bound is above zero the
  inflation is unlikely to be a sampling artefact. If the interval spans
  zero the inflation estimate is uncertain given the number of repeats.

- **Sign-flip p-value** (shown when `R_eff ≥ 5`): tests H₀: no leakage
  inflation (mean $`\Delta_r = 0`$) via a randomization test. A small
  p-value (e.g., \< 0.05) indicates the observed inflation is unlikely
  under the null of exchangeable pipeline labels.

- **Inference tier** and **`inference_ok`**: summarise how much
  statistical weight the design supports. Only tier A (`R_eff ≥ 20`,
  `inference_ok = TRUE`) provides a full combination of point estimate +
  CI + p-value. Lower tiers report what is reliable given the available
  repeats. The formal tier definitions and their rationale are given in
  the Statistical inference section below.

### Designing for a target inference tier

The formal tier thresholds have direct implications for how you set up
your repeated CV design. Plan the number of folds and repeats around the
tier you need:

- **Tier D** (`D_insufficient`, point estimates only): any
  $`R_{\text{eff}} < 5`$ — useful for exploratory comparisons where no
  formal test is required.
- **Tier C** (`C_signflip`, sign-flip p-value available): 5-fold
  $`\times`$ 5 repeats $`\to R_{\text{eff}} = 5`$. The exact two-sided
  sign-flip test has smallest achievable p-value
  $`2/2^{R_{\text{eff}}}`$, which is $`2/2^5 = 0.0625`$ at
  $`R_{\text{eff}} = 5`$: that design can never reach $`p < 0.05`$. Use
  at least 6 repeats ($`2/2^6 = 0.031`$) if you need a 5% test; the
  attainable floor is stored in `info$min_p_achievable`. For
  `exchangeability = "blocked_time"` the minimum p-value is
  $`2/2^{n_{\text{blocks}}}`$; at least 5 blocks are required for the
  test to be available.
- **Tier B** (`B_signflip_ci`, BCa CI additionally available): 5-fold
  $`\times`$ 10 repeats $`\to R_{\text{eff}} = 10`$. Both sign-flip
  p-value and BCa CI are reported.
- **Tier A** (`A_full_inference`, `inference_ok = TRUE`): 5-fold
  $`\times`$ 20 repeats $`\to R_{\text{eff}} = 20`$. Recommended for any
  published or pre-registered ΔLSI analysis.

For paired designs `R_eff` equals the number of matched repeats. For
unpaired designs `R_eff = 0` and only the raw point estimates are
available.

``` r

cat("Tier:          ", result_dlsi@tier, "\n")
#> Tier:           C_signflip
cat("R_eff:         ", result_dlsi@R_eff, "\n")
#> R_eff:          5
cat("inference_ok:  ", result_dlsi@inference_ok, "\n")
#> inference_ok:   FALSE
```

### Statistical inference

#### Sign-flip randomization test

The sign-flip test (Phipson & Smyth, 2010) assesses the null hypothesis
that there is no systematic performance inflation:

``` math
H_0: \mathbb{E}[\Delta_r] = 0 \quad \text{(no leakage-induced inflation)}
```

The test exploits sign-exchangeability: under $`H_0`$, each $`\Delta_r`$
is symmetrically distributed around zero, so randomly flipping the sign
of any $`\Delta_r`$ leaves the distribution of the test statistic
unchanged. The test statistic is the arithmetic mean of the inflation
scores:

``` math
T_{\text{obs}} = \frac{1}{R} \sum_{r=1}^{R} \Delta_r
```

For each sign vector
$`\boldsymbol{\epsilon} = (\epsilon_1, \ldots, \epsilon_R)`$ with
$`\epsilon_r \in \{-1, +1\}`$, a null-distribution realisation is:

``` math
T_{\boldsymbol{\epsilon}} = \frac{1}{R} \sum_{r=1}^{R} \epsilon_r \Delta_r
```

The two-sided p-value is the proportion of sign vectors yielding a
statistic at least as extreme as observed:

``` math
p = \frac{\#\bigl\{\boldsymbol{\epsilon} \in \mathcal{E} :
           |T_{\boldsymbol{\epsilon}}| \geq |T_{\text{obs}}|\bigr\}}{|\mathcal{E}|}
```

Two computation regimes are used:

- **Exact enumeration** ($`R_{\text{eff}} \leq 15`$): all $`2^R`$ sign
  vectors are evaluated — $`|\mathcal{E}| = 2^R`$. The p-value is
  $`\text{count}/2^R`$ with **no** continuity correction, because every
  sign combination including the observed one is explicitly enumerated.
  Because the test is two-sided, the observed signs and their mirror
  image always count, so the minimum achievable p-value is $`2/2^R`$.

- **Monte Carlo** ($`R_{\text{eff}} > 15`$): $`M_{\text{flip}}`$ sign
  vectors are drawn uniformly from $`\{-1,+1\}^R`$. The Phipson &
  Smyth (2010) correction prevents $`p = 0`$ in the extreme tail:
  ``` math
  p = \frac{B^{+} + 1}{M_{\text{flip}} + 1},
  \qquad B^{+} = \#\{b : |T_{\boldsymbol{\epsilon}_b}| \geq |T_{\text{obs}}|\}
  ```

The arithmetic mean $`\hat{\delta}_{\text{metric}}`$ is used as the test
statistic — not the Huber estimate $`\hat{\delta}_{\text{LSI}}`$. Under
a sign-flip null the mean has a symmetric, well-defined exchangeable
distribution; the Huber M-estimator’s nonlinear weighting changes with
each sign assignment and therefore does not preserve
sign-exchangeability. As a result, `@p_value` tests `@delta_metric`, not
`@delta_lsi`.

This distinction matters in practice. When the per-repeat inflation
vector $`\{\Delta_r\}`$ contains outliers, `delta_lsi` (Huber) and
`delta_metric` (mean) can diverge in magnitude while agreeing in sign. A
user can therefore see a significant `@p_value` with a BCa CI for
`delta_lsi` that spans zero — or vice versa — if a handful of repeats
pull the mean but not the robust estimate. The
[`summary()`](https://rdrr.io/r/base/summary.html) method flags this
disagreement explicitly when it occurs.

The sign-flip test is two-sided and, under exact enumeration, the
observed sign pattern and its mirror image are always “as extreme”, so
the minimum achievable p-value is $`2/2^{R}`$. With $`R = 5`$ this is
$`2/2^5 = 0.0625 > 0.05`$: tier C at $`R_{\text{eff}} = 5`$ reports a
p-value, but a significant result at the conventional level is
impossible regardless of the data. With $`R = 6`$,
$`2/2^6 = 0.03125 < 0.05`$, making it attainable. The floor for a given
fit is stored in `info$min_p_achievable`, and
[`summary()`](https://rdrr.io/r/base/summary.html) notes when it exceeds
0.05. For `exchangeability = "blocked_time"`, the relevant threshold is
$`n_{\text{blocks}} \geq 5`$ (not $`R_{\text{eff}} \geq 5`$), since
blocks are the independent units under the test; `@p_value` is set to
`NA` and a warning is issued when $`n_{\text{blocks}} < 5`$. For other
exchangeability values, `@p_value` is `NA` when $`R_{\text{eff}} < 5`$
or when the design is unpaired.

#### Huber M-estimator and robustness properties

The Huber M-estimator was described formally in the Mathematical
framework section above. Its inferential properties in this context are:

- **Breakdown point**: $`\hat{\delta}_{\text{LSI}}`$ remains a
  consistent estimate of the true location even when up to $`\sim 29\%`$
  of repeats are outliers ($`k = 1.345`$ yields 95% asymptotic
  efficiency relative to the mean under a Gaussian distribution of
  $`\{\Delta_r\}`$).
- **Fixed scale**: $`\hat{\sigma}`$ is set to
  $`1.4826 \times \operatorname{MAD}(\boldsymbol{\Delta})`$ before
  iteration and never updated. Fixing the scale prevents a single
  influential $`\Delta_r`$ from collapsing $`\hat{\sigma}`$ toward zero,
  which would force all weights to 1 and reduce the estimator to the
  ordinary mean.
- **Degenerate case**: when $`\hat{\sigma} = 0`$ (all repeat differences
  identical), the algorithm returns
  $`\operatorname{median}(\boldsymbol{\Delta})`$ directly, which equals
  the common value.

The `delta_lsi` and `delta_metric` estimates always have the same sign
because the Huber estimator and the arithmetic mean agree on direction.
They differ only in magnitude, and then only when one or more
$`\Delta_r`$ are extreme relative to the bulk of the distribution.

#### BCa bootstrap confidence interval

The BCa (bias-corrected and accelerated) interval (Efron, 1987) for
$`\hat{\delta}_{\text{LSI}}`$ corrects for both median bias and
distributional skewness in the bootstrap. Given $`M_{\text{boot}}`$
bootstrap resamples $`\boldsymbol{\Delta}_b^*`$ drawn with replacement
from $`\boldsymbol{\Delta}`$, let
$`\hat{\theta}_b^* = \hat{H}(\boldsymbol{\Delta}_b^*)`$ be the Huber
estimate on resample $`b`$ and
$`\hat{\theta} = \hat{H}(\boldsymbol{\Delta})`$ the observed estimate.

**Bias-correction constant** $`\hat{z}_0`$ accounts for the median bias
of the bootstrap distribution relative to the observed estimate:

``` math
\hat{z}_0 = \Phi^{-1}\!\left(
  \frac{\#\{b : \hat{\theta}_b^* < \hat{\theta}\}}{M_{\text{boot}}}
\right)
```

When $`\hat{\theta}`$ coincides with the exact median of the bootstrap
distribution, $`\hat{z}_0 = 0`$ and no bias correction is applied. A
nonzero $`\hat{z}_0`$ shifts both interval endpoints in the same
direction.

**Acceleration constant** $`\hat{a}`$ corrects for skewness via the
jackknife. Let
$`\hat{\theta}_{(-r)} = \hat{H}(\boldsymbol{\Delta}_{(-r)})`$ be the
Huber estimate with observation $`r`$ removed and
$`\bar{\theta}_{(\cdot)} = R^{-1}\sum_{r=1}^R \hat{\theta}_{(-r)}`$:

``` math
\hat{a} = \frac{\displaystyle\sum_{r=1}^{R}
                  \bigl(\bar{\theta}_{(\cdot)} - \hat{\theta}_{(-r)}\bigr)^3}
                 {6\,\left[\displaystyle\sum_{r=1}^{R}
                  \bigl(\bar{\theta}_{(\cdot)} - \hat{\theta}_{(-r)}\bigr)^2
                 \right]^{3/2}}
```

Here $`\hat{a} > 0`$ indicates a right-skewed sampling distribution;
$`\hat{a} < 0`$ indicates left skew.

**BCa-adjusted quantile levels**: for nominal coverage
$`1 - \alpha = 0.95`$, let $`z_L = \Phi^{-1}(\alpha/2)`$ and
$`z_U = \Phi^{-1}(1-\alpha/2)`$. The adjusted quantile levels are:

``` math
p_j = \Phi\!\left(\hat{z}_0 +
  \frac{\hat{z}_0 + z_j}{1 - \hat{a}\,(\hat{z}_0 + z_j)}
\right), \quad j \in \{L,\, U\}
```

The 95% BCa interval is the $`p_L`$ and $`p_U`$ empirical quantiles of
the bootstrap distribution $`\{\hat{\theta}_b^*\}`$:

``` math
\bigl[\hat{\theta}^*_{(p_L)},\; \hat{\theta}^*_{(p_U)}\bigr]
```

A separate BCa interval is computed for $`\hat{\delta}_{\text{metric}}`$
(using the arithmetic mean as the bootstrap statistic) and stored in
`@delta_metric_ci`. Both intervals are `NA` when
$`R_{\text{eff}} < 10`$, because fewer than 10 leave-one-out jackknife
samples do not provide a reliable estimate of $`\hat{a}`$.

#### Pairing condition and the inference tier system

Two `LeakFit` objects are treated as **paired** if and only if they have
the same number of repeats and every fold shares an identical test-index
composition:

``` math
\operatorname{sort}\!\bigl(\mathcal{T}_{r,k}^{\text{naive}}\bigr)
= \operatorname{sort}\!\bigl(\mathcal{T}_{r,k}^{\text{guarded}}\bigr)
\quad \forall\, r,\, k
```

This is verified by comparing the sorted integer index vectors for each
sequential fold position. When paired, $`R_{\text{eff}} = R`$; when
unpaired, $`R_{\text{eff}} = 0`$ and all inference is disabled. Pairing
is automatic when both fits are produced from the same `LeakSplits`
object.

All three inference procedures are gated by $`R_{\text{eff}}`$. The
four-tier system makes these dependencies explicit:

| Tier | Condition | `@delta_metric` | `@delta_lsi` | `@p_value` | `@delta_lsi_ci` | `@inference_ok` |
|----|----|:--:|:--:|:--:|:--:|:--:|
| `D_insufficient` | $`R_{\text{eff}} < 5`$ | ✓ | ✓ | `NA` | `NA` | `FALSE` |
| `C_signflip` | $`5 \leq R_{\text{eff}} < 10`$ | ✓ | ✓ | ✓ | `NA` | `FALSE` |
| `B_signflip_ci` | $`10 \leq R_{\text{eff}} < 20`$ | ✓ | ✓ | ✓ | ✓ | `FALSE` |
| `A_full_inference` | $`R_{\text{eff}} \geq 20`$ | ✓ | ✓ | ✓ | ✓ | `TRUE` |

Under `"iid"` exchangeability the exact two-sided sign-flip test has
minimum achievable p-value $`2/2^{R_{\text{eff}}}`$ ($`2/2^5 = 0.0625`$
at the tier C entry point; $`2/2^6 = 0.03125`$), so at least 6 paired
repeats are needed for $`p < 0.05`$ to be attainable. For
`exchangeability = "blocked_time"`, the independent units are blocks of
repeats, so the effective threshold is $`n_{\text{blocks}} \geq 5`$;
when fewer than five blocks are formed from the available repeats,
`@p_value` is set to `NA` even if $`R_{\text{eff}} \geq 5`$. The
threshold $`R_{\text{eff}} \geq 10`$ for BCa intervals ensures at least
10 jackknife leave-one-out estimates for a reliable $`\hat{a}`$. The
`inference_ok = TRUE` flag at tier A requires all three of
$`R_{\text{eff}} \geq 20`$, a finite `@p_value`, and finite
`@delta_lsi_ci` to hold simultaneously.

The code below shows direct slot access to confirm that the numbers
correspond exactly to the quantities defined above:

``` r

# delta_metric: arithmetic mean of {Δ_r}
cat("delta_metric:  ", round(result_dlsi@delta_metric, 4), "\n")
#> delta_metric:   0.2133
# delta_lsi: Huber M-estimate of {Δ_r} (k=1.345, fixed MAD scale)
cat("delta_lsi:     ", round(result_dlsi@delta_lsi,    4), "\n")
#> delta_lsi:      0.2159
# delta_lsi_ci: 95% BCa interval for delta_lsi (NA below tier B)
cat("95% BCa CI:    [",
    paste(round(result_dlsi@delta_lsi_ci, 4), collapse = ", "), "]\n")
#> 95% BCa CI:    [ NA, NA ]
# p_value: sign-flip randomization p-value (NA below tier C or if unpaired)
cat("p-value:       ",
    if (is.na(result_dlsi@p_value)) "NA (tier D or unpaired)"
    else format(result_dlsi@p_value, digits = 3), "\n")
#> p-value:        0.0625
```

### Paired vs. unpaired designs

**Paired** (recommended): both fits are produced from the same
`LeakSplits` object. The per-repeat inflation scores $`\Delta_r`$ are
computed as within-repeat differences, eliminating between-fold
variability. $`R_{\text{eff}} = R`$ and all inference tiers are
reachable. Equal repeat counts with differing index patterns produce a
warning and fall back to unpaired mode automatically.

**Unpaired**: the two fits have different fold structures or different
repeat counts. The point estimates use separately computed pipeline
means rather than paired differences; $`R_{\text{eff}} = 0`$ and all
inference is disabled. This mode is appropriate when the splitting
scheme itself is part of what the comparison is testing — for example,
evaluating whether subject-grouped splits alone (not a leaky feature)
account for a performance difference.

``` r

# Unpaired: naive uses sample-wise CV, guarded uses subject-grouped splits
splits_rand <- make_split_plan(
  df_dlsi_naive,
  outcome = "outcome",
  mode    = "subject_grouped",
  group   = "row_id",
  v       = 5L,
  repeats = 5L
)

fit_naive_rand <- fit_resample(
  df_dlsi_naive,
  outcome = "outcome",
  splits  = splits_rand,
  learner = spec,
  metrics = "auc",
  seed    = 1L
)
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%

r_unpaired <- suppressWarnings(
  delta_lsi(fit_naive_rand, fit_dlsi_guarded, metric = "auc", seed = 1L)
)

cat(sprintf("Paired:        %s\n", r_unpaired@info$paired))
#> Paired:        FALSE
cat(sprintf("R_eff:         %d  (0 = unpaired, inference disabled)\n", r_unpaired@R_eff))
#> R_eff:         0  (0 = unpaired, inference disabled)
cat(sprintf("delta_metric:  %.4f  (naive raw mean - guarded raw mean)\n",
            r_unpaired@delta_metric))
#> delta_metric:  -0.0585  (naive raw mean - guarded raw mean)
cat(sprintf("p_value:       %s\n",
            if (is.na(r_unpaired@p_value)) "NA" else r_unpaired@p_value))
#> p_value:       NA
```

### Metric direction with `higher_is_better`

For lower-is-better metrics (RMSE, MAE, log-loss, Brier score,
deviance), set `higher_is_better = FALSE` or let
[`delta_lsi()`](https://selcukorkmaz.github.io/bioLeak/reference/delta_lsi.md)
auto-detect from the metric name. Auto-detected lower-is-better metric
names: `rmse`, `mse`, `mae`, `log_loss`, `logloss`, `brier`, `error`,
`loss`, `deviance`.

When `higher_is_better = FALSE`, the sign of `delta_metric` is flipped
so that a **positive** value still means the naive pipeline is inflated
relative to guarded (i.e., naive has a *lower* error, which is
paradoxically *better* performance — indicating leakage):

``` r

r_hib_true  <- suppressWarnings(delta_lsi(
  fit_dlsi_naive, fit_dlsi_guarded,
  metric           = "auc",
  higher_is_better = TRUE,
  seed             = 1L
))
r_hib_false <- suppressWarnings(delta_lsi(
  fit_dlsi_naive, fit_dlsi_guarded,
  metric           = "auc",
  higher_is_better = FALSE,
  seed             = 1L
))

cat(sprintf("hib = TRUE:   delta_metric = %+.4f  (positive = naive inflated)\n",
            r_hib_true@delta_metric))
#> hib = TRUE:   delta_metric = +0.2133  (positive = naive inflated)
cat(sprintf("hib = FALSE:  delta_metric = %+.4f  (sign flipped)\n",
            r_hib_false@delta_metric))
#> hib = FALSE:  delta_metric = -0.2133  (sign flipped)
```

The `info$higher_is_better` slot records the value used (after
auto-detection), so analyses are fully reproducible:

``` r

cat("higher_is_better used:", result_dlsi@info$higher_is_better, "\n")
#> higher_is_better used: TRUE
```

### Accessing the `LeakDeltaLSI` object

[`delta_lsi()`](https://selcukorkmaz.github.io/bioLeak/reference/delta_lsi.md)
returns a `LeakDeltaLSI` S4 object. All slots are accessible directly
for programmatic downstream use:

| Slot | Type | Description |
|----|----|----|
| `@metric` | character | Metric name used (e.g. `"auc"`, `"rmse"`) |
| `@exchangeability` | character | Exchangeability assumption for the sign-flip test. `"blocked_time"` activates a block sign-flip procedure; `"by_group"` and `"within_batch"` are stored but inference still uses iid sign-flip (a warning is issued); `"iid"` (default) applies the standard independent sign-flip |
| `@tier` | character | Inference tier: `"A_full_inference"`, `"B_signflip_ci"`, `"C_signflip"`, or `"D_insufficient"` |
| `@R_eff` | integer | Effective paired repeats (0 when unpaired; all inference is based on this count) |
| `@delta_metric` | numeric | Arithmetic mean of per-repeat inflation scores $`\bar{\Delta}`$; interpret directly in metric units |
| `@delta_lsi` | numeric | Huber M-estimate of per-repeat inflation scores; more robust than `delta_metric` when repeat-level outliers are present |
| `@delta_metric_ci` | numeric\[2\] | 95% BCa bootstrap CI for `delta_metric` (both elements `NA` below tier B) |
| `@delta_lsi_ci` | numeric\[2\] | 95% BCa bootstrap CI for `delta_lsi` (both elements `NA` below tier B) |
| `@p_value` | numeric | Sign-flip randomization p-value testing H₀: mean inflation = 0 (`NA` below tier C or when unpaired) |
| `@inference_ok` | logical | `TRUE` only at tier A (`R_eff ≥ 20`), indicating full point + CI + p-value inference is available |
| `@folds_naive` | data.frame | Per-fold data for the naive pipeline: columns `fold`, `metric` (fold-level metric value), `repeat_id`, `n` (test-fold size) |
| `@folds_guarded` | data.frame | Per-fold data for the guarded pipeline: columns `fold`, `metric`, `repeat_id`, `n` — same structure as `@folds_naive` |
| `@repeats_naive` | data.frame | Per-repeat aggregates for the naive pipeline: columns `repeat_id`, `metric` (size-weighted mean across folds), `n_folds`, `total_n` |
| `@repeats_guarded` | data.frame | Per-repeat aggregates for the guarded pipeline: same structure as `@repeats_naive` |
| `@info` | list | Metadata: `paired` (logical), `R_naive`, `R_guarded` (number of repeats in each fit), `higher_is_better` (logical, after auto-detection), `metric_naive` and `metric_guarded` (Huber-robust overall mean for each pipeline), `block_size_used` and `n_blocks` (for `exchangeability = "blocked_time"`; `NA` otherwise), and optionally `delta_r` (per-repeat $`\Delta_r`$ vector) and the original fit objects when `return_details = TRUE` |

``` r

# Per-repeat summaries for custom plotting
head(result_dlsi@repeats_naive,   3)
#>   repeat_id    metric n_folds total_n
#> 1         1 0.7225532       5     120
#> 2         2 0.7018017       5     120
#> 3         3 0.7159788       5     120
head(result_dlsi@repeats_guarded, 3)
#>   repeat_id    metric n_folds total_n
#> 1         1 0.5078699       5     120
#> 2         2 0.4826270       5     120
#> 3         3 0.4834788       5     120

# Per-fold diagnostics
head(result_dlsi@folds_naive,   3)
#>   fold    metric repeat_id  n
#> 1    1 0.7132867         1 24
#> 2    2 0.6953125         1 24
#> 3    3 0.6888889         1 24
head(result_dlsi@folds_guarded, 3)
#>   fold    metric repeat_id  n
#> 1    1 0.5174825         1 24
#> 2    2 0.4765625         1 24
#> 3    3 0.4962963         1 24

# Metadata
str(result_dlsi@info)
#> List of 9
#>  $ paired          : logi TRUE
#>  $ R_naive         : int 5
#>  $ R_guarded       : int 5
#>  $ higher_is_better: logi TRUE
#>  $ metric_naive    : num 0.716
#>  $ metric_guarded  : num 0.502
#>  $ block_size_used : int NA
#>  $ n_blocks        : int NA
#>  $ min_p_achievable: num 0.0625
```

The `@repeats_naive` and `@repeats_guarded` data frames are the primary
substrate for diagnostic visualisation. Plotting `metric` by `repeat_id`
for both pipelines side-by-side reveals whether the performance gap is
consistent across all repeats — a clean, reproducible offset indicates
structural leakage — or concentrated in a handful of repeats, which may
indicate a data-quality issue or an unstable fold partition.

The `@folds_naive` and `@folds_guarded` data frames go one level deeper,
providing the per-fold metric value along with the test-fold size (`n`)
and the repeat it belongs to (`repeat_id`). These are useful for
identifying individual folds that behave as outliers within a repeat,
for example because a small test partition happened to produce extreme
metric values.

### Connection to `audit_leakage()`

[`delta_lsi()`](https://selcukorkmaz.github.io/bioLeak/reference/delta_lsi.md)
and
[`audit_leakage()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md)
answer complementary questions:

| Question | Tool |
|----|----|
| Does leakage produce non-random performance? | [`audit_leakage()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md) — permutation gap + p-value |
| What mechanisms drive the leakage? | [`audit_leakage()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md) — `batch_assoc`, `target_assoc`, `duplicates` |
| How large is the performance inflation? | [`delta_lsi()`](https://selcukorkmaz.github.io/bioLeak/reference/delta_lsi.md) — `delta_metric`, `delta_lsi` |
| Is the inflation statistically significant? | [`delta_lsi()`](https://selcukorkmaz.github.io/bioLeak/reference/delta_lsi.md) — sign-flip p-value |
| Do guards fully eliminate the inflation? | [`delta_lsi()`](https://selcukorkmaz.github.io/bioLeak/reference/delta_lsi.md) — compare before/after guarding |

A recommended two-step workflow:

1.  Run `audit_leakage(fit_leaky)` to detect leakage and diagnose which
    mechanism is responsible.
2.  Build a guarded pipeline that addresses those mechanisms, then run
    `delta_lsi(fit_leaky, fit_guarded)` to quantify how much inflation
    the leakage introduced and whether guarding eliminated it.

``` r

# Step 1: detect and characterise leakage
audit_naive_dlsi <- audit_leakage(fit_dlsi_naive, metric = "auc", B = 20L)
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |===                                                                   |   4%  |                                                                              |======                                                                |   8%  |                                                                              |========                                                              |  12%  |                                                                              |===========                                                           |  16%  |                                                                              |==============                                                        |  20%  |                                                                              |=================                                                     |  24%  |                                                                              |====================                                                  |  28%  |                                                                              |======================                                                |  32%  |                                                                              |=========================                                             |  36%  |                                                                              |============================                                          |  40%  |                                                                              |===============================                                       |  44%  |                                                                              |==================================                                    |  48%  |                                                                              |====================================                                  |  52%  |                                                                              |=======================================                               |  56%  |                                                                              |==========================================                            |  60%  |                                                                              |=============================================                         |  64%  |                                                                              |================================================                      |  68%  |                                                                              |==================================================                    |  72%  |                                                                              |=====================================================                 |  76%  |                                                                              |========================================================              |  80%  |                                                                              |===========================================================           |  84%  |                                                                              |==============================================================        |  88%  |                                                                              |================================================================      |  92%  |                                                                              |===================================================================   |  96%  |                                                                              |======================================================================| 100%
cat("Naive audit summary:\n")
#> Naive audit summary:
summary(audit_naive_dlsi)
#> 
#> ==============================
#>  bioLeak Leakage Audit Summary
#> ==============================
#> 
#> Task: binomial | Outcome: outcome | Splitting mode: subject_grouped | Positive class: control
#> Hash: h0cc183b7 | Folds: 25 | Repeats: 1
#> 
#> Label-Permutation Association Test:
#>   Method: refit per permutation (auto)
#>   Null: observed folds reused | Permutation: group_restricted | Summary: pooled
#>   Observed metric: 0.741
#>   Permuted mean ± SD: 0.446 ± 0.065
#>   Gap: 0.294 (larger gap = stronger non-random signal)
#>   This test does NOT diagnose information leakage. Use the Batch Association,
#>   Target Leakage Scan, and Duplicate Detection sections to check for leakage.
#> 
#> Batch / Study Association: none detected.
#> 
#> Target Leakage Scan: not available.
#> 
#> Multivariate Target Scan: not available.
#> 
#> No near-duplicates detected.
#> 
#> Mechanism Risk Assessment:
#>        mechanism_class flagged        evidence statistic  p_value
#>      non_random_signal    TRUE permutation_gap  0.294444 0.047619
#>  confounding_alignment   FALSE     batch_assoc        NA       NA
#>   proxy_target_leakage   FALSE    target_assoc        NA       NA
#>      duplicate_overlap   FALSE      duplicates        NA       NA
#>     temporal_lookahead   FALSE      duplicates        NA       NA
#> 
#> Interpretation:
#>   ✓ Strong non-random signal.

# Step 2: confirmed guards reduced inflation; check residual DLSI
cat("\ndelta_metric (naive vs guarded):", round(result_dlsi@delta_metric, 4), "\n")
#> 
#> delta_metric (naive vs guarded): 0.2133
if (!is.na(result_dlsi@p_value)) {
  cat("Sign-flip p-value:              ", format(result_dlsi@p_value, digits = 3), "\n")
}
#> Sign-flip p-value:               0.0625
```

A ΔLSI near zero after guarding and a non-significant p-value together
provide strong evidence that the leakage source has been eliminated. A
persistent positive ΔLSI after guarding indicates residual inflation —
revisit the preprocessing pipeline for remaining leakage pathways.

## Parallel Processing

`bioLeak` uses the `future` framework for parallelism (Windows, macOS,
Linux).
[`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md),
[`audit_leakage()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md),
[`audit_leakage_by_learner()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage_by_learner.md),
and
[`simulate_leakage_suite()`](https://selcukorkmaz.github.io/bioLeak/reference/simulate_leakage_suite.md)
honor the active plan when `parallel = TRUE` (or
`parallel_learners = TRUE`).

``` r

## future is a Suggests dependency. Install with install.packages("future")
## before running this code; otherwise library(future) errors.
library(future)

# Use multiple cores (works on all OS)
plan(multisession, workers = 4)

# Run a heavy simulation
sim <- simulate_leakage_suite(..., parallel = TRUE)

# Parallel folds or audits
# fit_resample(..., parallel = TRUE)
# audit_leakage(..., parallel = TRUE)
# audit_leakage_by_learner(..., parallel_learners = TRUE)

# Return to sequential processing
plan(sequential)
```

This chunk only configures the parallel plan, so it does not produce a
result object. Use it as a template before running compute-heavy
functions.

### Simulation suite

[`simulate_leakage_suite()`](https://selcukorkmaz.github.io/bioLeak/reference/simulate_leakage_suite.md)
runs Monte Carlo simulations that inject specific leakage mechanisms and
then evaluates detection with the audit pipeline. It is useful for
validating your leakage checks before applying them to real data.
Available leakage types include `none`, `subject_overlap`,
`batch_confounded`, `peek_norm`, and `lookahead`. Use `signal_strength`,
`prevalence`, `rho`, `K`, `repeats`, `horizon`, and `preprocess` to
control the simulation. The output is a `LeakSimResults` data frame with
one row per seed. Metrics are selected automatically based on outcome
type (AUC for binary classification).

``` r

if (requireNamespace("glmnet", quietly = TRUE)) {
  sim <- simulate_leakage_suite(
    n = 80,
    p = 6,
    mode = "subject_grouped",
    learner = "glmnet",
    leakage = "subject_overlap",
    seeds = 1:2,
    B = 20,
    parallel = FALSE,
    signal_strength = 1
  )
  
  # Simulation results (first 6 rows)
  head(sim)
} else {
  cat("glmnet not installed; skipping simulation suite example.\n")
}
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%
#>   seed metric_obs      gap  p_value         leakage            mode
#> 1    1   0.828914 0.333128 0.047619 subject_overlap subject_grouped
#> 2    2   0.825397 0.348524 0.047619 subject_overlap subject_grouped
```

Each row corresponds to one simulation seed.

- `metric_obs` is the observed performance (AUC for this simulation).
  Higher values suggest stronger apparent signal, which can be inflated
  by leakage.

- `gap` and `p_value` indicate that the leaked signal is statistically
  distinguishable from random noise.

Use `metric_obs` to gauge the magnitude of the leakage effect, and `gap`
to assess detection sensitivity.

### Objects and summaries

bioLeak uses S4 classes and list-like results to capture provenance and
diagnostics:

- `LeakSplits`: returned by
  [`make_split_plan()`](https://selcukorkmaz.github.io/bioLeak/reference/make_split_plan.md),
  printed with `show()`. Use
  [`check_split_overlap()`](https://selcukorkmaz.github.io/bioLeak/reference/check_split_overlap.md)
  to verify no-overlap invariants on grouping columns.
- `LeakFit`: returned by
  [`fit_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md),
  summarized with [`summary()`](https://rdrr.io/r/base/summary.html)
  includes fold diagnostics in `fit@info$fold_status`
- `LeakAudit`: returned by
  [`audit_leakage()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md),
  summarized with [`summary()`](https://rdrr.io/r/base/summary.html).
  Now includes `audit@info$mechanism_summary` — a data.frame classifying
  leakage evidence into mechanism classes (`non_random_signal`,
  `confounding_alignment`, `proxy_target_leakage`, `duplicate_overlap`,
  `temporal_lookahead`).
- `LeakAuditList`: returned by
  [`audit_leakage_by_learner()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage_by_learner.md),
  printed directly
- `GuardFit`: returned by
  [`guard_fit()`](https://selcukorkmaz.github.io/bioLeak/reference/guard_fit.md),
  printed and summarized. Use
  [`guard_to_recipe()`](https://selcukorkmaz.github.io/bioLeak/reference/guard_to_recipe.md)
  to convert guard preprocessing specs to a
  [`recipes::recipe`](https://recipes.tidymodels.org/reference/recipe.html)
  for tidymodels workflows.
- `LeakImpute`: returned by
  [`impute_guarded()`](https://selcukorkmaz.github.io/bioLeak/reference/impute_guarded.md)
  (guarded data + diagnostics)
- `LeakTune`: returned by
  [`tune_resample()`](https://selcukorkmaz.github.io/bioLeak/reference/tune_resample.md)
  (nested tuning summaries) includes `fold_status`, selected parameters,
  and optional `final_model` / `final_workflow` when `refit = TRUE`
- `LeakSimResults`: returned by
  [`simulate_leakage_suite()`](https://selcukorkmaz.github.io/bioLeak/reference/simulate_leakage_suite.md)
  (per-seed results)
- `LeakDeltaLSI`: returned by
  [`delta_lsi()`](https://selcukorkmaz.github.io/bioLeak/reference/delta_lsi.md)
  (paired pipeline comparison). Key slots: `@metric`,
  `@exchangeability`, `@tier`, `@R_eff`, `@delta_metric`,
  `@delta_metric_ci`, `@delta_lsi`, `@delta_lsi_ci`, `@p_value`,
  `@inference_ok`, `@folds_naive`, `@folds_guarded`, `@repeats_naive`,
  `@repeats_guarded`, `@info`. Summarise with
  [`summary()`](https://rdrr.io/r/base/summary.html).

Set `options(bioLeak.strict = TRUE)` to activate strict leakage mode,
which escalates trained-recipe/workflow warnings to errors, warns on
missing seeds, and auto-runs overlap checks after splitting.

These objects store hashes and metadata that help reproduce and audit
the workflow. Inspect their slots (`@metrics`, `@predictions`,
`@permutation_gap`, `@duplicates`) to drill into specific leakage
signals.
