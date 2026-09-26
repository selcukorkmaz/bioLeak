# Fit and evaluate with leakage guards over predefined splits

Performs cross-validated model training and evaluation using
leakage-protected preprocessing (guard_fit) and user-specified learners.

## Usage

``` r
fit_resample(
  x,
  outcome,
  splits,
  preprocess = list(impute = list(method = "median"), normalize = list(method =
    "zscore"), filter = list(var_thresh = 0, iqr_thresh = 0), fs = list(method = "none")),
  learner = c("glmnet", "ranger"),
  learner_args = list(),
  custom_learners = list(),
  metrics = c("auc", "pr_auc", "accuracy"),
  class_weights = NULL,
  positive_class = NULL,
  classification_threshold = 0.5,
  parallel = FALSE,
  refit = TRUE,
  seed = 1,
  split_cols = "auto",
  store_refit_data = TRUE,
  id_cols = NULL
)
```

## Arguments

- x:

  SummarizedExperiment or matrix/data.frame

- outcome:

  outcome column name (if x is SE or data.frame), or a length-2
  character vector of time/event column names for survival outcomes.

- splits:

  LeakSplits object from make_split_plan(), or an \`rsample\`
  rset/rsplit.

- preprocess:

  list(impute, normalize, filter=list(...), fs) or a \`recipes::recipe\`
  object. When a recipe is supplied, the guarded preprocessing pipeline
  is bypassed and the recipe is prepped on training data only.
  Recipe/workflow leakage guardrails run before fitting; configure
  policy via
  `options(bioLeak.validation_mode = "warn" | "error" | "off")`.

- learner:

  parsnip model_spec (or list of model_spec objects) describing the
  model(s) to fit, or a \`workflows::workflow\`. For legacy use, a
  character vector of learner names (e.g., "glmnet", "ranger") or custom
  learner IDs is still supported.

- learner_args:

  list of additional arguments passed to legacy learners (ignored when
  \`learner\` is a parsnip model_spec).

- custom_learners:

  named list of custom learner definitions used only with legacy
  character learners. Each entry must contain `fit` and `predict`
  functions. The `fit` function should accept `x`, `y`, `task`, and
  `weights`, and return a model object. The `predict` function should
  accept `object`, `newdata`, and `task`. For
  binomial/regression/survival tasks it should return a numeric vector;
  for multiclass tasks it should return either class labels or a
  matrix/data.frame of class probabilities.

- metrics:

  named list of metric functions, vector of metric names, or a
  \`yardstick::metric_set\`. When a yardstick metric set (or list of
  yardstick metric functions) is supplied, metrics are computed using
  yardstick with the positive class set to the second factor level.

- class_weights:

  optional named numeric vector of weights for binomial or multiclass
  outcomes

- positive_class:

  optional value indicating the positive class for binomial outcomes.
  When set, the outcome levels are reordered so that `positive_class` is
  treated as the positive class (level 2). If NULL, the second factor
  level is used.

- classification_threshold:

  Numeric threshold in `[0, 1]` used to convert binomial probabilities
  into class predictions for `pred_class` and accuracy metrics. Ignored
  for non-binomial tasks.

- parallel:

  logical, use future.apply for multicore execution

- refit:

  logical, if TRUE retrain final model on full data

- seed:

  integer, for reproducibility

- split_cols:

  Optional named list/character vector or \`"auto"\` (default)
  overriding group/batch/study/time column names when \`splits\` is an
  rsample object and its attributes are missing. \`"auto"\` falls back
  to common metadata column names (e.g., \`group\`, \`subject\`,
  \`batch\`, \`study\`, \`time\`). Supported names are \`group\`,
  \`batch\`, \`study\`, and \`time\`.

- store_refit_data:

  Logical; when TRUE (default), stores the original data and learner
  configuration inside the fit to enable refit-based permutation tests
  without manual \`perm_refit_spec\` setup.

- id_cols:

  Optional character vector of identifier or metadata columns (for
  example a `sample_id`) to exclude from the predictors. The outcome and
  the split's group/batch/study/time columns are always excluded; every
  other column of `x` is used as a predictor, and character/factor
  columns are one-hot encoded. With the built-in guarded preprocessing,
  a warning is raised when a character or factor column not listed here
  has nearly one unique value per row (at least 90% unique), because
  such columns are almost always identifiers.

## Value

A
[`LeakFit`](https://selcukorkmaz.github.io/bioLeak/reference/LeakClasses.md)
S4 object containing:

- `splits`:

  The `LeakSplits` object used for resampling.

- `metrics`:

  Data.frame of per-fold, per-learner performance metrics with columns
  `fold`, `learner`, and one column per requested metric.

- `metric_summary`:

  Data.frame summarizing metrics across folds for each learner with
  columns `learner`, and `<metric>_mean` and `<metric>_sd` for each
  requested metric.

- `audit`:

  Data.frame with per-fold audit information including `fold`,
  `n_train`, `n_test`, `learner`, and `features_final` (number of
  features after preprocessing).

- `predictions`:

  List of data.frames containing out-of-fold predictions with columns
  `id` (sample identifier), `truth` (true outcome), `pred` (predicted
  value or probability), `fold`, and `learner`. For classification
  tasks, includes `pred_class`. For multiclass, includes per-class
  probability columns.

- `preprocess`:

  List of preprocessing state objects from each fold, storing imputation
  parameters, normalization statistics, and feature selection results.

- `learners`:

  List of fitted model objects from each fold.

- `outcome`:

  Character string naming the outcome variable.

- `task`:

  Character string indicating the task type (`"binomial"`,
  `"multiclass"`, `"gaussian"`, or `"survival"`).

- `feature_names`:

  Character vector of feature names after preprocessing.

- `info`:

  List of additional metadata including `hash`, `metrics_used`,
  `class_weights`, `positive_class`, `sample_ids`, `fold_status`,
  `refit`, `final_model` (refitted model if `refit = TRUE`),
  `final_preprocess`, `learner_names`, and `perm_refit_spec` (for
  permutation-based audits).

Use [`summary()`](https://rdrr.io/r/base/summary.html) to print a
formatted report. For programmatic access to slot contents, the
recommended interface is the S4 accessor method registered for
`LeakFit`:

- [`fit_metrics`](https://selcukorkmaz.github.io/bioLeak/reference/fit_metrics.md)`(fit)`
  – per-fold metric data frame (the `metrics` slot).

Slots without dedicated accessors (`predictions`, `info`, `audit`,
`learners`, `preprocess`, `feature_names`, `outcome`, `task`, `splits`,
`metric_summary`) are read directly via the standard `@` operator when
needed; the list of all accessor methods for the class is available
through `methods(class = "LeakFit")`.

## Details

Preprocessing is fit on the training fold and applied to the test fold,
preventing leakage from global imputation, scaling, or feature
selection. When a \`recipes::recipe\` or \`workflows::workflow\` is
supplied, the recipe is prepped on the training fold and baked on the
test fold. For data.frame or matrix inputs, columns used to define
splits (outcome, group, batch, study, time) are excluded from the
predictor matrix. Use `learner_args` to pass model-specific arguments,
either as a named list keyed by learner or a single list applied to all
learners. For custom learners, `learner_args[[name]]` may be a list with
`fit` and `predict` sublists to pass distinct arguments to each stage.
For binomial tasks, predictions and metrics assume the positive class is
the second factor level; use `positive_class` to control this. Use
`classification_threshold` to change the probability cutoff used for
class labels and accuracy. Parsnip learners must support probability
predictions for binomial metrics (AUC/PR-AUC/accuracy) and multiclass
log-loss when requested.

## Examples

``` r
set.seed(1)
df <- data.frame(
  subject = rep(1:10, each = 2),
  outcome = rbinom(20, 1, 0.5),
  x1 = rnorm(20),
  x2 = rnorm(20)
)
splits <- make_split_plan(df, outcome = "outcome",
                      mode = "subject_grouped", group = "subject", v = 5)
#> subject_grouped: repeat 1/1 done.

# glmnet learner (requires glmnet package)
fit <- fit_resample(df, outcome = "outcome", splits = splits,
                      learner = "glmnet", metrics = "auc")
#>   |                                                                              |                                                                      |   0%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#>   |                                                                              |==============                                                        |  20%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#>   |                                                                              |============================                                          |  40%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#>   |                                                                              |==========================================                            |  60%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#>   |                                                                              |========================================================              |  80%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
#>   |                                                                              |======================================================================| 100%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> Warning: Option grouped=FALSE enforced in cv.glmnet, since < 3 observations per fold
summary(fit)
#> 
#> ===========================
#>  bioLeak Model Fit Summary
#> ===========================
#> 
#> Task: binomial
#> Outcome: outcome
#> Positive class: 1
#> Learners: glmnet
#> Total folds: 5
#> Fold status: 4 success, 1 skipped, 0 failed
#> Refit performed: Yes
#> Hash: hc667e323
#> 
#> Cross-validated metrics (mean ± SD):
#>   learner auc_mean auc_sd auc_ci_lo auc_ci_hi
#> 1  glmnet      0.5      0       0.5       0.5
#> 
#> Audit overview:
#>  fold n_train n_test learner features_final
#>     1      16      4  glmnet              2
#>     2      16      4  glmnet              2
#>     4      16      4  glmnet              2
#>     5      16      4  glmnet              2
#> 

# Custom learner (logistic regression) - no extra packages needed
custom <- list(
  glm = list(
    fit = function(x, y, task, weights, ...) {
      stats::glm(y ~ ., data = as.data.frame(x),
                 family = stats::binomial(), weights = weights)
    },
    predict = function(object, newdata, task, ...) {
      as.numeric(stats::predict(object, newdata = as.data.frame(newdata), type = "response"))
    }
  )
)
fit2 <- fit_resample(df, outcome = "outcome", splits = splits,
                     learner = "glm", custom_learners = custom,
                     metrics = "accuracy")
#>   |                                                                              |                                                                      |   0%  |                                                                              |==============                                                        |  20%  |                                                                              |============================                                          |  40%  |                                                                              |==========================================                            |  60%  |                                                                              |========================================================              |  80%  |                                                                              |======================================================================| 100%

summary(fit2)
#> 
#> ===========================
#>  bioLeak Model Fit Summary
#> ===========================
#> 
#> Task: binomial
#> Outcome: outcome
#> Positive class: 1
#> Learners: glm
#> Total folds: 5
#> Fold status: 5 success, 0 skipped, 0 failed
#> Refit performed: Yes
#> Hash: hc667e323
#> 
#> Cross-validated metrics (mean ± SD):
#>   learner accuracy_mean accuracy_sd accuracy_ci_lo accuracy_ci_hi
#> 1     glm           0.5        0.25          0.034          0.966
#> 
#> Audit overview:
#>  fold n_train n_test learner features_final
#>     1      16      4     glm              2
#>     2      16      4     glm              2
#>     3      16      4     glm              2
#>     4      16      4     glm              2
#>     5      16      4     glm              2
#> 
```
