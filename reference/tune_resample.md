# Leakage-aware nested tuning with tidymodels

Runs nested cross-validation for hyperparameter tuning using
leakage-aware splits. Inner resamples are constructed from each outer
training fold to avoid information leakage during tuning. Requires
tidymodels tuning packages and a workflow or recipe-based preprocessing.
Survival tasks are not yet supported.

## Usage

``` r
tune_resample(
  x,
  outcome,
  splits,
  learner,
  preprocess = NULL,
  grid = 10,
  metrics = NULL,
  positive_class = NULL,
  selection = c("best", "one_std_err"),
  selection_metric = NULL,
  inner_v = NULL,
  inner_repeats = 1,
  inner_seed = NULL,
  control = NULL,
  parallel = FALSE,
  refit = FALSE,
  seed = 1,
  split_cols = "auto",
  tune_threshold = FALSE,
  threshold_grid = seq(0.1, 0.9, by = 0.05),
  threshold_metric = "accuracy"
)
```

## Arguments

- x:

  SummarizedExperiment or matrix/data.frame.

- outcome:

  Outcome column name (if x is SE or data.frame).

- splits:

  LeakSplits object defining the outer resamples. If the splits do not
  already include inner folds, they are created from each outer training
  fold using the same split metadata. rsample splits must already
  include inner folds.

- learner:

  A parsnip model_spec with tunable parameters, or a workflows workflow.
  When a model_spec is provided, a workflow is built using
  \`preprocess\` or a formula.

- preprocess:

  Optional \`recipes::recipe\`. Required when you need preprocessing for
  tuning. Ignored when \`learner\` is already a workflow.
  Recipe/workflow leakage guardrails run before tuning; configure policy
  via `options(bioLeak.validation_mode = "warn" | "error" | "off")`.

- grid:

  Tuning grid passed to \`tune::tune_grid()\`. Can be a data.frame or an
  integer size.

- metrics:

  Character vector of metric names (\`auc\`, \`pr_auc\`, \`accuracy\`,
  \`macro_f1\`, \`log_loss\`, \`rmse\`) or a yardstick metric set/list.
  Metrics are computed with yardstick; unsupported metrics are dropped
  with a warning. For binomial tasks, if any inner assessment fold
  contains a single class, probability metrics (\`auc\`, \`roc_auc\`,
  \`pr_auc\`) are dropped for tuning with a warning.

- positive_class:

  Optional value indicating the positive class for binomial outcomes.
  When set, the outcome levels are reordered so the positive class is
  second.

- selection:

  Selection rule for tuning, either \`"best"\` or \`"one_std_err"\`.

- selection_metric:

  Metric name used for selecting hyperparameters. Defaults to the first
  metric in \`metrics\`. If the chosen metric yields no valid results,
  the first available metric is used with a warning.

- inner_v:

  Optional number of folds for inner CV when inner splits are not
  precomputed. Defaults to the outer \`v\`.

- inner_repeats:

  Optional number of repeats for inner CV when inner splits are not
  precomputed. Defaults to 1.

- inner_seed:

  Optional seed for inner split generation when inner splits are not
  precomputed. Defaults to the outer split seed.

- control:

  Optional \`tune::control_grid()\` settings for tuning.

- parallel:

  Logical; passed to \[fit_resample()\] when evaluating outer folds
  (single-fold, no refit).

- refit:

  Logical; if TRUE, refits a final tuned workflow on the full dataset
  using aggregated hyperparameters across all outer folds (median for
  numeric parameters, majority vote for categorical). This avoids
  nested-CV leakage that would occur from selecting a single fold's
  params.

- seed:

  Integer seed for reproducibility.

- split_cols:

  Optional named list/character vector or \`"auto"\` (default)
  overriding group/batch/study/time column names when \`splits\` is an
  rsample object and its attributes are missing. \`"auto"\` falls back
  to common metadata column names (e.g., \`group\`, \`subject\`,
  \`batch\`, \`study\`, \`time\`). Supported names are \`group\`,
  \`batch\`, \`study\`, and \`time\`.

- tune_threshold:

  Logical; when \`TRUE\` for binomial tasks, selects a probability
  threshold from inner-fold predictions and applies it only to the
  corresponding outer-fold evaluation.

- threshold_grid:

  Numeric vector of thresholds in \`\[0, 1\]\` considered when
  \`tune_threshold = TRUE\`.

- threshold_metric:

  Metric used to pick thresholds when \`tune_threshold = TRUE\`.
  Supported values are \`"accuracy"\`, \`"balanced_accuracy"\`, and
  \`"f1"\`, or a custom function with signature \`function(truth,
  pred_class, prob, threshold)\`.

## Value

A list of class \`"LeakTune"\` with components:

- metrics:

  Outer-fold metrics.

- metric_summary:

  Mean/SD metrics across outer folds with columns `learner`, and
  `<metric>_mean` and `<metric>_sd` for each metric.

- best_params:

  Best hyperparameters per outer fold.

- inner_results:

  List of inner tuning results.

- outer_fits:

  List of outer LeakFit objects.

- thresholds:

  Per-fold threshold choices when threshold tuning is enabled.

- fold_status:

  Outer-fold status log with stage, status, reason, and notes.

- final_model:

  Optional final workflow fit when \`refit = TRUE\`.

- info:

  Metadata about the tuning run.

## Examples

``` r
# \donttest{
  if (requireNamespace("tune", quietly = TRUE) &&
      requireNamespace("recipes", quietly = TRUE) &&
      requireNamespace("glmnet", quietly = TRUE) &&
      requireNamespace("rsample", quietly = TRUE) &&
      requireNamespace("workflows", quietly = TRUE) &&
      requireNamespace("yardstick", quietly = TRUE) &&
     requireNamespace("dials", quietly = TRUE)) {
    df <- data.frame(
      subject = rep(1:10, each = 2),
      outcome = factor(rep(c(0, 1), each = 10)),
      x1 = rnorm(20),
      x2 = rnorm(20)
    )
    splits <- make_split_plan(df, outcome = "outcome",
                         mode = "subject_grouped", group = "subject",
                         v = 3, nested = TRUE, stratify = TRUE)
    spec <- parsnip::logistic_reg(penalty = tune::tune(), mixture = 1) |>
      parsnip::set_engine("glmnet")
    rec <- recipes::recipe(outcome ~ x1 + x2, data = df)
    tuned <- tune_resample(df, outcome = "outcome", splits = splits,
                          learner = spec, preprocess = rec, grid = 5)
    tuned$metric_summary
  }
#> subject_grouped: repeat 1/1 done.
#> → A | warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#> There were issues with some computations   A: x1
#> There were issues with some computations   A: x3
#> 
#>   |                                                                              |                                                                      |   0%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |======================================================================| 100%
#> → A | warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |                                                                      |   0%
#> Warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |======================================================================| 100%
#> → A | warning: one multinomial or binomial class has fewer than 8  observations; dangerous ground
#>   |                                                                              |                                                                      |   0%  |                                                                              |======================================================================| 100%
#>               learner accuracy_mean accuracy_sd roc_auc_mean roc_auc_sd
#> 1 logistic_reg/glmnet           0.5           0          0.5          0
#>   pr_auc_mean pr_auc_sd accuracy_ci_lo accuracy_ci_hi roc_auc_ci_lo
#> 1        0.75         0            0.5            0.5           0.5
#>   roc_auc_ci_hi pr_auc_ci_lo pr_auc_ci_hi
#> 1           0.5         0.75         0.75
# }
```
