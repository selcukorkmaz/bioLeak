# Audit leakage per learner

Runs \[audit_leakage()\] separately for each learner recorded in a
\[LeakFit\] and returns a named list of \[LeakAudit\] objects. Use this
when a single fit contains predictions for multiple models and you want
model-specific audits. If predictions do not include learner IDs, only a
single audit can be run and requesting multiple learners is an error.

## Usage

``` r
audit_leakage_by_learner(
  fit,
  metric = c("auc", "pr_auc", "accuracy", "macro_f1", "log_loss", "rmse", "cindex"),
  learners = NULL,
  parallel_learners = FALSE,
  mc.cores = NULL,
  ...
)
```

## Arguments

- fit:

  A \[LeakFit\] object produced by \[fit_resample()\]. It must contain
  predictions and split metadata. Learner IDs must be present in
  predictions to audit multiple models.

- metric:

  Character scalar. One of \`"auc"\`, \`"pr_auc"\`, \`"accuracy"\`,
  \`"macro_f1"\`, \`"log_loss"\`, \`"rmse"\`, or \`"cindex"\`. Controls
  which metric is audited for each learner.

- learners:

  Character vector or NULL. If NULL (default), audits all learners found
  in predictions. If provided, must match learner IDs stored in the
  predictions. Supplying more than one learner requires learner IDs.

- parallel_learners:

  Logical scalar. If TRUE, runs per-learner audits in parallel using
  \`future.apply\` (if installed). This changes runtime but not the
  audit results.

- mc.cores:

  Integer scalar or NULL. Number of workers used when
  \`parallel_learners = TRUE\`. Defaults to the minimum of available
  cores and the number of learners.

- ...:

  Additional named arguments forwarded to \[audit_leakage()\] for each
  learner. These control the audit itself. Common options include: \`B\`
  (integer permutations), \`perm_stratify\` (logical or \`"auto"\`),
  \`perm_refit\` (logical), \`perm_refit_spec\` (list), \`time_block\`
  (character), \`block_len\` (integer or NULL), \`include_z\` (logical),
  \`ci_method\` (character), \`boot_B\` (integer), \`parallel\`
  (logical), \`seed\` (integer), \`return_perm\` (logical),
  \`batch_cols\` (character vector), \`coldata\` (data.frame), \`X_ref\`
  (matrix/data.frame), \`target_scan\` (logical), \`target_threshold\`
  (numeric), \`target_p_adjust\` (character), \`target_alpha\`
  (numeric), \`feature_space\` (character), \`sim_method\` (character),
  \`sim_threshold\` (numeric), \`nn_k\` (integer), \`max_pairs\`
  (integer), and \`duplicate_scope\` (character). See
  \[audit_leakage()\] for full definitions; changing these values
  changes each learner's audit.

## Value

A named list of
[`LeakAudit`](https://selcukorkmaz.github.io/bioLeak/reference/LeakClasses.md)
objects, where each element is keyed by the learner ID (character
string). Each `LeakAudit` object contains the same slots as described in
[`audit_leakage`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md):
`fit`, `permutation_gap`, `perm_values`, `batch_assoc`, `target_assoc`,
`duplicates`, `trail`, and `info`. Use
[`names()`](https://rdrr.io/r/base/names.html) to retrieve learner IDs,
and access individual audits with `[[learner_id]]` or `$learner_id`.
Each audit reflects the performance and diagnostics for that specific
learner's predictions.

## Examples

``` r
set.seed(1)
df <- data.frame(
  subject = rep(1:6, each = 2),
  outcome = factor(rep(c(0, 1), 6)),
  x1 = rnorm(12),
  x2 = rnorm(12)
)
splits <- make_split_plan(df, outcome = "outcome",
                      mode = "subject_grouped", group = "subject",
                      v = 3, progress = FALSE)
custom <- list(
  glm = list(
    fit = function(x, y, task, weights, ...) {
      stats::glm(y ~ ., data = data.frame(y = y, x),
                 family = stats::binomial(), weights = weights)
    },
    predict = function(object, newdata, task, ...) {
      as.numeric(stats::predict(object,
                                newdata = as.data.frame(newdata),
                                type = "response"))
    }
  )
)
custom$glm2 <- custom$glm
fit <- fit_resample(df, outcome = "outcome", splits = splits,
                    learner = c("glm", "glm2"), custom_learners = custom,
                    metrics = "auc", refit = FALSE, seed = 1)
#>   |                                                                              |                                                                      |   0%  |                                                                              |=======================                                               |  33%  |                                                                              |===============================================                       |  67%  |                                                                              |======================================================================| 100%
audits <- audit_leakage_by_learner(fit, metric = "auc", B = 10,
                                   perm_stratify = FALSE)
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
names(audits)
#> [1] "glm"  "glm2"
```
