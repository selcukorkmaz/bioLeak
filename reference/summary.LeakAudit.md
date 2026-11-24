# Summarize a leakage audit

Prints a concise, human-readable report for a \`LeakAudit\` object
produced by \[audit_leakage()\]. The summary surfaces four diagnostics
when available: label-permutation gap (prediction-label association by
default), batch/study association tests (batch-outcome and fold-batch),
target leakage scan (features strongly associated with the outcome), and
near-duplicate detection (high similarity in \`X_ref\`). The output
reflects the stored audit results only; it does not recompute any tests.

## Usage

``` r
# S3 method for class 'LeakAudit'
summary(object, digits = 3, ...)
```

## Arguments

- object:

  A \`LeakAudit\` object from \[audit_leakage()\]. The summary reads
  stored results from \`object\` and prints them to the console.

- digits:

  Integer number of digits to show when formatting numeric statistics in
  the console output. Defaults to \`3\`. Increasing \`digits\` shows
  more precision; decreasing it shortens the printout without changing
  the underlying values.

- ...:

  Unused. Included for S3 method compatibility; additional arguments are
  ignored.

## Value

Invisibly returns \`object\` after printing the summary.

## Details

The permutation test quantifies prediction-label association when using
fixed predictions; refit-based permutations require \`perm_refit =
TRUE\` (or \`"auto"\` with refit data). It does not by itself prove or
rule out leakage. The batch section reports, per batch/study column, the
batch-outcome test, the fold-batch test (how batches are spread over
folds) and the share of test rows whose batch also occurs in training.
Confounding is flagged in the mechanism assessment only when batch is
associated with the outcome and batches are shared by training and test
rows (see \[audit_leakage()\]). Target leakage scan uses univariate
feature-outcome associations and can miss multivariate proxies,
interaction leakage, or features not included in \`X_ref\`. The
multivariate scan (enabled by default for supported tasks) reports an
additional model-based score. Duplicate detection only considers the
provided \`X_ref\` features and the similarity threshold used during
\[audit_leakage()\]. By default, \`duplicate_scope = "train_test"\`
filters to pairs that cross train/test; set \`duplicate_scope = "all"\`
to include within-fold duplicates. Sections are reported as "not
available" when the corresponding audit component was not computed.

## See also

\[plot_perm_distribution()\], \[plot_fold_balance()\],
\[plot_overlap_checks()\]

## Examples

``` r
set.seed(1)
df <- data.frame(
  subject = rep(1:6, each = 2),
  outcome = rbinom(12, 1, 0.5),
  x1 = rnorm(12),
  x2 = rnorm(12)
)
splits <- make_split_plan(df, outcome = "outcome",
                      mode = "subject_grouped", group = "subject", v = 3)
#> subject_grouped: repeat 1/1 done.
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
fit <- fit_resample(df, outcome = "outcome", splits = splits,
                    learner = "glm", custom_learners = custom,
                    metrics = "auc", refit = FALSE, seed = 1)
#> 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |=======================                                               |  33%
#> Warning: glm.fit: algorithm did not converge
#> Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred
#> 
  |                                                                            
  |===============================================                       |  67%
  |                                                                            
  |======================================================================| 100%
audit <- audit_leakage(fit, metric = "auc", B = 5,
                       X_ref = df[, c("x1", "x2")], seed = 1)
#> audit_leakage(): perm_refit = "auto" refits the model on each permuted outcome: 5 permutations x 3 folds = 15 model fits, which can be slow. Set perm_refit = FALSE for the fixed-prediction test (no refitting).
#> 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |=======================                                               |  33%
  |                                                                            
  |===============================================                       |  67%
  |                                                                            
  |======================================================================| 100%
#> 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |=======================                                               |  33%
  |                                                                            
  |===============================================                       |  67%
  |                                                                            
  |======================================================================| 100%
#> 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |=======================                                               |  33%
  |                                                                            
  |===============================================                       |  67%
  |                                                                            
  |======================================================================| 100%
#> 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |=======================                                               |  33%
  |                                                                            
  |===============================================                       |  67%
  |                                                                            
  |======================================================================| 100%
#> 
  |                                                                            
  |                                                                      |   0%
  |                                                                            
  |=======================                                               |  33%
  |                                                                            
  |===============================================                       |  67%
  |                                                                            
  |======================================================================| 100%
summary(audit) # prints the audit report and returns `audit` invisibly
#> 
#> ==============================
#>  bioLeak Leakage Audit Summary
#> ==============================
#> 
#> Task: binomial | Outcome: outcome | Splitting mode: subject_grouped | Positive class: 1
#> Hash: h84f22b32 | Folds: 3 | Repeats: 1
#> 
#> Label-Permutation Association Test:
#>   Method: refit per permutation (auto)
#>   Null: observed folds reused | Permutation: group_restricted | Summary: pooled
#>   Observed metric: 0.611
#>   Permuted mean ± SD: 0.344 ± 0.150
#>   Gap: 0.267, permutation p = 0.167 (B = 5)
#>   This test does NOT diagnose information leakage. Use the Batch Association,
#>   Target Leakage Scan, and Duplicate Detection sections to check for leakage.
#> 
#> Batch / Study Association: not tested (no batch/study columns available).
#> 
#> Target Leakage Scan:
#>   Features checked: 2 | Flagged (score ≥ 0.900): 0
#>   No strong proxy features detected.
#> 
#> Multivariate Target Scan:
#>   Metric: auc | Score: 0.167 | p = 0.693
#>   Features: 2 | Components: 2 | Interactions: 1 | Permutations: 100
#> 
#> No near-duplicates detected.
#>   Expected by chance: 0.0 pairs
#> 
#> Mechanism Risk Assessment:
#>        mechanism_class flagged        evidence statistic  p_value
#>      non_random_signal   FALSE permutation_gap 0.2666670 0.166667
#>  confounding_alignment   FALSE     batch_assoc        NA       NA
#>   proxy_target_leakage   FALSE    target_assoc 0.4444444       NA
#>      duplicate_overlap   FALSE      duplicates        NA       NA
#>     temporal_lookahead   FALSE      duplicates        NA       NA
#> 
#> Interpretation:
#>   ⚠ No significant non-random signal (p = 0.167): performance is not distinguishable from the label-permutation null.
#>   A significant gap shows predictions carry signal; it does not show the
#>   absence of leakage, which also produces a large gap.
#>   No leakage mechanism flagged; these checks cannot rule leakage out.
```
