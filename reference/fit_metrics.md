# Per-fold metric data frame from a LeakFit

Returns the per-fold metric data frame stored in a \[\`LeakFit\`\]
object. Each row is one (fold, repeat, learner) combination; columns
include the requested metric values such as \`auc\` and any
task-specific performance scores.

## Usage

``` r
fit_metrics(fit)

# S4 method for class 'LeakFit'
fit_metrics(fit)
```

## Arguments

- fit:

  A \[\`LeakFit\`\] object returned by \[fit_resample()\].

## Value

A \`data.frame\` with one row per (fold, repeat, learner) combination.

## Details

Implemented as an S4 generic with a method for \[\`LeakFit\`\]; visible
via \`methods(class = "LeakFit")\`.

## See also

\[LeakClasses\], \[fit_resample()\], \[audit_perm_gap()\]

## Examples

``` r
set.seed(1)
df <- data.frame(
  subject = rep(1:6, each = 2),
  outcome = factor(rep(c("a","b"), 6), levels = c("a","b")),
  x1 = rnorm(12), x2 = rnorm(12)
)
splits <- make_split_plan(df, outcome = "outcome",
                          mode = "subject_grouped", group = "subject", v = 3)
#> subject_grouped: repeat 1/1 done.
if (FALSE) { # \dontrun{
  fit <- fit_resample(df, outcome = "outcome", splits = splits,
                      learner = parsnip::logistic_reg() |>
                                parsnip::set_engine("glm"),
                      metrics = "auc")
  fit_metrics(fit)
} # }
```
