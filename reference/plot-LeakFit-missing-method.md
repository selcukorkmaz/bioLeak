# Plot method for LeakFit

Diagnostic plot for a \[\`LeakFit\`\] object. The default diagnostic is
the fold-balance check produced by
[`plot_fold_balance`](https://selcukorkmaz.github.io/bioLeak/reference/plot_fold_balance.md),
which works for any classification task. Use the `which` argument to
switch to one of the other diagnostics in the package: `"overlap"`
([`plot_overlap_checks`](https://selcukorkmaz.github.io/bioLeak/reference/plot_overlap_checks.md)),
`"calibration"`
([`plot_calibration`](https://selcukorkmaz.github.io/bioLeak/reference/plot_calibration.md);
binary outcomes only), `"time_acf"`
([`plot_time_acf`](https://selcukorkmaz.github.io/bioLeak/reference/plot_time_acf.md);
time-ordered splits), or `"confounder_sensitivity"`
([`plot_confounder_sensitivity`](https://selcukorkmaz.github.io/bioLeak/reference/plot_confounder_sensitivity.md)).

## Usage

``` r
# S4 method for class 'LeakFit,missing'
plot(
  x,
  y,
  which = c("fold_balance", "overlap", "calibration", "time_acf",
    "confounder_sensitivity"),
  ...
)
```

## Arguments

- x:

  A \[\`LeakFit\`\] object.

- y:

  Unused; present for S4 compatibility with
  [`base::plot`](https://rdrr.io/r/base/plot.html).

- which:

  One of `"fold_balance"` (default), `"overlap"`, `"calibration"`,
  `"time_acf"`, `"confounder_sensitivity"`.

- ...:

  Additional arguments passed to the selected helper.

## Value

Invisibly returns the list produced by the selected helper.

## See also

[`plot_fold_balance`](https://selcukorkmaz.github.io/bioLeak/reference/plot_fold_balance.md),
[`plot_overlap_checks`](https://selcukorkmaz.github.io/bioLeak/reference/plot_overlap_checks.md),
[`plot_calibration`](https://selcukorkmaz.github.io/bioLeak/reference/plot_calibration.md),
[`plot_time_acf`](https://selcukorkmaz.github.io/bioLeak/reference/plot_time_acf.md),
[`plot_confounder_sensitivity`](https://selcukorkmaz.github.io/bioLeak/reference/plot_confounder_sensitivity.md),
[`LeakFit`](https://selcukorkmaz.github.io/bioLeak/reference/LeakClasses.md)
