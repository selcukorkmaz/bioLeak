# Plot method for LeakAudit

Diagnostic plot for a \[\`LeakAudit\`\] object. The default diagnostic
is the permutation-distribution histogram produced by
[`plot_perm_distribution`](https://selcukorkmaz.github.io/bioLeak/reference/plot_perm_distribution.md).

## Usage

``` r
# S4 method for class 'LeakAudit,missing'
plot(x, y, ...)
```

## Arguments

- x:

  A \[\`LeakAudit\`\] object.

- y:

  Unused; present for S4 compatibility with
  [`base::plot`](https://rdrr.io/r/base/plot.html).

- ...:

  Additional arguments passed to
  [`plot_perm_distribution`](https://selcukorkmaz.github.io/bioLeak/reference/plot_perm_distribution.md).

## Value

Invisibly returns the list produced by
[`plot_perm_distribution`](https://selcukorkmaz.github.io/bioLeak/reference/plot_perm_distribution.md)
(observed value, permuted mean, permutation values, ggplot object).

## See also

[`plot_perm_distribution`](https://selcukorkmaz.github.io/bioLeak/reference/plot_perm_distribution.md),
[`LeakAudit`](https://selcukorkmaz.github.io/bioLeak/reference/LeakClasses.md)
