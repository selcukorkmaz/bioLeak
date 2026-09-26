# Plot method for LeakDeltaLSI

Diagnostic plot for a \[\`LeakDeltaLSI\`\] object: per-repeat
\\\Delta_r\\ scatter with the Huber-robust point estimate, the
arithmetic mean, and the BCa bootstrap confidence interval band. This is
the diagnostic shown as Figure 4 panel (b) of the manuscript.

## Usage

``` r
# S4 method for class 'LeakDeltaLSI,missing'
plot(x, y, ...)
```

## Arguments

- x:

  A \[\`LeakDeltaLSI\`\] object.

- y:

  Unused; present for S4 compatibility with
  [`base::plot`](https://rdrr.io/r/base/plot.html).

- ...:

  Additional arguments (currently unused).

## Value

Invisibly returns the list produced by
[`plot_dlsi_repeats`](https://selcukorkmaz.github.io/bioLeak/reference/plot_dlsi_repeats.md):
per-repeat deltas, the robust and arithmetic-mean estimates, the BCa
interval, and the ggplot object.

## See also

[`plot_dlsi_repeats`](https://selcukorkmaz.github.io/bioLeak/reference/plot_dlsi_repeats.md),
[`LeakDeltaLSI`](https://selcukorkmaz.github.io/bioLeak/reference/LeakClasses.md)
