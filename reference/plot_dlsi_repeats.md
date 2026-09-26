# Plot per-repeat \\\Delta_r\\ values from a LeakDeltaLSI object

Visualises the per-repeat metric differences (leaky minus guarded) for a
[`LeakDeltaLSI`](https://selcukorkmaz.github.io/bioLeak/reference/LeakClasses.md)
object, overlaid with the robust Huber point estimate, the arithmetic
mean, and the BCa bootstrap confidence interval. This is the diagnostic
shown as Figure 4 panel (b) of the manuscript. Requires ggplot2.

## Usage

``` r
plot_dlsi_repeats(dlsi)
```

## Arguments

- dlsi:

  A
  [`LeakDeltaLSI`](https://selcukorkmaz.github.io/bioLeak/reference/LeakClasses.md)
  object produced by
  [`delta_lsi`](https://selcukorkmaz.github.io/bioLeak/reference/delta_lsi.md).

## Value

A list with the per-repeat deltas, the robust and arithmetic-mean
estimates, the BCa confidence interval, and the ggplot object.

## See also

[`delta_lsi`](https://selcukorkmaz.github.io/bioLeak/reference/delta_lsi.md),
[`dlsi_repeats`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_repeats.md)
