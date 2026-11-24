# Plot per-repeat \\\Delta_r\\ values from a LeakDeltaLSI object

Visualises the per-repeat inflation \\\Delta_r\\ (leaky minus guarded,
sign-adjusted so that positive always means the leaky pipeline looks
better, also for lower-is-better metrics such as RMSE) for a
[`LeakDeltaLSI`](https://selcukorkmaz.github.io/bioLeak/reference/LeakClasses.md)
object, overlaid with the robust Huber point estimate, the arithmetic
mean, and the corrected-t 95% confidence interval for the Huber
estimate. This is the diagnostic shown as Figure 4 panel (b) of the
manuscript. Requires ggplot2.

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
estimates, the confidence interval (`ci`; also returned as `bca_ci` for
backward compatibility), `n_repeats`, the object's `R_eff`, and the
ggplot object.

## See also

[`delta_lsi`](https://selcukorkmaz.github.io/bioLeak/reference/delta_lsi.md),
[`dlsi_repeats`](https://selcukorkmaz.github.io/bioLeak/reference/dlsi_repeats.md)
