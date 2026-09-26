# BCa confidence interval from a LeakDeltaLSI

Returns the bias-corrected and accelerated (BCa) bootstrap confidence
interval stored in a \[\`LeakDeltaLSI\`\] object. By default returns the
interval for the Huber-robust delta estimate; set \`which = "metric"\`
to return the interval for the raw metric difference instead.

## Usage

``` r
dlsi_ci(dlsi, which = c("robust", "metric"))

# S4 method for class 'LeakDeltaLSI'
dlsi_ci(dlsi, which = c("robust", "metric"))
```

## Arguments

- dlsi:

  A \[\`LeakDeltaLSI\`\] object returned by \[delta_lsi()\].

- which:

  Either \`"robust"\` (default) for the Huber estimate's confidence
  interval, or \`"metric"\` for the raw arithmetic mean's confidence
  interval.

## Value

A length-two numeric vector \`c(lower, upper)\`. Returns \`c(NA_real\_,
NA_real\_)\` when the interval is not computed (for example, when the
inference tier did not include CIs).

## Details

Implemented as an S4 generic with a method for \[\`LeakDeltaLSI\`\];
visible via \`methods(class = "LeakDeltaLSI")\`.

## See also

\[LeakClasses\], \[delta_lsi()\], \[dlsi_robust()\], \[dlsi_metric()\]
