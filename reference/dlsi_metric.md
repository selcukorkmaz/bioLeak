# Raw delta-metric estimate from a LeakDeltaLSI

Returns the arithmetic mean of the per-repeat raw metric differences
(leaky minus guarded) stored in a \[\`LeakDeltaLSI\`\] object.

## Usage

``` r
dlsi_metric(dlsi)

# S4 method for class 'LeakDeltaLSI'
dlsi_metric(dlsi)
```

## Arguments

- dlsi:

  A \[\`LeakDeltaLSI\`\] object returned by \[delta_lsi()\].

## Value

A length-one numeric scalar.

## Details

Implemented as an S4 generic with a method for \[\`LeakDeltaLSI\`\];
visible via \`methods(class = "LeakDeltaLSI")\`.

## See also

\[LeakClasses\], \[delta_lsi()\], \[dlsi_robust()\], \[dlsi_ci()\]
