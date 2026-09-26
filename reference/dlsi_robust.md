# Huber-robust delta_lsi point estimate from a LeakDeltaLSI

Returns the Huber-robust point estimate of the per-repeat delta values
stored in a \[\`LeakDeltaLSI\`\] object.

## Usage

``` r
dlsi_robust(dlsi)

# S4 method for class 'LeakDeltaLSI'
dlsi_robust(dlsi)
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

\[LeakClasses\], \[delta_lsi()\], \[dlsi_metric()\], \[dlsi_ci()\]
