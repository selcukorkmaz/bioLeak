# Per-repeat metric data frames from a LeakDeltaLSI

Returns the per-repeat metric data frame for one of the two pipelines
stored in a \[\`LeakDeltaLSI\`\] object. The naive (or leaky) pipeline's
repeats are returned by default; set \`which = "guarded"\` to return the
guarded pipeline's repeats.

## Usage

``` r
dlsi_repeats(dlsi, which = c("naive", "guarded"))

# S4 method for class 'LeakDeltaLSI'
dlsi_repeats(dlsi, which = c("naive", "guarded"))
```

## Arguments

- dlsi:

  A \[\`LeakDeltaLSI\`\] object returned by \[delta_lsi()\].

- which:

  Either \`"naive"\` (default) for the naive/leaky pipeline's per-repeat
  data frame, or \`"guarded"\` for the guarded pipeline's per-repeat
  data frame.

## Value

A \`data.frame\` with one row per repeat.

## Details

Implemented as an S4 generic with a method for \[\`LeakDeltaLSI\`\];
visible via \`methods(class = "LeakDeltaLSI")\`.

## See also

\[LeakClasses\], \[delta_lsi()\]
