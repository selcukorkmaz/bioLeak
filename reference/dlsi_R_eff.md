# Effective number of independent units from a LeakDeltaLSI

Returns \`R_eff\` stored in a \[\`LeakDeltaLSI\`\] object: the number of
folds per repeat (the number of studies for \`study_loocv\`) when the
comparison is paired, and 0 otherwise. Repeats re-partition the same
data, so they sharpen the estimate but do not add independent units.

## Usage

``` r
dlsi_R_eff(dlsi)

# S4 method for class 'LeakDeltaLSI'
dlsi_R_eff(dlsi)
```

## Arguments

- dlsi:

  A \[\`LeakDeltaLSI\`\] object returned by \[delta_lsi()\].

## Value

A length-one integer.

## Details

Implemented as an S4 generic with a method for \[\`LeakDeltaLSI\`\];
visible via \`methods(class = "LeakDeltaLSI")\`.

## See also

\[LeakClasses\], \[delta_lsi()\], \[dlsi_tier()\]
