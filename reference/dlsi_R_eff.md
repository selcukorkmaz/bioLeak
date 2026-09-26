# Effective number of paired repeats from a LeakDeltaLSI

Returns the effective number of paired repeats \`R_eff\` stored in a
\[\`LeakDeltaLSI\`\] object. This is the count of repeats that
contribute to the inference; it equals the smaller of the leaky and
guarded fits' repeat counts when the comparison is paired.

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
