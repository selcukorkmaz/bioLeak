# Sign-flip p-value from a LeakDeltaLSI

Returns the paired sign-flip randomization-test p-value stored in a
\[\`LeakDeltaLSI\`\] object.

## Usage

``` r
dlsi_p_value(dlsi)

# S4 method for class 'LeakDeltaLSI'
dlsi_p_value(dlsi)
```

## Arguments

- dlsi:

  A \[\`LeakDeltaLSI\`\] object returned by \[delta_lsi()\].

## Value

A length-one numeric scalar in \`\[0, 1\]\`, or \`NA_real\_\` when the
inference tier did not include hypothesis testing.

## Details

Implemented as an S4 generic with a method for \[\`LeakDeltaLSI\`\];
visible via \`methods(class = "LeakDeltaLSI")\`.

## See also

\[LeakClasses\], \[delta_lsi()\], \[dlsi_tier()\]
