# Inference tier from a LeakDeltaLSI

Returns the inference tier label stored in a \[\`LeakDeltaLSI\`\]
object. Possible values are \`"A_full_inference"\` (\`R_eff \>= 20\`),
\`"B_signflip_ci"\` (\`R_eff \>= 10\`), \`"C_signflip"\` (\`R_eff \>=
3\`), or \`"D_insufficient"\` (\`R_eff \< 3\` or unpaired). \`R_eff\` is
the number of folds per repeat; the labels are historical, and tiers A-C
all report the corrected repeated-CV t-test p-value and CIs.

## Usage

``` r
dlsi_tier(dlsi)

# S4 method for class 'LeakDeltaLSI'
dlsi_tier(dlsi)
```

## Arguments

- dlsi:

  A \[\`LeakDeltaLSI\`\] object returned by \[delta_lsi()\].

## Value

A length-one character string giving the tier label.

## Details

Implemented as an S4 generic with a method for \[\`LeakDeltaLSI\`\];
visible via \`methods(class = "LeakDeltaLSI")\`.

## See also

\[LeakClasses\], \[delta_lsi()\], \[dlsi_R_eff()\]
