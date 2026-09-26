# Target-leakage scan results from a LeakAudit

Returns the target-association data frame stored in a \[\`LeakAudit\`\]
object. Each row is one predictor with its association score (rescaled
AUC; \`\|AUC - 0.5\| \* 2\`), threshold-based flag, and (where
applicable) p-value.

## Usage

``` r
audit_target_assoc(audit)

# S4 method for class 'LeakAudit'
audit_target_assoc(audit)
```

## Arguments

- audit:

  A \[\`LeakAudit\`\] object returned by \[audit_leakage()\].

## Value

A \`data.frame\` with one row per predictor.

## Details

Implemented as an S4 generic with a method for \[\`LeakAudit\`\];
visible via \`methods(class = "LeakAudit")\`.

## See also

\[LeakClasses\], \[audit_leakage()\]
