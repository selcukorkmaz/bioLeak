# Batch / study association results from a LeakAudit

Returns the batch/study chi-squared association data frame stored in a
\[\`LeakAudit\`\] object. Columns include the metadata column, repeat,
chi-squared statistic, degrees of freedom, p-value, and Cramer's V.

## Usage

``` r
audit_batch_assoc(audit)

# S4 method for class 'LeakAudit'
audit_batch_assoc(audit)
```

## Arguments

- audit:

  A \[\`LeakAudit\`\] object returned by \[audit_leakage()\].

## Value

A \`data.frame\` with one row per (metadata column, repeat).

## Details

Implemented as an S4 generic with a method for \[\`LeakAudit\`\];
visible via \`methods(class = "LeakAudit")\`.

## See also

\[LeakClasses\], \[audit_leakage()\]
