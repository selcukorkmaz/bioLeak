# Permutation-gap test results from a LeakAudit

Returns the permutation-gap test data frame stored in a
\[\`LeakAudit\`\] object. Columns include the observed metric,
permuted-null mean and SD, gap, z-score, and permutation p-value.

## Usage

``` r
audit_perm_gap(audit)

# S4 method for class 'LeakAudit'
audit_perm_gap(audit)
```

## Arguments

- audit:

  A \[\`LeakAudit\`\] object returned by \[audit_leakage()\].

## Value

A \`data.frame\` with one row per (mechanism class, repeat), summarising
the permutation-gap test.

## Details

Implemented as an S4 generic with a method for \[\`LeakAudit\`\];
visible via \`methods(class = "LeakAudit")\`.

## See also

\[LeakClasses\], \[audit_leakage()\], \[audit_target_assoc()\]
