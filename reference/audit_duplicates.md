# Near-duplicate pairs from a LeakAudit

Returns the near-duplicate sample pairs data frame stored in a
\[\`LeakAudit\`\] object. Each row is a unique (row_a, row_b) pair above
the configured cosine-similarity threshold that crossed train/test
partitions in at least one fold.

## Usage

``` r
audit_duplicates(audit)

# S4 method for class 'LeakAudit'
audit_duplicates(audit)
```

## Arguments

- audit:

  A \[\`LeakAudit\`\] object returned by \[audit_leakage()\].

## Value

A \`data.frame\` with one row per detected near-duplicate pair.

## Details

Implemented as an S4 generic with a method for \[\`LeakAudit\`\];
visible via \`methods(class = "LeakAudit")\`.

## See also

\[LeakClasses\], \[audit_leakage()\]
