# Auxiliary information from a LeakAudit

Returns the auxiliary information list stored in a \[\`LeakAudit\`\]
object. The list typically contains multivariate-target-scan results,
configuration flags, permutation-test diagnostics, and provenance
metadata.

## Usage

``` r
audit_info(audit)

# S4 method for class 'LeakAudit'
audit_info(audit)
```

## Arguments

- audit:

  A \[\`LeakAudit\`\] object returned by \[audit_leakage()\].

## Value

A named \`list\`.

## Details

Implemented as an S4 generic with a method for \[\`LeakAudit\`\];
visible via \`methods(class = "LeakAudit")\`.

## See also

\[LeakClasses\], \[audit_leakage()\]
