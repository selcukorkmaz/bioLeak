# Quantile break cache for permutation stratification

Internal environment used to cache quantile breakpoints for numeric
outcomes during restricted permutation testing. This avoids recomputing
quantiles across repeated calls in
[`audit_leakage()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md).

## Usage

``` r
.quantile_break_cache
```

## Format

An environment used to cache quantile breakpoints.

## Value

An environment (internal data object, not a function).
