# Stationary bootstrap indices

Implements the stationary bootstrap of Politis & Romano (1994), which
resamples contiguous blocks of variable length to preserve weak temporal
dependence while maintaining ergodicity.

## Usage

``` r
.stationary_bootstrap(idx, mean_block)
```

## Arguments

- idx:

  Integer vector of ordered indices.

- mean_block:

  Positive numeric, expected block length.

## Value

Integer vector of permuted indices of the same length as \`idx\`.

## References

Politis, D. N., & Romano, J. P. (1994). \*The stationary bootstrap.\*
Journal of the American Statistical Association, 89(428), 1303-1313.
