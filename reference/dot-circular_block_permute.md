# Circular block permutation indices

Generates a permutation of time indices by concatenating random-length
blocks sampled circularly from the ordered sequence. Used for creating
block-permuted surrogates that preserve short-range temporal structure.

## Usage

``` r
.circular_block_permute(idx, block_len)
```

## Arguments

- idx:

  Integer vector of ordered indices.

- block_len:

  Positive integer block length (\>= 1).

## Value

Integer vector of permuted indices of the same length as \`idx\`.
