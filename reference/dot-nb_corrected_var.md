# Nadeau-Bengio corrected variance

Computes the corrected variance for repeated K-fold CV as described in
Nadeau & Bengio (2003). The correction factor accounts for the overlap
between training sets across folds.

## Usage

``` r
.nb_corrected_var(per_fold_vals, n_train = NULL, n_test = NULL)
```

## Arguments

- per_fold_vals:

  Numeric vector of per-fold metric values.

- n_train:

  Average number of training samples per fold. NULL to skip correction.

- n_test:

  Average number of test samples per fold. NULL to skip correction.

## Value

Scalar corrected variance, or `NA_real_` for single-fold input.
