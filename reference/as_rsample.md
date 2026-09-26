# Convert LeakSplits to an rsample resample set

Convert LeakSplits to an rsample resample set

## Usage

``` r
as_rsample(x, data = NULL, ...)
```

## Arguments

- x:

  LeakSplits object created by \[make_split_plan()\].

- data:

  Optional data.frame used to populate rsample splits. When NULL, the
  stored \`coldata\` from \`x\` is used (if available).

- ...:

  Additional arguments passed to methods (unused).

## Value

An rsample `rset` object compatible with tidymodels workflows. The
returned object is a tibble with class `rset` containing:

- `splits`:

  List-column of `rsplit` objects, each with `analysis` (training
  indices) and `assessment` (test indices).

- `id`:

  Character column with fold identifiers (e.g., "Fold1").

- `id2`:

  Character column with repeat identifiers (e.g., "Repeat1") when
  multiple repeats are present; otherwise absent.

The object also carries attributes for `group`, `batch`, `study`, `time`
(when available from the original `LeakSplits`), and `bioLeak_mode`
indicating the original splitting mode. This allows the splits to be
used with
[`tune::tune_grid()`](https://tune.tidymodels.org/reference/tune_grid.html),
`rsample::fit_resamples()`, and other tidymodels functions.

## Examples

``` r
if (requireNamespace("rsample", quietly = TRUE)) {
  df <- data.frame(
    subject = rep(1:10, each = 2),
    outcome = rbinom(20, 1, 0.5),
    x1 = rnorm(20),
    x2 = rnorm(20)
  )
  splits <- make_split_plan(df, outcome = "outcome",
                        mode = "subject_grouped", group = "subject", v = 5)
  rset <- as_rsample(splits, data = df)
}
#> subject_grouped: repeat 1/1 done.
```
