# Display summary for LeakSplits objects

Prints fold counts, sizes, and hash metadata for quick inspection.

## Usage

``` r
# S4 method for class 'LeakSplits'
show(object)
```

## Arguments

- object:

  LeakSplits object.

## Value

No return value, called for side effects (prints a summary to the
console showing mode, fold count, repeats, outcome, stratification
status, nested status, per-fold train/test sizes, and the
reproducibility hash).

## Examples

``` r
df <- data.frame(
  subject = rep(1:10, each = 2),
  outcome = rbinom(20, 1, 0.5),
  x1 = rnorm(20),
  x2 = rnorm(20)
)
splits <- make_split_plan(df, outcome = "outcome",
                      mode = "subject_grouped", group = "subject", v = 5)
#> subject_grouped: repeat 1/1 done.
show(splits)
#> LeakSplits object (mode = subject_grouped, v = 5, repeats = 1)
#> Outcome: outcome | Stratified: FALSE | Nested: FALSE
#> ------------------------------------------------------
#>   fold repeat_id train_n test_n
#> 1    1         1      16      4
#> 2    2         1      16      4
#> 3    3         1      16      4
#> 4    4         1      16      4
#> 5    5         1      16      4
#> ------------------------------------------------------
#> Total folds: 5 | Hash: 1afe5ebbc6914f8dad7d7006ea44a76f
```
