# Check split overlap invariants

Verifies that a
[`LeakSplits`](https://selcukorkmaz.github.io/bioLeak/reference/LeakClasses.md)
object satisfies the expected no-overlap constraints for one or more
grouping columns. For each fold, the function checks that no group-level
value appearing in the test set is also present in the training set.

## Usage

``` r
check_split_overlap(splits, coldata = NULL, cols = NULL, stop_on_fail = TRUE)
```

## Arguments

- splits:

  A `LeakSplits` object from
  [`make_split_plan`](https://selcukorkmaz.github.io/bioLeak/reference/make_split_plan.md).

- coldata:

  A data.frame of sample metadata. When `NULL` (default), the function
  uses `splits@info$coldata` if available.

- cols:

  Character vector of column names to check for overlap. When `NULL`
  (default), the function infers columns from the split mode (e.g.,
  `group` for `subject_grouped`, `batch` for `batch_blocked`, both axes
  for `combined`).

- stop_on_fail:

  Logical; if `TRUE` (default), raises an error when any overlap is
  detected.

## Value

A data.frame with one row per fold-by-column combination and columns
`fold`, `repeat_id`, `col`, `n_overlap` (number of overlapping group
values), and `pass` (logical). Invisible. Raises an error if any fold
fails and `stop_on_fail = TRUE`.
