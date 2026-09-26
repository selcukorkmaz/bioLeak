# Restricted permutation label factory

Builds a closure that generates permuted outcome vectors per fold while
respecting grouping/batch/study/time constraints used in
[`audit_leakage()`](https://selcukorkmaz.github.io/bioLeak/reference/audit_leakage.md).
Numeric outcomes can be stratified by quantiles to preserve outcome
structure under permutation.

## Usage

``` r
.permute_labels_factory(
  cd,
  outcome,
  mode,
  folds,
  perm_stratify,
  time_block,
  block_len,
  seed,
  group_col = NULL,
  batch_col = NULL,
  study_col = NULL,
  time_col = NULL,
  perm_refit = TRUE,
  verbose = FALSE
)
```

## Arguments

- cd:

  data.frame of sample metadata.

- outcome:

  outcome column name.

- mode:

  resampling mode (subject_grouped, batch_blocked, study_loocv,
  time_series).

- folds:

  list of fold descriptors from `LeakSplits`. When compact splits are
  used, fold assignments are read from the `fold_assignments` attribute.

- perm_stratify:

  logical or "auto"; if TRUE, permute within strata.

- time_block:

  time-series block permutation method.

- block_len:

  block length for time-series permutations.

- seed:

  integer seed.

- group_col, batch_col, study_col:

  optional metadata columns.

- time_col:

  optional metadata column name for time-series ordering.

- perm_refit:

  logical; if TRUE model is retrained on permuted labels (block
  permutation preserves subject structure); if FALSE predictions are
  fixed and simple label shuffle is used for `subject_grouped` mode.

- verbose:

  logical; print progress messages.

## Value

A function that returns a list of permuted outcome vectors, one per
fold.
