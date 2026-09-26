# Create leakage-resistant splits

Generates leakage-safe cross-validation splits for common biomedical
setups: subject-grouped, batch-blocked, study leave-one-out, and
time-series rolling-origin. Supports repeats, optional stratification,
nested inner CV, and optional prediction horizon/purge/embargo gaps for
time series. Note that splits store explicit indices, which can be
memory-intensive for large `n` and many repeats.

## Usage

``` r
make_split_plan(
  x,
  outcome = NULL,
  mode = c("subject_grouped", "batch_blocked", "study_loocv", "time_series", "combined"),
  group = NULL,
  batch = NULL,
  study = NULL,
  time = NULL,
  primary_axis = NULL,
  secondary_axis = NULL,
  constraints = NULL,
  v = 5,
  repeats = 1,
  stratify = FALSE,
  nested = FALSE,
  seed = 1,
  horizon = 0,
  purge = 0,
  embargo = 0,
  progress = TRUE,
  compact = FALSE,
  strict = TRUE
)
```

## Arguments

- x:

  SummarizedExperiment or data.frame/matrix (samples x features). If
  SummarizedExperiment, metadata are taken from colData(x). If
  data.frame, metadata are taken from x (columns referenced by `group`,
  `batch`, `study`, `time`, `outcome`).

- outcome:

  character, outcome column name (used for stratification).

- mode:

  one of
  "subject_grouped","batch_blocked","study_loocv","time_series","combined".

- group:

  subject/group id column (for subject_grouped). Required when mode is
  \`subject_grouped\`; use \`group = "row_id"\` to explicitly request
  sample-wise CV.

- batch:

  batch/plate/center column (for batch_blocked).

- study:

  study id column (for study_loocv).

- time:

  time column (numeric or POSIXct) for time_series.

- primary_axis:

  List with elements `type` (one of `"subject"`, `"batch"`, `"study"`)
  and `col` (column name). Used only when `mode = "combined"` to define
  the primary grouping axis. Deprecated in favor of `constraints`; still
  supported for backward compatibility.

- secondary_axis:

  List with elements `type` and `col`. Used only when
  `mode = "combined"` to define the secondary constraint axis. Training
  sets exclude samples whose secondary-axis levels appear in the test
  set. Deprecated in favor of `constraints`; still supported for
  backward compatibility.

- constraints:

  A list of constraint specifications for `mode = "combined"`. Each
  element is a list with `type` (one of `"subject"`, `"batch"`,
  `"study"`) and `col` (column name). The first element defines the
  primary grouping axis (fold driver); subsequent elements define
  exclusion constraints (training samples sharing constraint-axis levels
  with the test set are removed). Requires at least 2 elements. Cannot
  be used together with `primary_axis`/`secondary_axis`.

- v:

  integer, number of folds (k) or rolling partitions.

- repeats:

  integer, number of repeats (\>=1) for non-LOOCV modes.

- stratify:

  logical, keep outcome proportions similar across folds. For grouped
  modes, stratification is applied at the group level (by majority class
  per group) if `outcome` is provided; otherwise ignored.

- nested:

  logical, whether to attach inner CV splits (per outer fold) using the
  same `mode` on the outer training set (with `v` folds, 1 repeat).

- seed:

  integer seed.

- horizon:

  numeric (\>=0), minimal time gap for time_series so that the training
  set only contains samples with time \< min(test_time) when horizon =
  0, and time \<= min(test_time) - horizon otherwise.

- purge:

  numeric (\>=0), additional gap removed immediately before each
  time-series test block.

- embargo:

  numeric (\>=0), additional exclusion window anchored at the end of
  each time-series test block. Training rows with
  `time > max(test_time) - embargo` are removed.

- progress:

  logical, print progress for large jobs.

- compact:

  logical; store fold assignments instead of explicit train/test indices
  to reduce memory usage for large datasets. Not supported when
  `nested = TRUE`.

- strict:

  logical; deprecated and ignored. \`subject_grouped\` always requires a
  non-NULL \`group\`.

## Value

A
[`LeakSplits`](https://selcukorkmaz.github.io/bioLeak/reference/LeakClasses.md)
S4 object containing:

- `mode`:

  Character string indicating the splitting mode (`"subject_grouped"`,
  `"batch_blocked"`, `"study_loocv"`, or `"time_series"`).

- `indices`:

  List of fold descriptors, each containing `train` (integer vector of
  training indices), `test` (integer vector of test indices), `fold`
  (fold number), and `repeat_id` (repeat identifier). When
  `compact = TRUE`, indices are stored as fold assignments instead.

- `info`:

  List of metadata including `outcome`, `v`, `repeats`, `seed`, grouping
  columns (`group`, `batch`, `study`, `time`), `stratify`, `nested`,
  `horizon`, `purge`, `embargo`, `summary` (data.frame of fold sizes),
  `hash` (reproducibility checksum), `inner` (nested inner splits if
  `nested = TRUE`), and `coldata` (sample metadata).

Use the `show` method to print a summary; downstream access to the
indices and metadata is normally done through the functions that consume
a `LeakSplits` (for example
[`fit_resample`](https://selcukorkmaz.github.io/bioLeak/reference/fit_resample.md))
rather than by reading slots directly.

## Examples

``` r
set.seed(1)
df <- data.frame(
  subject = rep(1:10, each = 2),
  outcome = rbinom(20, 1, 0.5),
  x1 = rnorm(20),
  x2 = rnorm(20)
)
splits <- make_split_plan(df, outcome = "outcome",
                      mode = "subject_grouped", group = "subject", v = 5)
#> subject_grouped: repeat 1/1 done.
```
