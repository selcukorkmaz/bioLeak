# Convert a splitGraph split_spec into bioLeak splits

Consume a \`split_spec\` produced by splitGraph and build a
corresponding `LeakSplits` object via
[`make_split_plan`](https://selcukorkmaz.github.io/bioLeak/reference/make_split_plan.md).
The spec supplies the grouping/blocking/ordering assignments; the caller
supplies the observation frame (features + outcome), joined on sample
id.

## Usage

``` r
as_leaksplits(spec, data, outcome, sample_id_col = "sample_id", v = 5, ...)
```

## Arguments

- spec:

  A `split_spec` object from splitGraph.

- data:

  A data.frame (or SummarizedExperiment) containing at least one
  identifier column matching `sample_id_col` and an `outcome` column.

- outcome:

  Name of the outcome column in `data`.

- sample_id_col:

  Name of the sample-id column in `data` (default `"sample_id"`).

- v:

  Number of CV folds to request from
  [`make_split_plan()`](https://selcukorkmaz.github.io/bioLeak/reference/make_split_plan.md).

- ...:

  Additional arguments forwarded to
  [`make_split_plan`](https://selcukorkmaz.github.io/bioLeak/reference/make_split_plan.md)
  (e.g. `stratify`, `seed`, `horizon`, `purge`).

## Value

A `LeakSplits` object.

## Details

The mapping from `spec$constraint_mode` to `make_split_plan(mode=)` is:

- `"subject"` -\> `"subject_grouped"`

- `"batch"` -\> `"batch_blocked"`

- `"study"` -\> `"study_loocv"`

- `"time"` -\> `"time_series"`

- `"site"`, `"region"`, `"platform"`, `"assay"`, `"relatedness"`,
  `"spatial"`, `"composite"` -\> `"subject_grouped"` with
  `group = spec$group_var` (`"group_id"`)

For the grouping modes in the last row, splitGraph has already resolved
the dependency structure into one group id per sample (for
`"composite"`, the merged components of all relations in `via`), so
keeping each group intact within a fold honours the constraint.
Composite specs that require ordering (`spec$ordering_required`) cannot
be expressed as grouped folds and raise an error; derive a `"time"` spec
instead. Unrecognised future modes fall back to grouped folds on
`group_var` with a warning.

Blocking variables declared on the spec (`batch_group`, `study_group`)
and ordering (`order_rank`) are forwarded automatically when relevant.

## See also

[`make_split_plan`](https://selcukorkmaz.github.io/bioLeak/reference/make_split_plan.md),
[`as_rsample`](https://selcukorkmaz.github.io/bioLeak/reference/as_rsample.md)
