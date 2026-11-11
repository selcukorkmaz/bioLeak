# Convert a splitGraph split_spec into bioLeak splits

Consume a \`split_spec\` produced by splitGraph and build a
corresponding `LeakSplits` object via
[`make_split_plan`](https://selcukorkmaz.github.io/bioLeak/reference/make_split_plan.md).
The spec supplies the grouping/blocking/ordering assignments; the caller
supplies the observation frame (features + outcome), joined on sample
id.

## Usage

``` r
as_leaksplits(
  spec,
  data,
  outcome,
  sample_id_col = "sample_id",
  v = 5,
  use_block_constraints = FALSE,
  ...
)
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

- use_block_constraints:

  Logical (default `FALSE`). When `TRUE` and the spec uses a grouped
  mode with non-empty `spec$block_vars`, converts `group_var` and
  `block_vars` into multi-axis `constraints` in
  `make_split_plan(mode = "combined")`.

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
keeping each group intact within a fold honours the constraint. Specs
that require ordering (`spec$ordering_required = TRUE`) and supply a
valid ordering column (`"order_rank"` or `spec$time_var`) map to
`"time_series"` while also passing `group = spec$group_var` (for
non-`"time"` modes) so both temporal ordering and group disjointness are
enforced. Unrecognised future modes fall back to grouped folds on
`group_var` with a warning.

Blocking variables declared on the spec (`batch_group`, `study_group`,
`site_group`, `region_group`, `platform_group`, `assay_group`),
stratification (`stratum` or `spec$stratum_var`), and ordering
(`order_rank`) are forwarded automatically when relevant. When
`use_block_constraints = TRUE`, non-missing `spec$block_vars` on a
grouped spec are enforced as secondary exclusion axes via
`make_split_plan(mode = "combined")`.

## See also

[`make_split_plan`](https://selcukorkmaz.github.io/bioLeak/reference/make_split_plan.md),
[`as_rsample`](https://selcukorkmaz.github.io/bioLeak/reference/as_rsample.md)
