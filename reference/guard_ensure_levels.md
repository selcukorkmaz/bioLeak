# Ensure consistent categorical levels for guarded preprocessing

Converts character/logical columns to factors and aligns factor levels
with a training-time `levels_map`. Adds a dummy level when a column has
only one observed level so that downstream one-hot encoding retains a
column.

## Usage

``` r
guard_ensure_levels(df, levels_map = NULL, dummy_prefix = "__dummy__")
```

## Arguments

- df:

  data.frame to normalize factor levels.

- levels_map:

  optional named list of factor levels learned from training data.

- dummy_prefix:

  prefix used when adding a dummy level to single-level factors.

## Value

List with elements `data` (data.frame) and `levels` (named list of
levels).
