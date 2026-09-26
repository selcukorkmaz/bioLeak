# Summarize a nested tuning result

Prints a concise report for a \`LeakTune\` object produced by
\[tune_resample()\]. The report highlights the tuning strategy,
selection metric, and cross-validated performance across outer folds,
plus a glimpse of the selected hyperparameters.

## Usage

``` r
# S3 method for class 'LeakTune'
summary(object, digits = 3, ...)
```

## Arguments

- object:

  A \[LeakTune\] object returned by \[tune_resample()\].

- digits:

  Integer scalar. Number of decimal places to print in numeric summary
  tables. Defaults to 3.

- ...:

  Unused. Included for S3 method compatibility.

## Value

Invisibly returns \`object\$metric_summary\`, the data frame of
per-learner metric means and standard deviations computed across outer
folds.
