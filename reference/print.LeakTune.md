# Print a LeakTune object

Brief one-screen auto-print representation of a \`LeakTune\` result
returned by \[tune_resample()\]. Use \[summary()\] for the full
diagnostic report (outer-loop metrics, selected hyperparameters,
fold-by-fold detail, and refit summary).

## Usage

``` r
# S3 method for class 'LeakTune'
print(x, ...)
```

## Arguments

- x:

  A \`LeakTune\` object returned by \[tune_resample()\].

- ...:

  Ignored; present so that the S3 signature matches \[base::print()\].

## Value

Invisibly returns \`x\`.
