# Simulation benchmark matrix for leakage diagnostics

Runs a reproducible grid of simulation scenarios across modalities,
leakage mechanisms, and split modes using \[simulate_leakage_suite()\].
This function is designed as a benchmarking harness to quantify
detection rates and performance inflation under controlled settings.

## Usage

``` r
benchmark_leakage_suite(
  modalities = c("omics", "imaging_tabular", "ehr_tabular"),
  leakages = c("none", "subject_overlap", "batch_confounded", "peek_norm", "lookahead"),
  modes = c("subject_grouped", "batch_blocked", "time_series"),
  learner = c("glmnet", "ranger"),
  seeds = 1:5,
  B = 200,
  alpha = 0.05,
  parallel = FALSE
)
```

## Arguments

- modalities:

  Character vector selecting predefined modality profiles. Supported
  values: \`"omics"\`, \`"imaging_tabular"\`, \`"ehr_tabular"\`.

- leakages:

  Character vector of leakage mechanisms passed to
  \[simulate_leakage_suite()\].

- modes:

  Character vector of split modes passed to
  \[simulate_leakage_suite()\].

- learner:

  Character scalar. \`"glmnet"\` (default) or \`"ranger"\`.

- seeds:

  Integer vector of Monte Carlo seeds.

- B:

  Integer scalar. Number of permutations per scenario.

- alpha:

  Numeric scalar in (0, 1). Detection threshold applied to permutation
  p-values.

- parallel:

  Logical scalar. If TRUE, evaluates scenarios in parallel when
  \`future.apply\` is available.

## Value

A data.frame with one row per simulation seed/scenario and columns:
\`modality\`, \`leakage\`, \`mode\`, \`seed\`, observed metric, gap,
p-value, and a logical \`detected\` flag. A scenario-level summary is
attached as \`attr(x, "summary")\`.
