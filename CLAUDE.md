# bioLeak — rules for working on this package

bioLeak executes leakage-safe resampling, train-only preprocessing
(guards), leakage audits and randomization inference. Its sibling
package `../splitGraph` declares dataset dependency structure and emits
a `split_spec`; bioLeak consumes it through
[`as_leaksplits()`](https://selcukorkmaz.github.io/bioLeak/reference/as_leaksplits.md).

Open defects found in review live in `dev/OPEN-ISSUES.md`. Read it
before touching `delta_lsi.R`, `guards.R`, `make_split_plan.R` or
`splitgraph_adapter.R`, and update it when an issue is fixed or found.

## 1. Package boundary

- Fold execution, guards, audits and inference belong here. Graph
  modelling, dependency derivation and leakage *declaration* belong in
  splitGraph.
- Do not make splitGraph a hard dependency. It stays in `Suggests`;
  every use is behind
  [`requireNamespace("splitGraph", quietly = TRUE)`](https://rdrr.io/r/base/ns-load.html)
  or `skip_if_not_installed("splitGraph")`.
- A new `Imports` entry needs the user’s approval. Prefer `Suggests`
  plus a runtime check.

## 2. What “done” means

A milestone or fix is done only when all of these hold. Until then, call
it “implemented, unverified”, never “completed”.

1.  Unit tests cover the case where the feature must act, the case where
    it must not act, and at least one edge case (empty, single group,
    NA, unseen level).
2.  The relevant probe in `.claude/skills/verify-change/probes/` passes,
    or a new probe was added for behaviour the probes don’t yet cover.
3.  `NOT_CRAN=true Rscript -e 'devtools::test()'` reports 0 failures and
    no new skips.
4.  `Rscript -e 'roxygen2::roxygenise()'` leaves no diff in `man/` or
    `NAMESPACE`. Never hand-edit `man/*.Rd` or `NAMESPACE`.
5.  `NEWS.md` has an entry under the development version. Behaviour
    changes that can alter results go under “Bug fixes (results can
    change)” or a “Breaking changes” heading.
6.  The change is committed. Uncommitted work is not done, whatever a
    roadmap says.

### Iterative checking (required)

Verification runs in rounds, as defined in the `verify-change` skill: R1
claim-vs-diff review → R2 tests + probes → R3 adversarial probes (null
calibration, benign look-alikes, unseen levels, edge cases) → R4
independent review by a fresh subagent given the code, not the
conclusions. Any finding is fixed, then **every** round is re-run from
R2, because fixes break other things. Stop only after a full round finds
nothing new. After three fix cycles that still find new defects, stop
and report to the user instead of continuing. Each round’s result goes
in the report; a round is never skipped silently.

## 3. Statistical rules (these caused real defects; do not relax them)

- **Inference units must be independent.** Never assign one fold-level
  statistic to several units (groups, batches, cells) and then treat
  them as exchangeable. `R_eff` and every `n_units` must equal the
  number of *independent* units. If a per-unit statistic is needed,
  compute it from that unit’s own predictions.
- **Every new or changed p-value needs a null calibration check**:
  simulate under H0, require type-I error ≤ α + 3·sqrt(α(1−α)/B). Add a
  power sanity check under a clear alternative. Correct-looking code has
  shipped at 59% type-I error at α = 0.05.
- **No silent no-ops.** A learned transform must define what happens on
  levels unseen in training (batch, factor level, column). It either
  errors or warns and applies a documented fallback. The same holds for
  any algorithm that cannot reach its target: say so.
- **Train-only means train-only.** A guard’s
  [`transform()`](https://rdrr.io/r/base/transform.html) output for a
  row must not depend on other test rows. Fitted state comes only from
  training rows, and its hash must change when the configuration
  changes.
- **Removing a group effect can remove biology.** Batch/location-scale
  adjustment must accept protected covariates, or warn when the batch is
  associated with the outcome in the training fold.
- **Diagnostics need both error directions tested**: one fixture that
  must be flagged, and one benign look-alike that must not be flagged.
- Method names must match the method. Don’t call something “ComBat”
  unless it implements ComBat’s model; document deviations in
  `@details`.

## 4. The splitGraph ↔︎ bioLeak seam

- Fields read from `split_spec`: `sample_data` (`sample_id`, `group_id`,
  `batch_group`, `study_group`, `site_group`, `region_group`,
  `platform_group`, `assay_group`, `timepoint_id`, `order_rank`,
  `stratum`), `group_var`, `block_vars`, `time_var`, `stratum_var`,
  `ordering_required`, `constraint_mode`, `constraint_strategy`.
- Reading a new field, or changing how one is interpreted, means:
  1.  update splitGraph’s `tests/testthat/test-bioleak-contract.R` and
      `?as_split_spec` (“What downstream consumers read”);
  2.  check whether splitGraph’s `.depgraph_schema_version` must bump
      (rules in `../splitGraph/dev/release.md` §1);
  3.  run the `contract-check` skill.
- Contract tests use `skip_if_not_installed()`, which sees the
  **installed** copy, not the working tree. Install both development
  versions before trusting a contract test result (the `contract-check`
  skill does this).

## 5. Commands

    NOT_CRAN=true Rscript -e 'devtools::test()'                     # full suite
    Rscript -e 'devtools::test(filter = "delta-lsi")'               # one file
    Rscript -e 'roxygen2::roxygenise()'                             # docs + NAMESPACE
    Rscript .claude/skills/verify-change/probes/inference_calibration.R
    Rscript .claude/skills/verify-change/probes/guard_invariants.R
    R CMD build . && R CMD check --as-cran bioLeak_<ver>.tar.gz     # before release

Run probes from the package root. They exit non-zero when any FAIL line
is printed.

## 6. Working conventions

- Fix bugs test-first: write the failing test or probe, see it fail,
  then fix.
- Keep internal helpers dot-prefixed (`.dlsi_*`, `.guard_*`, `.bio_*`).
  Raise conditions through `.bio_stop()` / `.bio_warn()` with a
  condition class.
- Never change [`set.seed()`](https://rdrr.io/r/base/Random.html)
  behaviour of user-visible functions without a NEWS entry. New internal
  code that needs randomness saves and restores `.Random.seed` instead
  of leaving the user’s global RNG state changed.
- Commit or push only when the user asks. This repo often carries large
  uncommitted work; never `git checkout`, `reset` or `stash` it away.
