---
name: verify-change
description: Verify that a bioLeak change, fix, or roadmap milestone actually works, in iterative rounds (claim-vs-diff review, tests and probes, adversarial probes, independent review) until a full round finds nothing new. Use when asked to check, verify, audit or review bioLeak work, or before calling any milestone "completed".
---

# Verify a bioLeak change in rounds

Goal: decide, with evidence, whether each claimed behaviour works. Passing
tests are not evidence on their own: the 0.4.0 work passed every test while
`delta_lsi(exchangeability = "by_group")` ran at 59% type-I error.

Work from the package root. Keep a findings table as you go:

| ID | Round found | Claim / area | Evidence (command + output) | Severity | Status |

Severity: **critical** (wrong results or invalid inference), **major**
(silent wrong behaviour in a supported path), **minor** (edge case, naming,
docs). Check `dev/OPEN-ISSUES.md` first so known issues get their existing ID.

## Round 1: claims vs code

1. List every claim being verified (from the roadmap, NEWS, PR text, or the
   user). One row per claim.
2. Read the diff for each: `git diff -- R/<file>.R`, and `git status` to see
   what is uncommitted.
3. For each claim, write the precise condition it must satisfy before reading
   the code closely: inputs, expected output, and what must *not* happen.
4. Note anything the code does that differs from the claim, and any branch
   that silently returns its input, falls back, or swallows a condition.

## Round 2: tests and existing probes

```
NOT_CRAN=true Rscript -e 'devtools::test(reporter = "summary")'
Rscript .claude/skills/verify-change/probes/inference_calibration.R
Rscript .claude/skills/verify-change/probes/guard_invariants.R
```

- Record failures and **new skips** (a skip is not a pass).
- Check that tests exist for each claim from Round 1, covering both the case
  that must act and the case that must not. A missing test is a finding.
- If the change touches `as_leaksplits()` or anything `split_spec`-related,
  also run the `contract-check` skill.

## Round 3: adversarial probes

Write a probe for each claim the existing probes don't cover. Put durable
ones in `probes/` (same PASS/FAIL/WARN format, non-zero exit on FAIL) so the
next round and the next session re-run them. Required attacks by area:

- **p-values / inference**: simulate under H0 (B ≥ 300) and check type-I
  error ≤ α + 3·sqrt(α(1−α)/B); count the truly independent units and check
  `R_eff` / `n_units` equals that count; check the tier reported matches.
  Include a positive control that must pass, so a broken probe is visible.
- **splits**: no group, batch or study level in both train and test of any
  fold; time ordering (`max(train time) < min(test time)` when ordered);
  every sample tested at least once where the mode promises it; fold sizes
  not degenerate (train fraction, empty folds).
- **guards / preprocessing**: unseen batch or factor level at transform
  time; batch confounded with the outcome; `transform()` of a row does not
  depend on other test rows; state hash changes when config changes;
  one-column and all-NA inputs.
- **diagnostics / audits**: one leaky fixture that must be flagged and one
  benign look-alike that must not.
- **every new argument**: default value, `NULL`, wrong type, and an
  interaction with another option (e.g. `compact = TRUE`, `nested = TRUE`,
  `repeats > 1`).

## Round 4: independent review

Spawn a fresh subagent (`general-purpose`) and give it the diff and file
paths, **not** your findings or conclusions. Ask it to find cases where the
code gives wrong results, invalid inference or silent no-ops, each with a
concrete failing input. Verify every finding it returns with a probe before
adding it to the table; discard ones you cannot reproduce, and say so.

## Iterate

- If any round produced a new finding and the user asked for fixes: fix it
  test-first (failing test, then fix), then re-run **from Round 2**, all
  rounds. A fix that passes its own test can still break another claim.
- If the user asked only for verification, don't fix; report.
- Stop when one complete pass of Rounds 2–4 finds nothing new.
- After 3 fix cycles that still surface new defects, stop and report to the
  user with the table; the design probably needs a decision, not more
  patches.

## Report

1. One-line verdict per claim: ✅ works / ⚠️ works with caveats / ❌ broken.
2. The findings table, most severe first, each with the exact command that
   reproduces it.
3. Rounds actually run, and how many cycles. Say plainly if a round was
   skipped and why.
4. Update `dev/OPEN-ISSUES.md`: add new IDs, mark fixed ones with the fix
   location, and never delete an entry without it being fixed or disproven.
