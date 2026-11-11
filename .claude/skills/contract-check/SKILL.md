---
name: contract-check
description: Check the splitGraph -> bioLeak handoff end to end (split_spec -> as_leaksplits -> folds) against the development versions of both packages. Use whenever as_leaksplits(), split_spec fields, block_vars, ordering, stratification, or splitGraph constraint strategies change, and before releasing either package.
---

# Check the splitGraph ↔ bioLeak contract

Both packages live side by side: `../splitGraph` (relative to bioLeak) and
`../bioLeak` (relative to splitGraph). Run commands from the bioLeak root.

## Why installing matters

Contract tests on both sides use `skip_if_not_installed()`, which loads the
**installed** copy. A stale installed bioLeak or splitGraph makes contract
tests pass or skip while the working trees are broken. So:

1. Record what is installed now:
   `Rscript -e 'for (p in c("splitGraph","bioLeak")) cat(p, as.character(tryCatch(packageVersion(p), error=function(e) "none")), "\n")'`
2. Install both working trees into a temporary library, never over the
   user's library without asking:
   ```
   LIB=$(mktemp -d)
   R CMD INSTALL --no-docs --no-multiarch -l "$LIB" ../splitGraph
   R CMD INSTALL --no-docs --no-multiarch -l "$LIB" .
   ```
3. Run with that library first on the path:
   ```
   R_LIBS="$LIB" NOT_CRAN=true Rscript -e 'testthat::test_local("../splitGraph", filter = "bioleak-contract")'
   R_LIBS="$LIB" NOT_CRAN=true Rscript -e 'testthat::test_local(".", filter = "as-leaksplits")'
   R_LIBS="$LIB" Rscript .claude/skills/contract-check/probes/handoff_e2e.R
   ```
   Confirm from the output that the contract tests **ran** (no "skipped:
   bioLeak not installed").

## What the e2e probe asserts

For specs derived from one synthetic graph (subject, batch, study, composite
strict, composite min_cut, subject + time ordering), every fold produced by
`as_leaksplits()` must satisfy:

- no `group_id` in both train and test;
- with `use_block_constraints = TRUE`: no `block_vars` level in both;
- when `ordering_required`: `max(train order) < min(test order)`;
- `stratum_var` forwarded (`splits@info$stratify_by`);
- training fraction not degenerate (WARN below 40% of the available data).

## Rounds

This skill is Round 2 for seam changes in `verify-change`. If any assertion
fails: record it in `dev/OPEN-ISSUES.md` of the package that owns the fault
(adapter logic → bioLeak, spec content → splitGraph), fix test-first if asked,
and re-run this skill plus both packages' full suites. Stop after a clean run.

## Also check by reading

- splitGraph `?as_split_spec` "What downstream consumers read" matches the
  fields `R/splitgraph_adapter.R` reads.
- splitGraph's `.depgraph_schema_version` was bumped if a spec field or node
  attribute was added (rule in `../splitGraph/dev/release.md` §1).
- Both NEWS files mention a seam change, each from its own side.
