---
name: release-check
description: Pre-release gate for bioLeak (version bump, NEWS, docs, full tests, probes, contract check, R CMD check --as-cran). Use before tagging, submitting to CRAN, or calling a version "released", and when asked whether bioLeak is ready to ship.
---

# bioLeak release gate

Run from the bioLeak root. Report each item as ✅ / ❌ with the evidence.
Anything ❌ blocks the release; say so plainly rather than softening it.

## 0. Preconditions

- `dev/OPEN-ISSUES.md` has no open **critical** entry. Open major entries
  are listed in the report for the user to accept or fix.
- `verify-change` has been run on everything in this release, ending in a
  clean round, and `contract-check` passed against the splitGraph version
  this release targets.
- `git status` is clean apart from intended release edits; the release work
  is committed.

## 1. Version and changelog

- [ ] `DESCRIPTION` `Version:` is the release version (no `.9000`).
- [ ] `NEWS.md` top section is that version; every user-visible change in
      `git diff <last-tag>..HEAD -- R/` is mentioned. Result-changing fixes
      sit under "Bug fixes (results can change)". New arguments are named.
- [ ] Minimum splitGraph version the adapter needs is stated in NEWS and
      in `?as_leaksplits` if the seam changed.
- [ ] `README.md` version/year references are current (and `inst/CITATION`
      if one has been added).

## 2. Docs and code hygiene

- [ ] `Rscript -e 'roxygen2::roxygenise()'` leaves no diff.
- [ ] Every new exported argument has `@param`; every new method name
      describes what it computes (see CLAUDE.md §3).
- [ ] No `browser()`, `print(` debug lines, or `TODO`s in `R/` added since
      the last tag: `git diff <last-tag>..HEAD -- R/ | grep -nE '^\+.*(browser\(|TODO|FIXME)'`.

## 3. Tests and probes

- [ ] `NOT_CRAN=true Rscript -e 'devtools::test()'` → 0 failures; list skips
      and confirm each is intentional.
- [ ] `Rscript .claude/skills/verify-change/probes/inference_calibration.R`
      → 0 FAIL (use `PROBE_B=1000` for the release run).
- [ ] `Rscript .claude/skills/verify-change/probes/guard_invariants.R` → 0 FAIL.
- [ ] `contract-check` skill → 0 FAIL, contract tests ran (not skipped).
- [ ] Vignettes run: `Rscript -e 'for (v in list.files("vignettes", "[.]Rmd$", full.names=TRUE)) { f <- tempfile(fileext=".R"); knitr::purl(v, f, quiet=TRUE, documentation=0); source(f) }'`.

## 4. R CMD check

- [ ] `R CMD build .` then `R CMD check --as-cran bioLeak_<ver>.tar.gz`
      with all Suggests installed (including splitGraph), so nothing is
      skipped. 0 errors, 0 warnings; every NOTE explained in the report.
- [ ] The tarball excludes `CLAUDE.md`, `.claude/`, `dev/`, `paper/`,
      `docs/`: `tar -tzf bioLeak_<ver>.tar.gz | grep -E 'CLAUDE|\.claude|/dev/|/paper/'`
      prints nothing.
- [ ] `urlchecker::url_check()` clean.

## 5. After the gate

Do not tag, push, or submit. Report the checklist and ask the user; tagging
and CRAN submission are the user's actions.
