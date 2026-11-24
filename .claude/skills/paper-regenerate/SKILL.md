---
name: paper-regenerate
description: Regenerate every number, table, figure and listing output of the bioLeak JSS paper from a frozen package version, and diff the new values against the stored ones. Use after the package is frozen for the revision, or whenever the manuscript's numbers must be refreshed.
---

# Regenerate the paper's numbers

Every number in the manuscript is transcribed from the outputs this skill
produces, and from nothing else (CLAUDE.md 7.5).

## Preconditions

- Package frozen: version bumped and committed by the user; `release-check`
  clean. Record the version and commit hash.
- Dependencies available, including `curatedOvarianData` for the case study
  (Bioconductor; ask the user before installing anything).

## Steps

1. Install the frozen package into a temporary library and run against it,
   not the CRAN copy (the replication script calls `library(bioLeak)`):

   ```
   R CMD INSTALL -l <scratch>/lib .
   R_LIBS=<scratch>/lib Rscript -e 'packageVersion("bioLeak")'
   ```

2. Copy only the scripts of the replication working copy
   (`paper/submission/replication/`, which carries the revision's script
   edits) to an empty scratch directory: no `sim_results/`, no `.rds`
   outputs. The script resumes from per-task checkpoints keyed on the task
   settings, not the package version, so copied results would be reused
   silently ("All tasks already checkpointed; skipping compute"). Run it
   there with `R_LIBS` set, so outputs never land in the paper folder. Keep
   the full log, run long sections in the background, and check the runtime
   against the script's stated estimate.
3. Check the run: no errors; warnings reviewed (each either expected and
   explained, or a tracker item); session info saved.
4. Diff old against new values:

   ```r
   old <- readRDS("paper/submission/bioLeak_replication_v3/placeholder_values.rds")
   new <- readRDS("<scratch>/placeholder_values.rds")
   ```

   Produce a table of every placeholder: old, new, changed (yes/no). Add the
   changed ones to the tracker under their items; an unexpected change
   (not explained by a tracker item) is investigated before going on.
5. Regenerate the listings' printed output from the frozen version and diff
   it against the manuscript's listings; record the differences.
6. Hand the new `placeholder_values.rds`, the tables Section 6 prints, the
   figures and the listing outputs to the user for the manuscript and
   replication folders.

## Notes

- Statistical results are checked, not just copied: a type-I rate above
  alpha + 3*sqrt(alpha(1-alpha)/R) or a power drop is a finding, not a
  number to print.
- If the run cannot reproduce a stored number and no tracker item explains
  it, stop and report.
