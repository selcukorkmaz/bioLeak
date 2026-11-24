---
name: paper-claims-check
description: Check every statement the bioLeak JSS manuscript makes about the package (behaviour, defaults, statistical properties, outputs, listings, numbers) against the current code and regenerated outputs, at a high bar. Use before calling any manuscript change done and before the revision is submitted.
---

# Check the manuscript against the package

Reviewer K's point applies throughout: the paper must describe only what
exists. A wrong claim in the paper is as serious as a bug.

## Round structure

Same discipline as `verify-change`: findings are reproduced before they are
reported; the round is repeated after fixes; stop after a round with no new
findings or after three fix cycles (then report to the user).

## Steps

1. Extract the claims. Go through the manuscript and list every statement
   about what a bioLeak function does, its defaults, the conditions it
   handles, the statistics it computes, what it prints, and every number.
   Include captions, tables, the appendix and the code listings.
2. Check each claim against the code (`R/`, help pages, NEWS development
   section). For behaviour, run a minimal example; do not rely on reading
   alone.
3. Check numbers against the regenerated outputs only (`paper-regenerate`).
4. Re-run every listing with the frozen version and compare the printed
   output line by line with the manuscript.
5. Check the replication script itself: comments that describe results,
   calls to internal functions, arguments that are ignored or deprecated.

## High bar

Report only: a statement now false that would mislead a reader or reviewer,
a claim stronger than its evidence (e.g. "eliminates" where the results show
a reduction, or a general claim from one setting), a listing that errors or
prints something different, a number that does not reproduce, or a feature
presented as available that is not. No wording
style. For each finding: line, what it says, what is true, evidence
(command + output). Each finding becomes a tracker row.

## Independence

For the pre-submission round, give a fresh subagent the manuscript, the
code and this skill, not the tracker's conclusions.
