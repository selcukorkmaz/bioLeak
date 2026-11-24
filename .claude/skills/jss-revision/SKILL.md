---
name: jss-revision
description: Run the JSS revision of the bioLeak paper end to end - reviewer items, tracker, package changes, freeze, regeneration, manuscript edits and the response letter - under the rules in CLAUDE.md section 7. Use whenever work touches the reviewer comments, the manuscript, the replication script or the response letter.
---

# JSS revision workflow

Files:

| What | Where |
|---|---|
| Reviewer reports | `paper/submission/comments.rtf` (`textutil -convert txt -stdout`) |
| Manuscript | `paper/submission/manuscript/manuscript_jss.tex` (+ PDF) |
| Submitted replication (read-only baseline) | `paper/submission/bioLeak_replication_v3/replication.R`, `placeholder_values.rds` |
| Replication working copy (edits go here) | `paper/submission/replication/` (created from the v3 folder at the start of Phase 1) |
| Tracker (source of truth) | `dev/jss-revision/TRACKER.md` |
| Changes forced by the 2026-11 fixes | `dev/paper-updates.md` |
| Response letter draft | `dev/jss-revision/response-letter.md` |

The user may replace `paper/`; re-check the paths above at the start of each
session and update this table if they moved.

## Phase 0: tracker

1. Split every reviewer report into atomic requests. One request = one thing
   a reader could check was done. IDs: reviewer letter + number, sub-points
   with a dot (`L1.2`). Items from `dev/paper-updates.md` get `P` IDs.
2. One row per item:

   | ID | Reviewer's point (short quote) | Type | Proposed decision + reason | Action | Location | Status |

   Type: `package`, `replication`, `text`, `format`. Status lifecycle:
   `proposed` -> `agreed` (user confirmed) -> `in progress` -> `done`
   (verified, evidence linked) -> `in letter`.
3. Mark overlaps (e.g. K4 and P-ΔLSI) and merge their actions; keep both IDs.
4. Present the decisions to the user and record their answer before any
   work starts (CLAUDE.md 7.2).

## Phase 1: package items (type `package`)

1. Order by dependency: changes that alter outputs other items rely on go
   first (fold assignment, diagnostics), small items last.
2. For each item: write the precise condition the change must satisfy, the
   failing test first, then the fix. Follow CLAUDE.md sections 2-3 (null
   calibration + power for p-values; flagged fixture + benign look-alike for
   diagnostics; train-only; no silent no-ops).
3. Close the item with the `verify-change` rounds, a NEWS entry naming the
   item ID, and the evidence (test names, calibration numbers) in the tracker.
4. After all package items: one independent review round across them
   (fresh subagent, given the code and the tracker IDs, not the conclusions).

## Phase 2: freeze

Ask the user to bump the version and commit (CLAUDE.md 7.12). Run the
`release-check` skill. Nothing after this point changes package code; if a
later phase finds a package defect, it goes back to Phase 1 and the freeze
is repeated.

## Phase 3: replication edits, then regenerate

1. Apply all agreed `replication` items (new analyses, script fixes,
   comments) to the working copy `paper/submission/replication/`. These are
   script changes, not package changes, so they do not reopen the freeze.
2. Use the `paper-regenerate` skill. Record every changed placeholder
   (old -> new) in the tracker under the items it belongs to.

## Phase 4: manuscript

1. Go through the manuscript section by section. For each section, collect
   all tracker items whose Location points there; apply them in one pass.
2. Text about package behaviour is written from the code, not from memory;
   numbers are transcribed from the regenerated outputs and nowhere else
   (CLAUDE.md 7.5).
3. Write new text to CLAUDE.md 7.11: precise, compact, no overclaims;
   prefer replacing text over adding it.
4. Formatting items (decimal alignment, margins, figure/table choices) last.
5. Run `paper-claims-check`; every finding goes back into the tracker.

## Phase 5: response letter

For each item, in reviewer order: the point (quote or close paraphrase),
the response (agree / partly / respectfully disagree, with reasons), what
changed, where (section, page, line of the revised PDF), and the new text
or numbers quoted. Declines carry the argument and, where useful, evidence.
Start with a short summary of the main changes. Keep each response to what
was asked, in a few sentences (CLAUDE.md 7.8, 7.11). The user reviews and
owns the final letter.

## Before calling the revision ready

- Every tracker row is `in letter` (or explicitly deferred by the user).
- `paper-claims-check` round with no new findings.
- Manuscript compiles; listings and figures regenerated from the frozen
  version; the version in the manuscript matches the frozen package.
