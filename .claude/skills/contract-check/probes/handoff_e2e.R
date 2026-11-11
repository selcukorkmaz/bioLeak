# End-to-end splitGraph -> bioLeak handoff invariants.
# Run from the bioLeak package root (siblings: ../splitGraph):
#   Rscript .claude/skills/contract-check/probes/handoff_e2e.R
# Loads both working trees with pkgload, so it tests the code on disk.

source(".claude/skills/verify-change/probes/_probe_helpers.R")
suppressMessages({
  pkgload::load_all("../splitGraph", quiet = TRUE)
  pkgload::load_all(".", quiet = TRUE)
})

set.seed(4)
grid <- expand.grid(subj = 1:30, visit = 1:4)
meta <- data.frame(
  sample_id    = sprintf("S%03d", seq_len(nrow(grid))),
  subject_id   = paste0("P", grid$subj),
  timepoint_id = paste0("T", grid$visit),
  time_index   = grid$visit,
  batch_id     = paste0("B", (grid$subj %% 5) + 1),
  study_id     = paste0("ST", (grid$subj %% 3) + 1),
  stringsAsFactors = FALSE
)
g <- graph_from_metadata(meta, time_precedence = TRUE)
dat <- data.frame(sample_id = meta$sample_id,
                  y = factor(rbinom(nrow(meta), 1, 0.5)),
                  x1 = rnorm(nrow(meta)))

fold_checks <- function(label, spec, splits, block_cols = character()) {
  sd <- spec$sample_data
  idx <- match(dat$sample_id, sd$sample_id)
  gid <- sd$group_id[idx]
  ord <- if ("order_rank" %in% names(sd)) sd$order_rank[idx] else NULL
  folds <- splits@indices
  check("CTRL-2", length(folds) > 0L, sprintf("%s: %d folds produced", label, length(folds)))
  if (!length(folds) || is.null(folds[[1]]$train)) return(invisible())
  overlap <- vapply(folds, function(f) length(intersect(gid[f$train], gid[f$test])), integer(1))
  check("INV-4", all(overlap == 0L), sprintf("%s: group_id overlap per fold = %s", label, paste(overlap, collapse = ",")))
  for (bc in block_cols) {
    bv <- sd[[bc]][idx]
    ov <- vapply(folds, function(f) length(intersect(stats::na.omit(bv[f$train]), stats::na.omit(bv[f$test]))), integer(1))
    check("INV-5", all(ov == 0L), sprintf("%s: %s overlap per fold = %s", label, bc, paste(ov, collapse = ",")))
  }
  if (isTRUE(spec$ordering_required) && !is.null(ord)) {
    ok <- vapply(folds, function(f) max(ord[f$train]) < min(ord[f$test]), logical(1))
    check("INV-6", all(ok), sprintf("%s: train strictly precedes test in every fold", label))
  }
  frac <- vapply(folds, function(f) length(f$train) / (nrow(dat) - length(f$test)), numeric(1))
  check("BL-5", min(frac) >= 0.4,
        sprintf("%s: min training fraction of available rows = %.2f", label, min(frac)), level = "WARN")
}

run <- function(label, con, ..., block_cols = character()) {
  spec <- as_split_spec(con, graph = g)
  r <- capture(as_leaksplits(spec, dat, outcome = "y", v = 4, ...))
  if (!is.null(r$error)) {
    check("INV-4", FALSE, sprintf("%s: as_leaksplits() error: %s", label, r$error))
    return(invisible())
  }
  fold_checks(label, spec, r$value, block_cols)
}

run("subject", derive_split_constraints(g, mode = "subject"))
run("batch",   derive_split_constraints(g, mode = "batch"))
run("study",   derive_split_constraints(g, mode = "study"))
run("subject+time", derive_split_constraints(g, mode = "subject", time_ordering = TRUE))
run("composite strict", derive_split_constraints(g, mode = "composite", strategy = "strict",
                                                 via = c("subject", "batch")))
mc <- derive_split_constraints(g, mode = "composite", strategy = "min_cut",
                               via = c("subject", "batch", "study"), k = 4)
spec_mc <- as_split_spec(mc, graph = g)
run("composite min_cut", mc)

# SG-4: min_cut severs secondary bridges on purpose; the spec must not then
# declare those same axes as hard block_vars, or use_block_constraints = TRUE
# re-imposes what min_cut cut (observed: "No valid folds generated").
severed_axes <- unique(sub(":.*$", "", as.character(mc$metadata$severed_bridges)))
blocked_axes <- sub("_group$", "", as.character(spec_mc$block_vars))
clash <- intersect(severed_axes, blocked_axes)
check("SG-4", length(clash) == 0L,
      sprintf("min_cut severed axes (%s) vs spec block_vars (%s): conflicting = %s",
              paste(severed_axes, collapse = ","), paste(spec_mc$block_vars, collapse = ","),
              if (length(clash)) paste(clash, collapse = ",") else "none"))
r <- capture(as_leaksplits(spec_mc, dat, outcome = "y", v = 4, use_block_constraints = TRUE))
check("SG-4", is.null(r$error),
      sprintf("min_cut spec + use_block_constraints = TRUE: %s",
              if (is.null(r$error)) "folds produced" else paste("error:", r$error)), level = "WARN")

# stratum_var forwarding
spec_s <- as_split_spec(derive_split_constraints(g, mode = "subject"), graph = g)
if (!is.null(spec_s$stratum_var)) {
  r <- capture(as_leaksplits(spec_s, dat, outcome = "y", v = 4))
  check("INV-7", identical(r$value@info$stratify_by, spec_s$stratum_var),
        "stratum_var forwarded to splits@info$stratify_by")
} else {
  cat("SKIP [INV-7] spec has no stratum_var in this fixture\n")
}

probe_finish()
