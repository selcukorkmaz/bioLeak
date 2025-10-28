#' Summarize a leakage audit
#'
#' Provides a console summary of the leakage audit results,
#' including permutation gap, batch association, and duplicates.
#'
#' @seealso [plot.LeakAudit()], [bioLeak()], [make_splits()]
#'
#' @param object LeakAudit
#' @param digits number of digits for numeric display
#' @param ... ignored
#' @return Invisibly returns the input object.
#' @export
summary.LeakAudit <- function(object, digits = 3, ...) {
  if (!inherits(object, "LeakAudit"))
    stop("Object must be a 'LeakAudit'.", call. = FALSE)

  cat("\n==============================\n")
  cat(" bioLeak Leakage Audit Summary\n")
  cat("==============================\n\n")

  fit <- object@fit
  info <- fit@info
  warn_sym <- if (cli::is_utf8_output()) "⚠️" else "WARNING:"
  ok_sym   <- if (cli::is_utf8_output()) "✓" else "OK:"

  task <- fit@task %||% NA_character_
  outcome <- fit@outcome %||% NA_character_
  splits <- fit@splits
  mode <- if (!is.null(splits)) splits@mode %||% NA_character_ else NA_character_
  indices <- if (!is.null(splits)) splits@indices else NULL
  hash <- info$hash %||% "NA"

  cat(sprintf("Task: %s | Outcome: %s | Splitting mode: %s\n",
              task, outcome, mode))
  cat(sprintf("Hash: %s | Folds: %d | Repeats: %d\n\n",
              substr(hash, 1, 12),
              length(indices),
              info$repeats %||% 1))

  # --- Permutation gap ---
  if (!is.null(object@permutation_gap) && nrow(object@permutation_gap) > 0) {
    pg <- object@permutation_gap
    cat("Permutation Gap Test:\n")
    cat(sprintf("  Observed metric: %s\n",
                formatC(pg$metric, digits = digits, format = "f")))
    cat(sprintf("  Permuted mean ± SD: %s ± %s\n",
                formatC(pg$perm_mean, digits = digits, format = "f"),
                formatC(pg$perm_sd, digits = digits, format = "f")))
    cat(sprintf("  Gap: %s (larger gap → less leakage)\n\n",
                formatC(pg$gap, digits = digits, format = "f")))
  } else {
    cat("Permutation Gap Test: not available.\n\n")
  }

  # --- Batch association ---
  if (!is.null(object@batch_assoc) && nrow(object@batch_assoc) > 0) {
    ba <- object@batch_assoc
    cat("Batch / Study Association:\n")
    cat(sprintf("  χ² = %s (df = %s), p = %s\n\n",
                formatC(ba$stat, digits = digits, format = "f"),
                formatC(ba$df, digits = digits, format = "f"),
                formatC(ba$pval, digits = digits, format = "f")))
  } else {
    cat("Batch / Study Association: none detected.\n\n")
  }

  # --- Duplicate detection ---
  if (!is.null(object@duplicates) && nrow(object@duplicates) > 0) {
    dd <- object@duplicates
    cat("Near-Duplicate Samples:\n")
    cat(sprintf("  %d pairs detected above cosine ≥ %s\n",
                nrow(dd),
                formatC(object@info$duplicate_threshold, digits = digits, format = "f")))
    head_pairs <- utils::head(dd, 5)
    cat("  Example pairs:\n")
    print(head_pairs, row.names = FALSE)
    cat("\n")
  } else {
    cat("No near-duplicates detected.\n\n")
  }

  # --- Overall interpretation ---
  cat("Interpretation:\n")
  pg_available <- !is.null(object@permutation_gap) && nrow(object@permutation_gap) > 0

  if (!pg_available) {
    cat("  No permutation test results.\n")
  } else if (object@permutation_gap$gap < 0.01) {
    cat(sprintf("  %s Strong evidence of leakage (tiny permutation gap).\n", warn_sym))
  } else if (object@permutation_gap$gap < 0.05) {
    cat(sprintf("  %s Mild leakage signal detected.\n", warn_sym))
  } else {
    cat(sprintf("  %s No strong evidence of leakage.\n", ok_sym))
  }

  invisible(object)
}
