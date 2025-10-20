#' Summarize a leakage audit
#'
#' Provides a console summary of the leakage audit results,
#' including permutation gap, batch association, and duplicates.
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
  cat(sprintf("Task: %s | Outcome: %s | Splitting mode: %s\n",
              fit@task, fit@outcome, fit@splits@mode))
  cat(sprintf("Hash: %s | Folds: %d | Repeats: %d\n\n",
              substr(info$hash, 1, 12),
              length(fit@splits@indices),
              info$repeats %||% 1))

  # --- Permutation gap ---
  if (!is.null(object@permutation_gap) && nrow(object@permutation_gap) > 0) {
    pg <- object@permutation_gap
    cat("Permutation Gap Test:\n")
    cat(sprintf("  Observed metric: %.3f\n", pg$metric))
    cat(sprintf("  Permuted mean ± SD: %.3f ± %.3f\n", pg$perm_mean, pg$perm_sd))
    cat(sprintf("  Gap: %.3f (larger gap → less leakage)\n\n", pg$gap))
  } else {
    cat("Permutation Gap Test: not available.\n\n")
  }

  # --- Batch association ---
  if (!is.null(object@batch_assoc) && nrow(object@batch_assoc) > 0) {
    ba <- object@batch_assoc
    cat("Batch / Study Association:\n")
    cat(sprintf("  χ² = %.3f (df = %.0f), p = %.3g\n\n",
                ba$stat, ba$df, ba$pval))
  } else {
    cat("Batch / Study Association: none detected.\n\n")
  }

  # --- Duplicate detection ---
  if (!is.null(object@duplicates) && nrow(object@duplicates) > 0) {
    dd <- object@duplicates
    cat("Near-Duplicate Samples:\n")
    cat(sprintf("  %d pairs detected above cosine ≥ %.3f\n",
                nrow(dd), object@info$duplicate_threshold))
    head_pairs <- utils::head(dd, 5)
    cat("  Example pairs:\n")
    print(head_pairs, row.names = FALSE)
    cat("\n")
  } else {
    cat("No near-duplicates detected.\n\n")
  }

  # --- Overall interpretation ---
  cat("Interpretation:\n")
  cat(ifelse(
    !is.null(object@permutation_gap) && object@permutation_gap$gap < 0.01,
    "  ⚠️  Possible leakage (tiny permutation gap).\n",
    "  ✓  No strong evidence of leakage.\n"
  ))

  cat("\n----------------------------------------------\n")
  cat("Use `plot(object)` for visual summary (optional)\n")
  cat("----------------------------------------------\n\n")

  invisible(object)
}
