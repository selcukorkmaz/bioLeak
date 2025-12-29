#' Summarize a leakage audit
#'
#' Provides a console summary of the leakage audit results,
#' including permutation gap, batch association, target scan, and duplicates.
#'
#' @seealso [plot_perm_distribution()], [plot_fold_balance()], [plot_overlap_checks()]
#'
#' @param object LeakAudit
#' @param digits number of digits for numeric display
#' @param ... ignored
#' @return Invisibly returns the input object.
#' @examples
#' \dontrun{
#' # audit <- audit_leakage(fit, metric = "auc", B = 50)
#' summary(audit)
#' }
#' @export
summary.LeakAudit <- function(object, digits = 3, ...) {
  if (!inherits(object, "LeakAudit"))
    stop("Object must be a 'LeakAudit'.", call. = FALSE)

  cat("\n==============================\n")
  cat(" bioLeak Leakage Audit Summary\n")
  cat("==============================\n\n")

  fit <- object@fit
  info <- fit@info
  utf8_ok <- FALSE
  if (requireNamespace("cli", quietly = TRUE)) {
    utf8_ok <- isTRUE(cli::is_utf8_output())
  } else {
    loc <- tryCatch(l10n_info(), error = function(e) NULL)
    if (!is.null(loc) && isTRUE(loc$`UTF-8`)) utf8_ok <- TRUE
  }
  warn_sym <- if (utf8_ok) "⚠️" else "WARNING:"
  ok_sym   <- if (utf8_ok) "✓" else "OK:"

  task <- fit@task %||% NA_character_
  outcome <- fit@outcome %||% NA_character_
  splits <- fit@splits
  mode <- if (!is.null(splits)) splits@mode %||% NA_character_ else NA_character_
  indices <- if (!is.null(splits)) splits@indices else NULL
  hash <- info$hash %||% "NA"

  pos_class <- fit@info$positive_class %||% NA_character_
  pos_label <- ""
  if (!is.null(pos_class) && isTRUE(task == "binomial") &&
      !is.na(pos_class) && nzchar(as.character(pos_class))) {
    pos_label <- paste0(" | Positive class: ", as.character(pos_class))
  }
  cat(sprintf("Task: %s | Outcome: %s | Splitting mode: %s%s\n",
              task, outcome, mode, pos_label))
  cat(sprintf("Hash: %s | Folds: %d | Repeats: %d\n\n",
              substr(hash, 1, 12),
              length(indices),
              info$repeats %||% 1))

  # --- Permutation significance test ---
  if (!is.null(object@permutation_gap) && nrow(object@permutation_gap) > 0) {
    pg <- object@permutation_gap
    cat("Permutation Significance Test:\n")
    cat(sprintf("  Observed metric: %s\n",
                formatC(pg$metric_obs, digits = digits, format = "f")))
    cat(sprintf("  Permuted mean ± SD: %s ± %s\n",
                formatC(pg$perm_mean, digits = digits, format = "f"),
                formatC(pg$perm_sd, digits = digits, format = "f")))
    cat(sprintf("  Gap: %s (larger gap = stronger non-random signal)\n",
                formatC(pg$gap, digits = digits, format = "f")))
    cat("  Note: This tests if the model signal is non-random. It does NOT diagnose information leakage.\n")
    cat("  Use the Batch Association, Target Leakage Scan, and Duplicate Detection sections to check for leakage.\n\n")
  } else {
    cat("Permutation Significance Test: not available.\n\n")
  }

  # --- Batch association ---
  if (!is.null(object@batch_assoc) && nrow(object@batch_assoc) > 0) {
    ba <- object@batch_assoc
    cat("Batch / Study Association:\n")
    if ("batch_col" %in% names(ba)) {
      for (i in seq_len(nrow(ba))) {
        cat(sprintf("  %s: χ² = %s (df = %s), p = %s\n",
                    ba$batch_col[i],
                    formatC(ba$stat[i], digits = digits, format = "f"),
                    formatC(ba$df[i], digits = digits, format = "f"),
                    formatC(ba$pval[i], digits = digits, format = "f")))
      }
      cat("\n")
    } else {
      cat(sprintf("  χ² = %s (df = %s), p = %s\n\n",
                  formatC(ba$stat, digits = digits, format = "f"),
                  formatC(ba$df, digits = digits, format = "f"),
                  formatC(ba$pval, digits = digits, format = "f")))
    }
  } else {
    cat("Batch / Study Association: none detected.\n\n")
  }

  # --- Target leakage scan ---
  if (!is.null(object@target_assoc) && nrow(object@target_assoc) > 0) {
    ta <- object@target_assoc
    threshold <- object@info$target_threshold %||% 0.9
    flagged <- ta[!is.na(ta$score) & ta$score >= threshold, , drop = FALSE]
    cat("Target Leakage Scan:\n")
    cat(sprintf("  Features checked: %d | Flagged (score >= %s): %d\n",
                nrow(ta),
                formatC(threshold, digits = digits, format = "f"),
                nrow(flagged)))
    if (nrow(flagged) > 0) {
      flagged <- flagged[order(flagged$score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
      top <- utils::head(flagged, 5)
      cols <- intersect(c("feature", "metric", "value", "score", "p_value"), names(top))
      print(top[, cols, drop = FALSE], row.names = FALSE)
    } else {
      cat("  No strong proxy features detected.\n")
    }
    cat("\n")
  } else {
    cat("Target Leakage Scan: not available.\n\n")
  }

  # --- Duplicate detection ---
  if (!is.null(object@duplicates) && nrow(object@duplicates) > 0) {
    dd <- object@duplicates
    cat("Near-Duplicate Samples:\n")
    sim_label <- object@info$sim_method %||% "cosine"
    cat(sprintf("  %d pairs detected above %s ≥ %s\n",
                nrow(dd),
                sim_label,
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
    cat(sprintf("  %s Little non-random signal (gap near zero).\n", warn_sym))
  } else if (object@permutation_gap$gap < 0.05) {
    cat(sprintf("  %s Modest non-random signal.\n", ok_sym))
  } else {
    cat(sprintf("  %s Strong non-random signal.\n", ok_sym))
  }

  invisible(object)
}

#' Summarize a LeakFit object
#'
#' Provides a concise summary of resampled model performance,
#' including learners, folds, and mean ± SD of metrics.
#'
#' @param object A \code{LeakFit} object returned by [fit_resample()].
#' @param digits Number of digits to display.
#' @param ... Not used.
#' @return Invisibly returns the summary data frame.
#' @examples
#' \dontrun{
#' # fit <- fit_resample(...)
#' summary(fit)
#' }
#' @export
summary.LeakFit <- function(object, digits = 3, ...) {
  if (!inherits(object, "LeakFit"))
    stop("Object must be of class 'LeakFit'.")

  cat("\n===========================\n")
  cat(" bioLeak Model Fit Summary\n")
  cat("===========================\n\n")

  # Basic info
  info <- object@info
  cat(sprintf("Task: %s\n", object@task))
  cat(sprintf("Outcome: %s\n", object@outcome))
  if (identical(object@task, "binomial") && !is.null(info$positive_class) &&
      !is.na(info$positive_class) && nzchar(as.character(info$positive_class))) {
    cat(sprintf("Positive class: %s\n", as.character(info$positive_class)))
  }
  cat(sprintf("Learners: %s\n", paste(unique(object@metrics$learner), collapse = ", ")))
  cat(sprintf("Total folds: %d\n", length(object@splits@indices)))
  cat(sprintf("Refit performed: %s\n", if (isTRUE(info$refit)) "Yes" else "No"))
  cat(sprintf("Hash: %s\n\n", substr(info$hash, 1, 12)))

  # Metric summary
  if (nrow(object@metric_summary) > 0) {
    cat("Cross-validated metrics (mean ± SD):\n")
    ms <- object@metric_summary
    metrics_fmt <- as.data.frame(ms)
    num_cols <- vapply(metrics_fmt, is.numeric, logical(1))
    if (any(num_cols)) {
      metrics_fmt[num_cols] <- lapply(metrics_fmt[num_cols], round, digits = digits)
    }
    print(metrics_fmt)
    cat("\n")
  } else {
    cat("No metric summary available.\n\n")
  }

  # Audit information
  if (nrow(object@audit) > 0) {
    cat("Audit overview:\n")
    audit_df <- head(object@audit, 5)
    print(audit_df, row.names = FALSE)
    cat("\n")
  } else {
    cat("No audit information stored.\n\n")
  }

  invisible(object@metric_summary)
}
