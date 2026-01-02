#' Summarize a leakage audit
#'
#' Prints a concise, human-readable report for a `LeakAudit` object produced by
#' [audit_leakage()]. The summary surfaces four diagnostics when available:
#' label-permutation gap (prediction-label association by default), batch/study
#' association tests (metadata aligned with fold splits), target leakage scan
#' (features strongly associated with the outcome), and near-duplicate detection
#' (high similarity in `X_ref`). The output reflects the stored audit results
#' only; it does not recompute any tests.
#'
#' @details
#' The permutation test quantifies prediction-label association when using fixed
#' predictions; refit-based permutations require `perm_refit = TRUE`. It does
#' not by itself prove or rule out leakage.
#' Batch association flags metadata that align with fold assignment; this may
#' reflect study design rather than leakage.
#' Target leakage scan uses univariate feature-outcome associations and can miss
#' multivariate leakage or proxies not present in `X_ref`.
#' Duplicate detection only considers the provided `X_ref` features and the
#' similarity threshold used during [audit_leakage()]. By default,
#' `duplicate_scope = "train_test"` filters to pairs that cross train/test;
#' set `duplicate_scope = "all"` to include within-fold duplicates.
#' Sections are reported as "not available" when the corresponding audit
#' component was not computed.
#'
#' @seealso [plot_perm_distribution()], [plot_fold_balance()], [plot_overlap_checks()]
#'
#' @param object A `LeakAudit` object from [audit_leakage()]. The summary reads
#'   stored results from `object` and prints them to the console.
#' @param digits Integer number of digits to show when formatting numeric
#'   statistics in the console output. Defaults to `3`. Increasing `digits`
#'   shows more precision; decreasing it shortens the printout without changing
#'   the underlying values.
#' @param ... Unused. Included for S3 method compatibility; additional
#'   arguments are ignored.
#' @return Invisibly returns `object` after printing the summary.
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   subject = rep(1:6, each = 2),
#'   outcome = rbinom(12, 1, 0.5),
#'   x1 = rnorm(12),
#'   x2 = rnorm(12)
#' )
#' splits <- make_splits(df, outcome = "outcome",
#'                       mode = "subject_grouped", group = "subject", v = 3)
#' custom <- list(
#'   glm = list(
#'     fit = function(x, y, task, weights, ...) {
#'       stats::glm(y ~ ., data = as.data.frame(x),
#'                  family = stats::binomial(), weights = weights)
#'     },
#'     predict = function(object, newdata, task, ...) {
#'       as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
#'                                 type = "response"))
#'     }
#'   )
#' )
#' fit <- fit_resample(df, outcome = "outcome", splits = splits,
#'                     learner = "glm", custom_learners = custom,
#'                     metrics = "auc", refit = FALSE, seed = 1)
#' audit <- audit_leakage(fit, metric = "auc", B = 5,
#'                        X_ref = df[, c("x1", "x2")], seed = 1)
#' summary(audit) # prints the audit report and returns `audit` invisibly
#' @export
summary.LeakAudit <- function(object, digits = 3, ...) {
  if (!inherits(object, "LeakAudit"))
    stop("Object must be a 'LeakAudit'.", call. = FALSE)

  cat("\n==============================\n")
  cat(" bioLeak Leakage Audit Summary\n")
  cat("==============================\n\n")

  fit <- object@fit
  info <- fit@info
  warn_sym <- .bio_symbol("warn")
  ok_sym   <- .bio_symbol("check")
  sym_pm <- .bio_symbol("pm")
  sym_chi <- .bio_symbol("chi_sq")
  sym_ge <- .bio_symbol("ge")

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

  # --- Label-permutation association test ---
  if (!is.null(object@permutation_gap) && nrow(object@permutation_gap) > 0) {
    pg <- object@permutation_gap
    perm_method <- object@info$perm_method %||% "fixed"
    perm_label <- if (identical(perm_method, "refit")) {
      "refit per permutation"
    } else {
      "fixed predictions"
    }
    cat("Label-Permutation Association Test:\n")
    cat(sprintf("  Method: %s\n", perm_label))
    cat(sprintf("  Observed metric: %s\n",
                formatC(pg$metric_obs, digits = digits, format = "f")))
    cat(sprintf("  Permuted mean %s SD: %s %s %s\n",
                sym_pm,
                formatC(pg$perm_mean, digits = digits, format = "f"),
                sym_pm,
                formatC(pg$perm_sd, digits = digits, format = "f")))
    cat(sprintf("  Gap: %s (larger gap = stronger non-random signal)\n",
                formatC(pg$gap, digits = digits, format = "f")))
    if (!identical(perm_method, "refit")) {
      cat("  Note: Fixed-prediction label permutations quantify prediction-label association.\n")
      cat("  They do NOT refit models and are not a full null test of no signal.\n")
    }
    cat("  This test does NOT diagnose information leakage. Use the Batch Association,\n")
    cat("  Target Leakage Scan, and Duplicate Detection sections to check for leakage.\n\n")
  } else {
    cat("Label-Permutation Association Test: not available.\n\n")
  }

  # --- Batch association ---
  if (!is.null(object@batch_assoc) && nrow(object@batch_assoc) > 0) {
    ba <- object@batch_assoc
    cat("Batch / Study Association:\n")
    if ("batch_col" %in% names(ba)) {
      for (i in seq_len(nrow(ba))) {
        cat(sprintf("  %s: %s = %s (df = %s), p = %s\n",
                    ba$batch_col[i],
                    sym_chi,
                    formatC(ba$stat[i], digits = digits, format = "f"),
                    formatC(ba$df[i], digits = digits, format = "f"),
                    formatC(ba$pval[i], digits = digits, format = "f")))
      }
      cat("\n")
    } else {
      cat(sprintf("  %s = %s (df = %s), p = %s\n\n",
                  sym_chi,
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
    cat(sprintf("  Features checked: %d | Flagged (score %s %s): %d\n",
                nrow(ta),
                sym_ge,
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
    dup_scope <- object@info$duplicate_scope %||% "all"
    scope_label <- if (identical(dup_scope, "train_test")) "train/test only" else "all pairs"
    cat(sprintf("  Scope: %s\n", scope_label))
    cat(sprintf("  %d pairs detected above %s %s %s\n",
                nrow(dd),
                sim_label,
                sym_ge,
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
#' Prints a compact console report for a [LeakFit] object created by
#' [fit_resample()]. The report lists task/outcome metadata, learners,
#' total folds, and cross-validated metrics summarized as mean and standard
#' deviation across completed folds, plus a small audit table with per-fold
#' train/test sizes and retained feature counts.
#'
#' This summary is meant for quick sanity checks of the resampling setup and
#' performance. It does not run leakage diagnostics and will not detect target
#' leakage, duplicate samples, or batch/study confounding; use [audit_leakage()]
#' or `summary()` on a [LeakAudit] object for those checks.
#'
#' @param object A [LeakFit] object returned by [fit_resample()]. It should
#'   contain `metric_summary` and `audit` slots; missing entries result in empty
#'   sections in the printed report.
#' @param digits Integer scalar. Number of decimal places to print in numeric
#'   summary tables. Defaults to 3; affects printed output only, not the
#'   returned data.
#' @param ... Unused. Included for S3 method compatibility; changing these
#'   values has no effect.
#' @return Invisibly returns `object@metric_summary`, a data frame of per-learner
#'   metric means and standard deviations computed across folds. This function
#'   does not recompute metrics.
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   subject = rep(1:6, each = 2),
#'   outcome = factor(rep(c(0, 1), each = 6)),
#'   x1 = rnorm(12),
#'   x2 = rnorm(12)
#' )
#' splits <- make_splits(
#'   df,
#'   outcome = "outcome",
#'   mode = "subject_grouped",
#'   group = "subject",
#'   v = 3,
#'   stratify = TRUE,
#'   progress = FALSE
#' )
#' custom <- list(
#'   glm = list(
#'     fit = function(x, y, task, weights, ...) {
#'       stats::glm(y ~ ., data = data.frame(y = y, x),
#'                  family = stats::binomial(), weights = weights)
#'     },
#'     predict = function(object, newdata, task, ...) {
#'       as.numeric(stats::predict(object,
#'                                 newdata = as.data.frame(newdata),
#'                                 type = "response"))
#'     }
#'   )
#' )
#' fit <- fit_resample(df, outcome = "outcome", splits = splits,
#'                     learner = "glm", custom_learners = custom,
#'                     metrics = "auc", seed = 1)
#' summary_df <- summary(fit)
#' summary_df
#' @export
summary.LeakFit <- function(object, digits = 3, ...) {
  if (!inherits(object, "LeakFit"))
    stop("Object must be of class 'LeakFit'.")

  cat("\n===========================\n")
  cat(" bioLeak Model Fit Summary\n")
  cat("===========================\n\n")

  sym_pm <- .bio_symbol("pm")

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
    cat(sprintf("Cross-validated metrics (mean %s SD):\n", sym_pm))
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
