# Delta Leakage Sensitivity Index (ΔLSI) ────────────────────────────────────
# Compares a naive (potentially leaky) CV pipeline against a guarded pipeline
# and quantifies leakage-induced performance inflation.

# ── Internal helpers ──────────────────────────────────────────────────────────

# Resolve which learner to use from a LeakFit's @metrics.
.dlsi_resolve_learner <- function(fit, learner) {
  if (!is.null(learner)) return(as.character(learner[1L]))
  if ("learner" %in% names(fit@metrics)) {
    lrns <- unique(as.character(fit@metrics$learner))
    lrns <- lrns[nzchar(lrns)]
    if (!length(lrns)) return(NULL)
    if (length(lrns) > 1L)
      warning(sprintf(
        "[delta_lsi] Multiple learners in fit@metrics (%s); using '%s'.",
        paste(lrns, collapse = ", "), lrns[1L]
      ))
    return(lrns[1L])
  }
  NULL
}

# Extract per-fold metric values from fit@metrics.
# Returns data.frame(fold, metric).
.dlsi_fold_metrics <- function(fit, metric, learner) {
  mdf <- fit@metrics

  if (!is.null(learner) && "learner" %in% names(mdf))
    mdf <- mdf[as.character(mdf$learner) == learner, , drop = FALSE]

  if (metric %in% names(mdf) && "fold" %in% names(mdf)) {
    out <- mdf[, c("fold", metric), drop = FALSE]
    names(out)[2L] <- "metric"
    return(out)
  }

  # Fallback: recompute fold metrics from @predictions
  .dlsi_fold_metrics_from_preds(fit, metric, learner)
}

# Recompute fold metrics from @predictions when not in @metrics.
.dlsi_fold_metrics_from_preds <- function(fit, metric, learner) {
  task    <- if (length(fit@task)) fit@task else "binomial"
  pred_df <- .dlsi_combine_preds(fit, learner)
  if (is.null(pred_df) || !"fold" %in% names(pred_df))
    stop(sprintf(
      "[delta_lsi] Metric '%s' not in fit@metrics and no usable predictions.",
      metric
    ))
  folds <- sort(unique(pred_df$fold))
  vals  <- vapply(folds, function(f) {
    sub <- pred_df[pred_df$fold == f, , drop = FALSE]
    tryCatch(
      .metric_value(metric, task, sub$truth, sub$pred, pred_df = sub),
      error = function(e) NA_real_
    )
  }, numeric(1L))
  data.frame(fold = folds, metric = vals, stringsAsFactors = FALSE)
}

# Combine all per-fold prediction data frames into one.
.dlsi_combine_preds <- function(fit, learner = NULL) {
  preds <- fit@predictions
  if (!length(preds)) return(NULL)
  pred_df <- tryCatch(
    do.call(rbind, Filter(is.data.frame, preds)),
    error = function(e) NULL
  )
  if (is.null(pred_df) || !nrow(pred_df)) return(NULL)
  if (!is.null(learner) && "learner" %in% names(pred_df)) {
    pred_df <- pred_df[as.character(pred_df$learner) == learner, , drop = FALSE]
  } else if ("learner" %in% names(pred_df)) {
    first_lrn <- unique(as.character(pred_df$learner))[1L]
    pred_df   <- pred_df[as.character(pred_df$learner) == first_lrn, , drop = FALSE]
  }
  pred_df
}

# Get repeat_id for each sequential fold position (1..N).
.dlsi_repeat_ids <- function(fit) {
  vapply(fit@splits@indices, function(z) {
    as.integer(z$repeat_id %||% 1L)
  }, integer(1L))
}

# Get test-set size for each sequential fold position.
# Resolves compact splits through fold_assignments; falls back to the
# predictions-based count only when no test indices can be resolved.
.dlsi_fold_sizes <- function(fit, pred_df = NULL) {
  sizes <- vapply(seq_along(fit@splits@indices), function(i) {
    n <- length(.dlsi_fold_test_idx(fit@splits, i))
    if (n > 0L) as.integer(n) else NA_integer_
  }, integer(1L))

  if (any(is.na(sizes)) && !is.null(pred_df) && "fold" %in% names(pred_df)) {
    tab      <- table(pred_df$fold)
    fold_seq <- seq_along(fit@splits@indices)
    for (i in which(is.na(sizes))) {
      k <- as.character(fold_seq[i])
      if (k %in% names(tab)) sizes[i] <- as.integer(tab[[k]])
    }
  }
  sizes[is.na(sizes)] <- 1L
  sizes
}

# Resolve test-set indices for sequential fold position `i` (supports explicit
# and compact LeakSplits).
.dlsi_fold_test_idx <- function(splits, i) {
  idx_list <- splits@indices
  if (i < 1L || i > length(idx_list)) return(integer(0))
  z <- idx_list[[i]]
  if (!is.null(z$test)) return(as.integer(z$test))
  fa <- splits@info$fold_assignments
  if (isTRUE(splits@info$compact) && !is.null(fa) && length(fa)) {
    r <- as.integer(z$repeat_id %||% 1L)
    if (r >= 1L && r <= length(fa) && !is.null(fa[[r]])) {
      return(which(fa[[r]] == z$fold))
    }
  }
  integer(0)
}

# Check whether two fits share exactly the same fold-level test indices.
# Returns TRUE only when every fold's test set is identical (same observations,
# same order after sorting).  This guards against pairing fits that happen to
# have the same repeat count but were run on different data splits.
.dlsi_splits_match <- function(fit_naive, fit_guarded) {
  idx_n <- fit_naive@splits@indices
  idx_g <- fit_guarded@splits@indices
  if (length(idx_n) != length(idx_g)) return(FALSE)
  all(vapply(seq_along(idx_n), function(i) {
    te_n <- sort(.dlsi_fold_test_idx(fit_naive@splits, i))
    te_g <- sort(.dlsi_fold_test_idx(fit_guarded@splits, i))
    length(te_n) > 0L && identical(te_n, te_g)
  }, logical(1L)))
}

# Whether a higher value of `metric` is better, from its name. Lower-is-better
# patterns are checked first; an unrecognised name warns and assumes higher.
.dlsi_metric_direction <- function(metric) {
  m <- tolower(metric)
  # Skill scores and explained-deviance are higher-is-better even though their
  # names contain a loss term, so they are checked first.
  if (grepl("skill|explained", m))
    return(TRUE)
  if (grepl(paste0("rmse|mse|mae|mape|mase|rmsle|log_?loss|brier|error|loss|deviance|cost|",
                   "^mad$|^ece$|^mce$|^fpr$|^fnr$|^nll$|cross_?entropy|misclassification"), m))
    return(FALSE)
  if (grepl(paste0("auc|aun[pu]|accuracy|f1|f_meas|precision|recall|sens|spec|kappa|^kap$|",
                   "mcc|rsq|r2|cindex|c_index|concordance|ppv|npv|j_index|^ccc$|gain_capture"), m))
    return(TRUE)
  warning(sprintf(
    "[delta_lsi] Cannot tell whether a higher '%s' is better; assuming it is. Set `higher_is_better` explicitly.",
    metric), call. = FALSE)
  TRUE
}

# Huber M-estimator of location (iterative; Huber k = 1.345).
.dlsi_huber <- function(x, k = 1.345, tol = 1e-8, max_iter = 100L) {
  x <- as.numeric(x[is.finite(x)])
  n <- length(x)
  if (!n)    return(NA_real_)
  if (n == 1L) return(x[1L])
  mu    <- stats::median(x)
  sigma <- stats::mad(x, constant = 1.4826)
  if (!is.finite(sigma) || sigma < tol) return(mu)
  for (iter in seq_len(max_iter)) {
    r     <- (x - mu) / sigma
    w     <- ifelse(abs(r) <= k, 1.0, k / abs(r))
    mu_new <- sum(w * x) / sum(w)
    if (abs(mu_new - mu) < tol * (1.0 + abs(mu))) break
    mu <- mu_new
  }
  mu
}

# Training-set size for each sequential fold position. Explicit splits use
# length(train); compact splits rebuild the training set the way
# fit_resample() does (.bio_compact_train_indices()).
.dlsi_fold_train_sizes <- function(splits) {
  idx <- splits@indices
  if (!length(idx)) return(numeric(0))
  cd    <- splits@info$coldata
  n_all <- if (is.data.frame(cd)) nrow(cd) else NA_real_
  fa    <- splits@info$fold_assignments
  vapply(seq_along(idx), function(i) {
    z <- idx[[i]]
    if (!is.null(z$train)) return(as.numeric(length(z$train)))
    te    <- .dlsi_fold_test_idx(splits, i)
    r     <- as.integer(z$repeat_id %||% 1L)
    n_tot <- if (!is.null(fa) && length(fa) >= r && !is.null(fa[[r]])) length(fa[[r]]) else n_all
    if (!is.finite(n_tot) || !length(te)) return(NA_real_)
    n_tr  <- length(.bio_compact_train_indices(splits, te, as.integer(n_tot)))
    if (n_tr > 0L) as.numeric(n_tr) else NA_real_
  }, numeric(1L))
}

# Column whose levels act as clusters for exchangeability = "by_group" /
# "within_batch". The split mode's own blocking axis comes first, so a
# study_loocv plan clusters on its studies even when `group` is also set.
.dlsi_cluster_column <- function(info, mode, exchangeability) {
  cst <- info$constraints %||% list()
  cst_col <- function(types) {
    for (ax in cst) if (!is.null(ax$type) && ax$type %in% types) return(ax$col)
    NULL
  }
  if (identical(exchangeability, "by_group")) {
    mode_col <- switch(as.character(mode %||% ""),
      study_loocv     = info$study,
      batch_blocked   = info$batch,
      subject_grouped = info$group,
      time_series     = info$group,
      combined        = if (length(cst)) cst[[1L]]$col else NULL,
      NULL
    )
    mode_col %||% info$group %||% info$study %||% info$batch
  } else {
    # within_batch: batch-like axes only, never a subject/group column.
    info$batch %||% info$study %||% cst_col(c("batch", "study"))
  }
}

# ── Main function ─────────────────────────────────────────────────────────────

#' Delta Leakage Sensitivity Index (Delta LSI)
#'
#' Compares a naive (potentially leaky) cross-validation pipeline against a
#' guarded (leakage-protected) pipeline and quantifies leakage-induced
#' performance inflation using the Leakage Sensitivity Index (LSI).
#'
#' @details
#' \subsection{Method}{
#' For each fit, per-fold metric values are extracted from \code{fit@@metrics}
#' (or recomputed from \code{fit@@predictions} if necessary). Fold test-set
#' sizes are used as weights to aggregate fold metrics into per-repeat
#' estimates \eqn{\mu_r}. The repeat-level delta
#' \eqn{\Delta_r = s \cdot (\mu_r^{\text{naive}} - \mu_r^{\text{guarded}})}
#' captures leakage-induced performance inflation for each CV repeat, where
#' \eqn{s = +1} for higher-is-better metrics (e.g., AUC) and \eqn{s = -1}
#' for lower-is-better metrics (e.g., RMSE), so that \eqn{\Delta_r > 0}
#' always indicates the naive pipeline is more optimistic than the guarded one.
#'
#' The \strong{delta_lsi} point estimate is the Huber M-estimator (k = 1.345)
#' applied to \eqn{\{\Delta_r\}}, which is robust to occasional outlier
#' repeats. \strong{delta_metric} is the arithmetic mean of \eqn{\{\Delta_r\}}
#' for easy interpretation in the original metric's units.
#'
#' Pairing requires that \code{fit_leaky} and \code{fit_guarded} share
#' \emph{identical fold structures} (same test-set membership per fold).
#' When fold structures differ, a warning is issued and the fits are treated
#' as unpaired. Folds whose metric is missing in either fit are dropped from
#' both before pairing, so the paired delta always compares the same test sets.
#'
#' \strong{Inference} uses a corrected repeated cross-validation t-test on
#' the paired fold deltas \eqn{d_j}, with folds that re-use an identical
#' test set (e.g. repeated leave-one-batch-out) averaged into one, so \eqn{J}
#' counts distinct test sets. Repeated CV re-partitions the
#' same data, so fold deltas are correlated through overlapping training sets
#' and across repeats; treating them (or the repeats) as independent gives
#' badly anti-conservative tests. The variance of the mean is corrected as
#' \eqn{(1/J + n_{\text{test}}/n_{\text{train}})\, s_d^2} (Nadeau & Bengio,
#' 2003), and the t reference uses \eqn{k - 1} degrees of freedom, where
#' \eqn{k} is the number of folds per repeat (calibrated degrees of freedom,
#' Bouckaert, 2003): repeats sharpen the estimate but add no independent
#' information. The test statistic is \code{delta_metric / se}; the 95\% CIs
#' for \code{delta_metric} and \code{delta_lsi} are the estimate \eqn{\pm
#' t_{0.975, k-1}\, se}. In simulations with real model fits under
#' \eqn{H_0} (subject-grouped, batch-blocked and leave-one-study-out designs,
#' 1-20 repeats) the type-I error stayed at or below the nominal 5\%; with
#' 3-fold designs (2 df) the test is conservative.
#'
#' \code{R_eff} is the number of folds per repeat \eqn{k} (for
#' \code{study_loocv}, the number of studies), the number of effectively
#' independent units.
#' }
#'
#' \subsection{Inference tiers}{
#' Tier labels are kept from earlier versions for compatibility; every tier
#' except D reports a p-value and both CIs.
#' \describe{
#'   \item{\code{"A_full_inference"}}{R_eff >= 20}
#'   \item{\code{"B_signflip_ci"}}{10 <= R_eff < 20}
#'   \item{\code{"C_signflip"}}{3 <= R_eff < 10}
#'   \item{\code{"D_insufficient"}}{R_eff < 3 or unpaired: point estimate only}
#' }
#' \code{inference_ok} is \code{TRUE} whenever a finite p-value and CI were
#' computed.}
#'
#' @references
#' Nadeau, C. and Bengio, Y. (2003). Inference for the generalization error.
#' \emph{Machine Learning}, 52, 239-281.
#'
#' Bouckaert, R. R. (2003). Choosing between two learning algorithms based on
#' calibrated tests. \emph{Proceedings of ICML 2003}, 51-58.
#'
#' @param fit_leaky A \code{\linkS4class{LeakFit}} object from the leaky
#'   (unprotected) evaluation pipeline.
#' @param fit_guarded A \code{\linkS4class{LeakFit}} object from the guarded
#'   (leakage-protected) evaluation pipeline.
#' @param metric Character. Performance metric to compare. Must appear in
#'   \code{fit@@metrics} of both fits (e.g., \code{"auc"}, \code{"rmse"}).
#' @param exchangeability Character. One of \code{"iid"} (default),
#'   \code{"by_group"}, \code{"within_batch"}, \code{"blocked_time"}. Stored
#'   and reported; all settings use the same corrected t-test, because the
#'   correction already accounts for the dependence these settings used to
#'   address. For \code{"by_group"} / \code{"within_batch"} the number of
#'   group / batch levels is reported in \code{info$n_groups_used} /
#'   \code{info$n_batches_used}. Earlier versions ran sign-flip tests here
#'   that were anti-conservative with real model fits.
#' @param block_size Ignored; retained for compatibility (a warning is issued
#'   when supplied).
#' @param learner Optional character. Learner name to select from multi-learner
#'   fits. If \code{NULL}, the first learner found in \code{fit@@metrics} is used.
#' @param higher_is_better Logical or \code{NULL}. Whether a higher value of
#'   \code{metric} indicates better performance. When \code{NULL} (default),
#'   auto-detected from the metric name: names containing \code{"rmse"},
#'   \code{"mse"}, \code{"mae"}, \code{"mape"}, \code{"mase"},
#'   \code{"log_loss"}, \code{"brier"}, \code{"error"}, \code{"loss"},
#'   \code{"deviance"} or \code{"cost"} are lower-is-better; names such as
#'   \code{"auc"}, \code{"accuracy"}, \code{"f1"}, \code{"kappa"},
#'   \code{"rsq"} or \code{"cindex"} are higher-is-better; any other name
#'   warns and is treated as higher-is-better. Setting this correctly ensures
#'   that a positive \code{delta_lsi} always indicates leakage inflation
#'   (the naive pipeline is artificially more optimistic than the guarded one).
#' @param M_boot Ignored; retained for compatibility (no bootstrap is used).
#' @param M_flip Ignored; retained for compatibility (no sign-flip test is used).
#' @param strict Logical. If \code{TRUE}, error on insufficient R_eff instead
#'   of a warning.
#' @param return_details Logical. If \code{TRUE}, include the per-repeat
#'   \eqn{\Delta_r} vector, the paired fold deltas and the original fit objects
#'   in the \code{info} slot.
#' @param seed Ignored; retained for compatibility (the method is deterministic).
#' @param ... Unused. Reserved for deprecated aliases such as
#'   \code{fit_naive}.
#'
#' @return A \code{\linkS4class{LeakDeltaLSI}} object. Use
#'   \code{summary()} to print a formatted report. For programmatic
#'   access to inflation estimates, confidence intervals, p-values,
#'   tier labels, and the per-repeat fold-level deltas, use the S4
#'   accessor methods registered for \code{LeakDeltaLSI}:
#'   \itemize{
#'     \item \code{\link{dlsi_metric}(dlsi)} -- raw mean of
#'       per-repeat metric differences (the \code{delta_metric} slot).
#'     \item \code{\link{dlsi_robust}(dlsi)} -- Huber-robust point
#'       estimate (the \code{delta_lsi} slot).
#'     \item \code{\link{dlsi_ci}(dlsi, which = "robust" | "metric")} --
#'       corrected-t 95\% confidence interval for either the robust or
#'       the raw estimate.
#'     \item \code{\link{dlsi_p_value}(dlsi)} -- corrected repeated-CV
#'       t-test p-value.
#'     \item \code{\link{dlsi_tier}(dlsi)} -- inference-tier label
#'       (\code{"A_full_inference"}, \code{"B_signflip_ci"},
#'       \code{"C_signflip"}, or \code{"D_insufficient"}).
#'     \item \code{\link{dlsi_R_eff}(dlsi)} -- effective number of
#'       independent units (folds per repeat).
#'     \item \code{\link{dlsi_repeats}(dlsi, which = "naive" | "guarded")}
#'       -- per-repeat metric data frame for the requested pipeline.
#'   }
#'   The full list of accessor methods for the class is available
#'   through \code{methods(class = "LeakDeltaLSI")}.
#'
#' @seealso \code{\link{audit_leakage}}, \code{\link{fit_resample}},
#'   \code{\linkS4class{LeakDeltaLSI}}
#' @export
delta_lsi <- function(
    fit_leaky,
    fit_guarded,
    metric          = "auc",
    exchangeability = c("iid", "by_group", "within_batch", "blocked_time"),
    learner         = NULL,
    higher_is_better = NULL,
    block_size      = NULL,
    M_boot          = 2000L,
    M_flip          = 10000L,
    strict          = FALSE,
    return_details  = FALSE,
    seed            = 42L,
    ...
) {
  extras <- list(...)
  if ("fit_naive" %in% names(extras)) {
    if (!missing(fit_leaky)) {
      stop("Use either fit_leaky or the deprecated fit_naive argument, not both.")
    }
    fit_leaky <- extras$fit_naive
    extras$fit_naive <- NULL
    warning("delta_lsi(): fit_naive is deprecated; use fit_leaky instead.")
  }
  if (length(extras)) {
    stop(
      sprintf(
        "Unused argument(s): %s",
        paste(names(extras), collapse = ", ")
      )
    )
  }

  exchangeability <- match.arg(exchangeability)

  if (!inherits(fit_leaky,   "LeakFit")) stop("fit_leaky must be a LeakFit object.")
  if (!inherits(fit_guarded, "LeakFit")) stop("fit_guarded must be a LeakFit object.")

  # ── Resolve metric directionality ────────────────────────────────────────────
  if (is.null(higher_is_better)) {
    higher_is_better <- .dlsi_metric_direction(metric)
  }
  # sign_factor > 0 means "naive - guarded > 0 is leakage"; flip for lower-is-better
  sign_factor <- if (isTRUE(higher_is_better)) 1.0 else -1.0

  # Resolve learners (one per fit; may differ)
  lrn_n <- .dlsi_resolve_learner(fit_leaky,   learner)
  lrn_g <- .dlsi_resolve_learner(fit_guarded, learner)

  # Per-fold metric data frames: columns fold, metric
  fm_n <- .dlsi_fold_metrics(fit_leaky,   metric, lrn_n)
  fm_g <- .dlsi_fold_metrics(fit_guarded, metric, lrn_g)

  # Sequential fold positions (1..N) for each fit
  fold_seq_n <- seq_along(fit_leaky@splits@indices)
  fold_seq_g <- seq_along(fit_guarded@splits@indices)

  # Repeat IDs and fold sizes
  rep_ids_n <- .dlsi_repeat_ids(fit_leaky)
  rep_ids_g <- .dlsi_repeat_ids(fit_guarded)

  pred_n <- .dlsi_combine_preds(fit_leaky,   lrn_n)
  pred_g <- .dlsi_combine_preds(fit_guarded, lrn_g)

  sizes_n <- .dlsi_fold_sizes(fit_leaky,   pred_n)
  sizes_g <- .dlsi_fold_sizes(fit_guarded, pred_g)

  # Build fold-level data frames
  .build_folds <- function(fm, rep_ids, sizes, fold_seq) {
    idx <- match(fm$fold, fold_seq)
    data.frame(
      fold      = fm$fold,
      metric    = fm$metric,
      repeat_id = rep_ids[idx],
      n         = sizes[idx],
      stringsAsFactors = FALSE
    )
  }

  folds_n <- .build_folds(fm_n, rep_ids_n, sizes_n, fold_seq_n)
  folds_g <- .build_folds(fm_g, rep_ids_g, sizes_g, fold_seq_g)
  folds_n$n[is.na(folds_n$n)] <- 1L
  folds_g$n[is.na(folds_g$n)] <- 1L

  # Aggregate folds → per-repeat size-weighted mean metric
  .agg_repeats <- function(folds_df) {
    reps <- sort(unique(folds_df$repeat_id[is.finite(folds_df$metric)]))
    rows <- lapply(reps, function(r) {
      sub <- folds_df[folds_df$repeat_id == r & is.finite(folds_df$metric), ]
      if (!nrow(sub)) return(NULL)
      w  <- pmax(as.numeric(sub$n), 1.0)
      mu <- sum(w * sub$metric) / sum(w)
      data.frame(
        repeat_id = r, metric = mu,
        n_folds = nrow(sub), total_n = sum(sub$n),
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, Filter(Negate(is.null), rows))
    if (is.null(out)) {
      out <- data.frame(repeat_id = integer(0), metric = numeric(0),
                        n_folds = integer(0), total_n = integer(0),
                        stringsAsFactors = FALSE)
    }
    out
  }

  reps_n <- .agg_repeats(folds_n)
  reps_g <- .agg_repeats(folds_g)

  R_naive   <- nrow(reps_n)
  R_guarded <- nrow(reps_g)

  # ── Pairing: matching fold structures ───────────────────────────────────────
  # Fits are paired when every fold position has identical test members in
  # both. Repeat counts may still differ when one fit has an all-NA repeat;
  # the fold-level pairing below drops such folds from both fits.
  paired <- FALSE
  if (R_naive > 0L && R_guarded > 0L) {
    paired <- .dlsi_splits_match(fit_leaky, fit_guarded)
    if (!paired && R_naive == R_guarded) {
      warning(sprintf(
        paste0("[delta_lsi] fit_leaky and fit_guarded have R=%d repeats each ",
               "but different fold structures; treating as unpaired."),
        R_naive
      ))
    }
  }

  # ── Fold-level paired deltas (sign_factor: positive = leakage inflation) ────
  fd <- NULL
  if (paired) {
    # Every fold position of the (shared) plan, so a fold missing from either
    # fit's metrics is counted as dropped rather than silently ignored.
    pos <- seq_along(fit_guarded@splits@indices)
    fd <- data.frame(
      fold      = pos,
      repeat_id = rep_ids_n[pos],
      n_test    = as.numeric(pmax(sizes_n[pos], 1L)),
      d         = sign_factor * (folds_n$metric[match(pos, folds_n$fold)] -
                                 folds_g$metric[match(pos, folds_g$fold)]),
      stringsAsFactors = FALSE
    )
    tr_sizes   <- .dlsi_fold_train_sizes(fit_guarded@splits)
    fd$n_train <- as.numeric(tr_sizes[match(fd$fold, seq_along(tr_sizes))])
    fd$test_key <- vapply(fd$fold, function(i)
      paste(sort(.dlsi_fold_test_idx(fit_guarded@splits, i)), collapse = ","), character(1L))
    n_dropped  <- sum(!is.finite(fd$d))
    fd <- fd[is.finite(fd$d), , drop = FALSE]
    if (n_dropped > 0L) {
      warning(sprintf(
        paste0("[delta_lsi] %d fold(s) with a non-finite metric in either fit ",
               "were dropped from both fits before pairing."),
        n_dropped
      ), call. = FALSE)
    }
    if (!nrow(fd)) {
      paired <- FALSE
    } else {
      # Per-pipeline summaries use the same (paired) folds as the delta.
      folds_n$metric[!folds_n$fold %in% fd$fold] <- NA_real_
      folds_g$metric[!folds_g$fold %in% fd$fold] <- NA_real_
      reps_n <- .agg_repeats(folds_n)
      reps_g <- .agg_repeats(folds_g)
      R_naive   <- nrow(reps_n)
      R_guarded <- nrow(reps_g)
    }
  }

  # Per-repeat deltas: size-weighted over the folds valid in both fits.
  delta_r <- NULL
  if (paired) {
    delta_r <- unname(vapply(split(fd, fd$repeat_id), function(s) {
      w <- pmax(s$n_test, 1)
      sum(w * s$d) / sum(w)
    }, numeric(1L)))
  }

  # ── Effective units ──────────────────────────────────────────────────────────
  # Repeats re-partition the same data, so they sharpen the estimate but add no
  # independent information: the effective number of units is the number of
  # folds per repeat (Bouckaert, 2003), not folds x repeats.
  # Folds that re-use an identical test set (e.g. repeated leave-one-batch-out)
  # carry no new information about the variance; count each test set once.
  fu <- NULL
  if (paired) {
    first <- !duplicated(fd$test_key)
    fu <- data.frame(
      d       = as.numeric(tapply(fd$d, fd$test_key, mean)[fd$test_key[first]]),
      n_test  = fd$n_test[first],
      n_train = fd$n_train[first]
    )
  }
  J       <- if (paired) nrow(fu) else 0L
  k_folds <- if (paired) as.integer(floor(stats::median(table(fd$repeat_id)))) else 0L
  R_eff   <- as.integer(k_folds)

  tier <- if (!paired || R_eff < 3L) "D_insufficient" else
          if (R_eff >= 20L) "A_full_inference" else
          if (R_eff >= 10L) "B_signflip_ci"    else
                            "C_signflip"

  if (tier == "D_insufficient") {
    msg <- sprintf(
      paste0("[delta_lsi] R_eff = %d folds per repeat (R_leaky=%d, R_guarded=%d, paired=%s). ",
             "The corrected t-test needs paired fits with >= 3 folds per repeat; ",
             "reporting point estimate only."),
      R_eff, R_naive, R_guarded, paired
    )
    if (isTRUE(strict)) stop(msg) else warning(msg)
  }

  # Point estimates
  delta_lsi_est    <- if (!is.null(delta_r)) {
    .dlsi_huber(delta_r)
  } else {
    sign_factor * (.dlsi_huber(reps_n$metric) - .dlsi_huber(reps_g$metric))
  }
  delta_metric_est <- if (!is.null(delta_r)) {
    mean(delta_r, na.rm = TRUE)
  } else {
    sign_factor * (mean(reps_n$metric, na.rm = TRUE) - mean(reps_g$metric, na.rm = TRUE))
  }

  # ── Corrected repeated-CV t-test and CI ─────────────────────────────────────
  # Nadeau & Bengio (2003): Var(mean) = (1/J + n_test/n_train) * s^2 over the J
  # fold deltas, correcting for overlapping training sets; df = folds per
  # repeat - 1 (calibrated degrees of freedom). Calibrated against real model
  # fits under H0 across grouped, batch-blocked and leave-one-study-out designs.
  p_val  <- NA_real_
  se     <- NA_real_
  df_t   <- NA_integer_
  rho    <- NA_real_
  ci_lsi <- c(NA_real_, NA_real_)
  ci_dm  <- c(NA_real_, NA_real_)
  if (tier != "D_insufficient") {
    rho <- if (all(is.finite(fu$n_train)) && sum(fu$n_train) > 0) {
      mean(fu$n_test) / mean(fu$n_train)
    } else {
      1 / (k_folds - 1)   # k-fold approximation when training sizes are unknown
    }
    df_t <- as.integer(R_eff - 1L)
    se   <- sqrt(.nb_corrected_var(fu$d, n_train = 1, n_test = rho))  # (1/J + rho) * var
    if (is.finite(se)) {
      p_val <- if (se > 0) {
        2 * stats::pt(-abs(delta_metric_est / se), df_t)
      } else if (abs(delta_metric_est) > 0) 0 else 1
      q      <- stats::qt(0.975, df_t)
      ci_dm  <- delta_metric_est + c(-1, 1) * q * se
      ci_lsi <- delta_lsi_est    + c(-1, 1) * q * se
    }
  }

  inference_ok <- tier != "D_insufficient" && is.finite(p_val) && all(is.finite(ci_dm))

  # Informational: number of blocking-unit levels (does not change the test).
  n_groups_used  <- NA_integer_
  n_batches_used <- NA_integer_
  if (exchangeability %in% c("by_group", "within_batch")) {
    info_g <- fit_guarded@splits@info
    cd     <- info_g$coldata
    col    <- .dlsi_cluster_column(info_g, fit_guarded@splits@mode, exchangeability)
    if (is.data.frame(cd) && !is.null(col) && col %in% names(cd)) {
      n_lv <- as.integer(length(unique(stats::na.omit(cd[[col]]))))
      if (identical(exchangeability, "by_group")) n_groups_used <- n_lv else n_batches_used <- n_lv
    }
  }
  if (!is.null(block_size)) {
    warning("[delta_lsi] `block_size` is ignored: all exchangeability settings use the corrected t-test.",
            call. = FALSE)
  }

  # Metadata
  info_list <- list(
    paired           = paired,
    R_naive          = R_naive,
    R_guarded        = R_guarded,
    higher_is_better = higher_is_better,
    metric_naive     = .dlsi_huber(reps_n$metric),
    metric_guarded   = .dlsi_huber(reps_g$metric),
    test             = "corrected repeated-CV t (Nadeau & Bengio 2003; df = folds per repeat - 1)",
    se               = se,
    df               = df_t,
    rho              = rho,
    n_fold_deltas    = if (paired) nrow(fd) else 0L,
    n_test_sets      = as.integer(J),
    folds_per_repeat = as.integer(k_folds),
    n_groups_used    = n_groups_used,
    n_batches_used   = n_batches_used,
    block_size_used  = NA_integer_,
    n_blocks         = NA_integer_,
    min_p_achievable = NA_real_
  )
  if (return_details) {
    info_list$delta_r     <- delta_r
    info_list$fold_deltas <- fd
    info_list$fit_leaky   <- fit_leaky
    info_list$fit_naive   <- fit_leaky
    info_list$fit_guarded <- fit_guarded
  }

  methods::new("LeakDeltaLSI",
    metric          = metric,
    exchangeability = exchangeability,
    tier            = tier,
    strict          = strict,
    R_eff           = R_eff,
    delta_lsi       = delta_lsi_est,
    delta_lsi_ci    = as.numeric(ci_lsi),
    delta_metric    = delta_metric_est,
    delta_metric_ci = as.numeric(ci_dm),
    p_value         = p_val,
    inference_ok    = inference_ok,
    folds_naive     = folds_n,
    folds_guarded   = folds_g,
    repeats_naive   = reps_n,
    repeats_guarded = reps_g,
    info            = info_list
  )
}

# ── S4 show method ─────────────────────────────────────────────────────────────

#' @rdname LeakClasses
#' @param object A \code{\linkS4class{LeakDeltaLSI}} object.
setMethod("show", "LeakDeltaLSI", function(object) {
  hib <- object@info[["higher_is_better"]]
  hib_str <- if (is.null(hib)) "" else
             if (isTRUE(hib)) " [higher=better]" else " [lower=better]"
  cat("A LeakDeltaLSI object\n")
  cat(sprintf("  Metric:         %s%s\n", object@metric, hib_str))
  cat(sprintf("  Tier:           %s  (R_eff = %d)\n", object@tier, object@R_eff))
  fmt <- function(x) if (is.finite(x)) sprintf("%.4f", x) else "NA"
  cat(sprintf("  delta_metric:   %s", fmt(object@delta_metric)))
  if (all(is.finite(object@delta_metric_ci)))
    cat(sprintf("  [%.4f, %.4f] 95%% CI", object@delta_metric_ci[1L], object@delta_metric_ci[2L]))
  cat("\n")
  cat(sprintf("  delta_lsi:      %s", fmt(object@delta_lsi)))
  if (all(is.finite(object@delta_lsi_ci)))
    cat(sprintf("  [%.4f, %.4f] 95%% CI", object@delta_lsi_ci[1L], object@delta_lsi_ci[2L]))
  cat("\n")
  if (is.finite(object@p_value))
    cat(sprintf("  p-value:        %.4f  (corrected repeated-CV t; tests delta_metric)\n",
                object@p_value))
  invisible(object)
})

# ── summary method ─────────────────────────────────────────────────────────────

#' Summarize a LeakDeltaLSI object
#'
#' Prints a human-readable summary of the Delta LSI analysis comparing leaky vs
#' guarded evaluation pipelines.
#'
#' @param object A \code{LeakDeltaLSI} object from \code{\link{delta_lsi}}.
#' @param digits Integer. Number of decimal places to show (default 3).
#' @param ... Unused.
#' @return Invisibly returns \code{object}.
#' @export
summary.LeakDeltaLSI <- function(object, digits = 3L, ...) {
  fmt <- function(x) if (is.finite(x)) formatC(x, digits = digits, format = "f") else "NA"
  hib <- object@info[["higher_is_better"]] %||% TRUE

  cat("\n=====================================\n")
  cat(" bioLeak Delta LSI Summary\n")
  cat("=====================================\n\n")
  cat(sprintf("Metric:            %s\n", object@metric))
  cat(sprintf("Exchangeability:   %s  (reported only; same corrected test)\n",
              object@exchangeability))
  cat(sprintf("Inference tier:    %s%s\n", object@tier,
              if (identical(object@tier, "D_insufficient")) "" else
                "  (label kept for compatibility; test is the corrected t)"))
  cat(sprintf("R_eff:             %d folds per repeat  (R_leaky=%s, R_guarded=%s, paired=%s)\n",
              object@R_eff,
              object@info[["R_naive"]]   %||% "?",
              object@info[["R_guarded"]] %||% "?",
              object@info[["paired"]]    %||% "?"))

  cat("\nLeaky pipeline:    ", fmt(object@info[["metric_naive"]]   %||% NA_real_), "\n", sep = "")
  cat("Guarded pipeline:  ", fmt(object@info[["metric_guarded"]] %||% NA_real_), "\n", sep = "")

  direction <- if (isTRUE(hib)) "leaky - guarded" else "-(leaky - guarded)"
  cat(sprintf("\nPoint estimates (%s; positive = leakage inflation):\n", direction))
  cat(sprintf("  delta_metric:  %s  (raw metric difference)\n", fmt(object@delta_metric)))
  ci_dm <- object@delta_metric_ci
  if (all(is.finite(ci_dm)))
    cat(sprintf("  95%% CI:        [%s, %s]\n", fmt(ci_dm[1L]), fmt(ci_dm[2L])))
  cat(sprintf("  delta_lsi:     %s  (Huber-robust)\n", fmt(object@delta_lsi)))
  ci_lsi <- object@delta_lsi_ci
  if (all(is.finite(ci_lsi)))
    cat(sprintf("  95%% CI:        [%s, %s]\n", fmt(ci_lsi[1L]), fmt(ci_lsi[2L])))

  cat("\nHypothesis test  [H0: no systematic inflation]\n")
  cat("  Corrected repeated-CV t (Nadeau & Bengio 2003; df = folds per repeat - 1)\n")
  if (is.finite(object@p_value)) {
    p   <- object@p_value
    sig <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
    cat(sprintf("  t-test p:      %.4f %s  (df = %s)\n", p, sig,
                object@info[["df"]] %||% "?"))
  } else {
    cat("  t-test p:      NA  (needs paired fits with >= 3 folds per repeat)\n")
  }

  cat(sprintf("\nInference valid:   %s\n",
              if (object@inference_ok) sprintf("YES (%s)", object@tier) else
                sprintf("NO (%s)", object@tier)))
  cat("\n")
  invisible(object)
}
