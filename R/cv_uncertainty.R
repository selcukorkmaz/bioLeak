# Cross-validation uncertainty estimation -----------------------------------

#' Nadeau-Bengio corrected variance
#'
#' Computes the corrected variance for repeated K-fold CV as described in
#' Nadeau & Bengio (2003). The correction factor accounts for the overlap
#' between training sets across folds.
#'
#' @param per_fold_vals Numeric vector of per-fold metric values.
#' @param n_train Average number of training samples per fold. NULL to skip
#'   correction.
#' @param n_test Average number of test samples per fold. NULL to skip
#'   correction.
#' @return Scalar corrected variance, or \code{NA_real_} for single-fold input.
#' @keywords internal
.nb_corrected_var <- function(per_fold_vals, n_train = NULL, n_test = NULL) {
  per_fold_vals <- per_fold_vals[is.finite(per_fold_vals)]
  K <- length(per_fold_vals)
  if (K < 2L) return(NA_real_)
  v <- stats::var(per_fold_vals)
  if (!is.null(n_train) && !is.null(n_test) &&
      is.finite(n_train) && is.finite(n_test) && n_train > 0) {
    # Nadeau-Bengio correction: var_corrected = (1/K + n_test/n_train) * var
    corrected <- (1 / K + n_test / n_train) * v
  } else {
    # Fallback: simple var/K
    corrected <- v / K
  }
  corrected
}

# Per sequential fold of a LeakSplits: a key identifying its test set, its
# train/test sizes and repeat id. Repeats that re-use an identical test set
# (e.g. repeated leave-one-batch-out) add no information and are merged.
.bio_cv_fold_info <- function(splits) {
  idx <- splits@indices
  if (!length(idx)) return(NULL)
  keys <- vapply(seq_along(idx), function(i)
    paste(sort(.dlsi_fold_test_idx(splits, i)), collapse = ","), character(1))
  n_te <- vapply(seq_along(idx), function(i) length(.dlsi_fold_test_idx(splits, i)), numeric(1))
  data.frame(fold = seq_along(idx), test_set = keys, n_test = n_te,
             n_train = .dlsi_fold_train_sizes(splits),
             repeat_id = vapply(idx, function(z) as.integer(z$repeat_id %||% 1L), integer(1)),
             stringsAsFactors = FALSE)
}

# Fill n_train / n_test / folds_per_repeat from fold info (distinct test sets).
.bio_cv_defaults <- function(info, n_train, n_test, folds_per_repeat) {
  first <- !duplicated(info$test_set)
  if (is.null(n_train)) { v <- mean(info$n_train[first], na.rm = TRUE); if (is.finite(v) && v > 0) n_train <- v }
  if (is.null(n_test))  { v <- mean(info$n_test[first],  na.rm = TRUE); if (is.finite(v) && v > 0) n_test <- v }
  if (is.null(folds_per_repeat)) folds_per_repeat <- max(2L, as.integer(stats::median(table(info$repeat_id))))
  list(n_train = n_train, n_test = n_test, folds_per_repeat = folds_per_repeat)
}

#' Confidence intervals for cross-validated metrics
#'
#' Computes per-learner confidence intervals for each metric column in a
#' per-fold metrics data.frame. Supports the standard normal/t approach and the
#' Nadeau-Bengio (2003) corrected variance for repeated K-fold CV.
#'
#' Fold metrics from cross-validation are correlated: training sets overlap,
#' and repeats re-partition the same data. \code{"normal"} ignores this
#' (\code{sd / sqrt(K)} over all folds) and is too narrow for repeated CV.
#' \code{"nadeau_bengio"} inflates the variance to
#' \code{(1/K + n_test/n_train) * var} and uses a t reference on
#' \code{folds_per_repeat - 1} degrees of freedom (calibrated degrees of
#' freedom, Bouckaert 2003): repeats sharpen the estimate but add no
#' independent information. Without \code{folds_per_repeat} all folds are
#' treated as one repeat.
#'
#' @param metrics_df A \code{\linkS4class{LeakFit}} or a \code{LeakTune}
#'   (recommended: fold sizes, folds per repeat and repeated test sets are then
#'   taken from the object's splits), or a data.frame with columns \code{fold},
#'   \code{learner}, and one or more numeric metric columns, such as
#'   \code{fit@@metrics}. An optional character column \code{test_set}
#'   identifies each fold's test set: folds sharing a test set (repeats that
#'   re-use the same partition) are averaged and counted once.
#' @param level Confidence level (default 0.95).
#' @param method \code{"nadeau_bengio"} (default) or \code{"normal"}.
#'   \code{"normal"} treats folds as independent and is far too narrow for
#'   cross-validation (about 37\% coverage of a nominal 95\% interval for
#'   repeated 5-fold CV in simulation); it is kept for comparison.
#' @param n_train Average number of training samples per fold. Used only when
#'   \code{method = "nadeau_bengio"}. NULL to use fallback variance.
#' @param n_test Average number of test samples per fold. Used only when
#'   \code{method = "nadeau_bengio"}. Without \code{n_train} and
#'   \code{n_test} the correction cannot be applied: a warning is issued and
#'   the uncorrected variance \code{var / K} is used.
#' @param folds_per_repeat Number of folds in one CV repeat (e.g. 5 for
#'   5-fold CV repeated 10 times). Used only when
#'   \code{method = "nadeau_bengio"} to set the degrees of freedom. The
#'   resulting intervals are slightly conservative in simulation (about 98\%
#'   coverage at the 95\% level for repeated 5-fold CV). Intervals are not
#'   clipped to a metric's range (e.g. an AUC bound can exceed 1).
#' @return A data.frame with \code{learner} and, for each metric, columns
#'   \code{<metric>_mean}, \code{<metric>_sd}, \code{<metric>_ci_lo}, and
#'   \code{<metric>_ci_hi}.
#' @export
cv_ci <- function(metrics_df, level = 0.95,
                  method = c("nadeau_bengio", "normal"),
                  n_train = NULL, n_test = NULL, folds_per_repeat = NULL) {
  method <- match.arg(method)
  obj_splits <- NULL
  if (inherits(metrics_df, "LeakFit")) {
    obj_splits <- metrics_df@splits
    metrics_df <- metrics_df@metrics
  } else if (inherits(metrics_df, "LeakTune")) {
    obj_splits <- metrics_df$splits
    metrics_df <- metrics_df$metrics
  }
  if (inherits(obj_splits, "LeakSplits") && is.data.frame(metrics_df) &&
      "fold" %in% names(metrics_df)) {
    info <- .bio_cv_fold_info(obj_splits)
    if (!is.null(info)) {
      d <- .bio_cv_defaults(info, n_train, n_test, folds_per_repeat)
      n_train <- d$n_train; n_test <- d$n_test; folds_per_repeat <- d$folds_per_repeat
      if (!"test_set" %in% names(metrics_df)) {
        metrics_df$test_set <- info$test_set[match(metrics_df$fold, info$fold)]
      }
    }
  }
  scalar_pos <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0
  for (nm in c("n_train", "n_test")) {
    val <- get(nm)
    if (!is.null(val) && !scalar_pos(val)) {
      stop(sprintf("`%s` must be a single positive number (the average fold size).", nm), call. = FALSE)
    }
  }
  if (identical(method, "nadeau_bengio") && !is.null(folds_per_repeat) &&
      (!scalar_pos(folds_per_repeat) || folds_per_repeat < 2 ||
       folds_per_repeat != round(folds_per_repeat))) {
    stop("`folds_per_repeat` must be a single whole number >= 2.", call. = FALSE)
  }
  if (identical(method, "nadeau_bengio") &&
      (is.null(n_train) || is.null(n_test) || !is.finite(n_train) || !is.finite(n_test) ||
       n_train <= 0)) {
    .bio_warn(paste0(
      "cv_ci(method = \"nadeau_bengio\") needs n_train and n_test; without them the ",
      "correction cannot be applied and the uncorrected variance (var / K) is used, which is ",
      "too narrow for cross-validation."), "bioLeak_cv_ci_warning")
  }
  stopifnot(is.data.frame(metrics_df))
  stopifnot("fold" %in% names(metrics_df), "learner" %in% names(metrics_df))
  stopifnot(level > 0 && level < 1)

  metric_cols <- setdiff(names(metrics_df), c("fold", "learner", "test_set", "repeat_id"))
  if (!length(metric_cols)) {
    .bio_stop("No numeric metric columns found in metrics_df.",
              "bioLeak_input_error")
  }

  learners <- unique(metrics_df$learner)
  result_rows <- vector("list", length(learners))

  for (li in seq_along(learners)) {
    ln <- learners[[li]]
    sub <- metrics_df[metrics_df$learner == ln, , drop = FALSE]
    K <- nrow(sub)
    row <- list(learner = ln)

    for (mc in metric_cols) {
      vals <- sub[[mc]]
      if (!is.numeric(vals)) next
      if ("test_set" %in% names(sub) && !all(is.na(sub$test_set))) {
        # one value per distinct test set (repeats of the same split add nothing)
        key <- ifelse(is.na(sub$test_set), paste0("__fold_", sub$fold), sub$test_set)
        vals <- as.numeric(tapply(vals, factor(key, levels = unique(key)),
                                  function(v) if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)))
      }
      vals_finite <- vals[is.finite(vals)]
      K_eff <- length(vals_finite)
      m <- if (K_eff > 0L) mean(vals_finite) else NA_real_
      s <- if (K_eff > 1L) stats::sd(vals_finite) else NA_real_
      row[[paste0(mc, "_mean")]] <- m
      row[[paste0(mc, "_sd")]] <- s

      if (K_eff < 2L) {
        row[[paste0(mc, "_ci_lo")]] <- NA_real_
        row[[paste0(mc, "_ci_hi")]] <- NA_real_
        next
      }

      if (identical(method, "nadeau_bengio")) {
        corrected_var <- .nb_corrected_var(vals_finite, n_train = n_train, n_test = n_test)
        se <- sqrt(corrected_var)
      } else {
        se <- s / sqrt(K_eff)
      }

      alpha <- 1 - level
      if (identical(method, "nadeau_bengio")) {
        # Calibrated df: folds per repeat - 1 (repeats add no independent units).
        df_t <- min(K_eff, folds_per_repeat %||% K_eff) - 1
        q <- stats::qt(1 - alpha / 2, df = df_t)
      } else if (K_eff < 30) {
        q <- stats::qt(1 - alpha / 2, df = K_eff - 1)
      } else {
        q <- stats::qnorm(1 - alpha / 2)
      }

      row[[paste0(mc, "_ci_lo")]] <- m - q * se
      row[[paste0(mc, "_ci_hi")]] <- m + q * se
    }

    result_rows[[li]] <- as.data.frame(row, stringsAsFactors = FALSE)
  }

  do.call(rbind, result_rows)
}
