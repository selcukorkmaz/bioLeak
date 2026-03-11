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

#' Confidence intervals for cross-validated metrics
#'
#' Computes per-learner confidence intervals for each metric column in a
#' per-fold metrics data.frame. Supports the standard normal/t approach and the
#' Nadeau-Bengio (2003) corrected variance for repeated K-fold CV.
#'
#' @param metrics_df Data.frame with columns \code{fold}, \code{learner}, and
#'   one or more numeric metric columns.
#' @param level Confidence level (default 0.95).
#' @param method One of \code{"normal"} or \code{"nadeau_bengio"}.
#' @param n_train Average number of training samples per fold. Used only when
#'   \code{method = "nadeau_bengio"}. NULL to use fallback variance.
#' @param n_test Average number of test samples per fold. Used only when
#'   \code{method = "nadeau_bengio"}. NULL to use fallback variance.
#' @return A data.frame with \code{learner} and, for each metric, columns
#'   \code{<metric>_mean}, \code{<metric>_sd}, \code{<metric>_ci_lo}, and
#'   \code{<metric>_ci_hi}.
#' @export
cv_ci <- function(metrics_df, level = 0.95,
                  method = c("normal", "nadeau_bengio"),
                  n_train = NULL, n_test = NULL) {
  method <- match.arg(method)
  stopifnot(is.data.frame(metrics_df))
  stopifnot("fold" %in% names(metrics_df), "learner" %in% names(metrics_df))
  stopifnot(level > 0 && level < 1)

  metric_cols <- setdiff(names(metrics_df), c("fold", "learner"))
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

      # Use t-distribution for K_eff < 30, normal otherwise
      alpha <- 1 - level
      if (K_eff < 30) {
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
