# Additional diagnostic helpers ----------------------------------------------

.select_predictions_for_diagnostics <- function(fit, context, learner = NULL) {
  pred_df <- if (length(fit@predictions)) {
    do.call(rbind, lapply(fit@predictions, function(df) data.frame(df, stringsAsFactors = FALSE)))
  } else {
    NULL
  }
  if (is.null(pred_df) || !nrow(pred_df)) {
    stop(sprintf("No predictions available for %s.", context), call. = FALSE)
  }

  has_learner <- "learner" %in% names(pred_df)
  if (has_learner) {
    pred_df$learner <- as.character(pred_df$learner)
    learner_vals <- unique(pred_df$learner)
    if (is.null(learner)) {
      if (length(learner_vals) == 1L) {
        learner <- learner_vals[[1]]
      } else {
        stop("Multiple learners found in predictions; specify `learner` to select a single model.",
             call. = FALSE)
      }
    } else {
      if (length(learner) != 1L) stop("learner must be a single value.", call. = FALSE)
      if (!learner %in% learner_vals) {
        stop(sprintf("Learner '%s' not found in predictions. Available: %s",
                     learner, paste(learner_vals, collapse = ", ")),
             call. = FALSE)
      }
    }
    pred_df <- pred_df[pred_df$learner == learner, , drop = FALSE]
    if (!nrow(pred_df)) {
      stop(sprintf("No predictions available for learner '%s'.", learner), call. = FALSE)
    }
  } else {
    if (!is.null(learner)) {
      warning("`learner` ignored: predictions do not include learner IDs.", call. = FALSE)
    } else if (!is.null(fit@metrics) && "learner" %in% names(fit@metrics) &&
               length(unique(fit@metrics$learner)) > 1L) {
      warning("Multiple learners were fit but predictions lack learner IDs; diagnostics may mix learners. Refit with updated bioLeak.",
              call. = FALSE)
    }
  }

  pred_df
}

#' Calibration diagnostics for binomial predictions
#'
#' Computes reliability curve summaries and calibration metrics for a
#' binomial [LeakFit] using out-of-fold predictions.
#'
#' @param fit A [LeakFit] object from [fit_resample()].
#' @param bins Integer number of probability bins for the calibration curve.
#' @param min_bin_n Minimum samples per bin used in plotting; bins smaller than
#'   this are retained in the output but can be filtered by the caller.
#' @param learner Optional character scalar. When predictions include multiple
#'   learners, selects the learner to summarize.
#' @return A list with a `curve` data.frame and a one-row `metrics` data.frame
#'   containing ECE, MCE, and Brier score.
#' @examples
#' \donttest{
#' fit <- fit_resample(...)
#' cal <- calibration_summary(fit, bins = 10)
#' cal$metrics
#' }
#' @export
calibration_summary <- function(fit, bins = 10, min_bin_n = 5, learner = NULL) {
  stopifnot(inherits(fit, "LeakFit"))
  if (!identical(fit@task, "binomial")) {
    stop("calibration_summary is only available for binomial tasks.", call. = FALSE)
  }
  bins <- as.integer(bins)
  if (!is.finite(bins) || bins < 2L) {
    stop("bins must be an integer >= 2.", call. = FALSE)
  }
  min_bin_n <- as.integer(min_bin_n)
  if (!is.finite(min_bin_n) || min_bin_n < 1L) {
    stop("min_bin_n must be an integer >= 1.", call. = FALSE)
  }

  pred_df <- .select_predictions_for_diagnostics(fit, "calibration", learner = learner)
  if (!"pred" %in% names(pred_df)) {
    stop("Predictions missing 'pred' column for calibration.", call. = FALSE)
  }
  if (!"truth" %in% names(pred_df)) {
    stop("Predictions missing 'truth' values for calibration.", call. = FALSE)
  }

  pred <- as.numeric(pred_df$pred)
  truth <- pred_df$truth
  if (!is.factor(truth)) truth <- factor(truth)
  pos_class <- fit@info$positive_class
  if (is.null(pos_class) || !as.character(pos_class) %in% levels(truth)) {
    pos_class <- levels(truth)[2]
  }
  y <- as.integer(truth == pos_class)
  ok <- is.finite(pred) & !is.na(y)
  pred <- pred[ok]
  y <- y[ok]
  if (!length(pred)) {
    stop("No finite predictions available for calibration.", call. = FALSE)
  }

  if (any(pred < 0 | pred > 1, na.rm = TRUE)) {
    warning("Predictions outside [0,1] were clipped for calibration.", call. = FALSE)
    pred <- pmin(pmax(pred, 0), 1)
  }

  breaks <- seq(0, 1, length.out = bins + 1L)
  bin_id <- cut(pred, breaks = breaks, include.lowest = TRUE, right = TRUE)
  n_bin <- tapply(pred, bin_id, length)
  pred_mean <- tapply(pred, bin_id, mean)
  obs_rate <- tapply(y, bin_id, mean)
  pred_min <- tapply(pred, bin_id, min)
  pred_max <- tapply(pred, bin_id, max)

  curve <- data.frame(
    bin = factor(names(n_bin), levels = levels(bin_id), ordered = TRUE),
    n = as.integer(n_bin),
    pred_mean = as.numeric(pred_mean),
    obs_rate = as.numeric(obs_rate),
    pred_min = as.numeric(pred_min),
    pred_max = as.numeric(pred_max),
    stringsAsFactors = FALSE
  )

  valid <- is.finite(curve$pred_mean) & is.finite(curve$obs_rate) & curve$n > 0
  n_total <- sum(curve$n[valid])
  if (!n_total) {
    stop("No valid bins available for calibration metrics.", call. = FALSE)
  }
  ece <- sum((curve$n[valid] / n_total) * abs(curve$obs_rate[valid] - curve$pred_mean[valid]))
  mce <- max(abs(curve$obs_rate[valid] - curve$pred_mean[valid]))
  brier <- mean((pred - y)^2, na.rm = TRUE)

  metrics <- data.frame(
    n = n_total,
    bins = bins,
    min_bin_n = min_bin_n,
    ece = ece,
    mce = mce,
    brier = brier,
    stringsAsFactors = FALSE
  )

  list(curve = curve, metrics = metrics, positive_class = pos_class)
}

#' Confounder sensitivity summaries
#'
#' Computes performance metrics within confounder strata to surface potential
#' confounding. Requires aligned metadata in `coldata`.
#'
#' @param fit A [LeakFit] object from [fit_resample()].
#' @param confounders Character vector of columns in `coldata` to evaluate.
#'   Defaults to common batch/study identifiers when available.
#' @param metric Metric name to compute within each stratum. Defaults to the
#'   first metric used in the fit (or task defaults if unavailable).
#' @param min_n Minimum samples per stratum; smaller strata return NA metrics.
#' @param coldata Optional data.frame of sample metadata. Defaults to
#'   `fit@splits@info$coldata` when available.
#' @param numeric_bins Integer number of quantile bins for numeric confounders
#'   with many unique values.
#' @param learner Optional character scalar. When predictions include multiple
#'   learners, selects the learner to summarize.
#' @return A data.frame with per-confounder, per-level metrics and counts.
#' @examples
#' \donttest{
#' fit <- fit_resample(...)
#' confounder_sensitivity(fit, confounders = c("batch", "study"))
#' }
#' @export
confounder_sensitivity <- function(fit, confounders = NULL, metric = NULL,
                                   min_n = 10, coldata = NULL, numeric_bins = 4,
                                   learner = NULL) {
  stopifnot(inherits(fit, "LeakFit"))
  pred_df <- .select_predictions_for_diagnostics(fit, "confounder sensitivity", learner = learner)
  if (!"id" %in% names(pred_df)) {
    stop("Predictions are missing sample ids; confounder sensitivity unavailable.", call. = FALSE)
  }

  if (is.null(coldata)) coldata <- fit@splits@info$coldata
  if (is.null(coldata)) {
    stop("No coldata available for confounder sensitivity.", call. = FALSE)
  }
  coldata <- as.data.frame(coldata, check.names = FALSE)

  ids_chr <- as.character(pred_df$id)
  align_coldata <- function(cd, ids, sample_ids = NULL) {
    rn <- rownames(cd)
    if (!is.null(rn) && all(ids %in% rn)) {
      return(cd[ids, , drop = FALSE])
    }
    if (!is.null(sample_ids) && length(sample_ids) == nrow(cd)) {
      idx <- match(ids, as.character(sample_ids))
      if (all(!is.na(idx))) return(cd[idx, , drop = FALSE])
    }
    if ("row_id" %in% names(cd)) {
      rid <- as.character(cd[["row_id"]])
      if (!anyDuplicated(rid) && all(ids %in% rid)) {
        return(cd[match(ids, rid), , drop = FALSE])
      }
    }
    ids_int <- suppressWarnings(as.integer(ids))
    if (all(!is.na(ids_int)) && max(ids_int, na.rm = TRUE) <= nrow(cd)) {
      return(cd[ids_int, , drop = FALSE])
    }
    if (nrow(cd) == length(ids)) {
      warning("coldata rownames do not match prediction ids; assuming row order aligns to predictions.",
              call. = FALSE)
      return(cd)
    }
    warning("coldata not aligned to predictions; confounder sensitivity skipped.", call. = FALSE)
    NULL
  }

  aligned <- align_coldata(coldata, ids_chr, sample_ids = fit@info$sample_ids)
  if (is.null(aligned)) return(data.frame())

  if (is.null(confounders)) {
    confounders <- intersect(c("batch", "study", "plate", "center", "site", "group"),
                             names(aligned))
  }
  if (!length(confounders)) {
    stop("No confounder columns available in coldata.", call. = FALSE)
  }

  min_n <- as.integer(min_n)
  if (!is.finite(min_n) || min_n < 1L) {
    stop("min_n must be an integer >= 1.", call. = FALSE)
  }
  numeric_bins <- as.integer(numeric_bins)
  if (!is.finite(numeric_bins) || numeric_bins < 2L) {
    stop("numeric_bins must be an integer >= 2.", call. = FALSE)
  }

  default_metric <- fit@info$metrics_used %||% NULL
  if (is.null(metric)) {
    metric <- default_metric[1] %||% {
      if (identical(fit@task, "binomial")) "auc"
      else if (identical(fit@task, "multiclass")) "accuracy"
      else if (identical(fit@task, "survival")) "cindex"
      else "rmse"
    }
  }

  valid_metrics <- switch(fit@task,
                          binomial = c("auc", "roc_auc", "pr_auc", "accuracy", "log_loss"),
                          multiclass = c("accuracy", "macro_f1", "log_loss", "mn_log_loss"),
                          gaussian = c("rmse", "cindex"),
                          survival = c("cindex"),
                          c("auc", "rmse"))
  if (!metric %in% valid_metrics) {
    stop(sprintf("Metric '%s' is not supported for %s tasks.", metric, fit@task), call. = FALSE)
  }

  metric_internal <- metric
  if (identical(metric_internal, "roc_auc")) metric_internal <- "auc"
  if (identical(metric_internal, "mn_log_loss")) metric_internal <- "log_loss"

  direction <- if (metric %in% c("rmse", "log_loss", "mn_log_loss")) "lower" else "higher"

  pos_class <- fit@info$positive_class %||% NULL
  if (identical(fit@task, "binomial")) {
    if (!is.factor(pred_df$truth)) pred_df$truth <- factor(pred_df$truth)
    if (is.null(pos_class) || !as.character(pos_class) %in% levels(pred_df$truth)) {
      pos_class <- levels(pred_df$truth)[2]
    }
  }

  out <- list()
  for (conf in confounders) {
    if (!conf %in% names(aligned)) next
    vec <- aligned[[conf]]
    if (all(is.na(vec))) next
    if (is.numeric(vec)) {
      uniq <- unique(vec[is.finite(vec)])
      if (length(uniq) > numeric_bins) {
        probs <- seq(0, 1, length.out = numeric_bins + 1L)
        brks <- unique(stats::quantile(vec, probs = probs, na.rm = TRUE, names = FALSE))
        if (length(brks) >= 2L) {
          vec <- cut(vec, breaks = brks, include.lowest = TRUE, dig.lab = 10)
        } else {
          vec <- factor(vec)
        }
      } else {
        vec <- factor(vec)
      }
    } else {
      vec <- factor(vec)
    }

    df <- pred_df
    df$confounder_level <- vec
    lvls <- levels(df$confounder_level)
    for (lvl in lvls) {
      sub <- df[df$confounder_level == lvl & !is.na(df$confounder_level), , drop = FALSE]
      n_sub <- nrow(sub)
      value <- NA_real_
      if (n_sub >= min_n) {
        truth_val <- if ("truth" %in% names(sub)) sub$truth else sub$pred
        value <- .metric_value(metric_internal, fit@task, truth_val, sub$pred, pred_df = sub)
      }
      pos_rate <- NA_real_
      if (identical(fit@task, "binomial") && n_sub > 0L) {
        pos_rate <- mean(sub$truth == pos_class, na.rm = TRUE)
      }
      out[[length(out) + 1L]] <- data.frame(
        confounder = conf,
        level = as.character(lvl),
        metric = metric,
        direction = direction,
        value = value,
        n = n_sub,
        positive_rate = pos_rate,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(out)) return(data.frame())
  do.call(rbind, out)
}
