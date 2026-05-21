# Diagnostic plotting helpers --------------------------------------------------

#' Plot permutation distribution for a LeakAudit object
#'
#' Visualizes the label-permutation metric distribution and marks the observed
#' and permuted-mean values to help assess leakage signals. Requires ggplot2.
#'
#' @param audit LeakAudit.
#' @return A list containing the observed value, permuted mean, permutation values,
#'   and a ggplot object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(42)
#'   df <- data.frame(
#'     subject = rep(1:15, each = 2),
#'     outcome = factor(rep(c(0, 1), 15)),
#'     x1 = rnorm(30),
#'     x2 = rnorm(30)
#'   )
#'   splits <- make_split_plan(df, outcome = "outcome",
#'                             mode = "subject_grouped", group = "subject",
#'                             v = 3, progress = FALSE)
#'   custom <- list(
#'     glm = list(
#'       fit = function(x, y, task, weights, ...) {
#'         stats::glm(y ~ ., data = as.data.frame(x),
#'                    family = stats::binomial(), weights = weights)
#'       },
#'       predict = function(object, newdata, task, ...) {
#'         as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
#'                                   type = "response"))
#'       }
#'     )
#'   )
#'   fit <- fit_resample(df, outcome = "outcome", splits = splits,
#'                       learner = "glm", custom_learners = custom,
#'                       metrics = "auc", refit = FALSE, seed = 1)
#'   audit <- audit_leakage(fit, metric = "auc", B = 20)
#'   plot_perm_distribution(audit)
#' }
#'
#' @export
plot_perm_distribution <- function(audit) {
  stopifnot(inherits(audit, "LeakAudit"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install it to use plot_perm_distribution().",
         call. = FALSE)
  }
  perm <- audit@perm_values
  perm <- perm[is.finite(perm)]
  if (!length(perm)) {
    stop("No finite permutation values available for plotting.", call. = FALSE)
  }
  obs <- audit@permutation_gap$metric_obs
  perm_mean <- mean(perm, na.rm = TRUE)

  n_bins <- max(grDevices::nclass.FD(perm), 1L)
  brks <- seq(min(perm), max(perm), length.out = n_bins + 1L)
  hist_obj <- graphics::hist(perm, breaks = brks, plot = FALSE)
  df <- data.frame(mid = hist_obj$mids, count = hist_obj$counts)
  bin_width <- if (length(brks) > 1L) diff(brks)[1] else 1
  line_df <- data.frame(
    value = c(obs, perm_mean),
    type = factor(c("Observed", "Permuted mean"), levels = c("Observed", "Permuted mean"))
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = mid, y = count)) +
    ggplot2::geom_col(width = bin_width, fill = "grey80", color = "white") +
    ggplot2::geom_vline(data = line_df,
                        ggplot2::aes(xintercept = value, color = type, linetype = type),
                        linewidth = 1) +
    ggplot2::labs(title = "Permutation distribution", x = "Metric", y = "Count",
                  color = NULL, linetype = NULL) +
    ggplot2::scale_color_manual(values = c("Observed" = "red", "Permuted mean" = "blue")) +
    ggplot2::scale_linetype_manual(values = c("Observed" = "solid",
                                              "Permuted mean" = "dashed")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")

  if (interactive()) print(p)
  invisible(list(
    observed = obs,
    permuted_mean = perm_mean,
    perm_values = perm,
    plot = p
  ))
}

#' Plot fold balance of class counts per fold
#'
#' Displays a bar chart of class counts per fold. For binomial tasks, it also
#' overlays the positive proportion to diagnose stratification issues. The
#' positive class is taken from \code{fit@info$positive_class} when available;
#' otherwise the second factor level is used. For multiclass tasks, the plot
#' shows per-class counts without a proportion line. Only available for
#' classification tasks. Requires ggplot2.
#'
#' @param fit LeakFit.
#' @return A list containing the fold summary, positive class (if binomial),
#'   and a ggplot object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(42)
#'   df <- data.frame(
#'     subject = rep(1:15, each = 2),
#'     outcome = factor(rep(c(0, 1), 15)),
#'     x1 = rnorm(30),
#'     x2 = rnorm(30)
#'   )
#'   splits <- make_split_plan(df, outcome = "outcome",
#'                             mode = "subject_grouped", group = "subject",
#'                             v = 3, progress = FALSE)
#'   custom <- list(
#'     glm = list(
#'       fit = function(x, y, task, weights, ...) {
#'         stats::glm(y ~ ., data = as.data.frame(x),
#'                    family = stats::binomial(), weights = weights)
#'       },
#'       predict = function(object, newdata, task, ...) {
#'         as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
#'                                   type = "response"))
#'       }
#'     )
#'   )
#'   fit <- fit_resample(df, outcome = "outcome", splits = splits,
#'                       learner = "glm", custom_learners = custom,
#'                       metrics = "auc", refit = FALSE, seed = 1)
#'   plot_fold_balance(fit)
#' }
#'
#' @export
plot_fold_balance <- function(fit) {
  stopifnot(inherits(fit, "LeakFit"))
  if (!fit@task %in% c("binomial", "multiclass")) {
    stop("plot_fold_balance is only available for classification tasks.", call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install it to use plot_fold_balance().",
         call. = FALSE)
  }
  if (fit@task == "multiclass") {
    class_levels <- NULL
    for (df in fit@predictions) {
      if (is.factor(df$truth)) {
        class_levels <- levels(df$truth)
        break
      }
    }
    if (is.null(class_levels)) {
      class_levels <- sort(unique(unlist(lapply(fit@predictions, function(df) {
        as.character(df$truth)
      }))))
    }
    if (!length(class_levels)) {
      stop("No class labels available for plotting.", call. = FALSE)
    }
    tab <- lapply(seq_along(fit@predictions), function(i) {
      df <- fit@predictions[[i]]
      y <- factor(as.character(df$truth), levels = class_levels)
      counts <- table(y)
      data.frame(
        fold = i,
        class = factor(names(counts), levels = class_levels),
        count = as.numeric(counts),
        stringsAsFactors = FALSE
      )
    })
    tab <- do.call(rbind, tab)
    p <- ggplot2::ggplot(tab, ggplot2::aes(x = fold, y = count, fill = class)) +
      ggplot2::geom_col(width = 0.7, color = "white") +
      ggplot2::scale_x_continuous(breaks = unique(tab$fold)) +
      ggplot2::labs(title = "Fold class balance", x = "Fold", y = "Count", fill = "Class") +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = "top")
    if (interactive()) print(p)
    return(invisible(list(
      fold_summary = tab,
      positive_class = NA_character_,
      plot = p
    )))
  }
  pos_class <- fit@info$positive_class
  if (length(pos_class) != 1L) pos_class <- NULL
  pos_class <- if (!is.null(pos_class)) as.character(pos_class) else NULL
  if (!is.null(pos_class) && (is.na(pos_class) || !nzchar(pos_class))) {
    pos_class <- NULL
  }
  if (is.null(pos_class)) {
    for (df in fit@predictions) {
      if (is.factor(df$truth) && nlevels(df$truth) >= 2) {
        pos_class <- levels(df$truth)[2]
        break
      }
    }
  }
  resolve_pos_label <- function(y, pos_label) {
    if (!is.null(pos_label)) return(pos_label)
    if (is.factor(y) && nlevels(y) >= 2) return(levels(y)[2])
    "1"
  }
  tab <- lapply(seq_along(fit@predictions), function(i) {
    df <- fit@predictions[[i]]
    y <- df$truth
    pos_label <- resolve_pos_label(y, pos_class)
    y_chr <- as.character(y)
    is_pos <- y_chr == pos_label
    valid <- !is.na(is_pos)
    data.frame(fold = i,
               positives = sum(is_pos[valid]),
               negatives = sum(!is_pos[valid]))
  })
  tab <- do.call(rbind, tab)
  totals <- tab$positives + tab$negatives
  tab$prop_pos <- ifelse(totals > 0, tab$positives / totals, NA_real_)
  pos_legend <- if (!is.null(pos_class)) {
    paste0("Positives (", pos_class, ")")
  } else {
    "Positives"
  }
  df_counts <- data.frame(
    fold = rep(tab$fold, times = 2),
    class = factor(
      rep(c(pos_legend, "Negatives"), each = nrow(tab)),
      levels = c(pos_legend, "Negatives")
    ),
    count = c(tab$positives, tab$negatives)
  )
  max_count <- max(df_counts$count, na.rm = TRUE)
  if (!is.finite(max_count) || max_count <= 0) max_count <- 1
  prop_df <- data.frame(
    fold = tab$fold,
    prop_pos = tab$prop_pos,
    prop_scaled = tab$prop_pos * max_count
  )
  prop_df <- prop_df[is.finite(prop_df$prop_scaled), , drop = FALSE]
  prop_df$series <- "Positive proportion"

  p <- ggplot2::ggplot(df_counts, ggplot2::aes(x = fold, y = count, fill = class)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                      width = 0.7, color = "white") +
    ggplot2::scale_fill_manual(values = setNames(c("steelblue", "tan"),
                                                 c(pos_legend, "Negatives"))) +
    ggplot2::scale_x_continuous(breaks = tab$fold) +
    ggplot2::labs(title = "Fold class balance", x = "Fold", y = "Count", fill = NULL)
  if (nrow(prop_df) > 0) {
    p <- p +
      ggplot2::geom_line(data = prop_df,
                         ggplot2::aes(x = fold, y = prop_scaled,
                                      color = series, linetype = series),
                         inherit.aes = FALSE) +
      ggplot2::geom_point(data = prop_df,
                          ggplot2::aes(x = fold, y = prop_scaled, color = series, shape = series),
                          inherit.aes = FALSE) +
      ggplot2::scale_color_manual(values = c("Positive proportion" = "blue")) +
      ggplot2::scale_linetype_manual(values = c("Positive proportion" = "dashed")) +
      ggplot2::scale_shape_manual(values = c("Positive proportion" = 17)) +
      ggplot2::scale_y_continuous(
        sec.axis = ggplot2::sec_axis(~ . / max_count,
                                     name = "Positive proportion",
                                     labels = function(x) sprintf("%s%%", round(x * 100)))
      )
  }
  p <- p +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top") +
    ggplot2::guides(fill = ggplot2::guide_legend(order = 1),
                    color = ggplot2::guide_legend(order = 2),
                    linetype = ggplot2::guide_legend(order = 2),
                    shape = ggplot2::guide_legend(order = 2))
  if (interactive()) print(p)
  invisible(list(
    fold_summary = tab,
    positive_class = pos_class,
    plot = p
  ))
}

#' Plot overlap diagnostics between train/test groups
#'
#' @description
#' Checks whether the same group identifiers appear in both the training and
#' test partitions within each resample. This is designed to detect leakage
#' from grouped or repeated-measures data (for example, the same subject,
#' batch, plate, or study appearing on both sides of a fold) when group-wise
#' splitting is expected.
#'
#' @details
#' For each resample in `fit@splits@indices`, the function counts the number of
#' unique values of `column` in the train and test sets and the size of their
#' intersection. Any non-zero overlap indicates that at least one group appears
#' in both train and test for that resample. The check is metadata-based only:
#' it relies on exact matches of the supplied column and does not inspect
#' features or outcomes. It only checks train vs test within each resample, so
#' it will not detect overlaps across different resamples or other leakage
#' mechanisms. Inconsistent IDs or missing values in the metadata can hide or
#' inflate overlaps. `NA` values are treated as regular identifiers and will
#' count toward overlap if they appear in both partitions. Requires ggplot2.
#'
#' @param fit A `LeakFit` object produced by [fit_resample()]. It must contain
#'   the split indices and the associated metadata in `fit@splits@info$coldata`.
#'   The metadata rows must align with the data used to create the splits.
#' @param column Character scalar naming the metadata column to check (for
#'   example `"subject"` or `"batch"`). The function compares unique values of
#'   this column between train and test within each resample. There is no
#'   default: `NULL` or an unknown column triggers an error. Changing `column`
#'   changes which kind of leakage (subject-level, batch-level, etc.) is tested
#'   and therefore the overlap counts.
#' @return A list returned invisibly with:
#'   \itemize{
#'     \item `overlap_counts`: data.frame with one row per resample and columns
#'       `fold` (resample index in `fit@splits@indices`), `overlap` (unique IDs
#'       shared by train and test), `train` (unique IDs in train), and `test`
#'       (unique IDs in test).
#'     \item `column`: the metadata column name used for the check.
#'     \item `plot`: the ggplot object showing the three count series across folds.
#'   }
#'   The plot is also printed. When any overlap is detected, the plot adds a
#'   warning annotation.
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   subject = rep(1:6, each = 2),
#'   outcome = rbinom(12, 1, 0.5),
#'   x1 = rnorm(12),
#'   x2 = rnorm(12)
#' )
#' splits <- make_split_plan(df, outcome = "outcome",
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
#'                     metrics = "accuracy", refit = FALSE)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   out <- plot_overlap_checks(fit, column = "subject")
#'   out$overlap_counts
#' }
#'
#' @export
plot_overlap_checks <- function(fit, column = NULL) {
  stopifnot(inherits(fit, "LeakFit"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install it to use plot_overlap_checks().",
         call. = FALSE)
  }
  cd <- fit@splits@info$coldata
  if (is.null(cd) || is.null(column) || !column %in% names(cd)) {
    stop("Column not available in metadata.")
  }
  n <- nrow(cd)
  counts <- lapply(seq_along(fit@splits@indices), function(i) {
    idx <- fit@splits@indices[[i]]
    idx <- .bio_resolve_fold_indices(fit@splits, idx, n = n, data = cd)
    tr <- unique(cd[[column]][idx$train])
    te <- unique(cd[[column]][idx$test])
    data.frame(fold = i, overlap = length(intersect(tr, te)),
               train = length(tr), test = length(te))
  })
  counts <- do.call(rbind, counts)
  plot_df <- data.frame(
    fold = rep(counts$fold, times = 3),
    metric = factor(rep(c("Overlap", "Train unique", "Test unique"),
                        each = nrow(counts)),
                    levels = c("Overlap", "Train unique", "Test unique")),
    count = c(counts$overlap, counts$train, counts$test)
  )
  p <- ggplot2::ggplot(plot_df,
                       ggplot2::aes(x = fold, y = count, color = metric, linetype = metric)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = c("Overlap" = "red",
                                           "Train unique" = "grey40",
                                           "Test unique" = "grey70")) +
    ggplot2::scale_linetype_manual(values = c("Overlap" = "solid",
                                              "Train unique" = "dashed",
                                              "Test unique" = "dotted")) +
    ggplot2::scale_x_continuous(breaks = counts$fold) +
    ggplot2::labs(title = sprintf("Overlap diagnostics: %s", column),
                  x = "Fold", y = "Count", color = NULL, linetype = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")
  if (any(counts$overlap > 0)) {
    y_max <- max(plot_df$count, na.rm = TRUE)
    if (!is.finite(y_max)) y_max <- 1
    p <- p +
      ggplot2::expand_limits(y = y_max * 1.1) +
      ggplot2::annotate("text",
                        x = max(counts$fold),
                        y = y_max * 1.05,
                        label = "WARNING: Overlaps detected!",
                        hjust = 1, vjust = 0,
                        color = "red", fontface = "bold")
  }
  if (interactive()) print(p)
  invisible(list(
    overlap_counts = counts,
    column = column,
    plot = p
  ))
}

#' Plot ACF of test predictions for time-series leakage checks
#'
#' Uses the autocorrelation function of out-of-fold predictions to detect
#' temporal dependence that may indicate leakage. Predictions are ordered by
#' the split time column before computing the ACF. Requires numeric predictions
#' (regression or survival). Requires ggplot2.
#'
#' @param fit LeakFit.
#' @param lag.max maximum lag to show.
#' @return A list with the autocorrelation results, \code{lag.max}, and a ggplot object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(42)
#'   df <- data.frame(
#'     id = 1:30,
#'     time = seq.Date(as.Date("2020-01-01"), by = "day", length.out = 30),
#'     y = rnorm(30),
#'     x1 = rnorm(30),
#'     x2 = rnorm(30)
#'   )
#'   splits <- make_split_plan(df, outcome = "y", mode = "time_series",
#'                             time = "time", v = 3, progress = FALSE)
#'   custom <- list(
#'     lm = list(
#'       fit = function(x, y, task, weights, ...) {
#'         stats::lm(y ~ ., data = data.frame(y = y, x))
#'       },
#'       predict = function(object, newdata, task, ...) {
#'         as.numeric(stats::predict(object, newdata = as.data.frame(newdata)))
#'       }
#'     )
#'   )
#'   fit <- fit_resample(df, outcome = "y", splits = splits,
#'                       learner = "lm", custom_learners = custom,
#'                       metrics = "rmse", refit = FALSE, seed = 1)
#'   plot_time_acf(fit, lag.max = 10)
#' }
#'
#' @export
plot_time_acf <- function(fit, lag.max = 20) {
  stopifnot(inherits(fit, "LeakFit"))
  all_pred <- do.call(rbind, fit@predictions)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install it to use plot_time_acf().",
         call. = FALSE)
  }
  if (!"pred" %in% names(all_pred)) {
    stop("Predictions missing 'pred' column for ACF plotting.", call. = FALSE)
  }
  pred <- all_pred$pred
  if (!is.numeric(pred)) {
    stop("plot_time_acf requires numeric predictions (regression or survival).", call. = FALSE)
  }
  if (!"id" %in% names(all_pred)) {
    stop("Predictions are missing sample ids; time ordering unavailable.", call. = FALSE)
  }
  time_col <- fit@splits@info$time %||% NULL
  coldata <- fit@splits@info$coldata %||% NULL
  if (is.null(coldata)) {
    stop("plot_time_acf requires split metadata with time values.", call. = FALSE)
  }
  coldata <- as.data.frame(coldata, check.names = FALSE)
  if (is.null(time_col) || !time_col %in% names(coldata)) {
    stop("plot_time_acf requires a time column in split metadata.", call. = FALSE)
  }

  ids_chr <- as.character(all_pred$id)
  time_vals <- NULL
  rn <- rownames(coldata)
  if (!is.null(rn) && all(ids_chr %in% rn)) {
    time_vals <- coldata[match(ids_chr, rn), time_col]
  } else if ("row_id" %in% names(coldata)) {
    rid <- as.character(coldata[["row_id"]])
    if (!anyDuplicated(rid) && all(ids_chr %in% rid)) {
      time_vals <- coldata[match(ids_chr, rid), time_col]
    }
  } else if (!is.null(fit@info$sample_ids) &&
             length(fit@info$sample_ids) == nrow(coldata)) {
    sample_ids <- as.character(fit@info$sample_ids)
    idx <- match(ids_chr, sample_ids)
    if (all(!is.na(idx))) time_vals <- coldata[idx, time_col]
  } else {
    ids_int <- suppressWarnings(as.integer(ids_chr))
    if (all(!is.na(ids_int)) && max(ids_int, na.rm = TRUE) <= nrow(coldata)) {
      time_vals <- coldata[ids_int, time_col]
    }
  }
  if (is.null(time_vals) || length(time_vals) != length(pred)) {
    stop("plot_time_acf could not align time metadata to predictions.", call. = FALSE)
  }
  if (!is.numeric(time_vals) && !inherits(time_vals, c("POSIXct", "Date"))) {
    stop("plot_time_acf requires numeric, Date, or POSIXct time values.", call. = FALSE)
  }

  ok_time <- !is.na(time_vals)
  if (!all(ok_time)) {
    warning("plot_time_acf dropped predictions with missing time values.", call. = FALSE)
  }
  pred <- pred[ok_time]
  time_vals <- time_vals[ok_time]
  ok_pred <- is.finite(pred)
  pred <- pred[ok_pred]
  time_vals <- time_vals[ok_pred]
  if (length(pred) < 2) {
    stop("Not enough finite predictions for ACF plotting.", call. = FALSE)
  }
  ord <- order(time_vals, seq_along(time_vals), na.last = TRUE)
  pred <- pred[ord]
  acf_res <- stats::acf(pred, lag.max = lag.max, plot = FALSE)
  lag_vals <- as.numeric(acf_res$lag)
  acf_vals <- as.numeric(acf_res$acf)
  df <- data.frame(lag = lag_vals, acf = acf_vals)
  conf <- if (is.finite(acf_res$n.used) && acf_res$n.used > 0) {
    1.96 / sqrt(acf_res$n.used)
  } else {
    NA_real_
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = lag, y = acf)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey50") +
    ggplot2::geom_segment(ggplot2::aes(xend = lag, yend = 0), color = "steelblue") +
    ggplot2::geom_point(color = "steelblue") +
    ggplot2::labs(title = "Prediction autocorrelation", x = "Lag", y = "ACF") +
    ggplot2::theme_minimal()
  if (is.finite(conf)) {
    p <- p +
      ggplot2::geom_hline(yintercept = c(conf, -conf),
                          color = "red", linetype = "dashed")
  }
  if (interactive()) print(p)
  invisible(list(acf = acf_res, lag.max = lag.max, plot = p))
}

#' Plot calibration curve for binomial predictions
#'
#' Visualizes observed outcome rates versus predicted probabilities across
#' bins to diagnose calibration (binomial tasks only). Requires ggplot2.
#'
#' @param fit LeakFit.
#' @param bins Number of probability bins to use.
#' @param min_bin_n Minimum samples per bin shown in the plot.
#' @param learner Optional character scalar. When predictions include multiple
#'   learners, selects the learner to summarize.
#' @return A list containing the calibration curve, metrics, and a ggplot object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(42)
#'   df <- data.frame(
#'     subject = rep(1:15, each = 2),
#'     outcome = factor(rep(c(0, 1), 15)),
#'     x1 = rnorm(30),
#'     x2 = rnorm(30)
#'   )
#'   splits <- make_split_plan(df, outcome = "outcome",
#'                             mode = "subject_grouped", group = "subject",
#'                             v = 3, progress = FALSE)
#'   custom <- list(
#'     glm = list(
#'       fit = function(x, y, task, weights, ...) {
#'         stats::glm(y ~ ., data = as.data.frame(x),
#'                    family = stats::binomial(), weights = weights)
#'       },
#'       predict = function(object, newdata, task, ...) {
#'         as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
#'                                   type = "response"))
#'       }
#'     )
#'   )
#'   fit <- fit_resample(df, outcome = "outcome", splits = splits,
#'                       learner = "glm", custom_learners = custom,
#'                       metrics = "auc", refit = FALSE, seed = 1)
#'   plot_calibration(fit, bins = 5)
#' }
#'
#' @export
plot_calibration <- function(fit, bins = 10, min_bin_n = 5, learner = NULL) {
  stopifnot(inherits(fit, "LeakFit"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install it to use plot_calibration().",
         call. = FALSE)
  }
  cal <- calibration_summary(fit, bins = bins, min_bin_n = min_bin_n, learner = learner)
  df <- cal$curve
  df <- df[is.finite(df$pred_mean) & is.finite(df$obs_rate), , drop = FALSE]
  df_plot <- df[df$n >= min_bin_n, , drop = FALSE]
  if (!nrow(df_plot)) {
    stop("Not enough data for calibration plotting.", call. = FALSE)
  }
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = pred_mean, y = obs_rate)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                         color = "grey50") +
    ggplot2::geom_line(color = "steelblue") +
    ggplot2::geom_point(ggplot2::aes(size = n), color = "steelblue") +
    ggplot2::scale_size_continuous(range = c(2, 6)) +
    ggplot2::labs(title = "Calibration curve",
                  x = "Mean predicted probability",
                  y = "Observed event rate",
                  size = "Bin n") +
    ggplot2::theme_minimal()
  if (interactive()) print(p)
  invisible(list(curve = df, metrics = cal$metrics, plot = p))
}

#' Plot confounder sensitivity
#'
#' Shows performance metrics across confounder strata to assess sensitivity to
#' batch/study effects. Requires ggplot2.
#'
#' @param fit LeakFit.
#' @param confounders Character vector of columns in `coldata` to evaluate.
#' @param metric Metric name to compute within each stratum.
#' @param min_n Minimum samples per stratum to display.
#' @param coldata Optional data.frame of sample metadata.
#' @param numeric_bins Number of quantile bins for numeric confounders.
#' @param learner Optional character scalar. When predictions include multiple
#'   learners, selects the learner to summarize.
#' @return A list containing the sensitivity table and a ggplot object.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(42)
#'   df <- data.frame(
#'     subject = rep(1:15, each = 2),
#'     outcome = factor(rep(c(0, 1), 15)),
#'     batch = factor(rep(c("A", "B", "C"), 10)),
#'     x1 = rnorm(30),
#'     x2 = rnorm(30)
#'   )
#'   splits <- make_split_plan(df, outcome = "outcome",
#'                             mode = "subject_grouped", group = "subject",
#'                             v = 3, progress = FALSE)
#'   custom <- list(
#'     glm = list(
#'       fit = function(x, y, task, weights, ...) {
#'         stats::glm(y ~ ., data = as.data.frame(x),
#'                    family = stats::binomial(), weights = weights)
#'       },
#'       predict = function(object, newdata, task, ...) {
#'         as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
#'                                   type = "response"))
#'       }
#'     )
#'   )
#'   fit <- fit_resample(df, outcome = "outcome", splits = splits,
#'                       learner = "glm", custom_learners = custom,
#'                       metrics = "auc", refit = FALSE, seed = 1)
#'   plot_confounder_sensitivity(fit, confounders = "batch", coldata = df)
#' }
#'
#' @export
plot_confounder_sensitivity <- function(fit, confounders = NULL, metric = NULL,
                                        min_n = 10, coldata = NULL, numeric_bins = 4,
                                        learner = NULL) {
  stopifnot(inherits(fit, "LeakFit"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install it to use plot_confounder_sensitivity().",
         call. = FALSE)
  }
  df <- confounder_sensitivity(fit, confounders = confounders, metric = metric,
                               min_n = min_n, coldata = coldata,
                               numeric_bins = numeric_bins, learner = learner)
  if (is.null(df) || !nrow(df)) {
    stop("No confounder sensitivity results available for plotting.", call. = FALSE)
  }
  df_plot <- df[is.finite(df$value) & df$n >= min_n, , drop = FALSE]
  if (!nrow(df_plot)) {
    stop("No strata meet the minimum size for plotting.", call. = FALSE)
  }
  direction <- df_plot$direction[1] %||% "higher"
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = level, y = value)) +
    ggplot2::geom_col(fill = "steelblue") +
    ggplot2::facet_wrap(~ confounder, scales = "free_x") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Confounder sensitivity",
                  subtitle = sprintf("%s (better is %s)", df_plot$metric[1], direction),
                  x = NULL, y = "Metric value") +
    ggplot2::theme_minimal()
  if (interactive()) print(p)
  invisible(list(data = df, plot = p))
}


# =============================================================================
# Delta-LSI repeat-level diagnostic (matches manuscript Figure 4 panel b)
# =============================================================================

#' Plot per-repeat \eqn{\Delta_r} values from a LeakDeltaLSI object
#'
#' Visualises the per-repeat metric differences (leaky minus guarded) for a
#' \code{\linkS4class{LeakDeltaLSI}} object, overlaid with the robust Huber
#' point estimate, the arithmetic mean, and the BCa bootstrap confidence
#' interval. This is the diagnostic shown as Figure 4 panel (b) of the
#' manuscript. Requires ggplot2.
#'
#' @param dlsi A \code{\linkS4class{LeakDeltaLSI}} object produced by
#'   \code{\link{delta_lsi}}.
#' @return A list with the per-repeat deltas, the robust and arithmetic-mean
#'   estimates, the BCa confidence interval, and the ggplot object.
#' @seealso \code{\link{delta_lsi}}, \code{\link{dlsi_repeats}}
#' @export
plot_dlsi_repeats <- function(dlsi) {
  stopifnot(inherits(dlsi, "LeakDeltaLSI"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install it to use plot_dlsi_repeats().",
         call. = FALSE)
  }
  rn <- dlsi_repeats(dlsi, "naive")
  rg <- dlsi_repeats(dlsi, "guarded")
  if (is.null(rn) || is.null(rg) ||
      nrow(rn) == 0L || nrow(rg) == 0L ||
      !"metric" %in% names(rn) || !"metric" %in% names(rg)) {
    stop("No per-repeat metric data available for plotting.", call. = FALSE)
  }
  R_eff   <- min(nrow(rn), nrow(rg))
  deltas  <- rn$metric[seq_len(R_eff)] - rg$metric[seq_len(R_eff)]
  delta_robust <- dlsi_robust(dlsi)
  delta_mean   <- dlsi_metric(dlsi)
  ci           <- dlsi_ci(dlsi, which = "robust")

  df <- data.frame(repeat_idx = seq_len(R_eff), delta = deltas)
  ci_lo <- if (length(ci) == 2L && all(is.finite(ci))) ci[1] else NA_real_
  ci_hi <- if (length(ci) == 2L && all(is.finite(ci))) ci[2] else NA_real_

  p <- ggplot2::ggplot(df, ggplot2::aes(x = repeat_idx, y = delta))
  if (is.finite(ci_lo) && is.finite(ci_hi)) {
    p <- p + ggplot2::annotate("rect",
      xmin = 0.5, xmax = R_eff + 0.5,
      ymin = ci_lo, ymax = ci_hi,
      fill = "#0072B2", alpha = 0.12)
  }
  p <- p +
    ggplot2::geom_hline(yintercept = 0, color = "gray70", linetype = "dotted") +
    ggplot2::geom_hline(yintercept = delta_robust, color = "#0072B2",
                        linewidth = 0.9) +
    ggplot2::geom_hline(yintercept = delta_mean, color = "#0072B2",
                        linewidth = 0.6, linetype = "dashed") +
    ggplot2::geom_point(color = "gray40", size = 2) +
    ggplot2::labs(
      title = expression("Per-repeat" ~ Delta[italic(r)] ~
                         "(leaky " * minus * " guarded)"),
      subtitle = sprintf(
        "R_eff = %d   |   Delta_LSI = %.3f   |   Delta_metric = %.3f",
        R_eff, delta_robust, delta_mean),
      x = "Repeat", y = expression(Delta[italic(r)])) +
    ggplot2::theme_minimal()

  if (interactive()) print(p)
  invisible(list(
    deltas         = deltas,
    delta_robust   = delta_robust,
    delta_mean     = delta_mean,
    bca_ci         = c(ci_lo, ci_hi),
    R_eff          = R_eff,
    plot           = p
  ))
}


# =============================================================================
# S4 plot() methods for the bioLeak result classes
# =============================================================================
#
# Each method dispatches plot() on the class to the canonical diagnostic
# helper, so that users can call `plot(audit)`, `plot(fit)`, or `plot(dlsi)`
# without remembering the helper name. The helpers themselves are kept
# exported so power users can call them by name and pass additional
# arguments.
#
# Method signatures match base::plot(x, y, ...) so they integrate cleanly
# with the standard R plotting idiom.

#' @include classes.R accessors.R
NULL

#' Plot method for LeakAudit
#'
#' Diagnostic plot for a [`LeakAudit`] object. The default diagnostic is the
#' permutation-distribution histogram produced by
#' \code{\link{plot_perm_distribution}}.
#'
#' @param x A [`LeakAudit`] object.
#' @param y Unused; present for S4 compatibility with \code{base::plot}.
#' @param ... Additional arguments passed to \code{\link{plot_perm_distribution}}.
#' @return Invisibly returns the list produced by
#'   \code{\link{plot_perm_distribution}} (observed value, permuted mean,
#'   permutation values, ggplot object).
#' @seealso \code{\link{plot_perm_distribution}},
#'   \code{\linkS4class{LeakAudit}}
#' @docType methods
#' @aliases plot,LeakAudit-method
#' @export
setMethod("plot", signature(x = "LeakAudit", y = "missing"),
          function(x, y, ...) {
            plot_perm_distribution(x, ...)
          })

#' Plot method for LeakFit
#'
#' Diagnostic plot for a [`LeakFit`] object. The default diagnostic is the
#' fold-balance check produced by \code{\link{plot_fold_balance}}, which
#' works for any classification task. Use the \code{which} argument to
#' switch to one of the other diagnostics in the package:
#' \code{"overlap"} (\code{\link{plot_overlap_checks}}),
#' \code{"calibration"} (\code{\link{plot_calibration}}; binary outcomes
#' only), \code{"time_acf"} (\code{\link{plot_time_acf}}; time-ordered
#' splits), or \code{"confounder_sensitivity"}
#' (\code{\link{plot_confounder_sensitivity}}).
#'
#' @param x A [`LeakFit`] object.
#' @param y Unused; present for S4 compatibility with \code{base::plot}.
#' @param which One of \code{"fold_balance"} (default),
#'   \code{"overlap"}, \code{"calibration"}, \code{"time_acf"},
#'   \code{"confounder_sensitivity"}.
#' @param ... Additional arguments passed to the selected helper.
#' @return Invisibly returns the list produced by the selected helper.
#' @seealso \code{\link{plot_fold_balance}}, \code{\link{plot_overlap_checks}},
#'   \code{\link{plot_calibration}}, \code{\link{plot_time_acf}},
#'   \code{\link{plot_confounder_sensitivity}},
#'   \code{\linkS4class{LeakFit}}
#' @docType methods
#' @aliases plot,LeakFit-method
#' @export
setMethod("plot", signature(x = "LeakFit", y = "missing"),
          function(x, y,
                   which = c("fold_balance", "overlap",
                             "calibration", "time_acf",
                             "confounder_sensitivity"),
                   ...) {
            which <- match.arg(which)
            switch(which,
              fold_balance           = plot_fold_balance(x, ...),
              overlap                = plot_overlap_checks(x, ...),
              calibration            = plot_calibration(x, ...),
              time_acf               = plot_time_acf(x, ...),
              confounder_sensitivity = plot_confounder_sensitivity(x, ...)
            )
          })

#' Plot method for LeakDeltaLSI
#'
#' Diagnostic plot for a [`LeakDeltaLSI`] object: per-repeat \eqn{\Delta_r}
#' scatter with the Huber-robust point estimate, the arithmetic mean,
#' and the BCa bootstrap confidence interval band. This is the diagnostic
#' shown as Figure 4 panel (b) of the manuscript.
#'
#' @param x A [`LeakDeltaLSI`] object.
#' @param y Unused; present for S4 compatibility with \code{base::plot}.
#' @param ... Additional arguments (currently unused).
#' @return Invisibly returns the list produced by
#'   \code{\link{plot_dlsi_repeats}}: per-repeat deltas, the robust and
#'   arithmetic-mean estimates, the BCa interval, and the ggplot object.
#' @seealso \code{\link{plot_dlsi_repeats}},
#'   \code{\linkS4class{LeakDeltaLSI}}
#' @docType methods
#' @aliases plot,LeakDeltaLSI-method
#' @export
setMethod("plot", signature(x = "LeakDeltaLSI", y = "missing"),
          function(x, y, ...) {
            plot_dlsi_repeats(x, ...)
          })
