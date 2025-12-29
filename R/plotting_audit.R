# Diagnostic plotting helpers --------------------------------------------------

#' Plot permutation distribution for a LeakAudit object
#'
#' Visualizes the permutation metric distribution and marks the observed and
#' permuted-mean values to help assess leakage signals. Requires ggplot2.
#'
#' @param audit LeakAudit.
#' @return A list containing the observed value, permuted mean, permutation values,
#'   and a ggplot object.
#' @examples
#' \dontrun{
#' # audit <- audit_leakage(fit, metric = "auc", B = 100)
#' plot_perm_distribution(audit)
#' }
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

  hist_obj <- graphics::hist(perm, breaks = "FD", plot = FALSE)
  df <- data.frame(mid = hist_obj$mids, count = hist_obj$counts)
  bin_width <- if (length(hist_obj$breaks) > 1L) diff(hist_obj$breaks)[1] else 1
  line_df <- data.frame(
    value = c(obs, perm_mean),
    type = factor(c("Observed", "Permuted mean"), levels = c("Observed", "Permuted mean"))
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = mid, y = count)) +
    ggplot2::geom_col(width = bin_width, fill = "grey80", color = "white") +
    ggplot2::geom_vline(data = line_df,
                        ggplot2::aes(xintercept = value, color = type, linetype = type),
                        size = 1) +
    ggplot2::labs(title = "Permutation distribution", x = "Metric", y = "Count",
                  color = NULL, linetype = NULL) +
    ggplot2::scale_color_manual(values = c("Observed" = "red", "Permuted mean" = "blue")) +
    ggplot2::scale_linetype_manual(values = c("Observed" = "solid",
                                              "Permuted mean" = "dashed")) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")

  print(p)
  invisible(list(
    observed = obs,
    permuted_mean = perm_mean,
    perm_values = perm,
    plot = p
  ))
}

#' Plot fold balance of positives/negatives per fold
#'
#' Displays a bar chart of class counts per fold and overlays the positive
#' proportion to diagnose stratification issues. The positive class is taken
#' from \code{fit@info$positive_class} when available; otherwise the second
#' factor level is used. Requires ggplot2.
#'
#' @param fit LeakFit.
#' @return A list containing the fold summary, positive class, and a ggplot object.
#' @examples
#' \dontrun{
#' # fit <- fit_resample(...)
#' plot_fold_balance(fit)
#' }
#' @export
plot_fold_balance <- function(fit) {
  stopifnot(inherits(fit, "LeakFit"))
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install it to use plot_fold_balance().",
         call. = FALSE)
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
  print(p)
  invisible(list(
    fold_summary = tab,
    positive_class = pos_class,
    plot = p
  ))
}

#' Plot overlap diagnostics between train/test groups
#'
#' Compares unique group counts per fold and highlights overlaps between train
#' and test groups (e.g., subject or batch IDs). Requires ggplot2.
#'
#' @param fit LeakFit.
#' @param column metadata column name (e.g., subject, batch).
#' @return A list containing overlap counts, the requested column, and a ggplot object.
#' @examples
#' \dontrun{
#' # fit <- fit_resample(...)
#' plot_overlap_checks(fit, column = "subject")
#' }
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
  counts <- lapply(seq_along(fit@splits@indices), function(i) {
    idx <- fit@splits@indices[[i]]
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
  print(p)
  invisible(list(
    overlap_counts = counts,
    column = column,
    plot = p
  ))
}

#' Plot ACF of test predictions for time-series leakage checks
#'
#' Uses the autocorrelation function of out-of-fold predictions to detect
#' temporal dependence that may indicate leakage. Requires ggplot2.
#'
#' @param fit LeakFit.
#' @param lag.max maximum lag to show.
#' @return A list with the autocorrelation results, \code{lag.max}, and a ggplot object.
#' @examples
#' \dontrun{
#' # fit <- fit_resample(...)
#' plot_time_acf(fit, lag.max = 20)
#' }
#' @export
plot_time_acf <- function(fit, lag.max = 20) {
  stopifnot(inherits(fit, "LeakFit"))
  all_pred <- do.call(rbind, fit@predictions)
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting. Install it to use plot_time_acf().",
         call. = FALSE)
  }
  pred <- all_pred$pred
  pred <- pred[is.finite(pred)]
  if (length(pred) < 2) {
    stop("Not enough finite predictions for ACF plotting.", call. = FALSE)
  }
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
  print(p)
  invisible(list(acf = acf_res, lag.max = lag.max, plot = p))
}
