# Diagnostic plotting helpers --------------------------------------------------

#' Plot permutation distribution for a LeakAudit object
#' @param audit LeakAudit
#' @export
plot_perm_distribution <- function(audit) {
  stopifnot(inherits(audit, "LeakAudit"))
  perm <- audit@perm_distribution
  obs <- audit@permutation_gap$metric_obs
  graphics::hist(perm, breaks = "FD", col = "grey80", border = "white",
                 main = "Permutation distribution", xlab = "Metric")
  graphics::abline(v = obs, col = "red", lwd = 2)
  perm_mean <- mean(perm, na.rm = TRUE)
  graphics::abline(v = perm_mean, col = "blue", lty = 2, lwd = 2)
  graphics::mtext(sprintf("obs = %.3f", obs), side = 3, col = "red")
  graphics::legend("topright",
                   legend = c("Observed", "Permuted mean"),
                   col = c("red", "blue"),
                   lwd = 2,
                   lty = c(1, 2),
                   bg = "white")
}

#' Plot fold balance of positives/negatives per fold
#' @param fit LeakFit
#' @export
plot_fold_balance <- function(fit) {
  stopifnot(inherits(fit, "LeakFit"))
  tab <- lapply(seq_along(fit@predictions), function(i) {
    df <- fit@predictions[[i]]
    y <- if (is.factor(df$truth)) as.numeric(df$truth) - 1 else as.numeric(df$truth)
    data.frame(fold = i, positives = sum(y == 1, na.rm = TRUE),
               negatives = sum(y == 0, na.rm = TRUE))
  })
  tab <- do.call(rbind, tab)
  mat <- t(as.matrix(tab[, c("positives", "negatives")]))
  bp <- graphics::barplot(mat,
                          beside = TRUE,
                          col = c("steelblue", "tan"),
                          border = "white",
                          ylim = c(0, max(mat, na.rm = TRUE) * 1.1),
                          main = "Fold class balance",
                          xlab = "Fold",
                          ylab = "Count")
  totals <- tab$positives + tab$negatives
  tab$prop_pos <- ifelse(totals > 0, tab$positives / totals, NA_real_)
  legend_entries <- c("Positives", "Negatives")
  legend_fill <- c("steelblue", "tan")
  legend_border <- c("white", "white")
  legend_lty <- c(NA, NA)
  legend_pch <- c(NA, NA)
  legend_col <- c("steelblue", "tan")
  if (any(!is.na(tab$prop_pos))) {
    fold_positions <- colMeans(bp, na.rm = TRUE)
    graphics::par(new = TRUE)
    graphics::plot(fold_positions, tab$prop_pos,
                   type = "b", axes = FALSE, xlab = "", ylab = "",
                   col = "blue", pch = 17, lty = 2)
    axis_vals <- pretty(c(0, 1, tab$prop_pos), n = 5)
    graphics::axis(4, at = axis_vals,
                   labels = sprintf("%s%%", round(axis_vals * 100)))
    graphics::mtext("Positive proportion", side = 4, line = 3)
    legend_entries <- c(legend_entries, "Positive proportion")
    legend_fill <- c(legend_fill, NA)
    legend_border <- c(legend_border, NA)
    legend_lty <- c(legend_lty, 2)
    legend_pch <- c(legend_pch, 17)
    legend_col <- c(legend_col, "blue")
  }
  graphics::legend("topright",
                   legend = legend_entries,
                   fill = legend_fill,
                   border = legend_border,
                   lty = legend_lty,
                   pch = legend_pch,
                   col = legend_col,
                   pt.bg = rep(NA, length(legend_entries)),
                   bg = "white")
}

#' Plot overlap diagnostics between train/test groups
#' @param fit LeakFit
#' @param column metadata column name (e.g., subject, batch)
#' @export
plot_overlap_checks <- function(fit, column = NULL) {
  stopifnot(inherits(fit, "LeakFit"))
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
  graphics::matplot(counts$fold, counts[, c("overlap", "train", "test")],
                    type = "b", pch = 16,
                    lty = c(1, 2, 3), col = c("red", "grey40", "grey70"),
                    xlab = "Fold", ylab = "Count",
                    main = sprintf("Overlap diagnostics: %s", column))
  graphics::legend("topright", legend = c("Overlap", "Train unique", "Test unique"),
                   col = c("red", "grey40", "grey70"), lty = c(1, 2, 3), pch = 16)
  if (any(counts$overlap > 0)) {
    graphics::mtext("WARNING: Overlaps detected!", side = 3, line = 0.5,
                    col = "red", font = 2)
  }
}

#' Plot ACF of test predictions for time-series leakage checks
#' @param fit LeakFit
#' @export
plot_time_acf <- function(fit, lag.max = NULL) {
  stopifnot(inherits(fit, "LeakFit"))
  all_pred <- do.call(rbind, fit@predictions)
  graphics::acf(all_pred$pred, main = "Prediction autocorrelation",
                lag.max = lag.max)
}
