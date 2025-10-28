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
  graphics::mtext(sprintf("obs = %.3f", obs), side = 3, col = "red")
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
  graphics::barplot(mat, beside = TRUE, legend.text = TRUE,
                    args.legend = list(x = "topright"),
                    main = "Fold class balance")
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
}

#' Plot ACF of test predictions for time-series leakage checks
#' @param fit LeakFit
#' @export
plot_time_acf <- function(fit) {
  stopifnot(inherits(fit, "LeakFit"))
  all_pred <- do.call(rbind, fit@predictions)
  graphics::acf(all_pred$pred, main = "Prediction autocorrelation")
}
