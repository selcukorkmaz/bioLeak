# Metric utilities -----------------------------------------------------------

.cindex_pairwise <- function(pred, truth) {
  pred <- as.numeric(pred)
  truth <- as.numeric(truth)
  ok <- is.finite(pred) & is.finite(truth)
  pred <- pred[ok]
  truth <- truth[ok]
  n <- length(pred)
  if (n < 2L) return(NA_real_)
  if (length(unique(truth)) < 2L) return(NA_real_)

  conc <- 0L
  ties <- 0L
  total <- 0L

  for (i in seq_len(n - 1L)) {
    yi <- truth[i]
    pi <- pred[i]
    dy <- yi - truth[(i + 1L):n]
    valid <- dy != 0
    if (!any(valid)) next
    pj <- pred[(i + 1L):n][valid]
    dy <- dy[valid]
    dp <- pi - pj
    total <- total + length(dp)
    prod <- dp * dy
    conc <- conc + sum(prod > 0)
    ties <- ties + sum(dp == 0)
  }

  if (!total) return(NA_real_)
  (conc + 0.5 * ties) / total
}

.cindex_survival <- function(pred, truth) {
  if (!inherits(truth, "Surv")) return(NA_real_)
  if (!requireNamespace("survival", quietly = TRUE)) return(NA_real_)
  df <- data.frame(pred = as.numeric(pred))
  concord <- try(survival::concordance(truth ~ pred, data = df), silent = TRUE)
  if (inherits(concord, "try-error")) return(NA_real_)
  as.numeric(concord$concordance)
}

.multiclass_accuracy <- function(truth, pred_class) {
  if (is.null(pred_class)) return(NA_real_)
  mean(pred_class == truth, na.rm = TRUE)
}

.multiclass_macro_f1 <- function(truth, pred_class) {
  if (is.null(pred_class)) return(NA_real_)
  truth <- factor(truth)
  pred_class <- factor(pred_class, levels = levels(truth))
  lvls <- levels(truth)
  f1_vals <- vapply(lvls, function(lbl) {
    tp <- sum(pred_class == lbl & truth == lbl, na.rm = TRUE)
    fp <- sum(pred_class == lbl & truth != lbl, na.rm = TRUE)
    fn <- sum(pred_class != lbl & truth == lbl, na.rm = TRUE)
    prec <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    rec <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    if (is.na(prec) && is.na(rec)) return(NA_real_)
    if (is.na(prec) || is.na(rec) || (prec + rec) == 0) return(0)
    2 * prec * rec / (prec + rec)
  }, numeric(1))
  if (all(is.na(f1_vals))) return(NA_real_)
  mean(f1_vals, na.rm = TRUE)
}

.multiclass_log_loss <- function(truth, prob, eps = 1e-15) {
  if (is.null(prob)) return(NA_real_)
  truth <- factor(truth)
  levels_truth <- levels(truth)
  if (is.data.frame(prob)) prob <- as.matrix(prob)
  if (ncol(prob) != length(levels_truth)) return(NA_real_)
  prob <- pmin(pmax(prob, eps), 1 - eps)
  idx <- cbind(seq_len(nrow(prob)), match(truth, levels_truth))
  vals <- prob[idx]
  -mean(log(vals), na.rm = TRUE)
}
