# Metric utilities -----------------------------------------------------------

# Resolve the (negative, positive) class levels of a binary outcome. For a
# factor, the second level is the positive class (fit_resample() reorders
# levels so that `positive_class` is second); `positive_class` overrides this.
.binary_levels <- function(truth, positive_class = NULL) {
  if (is.factor(truth)) {
    lev <- levels(truth)
  } else if (is.logical(truth)) {
    lev <- c("FALSE", "TRUE")
  } else {
    u <- sort(unique(truth[!is.na(truth)]))
    lev <- if (all(u %in% c(0, 1))) c("0", "1") else as.character(u)
  }
  if (!is.null(positive_class)) {
    pos <- as.character(positive_class)
    if (pos %in% lev) lev <- c(setdiff(lev, pos), pos)
  }
  lev
}

# ROC AUC with a fixed orientation: higher `pred` means more likely to be the
# positive class. Anti-correlated predictions give AUC < 0.5, so permutation
# nulls centre at 0.5. Equivalent to pROC::roc(direction = "<") with explicit
# levels; never uses pROC's default direction = "auto".
.auc_binary <- function(truth, pred, positive_class = NULL) {
  pred <- as.numeric(pred)
  lev <- .binary_levels(truth, positive_class)
  if (length(lev) != 2L) return(NA_real_)
  truth_chr <- as.character(truth)
  ok <- is.finite(pred) & !is.na(truth_chr) & truth_chr %in% lev
  pred <- pred[ok]
  truth_chr <- truth_chr[ok]
  n_pos <- sum(truth_chr == lev[2])
  n_neg <- sum(truth_chr == lev[1])
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  if (requireNamespace("pROC", quietly = TRUE)) {
    roc <- pROC::roc(response = truth_chr, predictor = pred,
                     levels = lev, direction = "<", quiet = TRUE)
    return(as.numeric(pROC::auc(roc)))
  }
  r <- rank(pred, ties.method = "average")
  (sum(r[truth_chr == lev[2]]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

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
    if (is.na(prec) || is.na(rec)) return(NA_real_)
    if ((prec + rec) == 0) return(0)
    2 * prec * rec / (prec + rec)
  }, numeric(1))
  if (all(is.na(f1_vals))) return(NA_real_)
  f1_vals[is.na(f1_vals)] <- 0
  mean(f1_vals)
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
