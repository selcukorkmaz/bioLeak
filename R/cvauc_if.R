# Influence-function based cvAUC standard errors --------------------------------

.psi_auc <- function(pos, neg) {
  outer(pos, neg, function(a, b) {
    (a > b) + 0.5 * (a == b)
  })
}

.cvauc_if <- function(pred, truth, weights = NULL) {
  stopifnot(length(pred) == length(truth))
  m <- length(pred)
  fold_auc <- numeric(m)
  fold_var <- numeric(m)
  fold_n <- integer(m)
  fold_details <- vector("list", m)
  for (i in seq_len(m)) {
    pr <- pred[[i]]
    tr <- truth[[i]]
    if (length(pr) != length(tr)) stop("pred and truth length mismatch per fold")
    y <- if (is.factor(tr)) as.numeric(tr) - 1 else as.numeric(tr)
    pos <- pr[y == 1]
    neg <- pr[y == 0]
    n1 <- length(pos); n0 <- length(neg)
    if (n1 == 0L || n0 == 0L) {
      fold_auc[i] <- NA_real_
      fold_var[i] <- NA_real_
      fold_n[i] <- n1 + n0
      fold_details[[i]] <- data.frame(
        fold = i, auc = NA_real_, se = NA_real_, n_pos = n1, n_neg = n0
      )
      next
    }
    psi <- .psi_auc(pos, neg)
    v_i <- rowMeans(psi)
    w_j <- colMeans(psi)
    auc <- mean(v_i)
    var_pos <- if (n1 > 1L) stats::var(v_i) else 0
    var_neg <- if (n0 > 1L) stats::var(w_j) else 0
    var_fold <- var_pos / n1 + var_neg / n0
    fold_auc[i] <- auc
    fold_var[i] <- var_fold
    fold_n[i] <- n1 + n0
    fold_details[[i]] <- data.frame(
      fold = i, auc = auc, se = sqrt(max(var_fold, 0)), n_pos = n1, n_neg = n0
    )
  }
  if (is.null(weights)) {
    weights <- fold_n / sum(fold_n)
  }
  weights <- weights / sum(weights)
  auc_mean <- sum(weights * fold_auc, na.rm = TRUE)
  var_mean <- sum((weights^2) * fold_var, na.rm = TRUE)
  fold_stats <- fold_details[!vapply(fold_details, is.null, logical(1))]
  list(
    mean_auc = auc_mean,
    var_auc = var_mean,
    se_auc = sqrt(max(var_mean, 0)),
    fold_stats = if (length(fold_stats)) do.call(rbind, fold_stats) else data.frame()
  )
}
