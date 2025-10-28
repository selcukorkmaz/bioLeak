# Influence-function based cvAUC standard errors --------------------------------
.psi_auc <- function(pos, neg) {
  outer(pos, neg, function(a, b) (a > b) + 0.5 * (a == b))
}

.cvauc_if <- function(pred, truth, weights = NULL) {
  stopifnot(length(pred) == length(truth))
  m <- length(pred)

  fold_auc <- numeric(m)
  fold_var <- numeric(m)
  fold_n   <- integer(m)
  fold_details <- vector("list", m)

  for (i in seq_len(m)) {
    pr <- pred[[i]]
    tr <- truth[[i]]
    if (length(pr) != length(tr)) stop("pred and truth length mismatch per fold")

    y <- if (is.factor(tr)) as.numeric(tr) - 1 else as.numeric(tr)
    pos <- pr[y == 1]
    neg <- pr[y == 0]
    n1 <- length(pos); n0 <- length(neg)
    fold_n[i] <- n1 + n0

    if (n1 == 0L || n0 == 0L) {
      fold_auc[i] <- NA_real_
      fold_var[i] <- NA_real_
      fold_details[[i]] <- data.frame(
        fold = i, auc = NA_real_, se = NA_real_, n_pos = n1, n_neg = n0
      )
      next
    }

    psi <- .psi_auc(pos, neg)
    v_i <- rowMeans(psi)
    w_j <- colMeans(psi)
    auc <- mean(v_i)
    var_pos <- if (n1 > 1L) stats::var(v_i, na.rm = TRUE) else 0
    var_neg <- if (n0 > 1L) stats::var(w_j, na.rm = TRUE) else 0
    var_fold <- var_pos / n1 + var_neg / n0
    fold_auc[i] <- auc
    fold_var[i] <- var_fold
    fold_details[[i]] <- data.frame(
      fold = i, auc = auc, se = sqrt(max(var_fold, 0)),
      n_pos = n1, n_neg = n0
    )
  }

  if (sum(fold_n, na.rm = TRUE) == 0L) {
    return(list(mean_auc = NA_real_, var_auc = NA_real_, se_auc = NA_real_,
                fold_stats = data.frame()))
  }

  if (is.null(weights)) weights <- fold_n / sum(fold_n)
  valid <- is.finite(fold_auc) & is.finite(fold_var)
  if (!any(valid)) {
    return(list(mean_auc = NA_real_, var_auc = NA_real_, se_auc = NA_real_,
                fold_stats = data.frame()))
  }
  weights <- weights[valid] / sum(weights[valid])
  fold_auc <- fold_auc[valid]
  fold_var <- fold_var[valid]

  auc_mean <- sum(weights * fold_auc)
  var_mean <- sum((weights^2) * fold_var)
  fold_stats <- do.call(rbind, fold_details[valid])

  list(
    mean_auc = auc_mean,
    var_auc  = var_mean,
    se_auc   = sqrt(max(var_mean, 0)),
    fold_stats = fold_stats
  )
}
