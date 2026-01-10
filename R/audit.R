# audit_leakage(): permutation gap, batch association, target scan, duplicates, trail ----

.cosine_sim_block <- function(A, B = NULL) {
  if (is.null(B)) B <- A
  A <- as.matrix(A); B <- as.matrix(B)
  An <- sqrt(rowSums(A * A)); Bn <- sqrt(rowSums(B * B))
  S <- A %*% t(B)
  S / (An %o% Bn + 1e-12)
}

# Row-normalize (L2) to allow cosine distance via Euclidean NN
.row_l2_normalize <- function(M) {
  M <- as.matrix(M)
  nrm <- sqrt(rowSums(M * M))
  nrm[nrm == 0] <- 1
  M / nrm
}

.row_center <- function(M) {
  M <- as.matrix(M)
  mu <- rowMeans(M, na.rm = TRUE)
  sweep(M, 1, mu, "-")
}

# Safe chi^2 with Cramer's V
.chisq_assoc <- function(tab) {
  if (is.null(tab) || length(tab) == 0) {
    return(list(stat = NA_real_, df = NA_integer_, pval = NA_real_, cramer_v = NA_real_))
  }
  if (!is.null(dim(tab))) {
    row_sums <- rowSums(tab)
    col_sums <- colSums(tab)
    if (any(row_sums == 0) || any(col_sums == 0)) {
      tab <- tab[row_sums > 0, col_sums > 0, drop = FALSE]
    }
  }
  if (any(dim(tab) < 2) || sum(tab) == 0) {
    return(list(stat = NA_real_, df = NA_integer_, pval = NA_real_, cramer_v = NA_real_))
  }
  cs <- suppressWarnings(stats::chisq.test(tab))
  n <- sum(tab)
  v <- sqrt(unname(cs$statistic) / (n * (min(dim(tab)) - 1)))
  list(stat = unname(cs$statistic), df = unname(cs$parameter),
       pval = unname(cs$p.value), cramer_v = as.numeric(v))
}

.chisq_assoc_soft <- function(tab) {
  if (any(dim(tab) < 2) || sum(tab) == 0) {
    return(list(stat = NA_real_, df = NA_integer_, pval = NA_real_, cramer_v = NA_real_))
  }
  cs <- suppressWarnings(stats::chisq.test(tab))
  n <- sum(tab)
  v <- sqrt(unname(cs$statistic) / (n * (min(dim(tab)) - 1)))
  list(stat = unname(cs$statistic), df = unname(cs$parameter),
       pval = unname(cs$p.value), cramer_v = as.numeric(v))
}

.auc_rank <- function(x, y01) {
  x <- as.numeric(x)
  y01 <- as.numeric(y01)
  ok <- is.finite(x) & !is.na(y01)
  x <- x[ok]
  y01 <- y01[ok]
  if (!length(x)) return(NA_real_)
  n_pos <- sum(y01 == 1)
  n_neg <- sum(y01 == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(x, ties.method = "average")
  (sum(r[y01 == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

.cor_pval <- function(r, n) {
  if (!is.finite(r) || is.na(r) || n < 3L) return(NA_real_)
  if (abs(r) >= 1) return(0)
  tval <- r * sqrt((n - 2) / (1 - r^2))
  2 * stats::pt(-abs(tval), df = n - 2)
}

.align_by_ids <- function(X, ids_chr, sample_ids = NULL, warn = TRUE) {
  if (is.null(X) || (!is.data.frame(X) && !is.matrix(X))) return(NULL)
  rn <- rownames(X)
  if (!is.null(rn) && all(ids_chr %in% rn)) {
    return(X[ids_chr, , drop = FALSE])
  }
  if (!is.null(sample_ids) && length(sample_ids) == nrow(X)) {
    idx <- match(ids_chr, as.character(sample_ids))
    if (all(!is.na(idx))) return(X[idx, , drop = FALSE])
  }
  ids_int <- suppressWarnings(as.integer(ids_chr))
  if (all(!is.na(ids_int)) && max(ids_int, na.rm = TRUE) <= nrow(X)) {
    return(X[ids_int, , drop = FALSE])
  }
  if (nrow(X) == length(ids_chr)) {
    if (warn) {
      warning("X_ref rownames do not match prediction ids; assuming row order aligns to predictions.",
              call. = FALSE)
    }
    return(X)
  }
  if (warn) warning("X_ref not aligned to predictions; target leakage scan skipped.", call. = FALSE)
  NULL
}

.target_assoc_scan <- function(X, y, task, positive_class = NULL, threshold = 0.9) {
  if (is.null(X) || is.null(y)) return(data.frame())
  if (!is.data.frame(X) && !is.matrix(X)) return(data.frame())
  is_df <- is.data.frame(X)
  if (is_df) {
    n_rows <- nrow(X)
    n_cols <- ncol(X)
    col_names <- colnames(X)
    get_col <- function(j) X[[j]]
  } else {
    X_mat <- as.matrix(X)
    n_rows <- nrow(X_mat)
    n_cols <- ncol(X_mat)
    col_names <- colnames(X_mat)
    get_col <- function(j) X_mat[, j]
  }
  if (!n_rows || !n_cols) return(data.frame())
  if (length(y) != n_rows) {
    warning("X_ref rows do not match outcome length; target leakage scan skipped.", call. = FALSE)
    return(data.frame())
  }

  if (is.null(col_names)) {
    col_names <- paste0("feature_", seq_len(n_cols))
  }

  y01 <- NULL
  y_num <- NULL
  y_fac <- NULL
  if (task == "survival") {
    warning("Target leakage scan is not available for survival outcomes.", call. = FALSE)
    return(data.frame())
  }
  if (task == "binomial") {
    y01 <- y
    if (is.factor(y01)) {
      y01 <- droplevels(y01)
      if (!is.null(positive_class)) {
        pos_chr <- as.character(positive_class)
        if (pos_chr %in% levels(y01) && !identical(levels(y01)[2], pos_chr)) {
          y01 <- factor(y01, levels = c(setdiff(levels(y01), pos_chr), pos_chr))
        }
      }
      y01 <- as.numeric(y01) - 1
    } else if (is.logical(y01)) {
      y01 <- as.numeric(y01)
    } else {
      y01 <- as.numeric(y01)
      uniq <- sort(unique(y01[is.finite(y01)]))
      if (length(uniq) == 2 && !all(uniq %in% c(0, 1))) {
        y01 <- as.numeric(factor(y01, levels = uniq)) - 1
      }
    }
  } else if (task == "multiclass") {
    y_fac <- as.factor(y)
  } else {
    y_num <- as.numeric(y)
  }

  results <- lapply(seq_len(n_cols), function(j) {
    x <- get_col(j)
    if (inherits(x, c("Date", "POSIXct", "POSIXt"))) {
      x <- as.numeric(x)
    }
    is_cat <- is.factor(x) || is.character(x)

    if (is_cat) {
      x_fac <- as.factor(x)
      if (task == "binomial") {
        ok <- !is.na(x_fac) & !is.na(y01)
        x_ok <- droplevels(x_fac[ok])
        y_ok <- y01[ok]
        if (length(x_ok) < 2L || length(unique(x_ok)) < 2L || length(unique(y_ok)) < 2L) {
          return(NULL)
        }
        tab <- table(x_ok, y_ok)
        assoc <- .chisq_assoc_soft(tab)
        data.frame(
          feature = col_names[j],
          type = "categorical",
          metric = "cramer_v",
          value = assoc$cramer_v,
          score = assoc$cramer_v,
          p_value = assoc$pval,
          n = length(y_ok),
          n_levels = nlevels(x_ok),
          stringsAsFactors = FALSE
        )
      } else if (task == "multiclass") {
        ok <- !is.na(x_fac) & !is.na(y_fac)
        x_ok <- droplevels(x_fac[ok])
        y_ok <- droplevels(y_fac[ok])
        if (length(x_ok) < 2L || length(unique(x_ok)) < 2L || length(unique(y_ok)) < 2L) {
          return(NULL)
        }
        tab <- table(x_ok, y_ok)
        assoc <- .chisq_assoc_soft(tab)
        data.frame(
          feature = col_names[j],
          type = "categorical",
          metric = "cramer_v",
          value = assoc$cramer_v,
          score = assoc$cramer_v,
          p_value = assoc$pval,
          n = length(y_ok),
          n_levels = nlevels(x_ok),
          stringsAsFactors = FALSE
        )
      } else {
        ok <- !is.na(x_fac) & is.finite(y_num)
        x_ok <- droplevels(x_fac[ok])
        y_ok <- y_num[ok]
        if (length(y_ok) < 3L || nlevels(x_ok) < 2L) return(NULL)
        fit <- try(stats::lm(y_ok ~ x_ok), silent = TRUE)
        if (inherits(fit, "try-error")) return(NULL)
        an <- stats::anova(fit)
        ss_total <- sum(an$`Sum Sq`, na.rm = TRUE)
        eta_sq <- if (is.finite(ss_total) && ss_total > 0) an$`Sum Sq`[1] / ss_total else NA_real_
        p_val <- an$`Pr(>F)`[1]
        data.frame(
          feature = col_names[j],
          type = "categorical",
          metric = "eta_sq",
          value = eta_sq,
          score = eta_sq,
          p_value = p_val,
          n = length(y_ok),
          n_levels = nlevels(x_ok),
          stringsAsFactors = FALSE
        )
      }
    } else {
      x_num <- as.numeric(x)
      if (task == "binomial") {
        ok <- is.finite(x_num) & !is.na(y01)
        x_ok <- x_num[ok]
        y_ok <- y01[ok]
        if (length(x_ok) < 3L || length(unique(y_ok)) < 2L) return(NULL)
        if (stats::sd(x_ok) == 0) return(NULL)
        auc <- .auc_rank(x_ok, y_ok)
        score <- if (is.na(auc)) NA_real_ else abs(auc - 0.5) * 2
        data.frame(
          feature = col_names[j],
          type = "numeric",
          metric = "auc",
          value = auc,
          score = score,
          p_value = NA_real_,
          n = length(y_ok),
          n_levels = NA_integer_,
          stringsAsFactors = FALSE
        )
      } else if (task == "multiclass") {
        ok <- is.finite(x_num) & !is.na(y_fac)
        x_ok <- x_num[ok]
        y_ok <- droplevels(y_fac[ok])
        if (length(x_ok) < 3L || length(unique(y_ok)) < 2L) return(NULL)
        if (stats::sd(x_ok) == 0) return(NULL)
        fit <- try(stats::lm(x_ok ~ y_ok), silent = TRUE)
        if (inherits(fit, "try-error")) return(NULL)
        an <- stats::anova(fit)
        ss_total <- sum(an$`Sum Sq`, na.rm = TRUE)
        eta_sq <- if (is.finite(ss_total) && ss_total > 0) an$`Sum Sq`[1] / ss_total else NA_real_
        p_val <- an$`Pr(>F)`[1]
        data.frame(
          feature = col_names[j],
          type = "numeric",
          metric = "eta_sq",
          value = eta_sq,
          score = eta_sq,
          p_value = p_val,
          n = length(y_ok),
          n_levels = nlevels(y_ok),
          stringsAsFactors = FALSE
        )
      } else {
        ok <- is.finite(x_num) & is.finite(y_num)
        x_ok <- x_num[ok]
        y_ok <- y_num[ok]
        if (length(x_ok) < 3L) return(NULL)
        if (stats::sd(x_ok) == 0 || stats::sd(y_ok) == 0) return(NULL)
        r <- stats::cor(x_ok, y_ok)
        p_val <- .cor_pval(r, length(x_ok))
        data.frame(
          feature = col_names[j],
          type = "numeric",
          metric = "cor",
          value = r,
          score = abs(r),
          p_value = p_val,
          n = length(y_ok),
          n_levels = NA_integer_,
          stringsAsFactors = FALSE
        )
      }
    }
  })

  results <- results[!vapply(results, is.null, logical(1))]
  if (!length(results)) return(data.frame())
  out <- do.call(rbind, results)
  out$flag <- !is.na(out$score) & out$score >= threshold
  out <- out[order(out$score, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.target_scan_prepare_matrix <- function(X) {
  if (is.null(X) || (!is.data.frame(X) && !is.matrix(X))) return(NULL)
  X_mat <- if (is.data.frame(X)) {
    mm <- try(stats::model.matrix(~ . - 1, data = X), silent = TRUE)
    if (inherits(mm, "try-error")) return(NULL)
    mm
  } else {
    as.matrix(X)
  }
  if (!nrow(X_mat) || !ncol(X_mat)) return(NULL)
  storage.mode(X_mat) <- "numeric"
  for (j in seq_len(ncol(X_mat))) {
    col <- X_mat[, j]
    med <- stats::median(col, na.rm = TRUE)
    if (!is.finite(med)) med <- 0
    col[is.na(col)] <- med
    X_mat[, j] <- col
  }
  sdv <- vapply(seq_len(ncol(X_mat)), function(j) {
    stats::sd(X_mat[, j], na.rm = TRUE)
  }, numeric(1))
  keep <- is.finite(sdv) & sdv > 0
  X_mat <- X_mat[, keep, drop = FALSE]
  if (!ncol(X_mat)) return(NULL)
  X_mat
}

.target_scan_multivariate <- function(X, y, task, splits, folds, coldata, seed, B,
                                      max_components = 10L, interactions = TRUE,
                                      positive_class = NULL, permute_outcome = NULL) {
  if (is.null(X) || is.null(y)) return(data.frame())
  if (!task %in% c("binomial", "gaussian")) return(data.frame())
  X_mat <- .target_scan_prepare_matrix(X)
  if (is.null(X_mat)) return(data.frame())
  if (length(y) != nrow(X_mat)) return(data.frame())

  if (isTRUE(splits@info$compact) && identical(splits@mode, "time_series") && is.null(coldata)) {
    warning("time_series compact splits require coldata for multivariate target scan; skipping.",
            call. = FALSE)
    return(data.frame())
  }

  max_components <- as.integer(max_components)
  if (!is.finite(max_components) || max_components < 1L) return(data.frame())
  n <- nrow(X_mat)
  k <- min(max_components, n - 1L, ncol(X_mat))
  if (!is.finite(k) || k < 1L) return(data.frame())

  X_scaled <- scale(X_mat, center = TRUE, scale = TRUE)
  pc <- try(stats::prcomp(X_scaled, center = FALSE, scale. = FALSE, rank. = k), silent = TRUE)
  if (inherits(pc, "try-error") || is.null(pc$x)) return(data.frame())
  k_eff <- ncol(pc$x)
  if (!is.finite(k_eff) || k_eff < 1L) return(data.frame())
  X_design <- pc$x[, seq_len(k_eff), drop = FALSE]
  colnames(X_design) <- paste0("PC", seq_len(ncol(X_design)))

  n_interactions <- 0L
  if (isTRUE(interactions) && ncol(X_design) >= 2L) {
    k_int <- min(5L, ncol(X_design))
    comb <- utils::combn(k_int, 2)
    int_mat <- matrix(NA_real_, nrow = n, ncol = ncol(comb))
    colnames(int_mat) <- apply(comb, 2, function(idx) {
      paste0("PC", idx[1], "_x_PC", idx[2])
    })
    for (j in seq_len(ncol(comb))) {
      int_mat[, j] <- X_design[, comb[1, j]] * X_design[, comb[2, j]]
    }
    X_design <- cbind(X_design, int_mat)
    n_interactions <- ncol(int_mat)
  }

  score_for_y <- function(y_vec) {
    if (identical(task, "binomial")) {
      y_fac <- if (is.factor(y_vec)) droplevels(y_vec) else factor(y_vec)
      if (nlevels(y_fac) < 2L) {
        return(list(score = NA_real_, value = NA_real_, n_obs = 0L, metric = "auc"))
      }
      pos <- positive_class
      if (is.null(pos) || !as.character(pos) %in% levels(y_fac)) {
        pos <- levels(y_fac)[2]
      }
      y01 <- as.numeric(y_fac == pos)
    } else {
      y_num <- as.numeric(y_vec)
    }

    preds <- rep(NA_real_, n)
    for (i in seq_along(folds)) {
      fold_full <- try(.bio_resolve_fold_indices(splits, folds[[i]], n = n, data = coldata),
                       silent = TRUE)
      if (inherits(fold_full, "try-error")) next
      tr <- fold_full$train
      te <- fold_full$test
      if (!length(tr) || !length(te)) next
      if (identical(task, "binomial")) {
        y_tr <- y01[tr]
        ok_tr <- which(!is.na(y_tr))
        if (length(ok_tr) < 2L) next
        y_tr <- y_tr[ok_tr]
        if (length(unique(y_tr)) < 2L) next
        df_tr <- data.frame(y = y_tr, X_design[tr[ok_tr], , drop = FALSE], check.names = FALSE)
        fit <- try(stats::glm(y ~ ., data = df_tr, family = stats::binomial()), silent = TRUE)
        if (inherits(fit, "try-error")) next
        df_te <- as.data.frame(X_design[te, , drop = FALSE], check.names = FALSE)
        pr <- try(stats::predict(fit, newdata = df_te, type = "response"), silent = TRUE)
        if (inherits(pr, "try-error")) next
        preds[te] <- as.numeric(pr)
      } else {
        y_tr <- y_num[tr]
        ok_tr <- which(is.finite(y_tr))
        if (length(ok_tr) < 2L) next
        df_tr <- data.frame(y = y_tr[ok_tr], X_design[tr[ok_tr], , drop = FALSE], check.names = FALSE)
        fit <- try(stats::lm(y ~ ., data = df_tr), silent = TRUE)
        if (inherits(fit, "try-error")) next
        df_te <- as.data.frame(X_design[te, , drop = FALSE], check.names = FALSE)
        pr <- try(stats::predict(fit, newdata = df_te), silent = TRUE)
        if (inherits(pr, "try-error")) next
        preds[te] <- as.numeric(pr)
      }
    }

    if (identical(task, "binomial")) {
      ok <- is.finite(preds) & !is.na(y01)
      if (sum(ok) < 3L) {
        return(list(score = NA_real_, value = NA_real_, n_obs = sum(ok), metric = "auc"))
      }
      auc <- .auc_rank(preds[ok], y01[ok])
      return(list(score = auc, value = auc, n_obs = sum(ok), metric = "auc"))
    }
    ok <- is.finite(preds) & is.finite(y_num)
    if (sum(ok) < 3L) {
      return(list(score = NA_real_, value = NA_real_, n_obs = sum(ok), metric = "cor"))
    }
    r <- stats::cor(preds[ok], y_num[ok])
    list(score = abs(r), value = r, n_obs = sum(ok), metric = "cor")
  }

  obs <- score_for_y(y)
  perm_scores <- numeric(0)
  if (is.finite(B) && B >= 1L) {
    perm_scores <- vapply(seq_len(B), function(b) {
      set.seed(seed + b)
      y_perm <- if (!is.null(permute_outcome)) permute_outcome(b) else sample(y)
      if (length(y_perm) != length(y)) return(NA_real_)
      score_for_y(y_perm)$score
    }, numeric(1))
  }
  pval <- NA_real_
  if (is.finite(obs$score) && length(perm_scores) && any(is.finite(perm_scores))) {
    finite_perm <- is.finite(perm_scores)
    pval <- (1 + sum(perm_scores[finite_perm] >= obs$score)) / (1 + sum(finite_perm))
  }

  data.frame(
    scan = "multivariate",
    metric = obs$metric,
    value = obs$value,
    score = obs$score,
    p_value = pval,
    n = obs$n_obs,
    n_features = ncol(X_mat),
    n_components = k_eff,
    n_interactions = n_interactions,
    n_perm = length(perm_scores),
    stringsAsFactors = FALSE
  )
}

# Compute metric from predictions
.metric_value <- function(metric, task, truth, pred, pred_df = NULL) {
  if (task == "survival") {
    if (!inherits(truth, "Surv") && !is.null(pred_df)) {
      if (!requireNamespace("survival", quietly = TRUE)) return(NA_real_)
      if (all(c("truth_time", "truth_event") %in% names(pred_df))) {
        truth <- survival::Surv(pred_df$truth_time, pred_df$truth_event)
      } else if (all(c("time", "status") %in% names(pred_df))) {
        truth <- survival::Surv(pred_df$time, pred_df$status)
      }
    }
  }
  if (metric == "rmse") {
    return(sqrt(mean((as.numeric(truth) - pred)^2)))
  }
  if (metric == "accuracy" && task %in% c("binomial", "multiclass")) {
    pred_class <- if (!is.null(pred_df) && "pred_class" %in% names(pred_df)) {
      pred_df$pred_class
    } else if (is.factor(pred) || is.character(pred)) {
      pred
    } else {
      if (task == "binomial") {
        if (is.factor(truth)) {
          factor(ifelse(pred >= 0.5, levels(truth)[2], levels(truth)[1]),
                 levels = levels(truth))
        } else {
          ifelse(pred >= 0.5, 1, 0)
        }
      } else {
        pred
      }
    }
    return(.multiclass_accuracy(truth, pred_class))
  }
  if (metric == "macro_f1" && task == "multiclass") {
    pred_class <- if (!is.null(pred_df) && "pred_class" %in% names(pred_df)) {
      pred_df$pred_class
    } else {
      pred
    }
    return(.multiclass_macro_f1(truth, pred_class))
  }
  if (metric == "log_loss" && task == "multiclass") {
    prob <- NULL
    if (!is.null(pred_df)) {
      prob_cols <- grep("^\\.pred_", names(pred_df), value = TRUE)
      if (length(prob_cols)) {
        prob <- as.matrix(pred_df[, prob_cols, drop = FALSE])
      }
    }
    return(.multiclass_log_loss(truth, prob))
  }
  if (metric == "log_loss" && task == "binomial") {
    yb <- if (is.factor(truth)) as.numeric(truth) - 1 else as.numeric(truth)
    p <- pmin(pmax(as.numeric(pred), 1e-15), 1 - 1e-15)
    return(-mean(yb * log(p) + (1 - yb) * log(1 - p), na.rm = TRUE))
  }
  if (metric == "auc") {
    if (requireNamespace("pROC", quietly = TRUE)) {
      roc <- pROC::roc(truth, pred, quiet = TRUE)
      return(as.numeric(pROC::auc(roc)))
    }
    yb <- if (is.factor(truth)) as.numeric(truth) - 1 else truth
    pos <- pred[yb == 1]; neg <- pred[yb == 0]
    if (length(pos) && length(neg)) {
      comp <- outer(pos, neg, function(a, b) (a > b) + 0.5 * (a == b))
      return(mean(comp))
    }
    return(NA_real_)
  }
  if (metric == "pr_auc") {
    if (requireNamespace("PRROC", quietly = TRUE)) {
      yb <- if (is.factor(truth)) as.numeric(truth) - 1 else truth
      pr <- PRROC::pr.curve(scores.class0 = pred[yb == 1],
                            scores.class1 = pred[yb == 0], curve = FALSE)
      return(pr$auc.integral)
    }
    return(NA_real_)
  }
  if (metric == "cindex") {
    if (task == "survival") {
      return(.cindex_survival(pred, truth))
    }
    return(.cindex_pairwise(pred, truth))
  }
  NA_real_
}

#' Audit leakage and confounding
#'
#' @description
#' Computes a post-hoc leakage audit for a resampled model fit. The audit
#' (1) compares observed cross-validated performance to a label-permutation
#' null (by default refitting when data are available; otherwise using fixed
#' predictions), (2) tests whether fold assignments
#' are associated with batch or study metadata (confounding by design),
#' (3) scans features for unusually strong outcome proxies, and (4) flags
#' duplicate or near-duplicate samples in a reference feature matrix.
#'
#' The returned [LeakAudit] summarizes these diagnostics. It relies on the
#' stored predictions, splits, and optional metadata; it does not refit models
#' unless `perm_refit = TRUE` (or `perm_refit = "auto"` with a valid
#' `perm_refit_spec`). Results are conditional on the chosen metric
#' and supplied metadata/features and should be interpreted as diagnostics,
#' not proof of leakage or its absence.
#'
#' @param fit A [LeakFit] object from [fit_resample()] containing cross-validated
#'   predictions and split metadata. If predictions include learner IDs for
#'   multiple models, you must supply `learner` to select one; if learner IDs are
#'   absent, the audit uses all predictions and may mix learners.
#' @param metric Character scalar. One of `"auc"`, `"pr_auc"`, `"accuracy"`,
#'   `"macro_f1"`, `"log_loss"`, `"rmse"`, or `"cindex"`. Defaults to `"auc"`.
#'   This controls the observed performance
#'   statistic, the permutation null, and the sign of the reported gap.
#' @param B Integer scalar. Number of permutations used to build the null
#'   distribution (default 200). Larger values reduce Monte Carlo error but
#'   increase runtime.
#' @param perm_stratify Logical scalar or `"auto"`. If TRUE (default), permutations
#'   are stratified within each fold (factor levels; numeric outcomes are binned
#'   into quantiles when enough non-missing values are available). If FALSE, no
#'   stratification is used. Stratification only applies when `coldata` supplies
#'   the outcome; otherwise labels are shuffled within each fold.
#' @param perm_refit Logical scalar or `"auto"`. If FALSE, permutations keep
#'   predictions fixed and shuffle labels (association test). If TRUE, each
#'   permutation refits the model on permuted outcomes using `perm_refit_spec`.
#'   Refit-based permutations are slower but better approximate a full null
#'   distribution. The default is `"auto"`, which refits only when
#'   `perm_refit_spec` is provided and `B` is less than or equal to
#'   `perm_refit_auto_max`; otherwise it falls back to fixed-prediction
#'   permutations.
#' @param perm_refit_auto_max Integer scalar. Maximum `B` allowed for
#'   `perm_refit = "auto"` to trigger refitting. Defaults to 200.
#' @param perm_refit_spec List of inputs used when `perm_refit = TRUE`.
#'   Required elements: `x` (data used for fitting) and `learner` (parsnip
#'   model_spec, workflow, or legacy learner). Optional elements: `outcome`
#'   (defaults to `fit@outcome`), `preprocess`, `learner_args`,
#'   `custom_learners`, `class_weights`, `positive_class`, and `parallel`.
#'   Survival outcomes are not supported for refit-based permutations.
#' @param perm_mode Optional character scalar to override the permutation mode
#'   used for restricted shuffles. One of `"subject_grouped"`, `"batch_blocked"`,
#'   `"study_loocv"`, or `"time_series"`. Defaults to the split metadata when
#'   available (including rsample-derived modes).
#' @param time_block Character scalar, `"circular"` or `"stationary"`. Controls
#'   block permutation for `time_series` splits; ignored for other split modes.
#'   Default is `"circular"`.
#' @param block_len Integer scalar or NULL. Block length for time-series
#'   permutations. NULL selects `max(5, floor(0.1 * fold_size))`. Larger values
#'   preserve more temporal structure and yield a more conservative null.
#' @param include_z Logical scalar. If TRUE (default), include the z-score for the
#'   permutation gap when a standard error is available; if FALSE, `z` is NA.
#' @param ci_method Character scalar, `"if"` or `"bootstrap"`. Controls how the
#'   standard error and confidence interval for the permutation gap are estimated.
#'   Default is `"if"`. `"if"` uses an influence-function estimate when available; `"bootstrap"`
#'   resamples permutation values `boot_B` times. Failed estimates yield NA.
#' @param boot_B Integer scalar. Number of bootstrap resamples when
#'   `ci_method = "bootstrap"` (default 400). Larger values are more stable but slower.
#' @param parallel Logical scalar. If TRUE and `future.apply` is available,
#'   permutations run in parallel. Results should match sequential execution.
#'   Default is FALSE.
#' @param seed Integer scalar. Random seed used for permutations and bootstrap
#'   resampling; changing it changes the randomization but not the observed metric.
#'   Default is 1.
#' @param return_perm Logical scalar. If TRUE (default), stores the permutation
#'   distribution in `audit@perm_values`. Set FALSE to reduce memory use.
#' @param batch_cols Character vector. Names of `coldata` columns to test for
#'   association with fold assignment. If NULL, defaults to any of
#'   `"batch"`, `"plate"`, `"center"`, `"site"`, `"study"` found in `coldata`.
#'   Changing this controls which batch tests appear in `batch_assoc`.
#' @param coldata Optional data.frame of sample-level metadata. Rows must align
#'   to prediction ids via row names, a `row_id` column, or row order. Used to
#'   build restricted permutations (when the outcome column is present), compute
#'   batch associations, and supply outcomes for target scans. If NULL, uses
#'   `fit@splits@info$coldata` when available. If alignment fails, restricted
#'   permutations are disabled with a warning.
#' @param X_ref Optional numeric matrix/data.frame (samples x features). Used for
#'   duplicate detection and the target leakage scan. If NULL, uses
#'   `fit@info$X_ref` when available. Rows must align to sample ids (split order)
#'   via row names, a `row_id` column, or row order; misalignment disables these
#'   checks.
#' @param target_scan Logical scalar. If TRUE (default), computes per-feature
#'   outcome associations on `X_ref` and flags proxy features; if FALSE, or if
#'   `X_ref`/outcomes are unavailable, `target_assoc` is empty. Not available
#'   for survival outcomes.
#' @param target_scan_multivariate Logical scalar. If TRUE (default), fits a simple
#'   multivariate/interaction model on `X_ref` using the stored splits and
#'   reports a permutation-based score/p-value. This is slower and only
#'   implemented for binomial and gaussian tasks.
#' @param target_scan_multivariate_B Integer scalar. Number of permutations for
#'   the multivariate scan (default 100). Larger values stabilize the p-value.
#' @param target_scan_multivariate_components Integer scalar. Maximum number of
#'   principal components used in the multivariate scan (default 10).
#' @param target_scan_multivariate_interactions Logical scalar. If TRUE (default),
#'   adds pairwise interactions among the top components in the multivariate scan.
#' @param target_threshold Numeric scalar in (0,1). Threshold applied to the
#'   association score used to flag proxy features. Higher values are stricter.
#'   Default is 0.9.
#' @param feature_space Character scalar, `"raw"` or `"rank"`. If `"rank"`,
#'   each row of `X_ref` is rank-transformed before similarity calculations.
#'   This affects duplicate detection only. Default is `"raw"`.
#' @param sim_method Character scalar, `"cosine"` or `"pearson"`. Similarity
#'   metric for duplicate detection. `"pearson"` row-centers before cosine.
#'   Default is `"cosine"`.
#' @param sim_threshold Numeric scalar in (0,1). Similarity cutoff for reporting
#'   duplicate pairs (default 0.995). Higher values yield fewer pairs.
#' @param nn_k Integer scalar. For large datasets (`n > 3000`) with `RANN`
#'   installed, checks only the nearest `nn_k` neighbors per row. Larger values
#'   increase sensitivity but slow the search. Ignored when full comparisons are used.
#'   Default is 50.
#' @param max_pairs Integer scalar. Maximum number of duplicate pairs returned.
#'   If more pairs are found, only the most similar are kept. This does not
#'   affect permutation results. Default is 5000.
#' @param duplicate_scope Character scalar. One of `"train_test"` (default) or
#'   `"all"`. `"train_test"` retains only near-duplicate pairs that appear in
#'   train vs test in at least one repeat; `"all"` reports all near-duplicate
#'   pairs in `X_ref` regardless of fold assignment.
#' @param learner Optional character scalar. When predictions include multiple
#'   learner IDs, selects the learner to audit. If NULL and multiple learners
#'   are present, the function errors; if predictions lack learner IDs, this
#'   argument is ignored with a warning. Default is NULL.
#' @return A [LeakAudit] object with slots:
#'   `permutation_gap` (one-row data.frame), `perm_values` (optional numeric
#'   vector), `batch_assoc`, `target_assoc`, `duplicates`, and `trail`.
#' @details
#' The `permutation_gap` slot reports `metric_obs`, `perm_mean`, `perm_sd`,
#' `gap`, `z`, `p_value`, and `n_perm`. The gap is defined as
#' `metric_obs - perm_mean` for metrics where higher is better (AUC, PR-AUC,
#' accuracy, macro-F1, C-index) and `perm_mean - metric_obs` for RMSE/log-loss.
#
#' By default, `perm_refit = "auto"` refits models when refit data are available
#' and `B` is not too large; otherwise it keeps predictions fixed and shuffles
#' labels. Fixed-prediction permutations quantify prediction-label association
#' rather than a full refit null. Set `perm_refit = FALSE` to force fixed
#' predictions, or `perm_refit = TRUE` (with `perm_refit_spec`) to always refit.
#'
#' `batch_assoc` contains chi-square tests between fold assignment and each
#' `batch_cols` variable (`stat`, `df`, `pval`, `cramer_v`). `target_assoc`
#' reports feature-wise outcome associations on `X_ref`; numeric features use
#' AUC (binomial), `eta_sq` (multiclass), or correlation (gaussian), while
#' categorical features use Cramer's V (binomial/multiclass) or `eta_sq` from a
#' one-way ANOVA (gaussian). The
#' `score` column is the scaled effect size used for flagging
#' (`flag = score >= target_threshold`).
#' The univariate target leakage scan can miss multivariate proxies, interaction
#' leakage, or features not included in `X_ref`. The multivariate scan (enabled
#' by default for supported tasks) adds a model-based proxy check but still only
#' covers features present in `X_ref`.
#'
#' Duplicate detection compares rows of `X_ref` using the chosen `sim_method`
#' (cosine on L2-normalized rows, or Pearson via row-centering), optionally after
#' rank transformation (`feature_space = "rank"`). By default,
#' `duplicate_scope = "train_test"` filters to pairs that appear in train vs test
#' in at least one repeat; set `duplicate_scope = "all"` to include within-fold
#' duplicates. The `duplicates` slot returns index pairs and similarity values
#' for near-duplicate samples. Only duplicates present in `X_ref` can be
#' detected, and checks are skipped if inputs cannot be aligned to splits.
#' @examples
#' \dontrun{
#' set.seed(1)
#' df <- data.frame(
#'   subject = rep(1:6, each = 2),
#'   outcome = rbinom(12, 1, 0.5),
#'   x1 = rnorm(12),
#'   x2 = rnorm(12)
#' )
#'
#' splits <- make_split_plan(df, outcome = "outcome",
#'                       mode = "subject_grouped", group = "subject", v = 3,
#'                       progress = FALSE)
#'
#' custom <- list(
#'   glm = list(
#'     fit = function(x, y, task, weights, ...) {
#'       stats::glm(y ~ ., data = as.data.frame(x),
#'                  family = stats::binomial(), weights = weights)
#'     },
#'     predict = function(object, newdata, task, ...) {
#'       as.numeric(stats::predict(object,
#'                                 newdata = as.data.frame(newdata),
#'                                 type = "response"))
#'     }
#'   )
#' )
#'
#' fit <- fit_resample(df, outcome = "outcome", splits = splits,
#'                     learner = "glm", custom_learners = custom,
#'                     metrics = "auc", refit = FALSE, seed = 1)
#'
#' audit <- audit_leakage(fit, metric = "auc", B = 10,
#'                        X_ref = df[, c("x1", "x2")])
#'
#' }
#' @export
audit_leakage <- function(fit,
                          metric = c("auc", "pr_auc", "accuracy", "macro_f1", "log_loss", "rmse", "cindex"),
                          B = 200,
                          perm_stratify = FALSE,
                          perm_refit = "auto",
                          perm_refit_auto_max = 200,
                          perm_refit_spec = NULL,
                          perm_mode = NULL,
                          time_block = c("circular", "stationary"),
                          block_len = NULL,
                          include_z = TRUE,
                          ci_method = c("if", "bootstrap"),
                          boot_B = 400,
                          parallel = FALSE,
                          seed = 1,
                          return_perm = TRUE,
                          batch_cols = NULL,
                          coldata = NULL,
                          X_ref = NULL,
                          target_scan = TRUE,
                          target_scan_multivariate = TRUE,
                          target_scan_multivariate_B = 100,
                          target_scan_multivariate_components = 10,
                          target_scan_multivariate_interactions = TRUE,
                          target_threshold = 0.9,
                          feature_space = c("raw", "rank"),
                          sim_method = c("cosine", "pearson"),
                          sim_threshold = 0.995,
                          nn_k = 50,
                          max_pairs = 5000,
                          duplicate_scope = c("train_test", "all"),
                          learner = NULL) {

  # --- CRITICAL PATCH: Support LeakTune objects (via tune_resample) ---
  if (inherits(fit, "LeakTune")) {
    # If using 'tune_resample', we construct a proxy LeakFit object
    # by aggregating predictions from the outer resampling folds.

    # 1. Collect all predictions
    all_preds_list <- list()
    for (i in seq_along(fit$outer_fits)) {
      of <- fit$outer_fits[[i]]
      # LeakFit predictions are lists of dataframes
      if (length(of@predictions) > 0) {
        # Combine inner predictions (usually just one DF for outer fit)
        p_df <- do.call(rbind, lapply(of@predictions, as.data.frame))

        # CRITICAL FIX 1: Overwrite 'fold' index to match the outer loop index 'i'.
        # Individual outer_fits are trained as single-fold objects (fold=1),
        # so we must re-index them to match the original splits (1..5).
        if (nrow(p_df) > 0) {
          p_df$fold <- i
        }

        all_preds_list[[i]] <- p_df
      }
    }

    if (length(all_preds_list) == 0) {
      stop("LeakTune object contains no outer loop predictions to audit.")
    }

    # 2. Reconstruct a minimal LeakFit for auditing
    first_fit <- fit$outer_fits[[1]]

    # CRITICAL FIX 2: Safely extract feature names to satisfy S4 validation
    # (must be character, not NULL)
    fnames <- tryCatch(first_fit@feature_names, error = function(e) character(0))
    if (is.null(fnames)) fnames <- character(0)

    # Create the proxy S4 object using new() to bypass validation if needed
    fit <- new("LeakFit",
               splits = fit$splits,
               metrics = fit$metrics,
               metric_summary = fit$metric_summary,
               audit = data.frame(),
               predictions = all_preds_list,
               preprocess = list(),
               learners = list(),
               outcome = fit$info$outcome %||% first_fit@outcome,
               task = fit$info$task %||% first_fit@task,
               feature_names = fnames,
               info = list(
                 hash = fit$splits@info$hash,
                 metrics_requested = fit$info$metrics_requested,
                 metrics_used = fit$info$metrics_used,
                 class_weights = NULL,
                 positive_class = fit$info$positive_class,
                 # CRITICAL FIX 3: Use sample_ids from the first outer fit to enable valid permutation
                 sample_ids = first_fit@info$sample_ids,
                 fold_seeds = NULL,
                 refit = FALSE,
                 final_model = NULL,
                 final_preprocess = NULL,
                 learner_names = unique(fit$metrics$learner),
                 perm_refit_spec = NULL
               )
    )
  }
  # --- END PATCH ---

  metric <- match.arg(metric)
  feature_space <- match.arg(feature_space)
  sim_method    <- match.arg(sim_method)
  duplicate_scope <- match.arg(duplicate_scope)
  time_block <- match.arg(time_block)
  ci_method <- match.arg(ci_method)

  if (!is.numeric(target_threshold) || length(target_threshold) != 1L ||
      !is.finite(target_threshold) || target_threshold <= 0 || target_threshold >= 1) {
    stop("target_threshold must be a single numeric value in (0,1).")
  }
  if (!is.logical(target_scan_multivariate) || length(target_scan_multivariate) != 1L) {
    stop("target_scan_multivariate must be TRUE or FALSE.")
  }
  target_scan_multivariate <- isTRUE(target_scan_multivariate)
  target_scan_multivariate_B <- as.integer(target_scan_multivariate_B)
  if (isTRUE(target_scan_multivariate)) {
    if (!is.finite(target_scan_multivariate_B) || target_scan_multivariate_B < 1L) {
      stop("target_scan_multivariate_B must be an integer >= 1.")
    }
    target_scan_multivariate_components <- as.integer(target_scan_multivariate_components)
    if (!is.finite(target_scan_multivariate_components) ||
        target_scan_multivariate_components < 1L) {
      stop("target_scan_multivariate_components must be an integer >= 1.")
    }
    if (!is.logical(target_scan_multivariate_interactions) ||
        length(target_scan_multivariate_interactions) != 1L) {
      stop("target_scan_multivariate_interactions must be TRUE or FALSE.")
    }
  }

  set.seed(seed)
  if (is.null(perm_refit_spec) && !is.null(fit@info$perm_refit_spec)) {
    perm_refit_spec <- fit@info$perm_refit_spec
  }


  select_refit_learner <- function(spec_learner, learner, learner_names = NULL) {
    if (is.null(spec_learner) || is.null(learner)) return(spec_learner)
    if (inherits(spec_learner, "workflow") || inherits(spec_learner, "model_spec")) {
      return(spec_learner)
    }
    if (is.list(spec_learner)) {
      nm <- names(spec_learner)
      if (!is.null(nm) && learner %in% nm) {
        return(spec_learner[[learner]])
      }
      if (!is.null(learner_names)) {
        idx <- match(learner, learner_names)
        if (is.finite(idx) && idx >= 1L && idx <= length(spec_learner)) {
          return(spec_learner[[idx]])
        }
      }
    }
    if (is.character(spec_learner) && length(spec_learner) > 1L) {
      nm <- names(spec_learner)
      if (!is.null(nm) && learner %in% nm) {
        return(spec_learner[[learner]])
      }
      if (!is.null(learner_names)) {
        idx <- match(learner, learner_names)
        if (is.finite(idx) && idx >= 1L && idx <= length(spec_learner)) {
          return(spec_learner[[idx]])
        }
      }
    }
    spec_learner
  }

  if (!is.null(perm_refit_spec) && is.list(perm_refit_spec) && !is.null(learner)) {
    perm_refit_spec$learner <- select_refit_learner(
      perm_refit_spec$learner %||% NULL,
      learner,
      fit@info$learner_names %||% NULL
    )
  }
  perm_refit_mode <- "fixed"
  perm_refit_reason <- NULL
  perm_refit_raw <- perm_refit
  if (is.character(perm_refit_raw)) {
    if (!identical(perm_refit_raw, "auto")) {
      stop("perm_refit must be TRUE, FALSE, or \"auto\".", call. = FALSE)
    }
    perm_refit_mode <- "auto-fixed"
    if (!is.numeric(perm_refit_auto_max) || length(perm_refit_auto_max) != 1L ||
        !is.finite(perm_refit_auto_max) || perm_refit_auto_max < 1L) {
      stop("perm_refit_auto_max must be a single integer >= 1.", call. = FALSE)
    }
    if (is.null(perm_refit_spec) || !is.list(perm_refit_spec)) {
      perm_refit_reason <- "perm_refit_spec missing; using fixed predictions."
      perm_refit <- FALSE
    } else if (!is.finite(B) || B > perm_refit_auto_max) {
      perm_refit_reason <- sprintf("B=%s exceeds perm_refit_auto_max=%s; using fixed predictions.",
                                   as.character(B), as.character(perm_refit_auto_max))
      perm_refit <- FALSE
    } else {
      perm_refit <- TRUE
      perm_refit_mode <- "auto-refit"
    }
  } else if (isTRUE(perm_refit_raw)) {
    perm_refit <- TRUE
    perm_refit_mode <- "refit"
  } else if (isFALSE(perm_refit_raw)) {
    perm_refit <- FALSE
    perm_refit_mode <- "fixed"
  } else {
    stop("perm_refit must be TRUE, FALSE, or \"auto\".", call. = FALSE)
  }
  perm_method <- if (isTRUE(perm_refit)) "refit" else "fixed"

  refit_x <- NULL
  refit_outcome <- NULL
  refit_learner <- NULL
  refit_preprocess <- NULL
  refit_learner_args <- list()
  refit_custom_learners <- list()
  refit_class_weights <- NULL
  refit_positive_class <- NULL
  refit_parallel <- FALSE
  refit_coldata <- NULL
  refit_coldata_supplied <- FALSE
  default_preprocess <- list(
    impute = list(method = "median"),
    normalize = list(method = "zscore"),
    filter = list(var_thresh = 0, iqr_thresh = 0),
    fs = list(method = "none")
  )
  if (perm_refit) {
    if (is.null(perm_refit_spec) || !is.list(perm_refit_spec)) {
      stop("perm_refit=TRUE requires perm_refit_spec as a list with x and learner.")
    }
    refit_x <- perm_refit_spec$x %||% NULL
    if (is.null(refit_x)) {
      stop("perm_refit_spec$x is required when perm_refit=TRUE.")
    }
    refit_outcome <- perm_refit_spec$outcome %||% fit@outcome
    if (is.null(refit_outcome) || !length(refit_outcome)) {
      stop("perm_refit_spec$outcome (or fit@outcome) is required when perm_refit=TRUE.")
    }
    refit_learner <- perm_refit_spec$learner %||% NULL
    if (is.null(refit_learner)) {
      stop("perm_refit_spec$learner is required when perm_refit=TRUE.")
    }
    if (!inherits(refit_learner, c("workflow", "model_spec")) &&
        ((is.list(refit_learner) && length(refit_learner) != 1L) ||
         (is.character(refit_learner) && length(refit_learner) != 1L))) {
      stop("perm_refit requires a single learner; supply learner= or a single perm_refit_spec$learner.",
           call. = FALSE)
    }
    refit_preprocess <- perm_refit_spec$preprocess %||% default_preprocess
    refit_learner_args <- perm_refit_spec$learner_args %||% list()
    refit_custom_learners <- perm_refit_spec$custom_learners %||% list()
    refit_class_weights <- perm_refit_spec$class_weights %||% fit@info$class_weights %||% NULL
    refit_positive_class <- perm_refit_spec$positive_class %||% fit@info$positive_class %||% NULL
    refit_parallel <- isTRUE(perm_refit_spec$parallel)
    refit_coldata_supplied <- !is.null(perm_refit_spec$coldata)
    refit_coldata <- perm_refit_spec$coldata %||% NULL
  }

  # Trail / provenance
  trail <- list(
    indices_hash = .bio_hash_indices(fit@splits@indices),
    mode = fit@splits@mode,
    info = fit@splits@info,
    seed = seed
  )
  trail$metric <- metric
  trail$perm_method <- perm_method
  if (!is.null(fit@info$learner)) trail$learner <- fit@info$learner

  # --- Reconstruct main metric from CV predictions --------------------------
  task <- fit@task
  valid_metrics <- switch(task,
                          binomial = c("auc", "pr_auc", "accuracy", "log_loss"),
                          multiclass = c("accuracy", "macro_f1", "log_loss"),
                          gaussian = c("rmse", "cindex"),
                          survival = c("cindex"),
                          c("auc", "pr_auc", "rmse", "cindex"))
  if (!metric %in% valid_metrics) {
    stop(sprintf("Metric '%s' is not supported for %s tasks.", metric, task))
  }
  folds <- fit@splits@indices
  compact <- isTRUE(fit@splits@info$compact)
  fold_assignments <- fit@splits@info$fold_assignments
  perm_mode_use <- .bio_perm_mode(fit@splits)
  if (!is.null(perm_mode)) {
    valid_modes <- c("subject_grouped", "batch_blocked", "study_loocv", "time_series")
    perm_mode <- as.character(perm_mode)
    perm_mode <- perm_mode[!is.na(perm_mode) & nzchar(perm_mode)]
    if (length(perm_mode) != 1L || !perm_mode %in% valid_modes) {
      stop("perm_mode must be one of: subject_grouped, batch_blocked, study_loocv, time_series.",
           call. = FALSE)
    }
    perm_mode_use <- perm_mode
  }
  if (identical(fit@splits@mode, "rsample") && identical(perm_mode_use, "rsample")) {
    stop(paste0("rsample splits require an explicit perm_mode; pass perm_mode= to ",
                "audit_leakage() or set split_cols/attr(splits, 'bioLeak_perm_mode') when fitting."),
         call. = FALSE)
  }
  trail$perm_mode <- perm_mode_use
  resolve_test_idx <- function(fold) {
    if (!isTRUE(compact) && !is.null(fold$test)) return(fold$test)
    if (is.null(fold_assignments) || !length(fold_assignments)) {
      stop("Compact splits require fold assignments to resolve test indices.")
    }
    r <- fold$repeat_id
    if (is.null(r) || !is.finite(r)) r <- 1L
    assign_vec <- fold_assignments[[r]]
    if (is.null(assign_vec)) {
      stop(sprintf("Missing fold assignments for repeat %s.", r))
    }
    which(assign_vec == fold$fold)
  }
  if (is.null(coldata) && !is.null(fit@splits@info$coldata)) {
    coldata <- fit@splits@info$coldata
  }
  outcome_col <- fit@splits@info$outcome
  if (is.null(outcome_col)) outcome_col <- fit@outcome
  if (length(outcome_col) != 1L || identical(task, "survival")) outcome_col <- NULL
  if (isTRUE(target_scan) && identical(task, "survival")) {
    warning("Target leakage scan is not available for survival outcomes.", call. = FALSE)
    target_scan <- FALSE
  }
  if (isTRUE(target_scan_multivariate) && !task %in% c("binomial", "gaussian")) {
    warning("Multivariate target scan is only available for binomial or gaussian tasks.",
            call. = FALSE)
    target_scan_multivariate <- FALSE
  }

  pred_list_raw <- lapply(fit@predictions, function(df) data.frame(df, stringsAsFactors = FALSE))
  pred_df <- if (length(pred_list_raw)) do.call(rbind, pred_list_raw) else NULL
  if (is.null(pred_df) || !nrow(pred_df)) {
    stop("No predictions available in LeakFit object.")
  }

  has_learner <- "learner" %in% names(pred_df)
  if (has_learner) {
    pred_df$learner <- as.character(pred_df$learner)
    learner_vals <- unique(pred_df$learner)
    if (is.null(learner)) {
      if (length(learner_vals) == 1L) {
        learner <- learner_vals[[1]]
      } else {
        stop("Multiple learners found in predictions; specify `learner` to audit a single model.")
      }
    } else {
      if (length(learner) != 1L) stop("learner must be a single value.")
      if (!learner %in% learner_vals) {
        stop(sprintf("Learner '%s' not found in predictions. Available: %s",
                     learner, paste(learner_vals, collapse = ", ")))
      }
    }
    pred_df <- pred_df[pred_df$learner == learner, , drop = FALSE]
    if (!nrow(pred_df)) {
      stop(sprintf("No predictions available for learner '%s'.", learner))
    }
  } else {
    if (!is.null(learner)) {
      warning("`learner` ignored: predictions do not include learner IDs.")
    } else if (!is.null(fit@metrics) && length(unique(fit@metrics$learner)) > 1L) {
      warning("Multiple learners were fit but predictions lack learner IDs; audit may mix learners. Refit with updated bioLeak.")
    }
  }
  if (!is.null(learner)) {
    trail$learner <- learner
  }

  resolve_sample_ids <- function(fit_obj, fallback_n = NULL) {
    ids <- fit_obj@info$sample_ids %||% NULL
    if (!is.null(ids) && length(ids)) return(as.character(ids))
    cd_lookup <- fit_obj@splits@info$coldata %||% NULL
    if (!is.null(cd_lookup)) {
      rn_cd <- rownames(cd_lookup)
      if (!is.null(rn_cd) && !anyNA(rn_cd) && all(nzchar(rn_cd)) && !anyDuplicated(rn_cd)) {
        return(as.character(rn_cd))
      }
      if ("row_id" %in% names(cd_lookup)) {
        rid <- as.character(cd_lookup[["row_id"]])
        if (length(rid) == nrow(cd_lookup) && !anyNA(rid) &&
            !anyDuplicated(rid) && all(nzchar(rid))) {
          return(rid)
        }
      }
    }
    if (!is.null(fallback_n) && is.finite(fallback_n)) {
      return(as.character(seq_len(fallback_n)))
    }
    NULL
  }

  align_coldata_for_perm <- function(cd, sample_ids, context) {
    if (is.null(cd)) return(NULL)
    if (is.null(sample_ids) || !length(sample_ids)) {
      warning(sprintf("%s coldata alignment failed (missing sample ids); restricted permutations disabled.", context),
              call. = FALSE)
      return(NULL)
    }
    cd <- as.data.frame(cd, check.names = FALSE)
    sample_ids <- as.character(sample_ids)
    rn <- rownames(cd)
    if (!is.null(rn) && length(rn)) {
      if (anyDuplicated(rn)) {
        warning(sprintf("%s coldata rownames are duplicated; restricted permutations disabled.", context),
                call. = FALSE)
        return(NULL)
      }
      if (all(sample_ids %in% rn)) {
        return(cd[match(sample_ids, rn), , drop = FALSE])
      }
      warning(sprintf("%s coldata rownames do not match sample ids; restricted permutations disabled.",
                      context), call. = FALSE)
      return(NULL)
    }
    if ("row_id" %in% names(cd)) {
      rid <- as.character(cd[["row_id"]])
      if (anyDuplicated(rid)) {
        warning(sprintf("%s coldata row_id values are duplicated; restricted permutations disabled.",
                        context), call. = FALSE)
        return(NULL)
      }
      if (all(sample_ids %in% rid)) {
        return(cd[match(sample_ids, rid), , drop = FALSE])
      }
      warning(sprintf("%s coldata row_id values do not match sample ids; restricted permutations disabled.",
                      context), call. = FALSE)
      return(NULL)
    }
    if (nrow(cd) == length(sample_ids)) {
      warning(sprintf("%s coldata has no ids; assuming row order aligns to splits for permutations.",
                      context), call. = FALSE)
      return(cd)
    }
    warning(sprintf("%s coldata not aligned to splits; restricted permutations disabled.", context),
            call. = FALSE)
    NULL
  }

  if ("fold" %in% names(pred_df)) {
    pred_list <- lapply(seq_along(folds), function(i) {
      pred_df[pred_df$fold == i, , drop = FALSE]
    })
  } else if (has_learner) {
    pred_list <- lapply(pred_list_raw, function(df) {
      if ("learner" %in% names(df)) df[df$learner == learner, , drop = FALSE] else df
    })
  } else {
    pred_list <- pred_list_raw
  }

  all_pred <- do.call(rbind, pred_list)

  metric_obs <- .metric_value(metric, task, all_pred$truth, all_pred$pred, pred_df = all_pred)

  delta <- NA_real_
  perm_mean <- NA_real_
  perm_sd <- NA_real_
  pval <- NA_real_
  p_se <- NA_real_

  if (is.na(metric_obs) || !is.finite(metric_obs)) {
    warning("Observed metric is NA or non-finite; skipping permutation gap calculation.")
    metric_obs <- NA_real_
  }

  # Weights per fold for IF aggregation
  weights <- vapply(folds, function(f) length(resolve_test_idx(f)), integer(1))
  weights <- weights / sum(weights)

  se_obs <- NA_real_
  if (metric == "auc" && exists(".cvauc_if", mode = "function")) {
    obs_truth <- lapply(pred_list, `[[`, "truth")
    obs_pred <- lapply(pred_list, `[[`, "pred")
    if_stat <- try(.cvauc_if(pred = obs_pred, truth = obs_truth, weights = weights), silent = TRUE)
    if (!inherits(if_stat, "try-error") && !is.null(if_stat$se_auc)) {
      se_obs <- if_stat$se_auc
    }
  }

  higher_better <- metric %in% c("auc", "pr_auc", "cindex", "accuracy", "macro_f1")

  # --- Permutations ----------------------------------------------------------
  perm_vals <- numeric(0)
  if (perm_refit) {
    if (identical(task, "survival")) {
      stop("perm_refit=TRUE is not supported for survival tasks.", call. = FALSE)
    }
    if (length(refit_outcome) != 1L) {
      stop("perm_refit=TRUE requires a single outcome column name.", call. = FALSE)
    }

    refit_x_mat <- .bio_get_x(refit_x)
    n_refit <- nrow(refit_x_mat)
    if (!is.finite(n_refit) || n_refit < 1L) {
      stop("perm_refit_spec$x has no rows.", call. = FALSE)
    }
    if (isTRUE(compact) && length(fold_assignments)) {
      expected_n <- max(vapply(fold_assignments, length, integer(1)), na.rm = TRUE)
      if (is.finite(expected_n) && n_refit != expected_n) {
        stop("perm_refit_spec$x row count does not match compact split assignments.", call. = FALSE)
      }
    } else {
      max_idx <- suppressWarnings(max(vapply(folds, function(z) {
        idx <- c(z$train, z$test)
        if (length(idx)) max(idx) else 0L
      }, integer(1)), na.rm = TRUE))
      if (is.finite(max_idx) && n_refit < max_idx) {
        stop("perm_refit_spec$x has fewer rows than required by splits.", call. = FALSE)
      }
    }

    if (is.null(refit_coldata)) {
      if (.bio_is_se(refit_x)) {
        refit_coldata <- as.data.frame(SummarizedExperiment::colData(refit_x))
      } else if (is.data.frame(refit_x)) {
        refit_coldata <- refit_x
      }
    } else if (isTRUE(refit_coldata_supplied)) {
      refit_ids <- NULL
      if (.bio_is_se(refit_x)) {
        refit_ids <- rownames(SummarizedExperiment::colData(refit_x))
      } else if (is.data.frame(refit_x)) {
        refit_ids <- rownames(refit_x)
        if ((is.null(refit_ids) || !all(nzchar(refit_ids))) &&
            "row_id" %in% names(refit_x)) {
          refit_ids <- as.character(refit_x[["row_id"]])
        }
      } else if (is.matrix(refit_x)) {
        refit_ids <- rownames(refit_x)
      }
      if (!is.null(refit_ids) && all(nzchar(refit_ids))) {
        refit_coldata <- align_coldata_for_perm(refit_coldata, refit_ids,
                                                context = "perm_refit")
      } else if (is.null(refit_ids) && nrow(refit_coldata) == n_refit) {
        warning("perm_refit coldata has no ids; assuming row order aligns to refit data.",
                call. = FALSE)
      } else if (is.null(refit_ids)) {
        warning("perm_refit coldata not aligned to refit data; restricted permutations disabled.",
                call. = FALSE)
        refit_coldata <- NULL
      }
    }

    y_base <- NULL
    if (!is.null(refit_coldata) && refit_outcome %in% names(refit_coldata)) {
      y_base <- refit_coldata[[refit_outcome]]
    } else {
      y_base <- .bio_get_y(refit_x, refit_outcome)
    }
    if (length(y_base) != n_refit) {
      stop("Outcome length does not match rows in perm_refit_spec$x.", call. = FALSE)
    }

    perm_full_source <- NULL
    if (!is.null(refit_coldata) &&
        refit_outcome %in% names(refit_coldata) &&
        perm_mode_use %in% c("subject_grouped", "batch_blocked", "study_loocv", "time_series")) {
      folds_perm <- list(list(test = seq_len(n_refit), fold = 1L, repeat_id = 1L))
      perm_full_source <- .permute_labels_factory(
        cd = refit_coldata, outcome = refit_outcome, mode = perm_mode_use,
        folds = folds_perm, perm_stratify = perm_stratify, time_block = time_block,
        block_len = block_len, seed = seed,
        group_col = fit@splits@info$group, batch_col = fit@splits@info$batch,
        study_col = fit@splits@info$study, time_col = fit@splits@info$time
      )
    }

    permute_outcome <- function(b) {
      if (!is.null(perm_full_source)) {
        perm_full_source(b)[[1]]
      } else {
        sample(y_base)
      }
    }

    apply_outcome <- function(x, outcome, y_perm) {
      y_perm <- .coerce_truth_like(y_base, y_perm)
      if (.bio_is_se(x)) {
        cd <- SummarizedExperiment::colData(x)
        if (!outcome %in% colnames(cd)) {
          stop("Outcome column not found in SummarizedExperiment colData.", call. = FALSE)
        }
        cd[[outcome]] <- y_perm
        SummarizedExperiment::colData(x) <- cd
        return(x)
      }
      if (is.data.frame(x)) {
        if (!outcome %in% names(x)) {
          stop("Outcome column not found in data.frame.", call. = FALSE)
        }
        x[[outcome]] <- y_perm
        return(x)
      }
      if (is.matrix(x)) {
        if (is.character(outcome) && outcome %in% colnames(x)) {
          x[, outcome] <- y_perm
          return(x)
        }
        stop("Outcome column not found in matrix input.", call. = FALSE)
      }
      stop("perm_refit_spec$x must be a SummarizedExperiment, data.frame, or matrix.", call. = FALSE)
    }

    perm_eval_refit <- function(b) {
      set.seed(seed + b)
      y_perm <- permute_outcome(b)
      x_perm <- apply_outcome(refit_x, refit_outcome, y_perm)
      fit_perm <- try(fit_resample(
        x_perm,
        outcome = refit_outcome,
        splits = fit@splits,
        preprocess = refit_preprocess,
        learner = refit_learner,
        learner_args = refit_learner_args,
        custom_learners = refit_custom_learners,
        metrics = metric,
        class_weights = refit_class_weights,
        positive_class = refit_positive_class,
        parallel = refit_parallel,
        refit = FALSE,
        seed = seed + b
      ), silent = TRUE)
      if (inherits(fit_perm, "try-error")) {
        warning(sprintf("Permutation refit %d failed: %s", b, attr(fit_perm, "condition")$message),
                call. = FALSE)
        return(NA_real_)
      }
      perm_pred <- if (length(fit_perm@predictions)) {
        do.call(rbind, lapply(fit_perm@predictions, function(df) data.frame(df, stringsAsFactors = FALSE)))
      } else {
        NULL
      }
      if (is.null(perm_pred) || !nrow(perm_pred)) return(NA_real_)
      .metric_value(metric, task, perm_pred$truth, perm_pred$pred, pred_df = perm_pred)
    }

    if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
      perm_vals <- future.apply::future_sapply(seq_len(B), perm_eval_refit, future.seed = TRUE)
    } else {
      perm_vals <- sapply(seq_len(B), perm_eval_refit)
    }
  } else {
    perm_source <- NULL
    perm_coldata <- NULL
    if (!is.null(coldata)) {
      sample_ids <- resolve_sample_ids(fit, fallback_n = nrow(fit@splits@info$coldata %||% data.frame()))
      perm_coldata <- align_coldata_for_perm(coldata, sample_ids, context = "Permutation")
    }
    if (!is.null(perm_coldata) && !is.null(outcome_col) && outcome_col %in% names(perm_coldata)) {
      folds_perm <- folds
      if (isTRUE(compact)) {
        attr(folds_perm, "fold_assignments") <- fold_assignments
      }
      perm_source <- .permute_labels_factory(
        cd = perm_coldata, outcome = outcome_col, mode = perm_mode_use,
        folds = folds_perm, perm_stratify = perm_stratify, time_block = time_block,
        block_len = block_len, seed = seed,
        group_col = fit@splits@info$group, batch_col = fit@splits@info$batch,
        study_col = fit@splits@info$study, time_col = fit@splits@info$time
      )
    }
    if (is.null(perm_source)) {
      perm_source <- function(b) {
        lapply(pred_list, function(df) {
          if (identical(task, "survival")) {
            if (!requireNamespace("survival", quietly = TRUE)) {
              stop("Package 'survival' is required for survival permutations.")
            }
            if (all(c("truth_time", "truth_event") %in% names(df))) {
              idx <- sample(seq_len(nrow(df)))
              return(survival::Surv(df$truth_time[idx], df$truth_event[idx]))
            }
          }
          sample(df$truth)
        })
      }
    }

    perm_eval <- function(b) {
      truths <- perm_source(b)
      if (length(truths) != length(pred_list)) {
        stop(
          "Permutation source returned ", length(truths),
          " truth sets for ", length(pred_list), " prediction tables"
        )
      }
      new_preds <- Map(function(df, tr, fold_idx) {
        if (!nrow(df)) return(df)
        tr_len <- if (inherits(tr, "Surv")) nrow(tr) else length(tr)
        if (tr_len != nrow(df)) {
          stop(
            "Permutation truth length (", tr_len,
            ") does not match predictions (", nrow(df),
            ") for fold ", fold_idx, "."
          )
        }
        if (identical(task, "survival") &&
            all(c("truth_time", "truth_event") %in% names(df))) {
          tr_mat <- as.matrix(tr)
          time_col <- if ("time" %in% colnames(tr_mat)) "time" else colnames(tr_mat)[1]
          status_col <- if ("status" %in% colnames(tr_mat)) "status" else colnames(tr_mat)[ncol(tr_mat)]
          df$truth_time <- tr_mat[, time_col]
          df$truth_event <- tr_mat[, status_col]
        } else {
          df$truth <- .coerce_truth_like(df$truth, tr)
        }
        df
      }, pred_list, truths, seq_along(pred_list))
      agg <- do.call(rbind, new_preds)
      .metric_value(metric, task, agg$truth, agg$pred, pred_df = agg)
    }

    if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
      perm_vals <- future.apply::future_sapply(seq_len(B), function(i) {
        set.seed(seed + i); perm_eval(i)
      }, future.seed = TRUE)
    } else {
      perm_vals <- sapply(seq_len(B), function(i) { set.seed(seed + i); perm_eval(i) })
    }
  }

  if (!is.na(metric_obs) && is.finite(metric_obs)) {
    if (!any(is.finite(perm_vals))) {
      perm_mean <- NA_real_
      perm_sd <- NA_real_
      pval <- NA_real_
      p_se <- NA_real_
      delta <- NA_real_
    } else {
      perm_mean <- mean(perm_vals, na.rm = TRUE)
      perm_sd <- stats::sd(perm_vals, na.rm = TRUE)
      finite_perm <- is.finite(perm_vals)
      pval <- if (higher_better) {
        (1 + sum(perm_vals[finite_perm] >= metric_obs, na.rm = TRUE)) / (1 + sum(finite_perm))
      } else {
        (1 + sum(perm_vals[finite_perm] <= metric_obs, na.rm = TRUE)) / (1 + sum(finite_perm))
      }
      p_se <- sqrt(pval * (1 - pval) / (sum(finite_perm) + 1))
      delta <- if (higher_better) metric_obs - perm_mean else perm_mean - metric_obs
    }
  }

  seci <- list(se = NA_real_, ci = c(NA_real_, NA_real_), z = NA_real_)
  if (ci_method == "if" && exists(".se_ci_delta", mode = "function") && is.finite(se_obs) && is.finite(delta)) {
    seci_try <- try(.se_ci_delta(delta, se_obs, perm_vals), silent = TRUE)
    if (!inherits(seci_try, "try-error")) {
      seci <- seci_try
    } else {
      seci <- list(se = NA_real_, ci = c(NA_real_, NA_real_), z = NA_real_)
    }
  } else if (ci_method == "bootstrap" && any(is.finite(perm_vals)) && is.finite(metric_obs)) {
    perm_vals_finite <- perm_vals[is.finite(perm_vals)]
    boot_vals <- replicate(boot_B, {
      sample_vals <- sample(perm_vals_finite, replace = TRUE)
      if (higher_better) metric_obs - mean(sample_vals) else mean(sample_vals) - metric_obs
    })
    alpha <- 0.025
    se_boot <- stats::sd(boot_vals)
    seci <- list(
      se = se_boot,
      ci = as.numeric(stats::quantile(boot_vals, probs = c(alpha, 1 - alpha))),
      z = ifelse(se_boot > 0, delta / se_boot, NA_real_)
    )
  }

  z_val <- if (isTRUE(include_z)) seci$z else NA_real_

  perm_df <- data.frame(
    metric_obs = metric_obs,
    perm_mean  = perm_mean,
    perm_sd    = perm_sd,
    gap        = delta,
    z          = z_val,
    p_value    = pval,
    n_perm     = length(perm_vals)
  )

  perm_df[] <- lapply(perm_df, function(x)
    if (is.numeric(x)) round(x, 6) else x)

  # --- Batch / study association with folds ---------------------------------
  ids_pred <- all_pred$id
  ids_pred_chr <- as.character(ids_pred)
  ids_all <- unique(ids_pred)
  ids_all_chr <- as.character(ids_all)
  fold_id <- if ("fold" %in% names(all_pred)) all_pred$fold else NULL
  repeat_vals <- if (length(fit@splits@indices)) {
    vapply(fit@splits@indices, function(z) as.integer(z$repeat_id %||% 1L), integer(1))
  } else {
    1L
  }
  repeat_vals <- repeat_vals[is.finite(repeat_vals)]
  has_repeats <- length(unique(repeat_vals)) > 1L
  skip_batch_assoc <- FALSE
  if ((is.null(fold_id) || !length(fold_id) || all(is.na(fold_id))) && isTRUE(has_repeats)) {
    warning("Predictions lack fold identifiers for repeated CV; batch association skipped to avoid pooling repeats.",
            call. = FALSE)
    skip_batch_assoc <- TRUE
  }

  sample_ids <- fit@info$sample_ids %||% NULL
  if (is.null(sample_ids) || !length(sample_ids)) {
    if (!is.null(fit@splits@info$coldata)) {
      cd_lookup <- fit@splits@info$coldata
      rn_cd <- rownames(cd_lookup)
      if (!is.null(rn_cd) && !anyNA(rn_cd) && all(nzchar(rn_cd)) && !anyDuplicated(rn_cd)) {
        sample_ids <- rn_cd
      } else if ("row_id" %in% names(cd_lookup)) {
        rid <- as.character(cd_lookup[["row_id"]])
        if (length(rid) == nrow(cd_lookup) && !anyNA(rid) &&
            !anyDuplicated(rid) && all(nzchar(rid))) {
          sample_ids <- rid
        }
      }
    }
  }
  if (!skip_batch_assoc && (is.null(fold_id) || !length(fold_id) || all(is.na(fold_id)))) {
    fold_id_map <- rep(NA_integer_, length(ids_all_chr))
    names(fold_id_map) <- ids_all_chr
    if (!is.null(sample_ids) && length(sample_ids)) {
      sample_ids <- as.character(sample_ids)
      max_idx <- max(vapply(fit@splits@indices, function(z) {
        te_idx <- resolve_test_idx(z)
        if (length(te_idx)) max(te_idx) else 0L
      }, integer(1)), na.rm = TRUE)
      if (length(sample_ids) < max_idx) {
        warning("Sample ID mapping is shorter than split indices; fold mapping may be incomplete.",
                call. = FALSE)
      } else {
        for (i in seq_along(fit@splits@indices)) {
          te_idx <- resolve_test_idx(fit@splits@indices[[i]])
          te_ids <- sample_ids[te_idx]
          fold_id_map[as.character(te_ids)] <- i
        }
      }
    } else if (is.numeric(ids_all)) {
      for (i in seq_along(fit@splits@indices)) {
        te_idx <- resolve_test_idx(fit@splits@indices[[i]])
        te_ids <- intersect(te_idx, ids_all)
        fold_id_map[as.character(te_ids)] <- i
      }
    }
    fold_id <- fold_id_map[ids_pred_chr]
  }

  fold_label <- fold_id
  repeat_id <- NULL
  if (!skip_batch_assoc && !is.null(fold_id) && length(fold_id) && length(fit@splits@indices)) {
    fold_meta <- data.frame(
      fold_seq = seq_along(fit@splits@indices),
      fold = vapply(fit@splits@indices, function(z) {
        as.integer(z$fold %||% NA_integer_)
      }, integer(1)),
      repeat_id = vapply(fit@splits@indices, function(z) {
        as.integer(z$repeat_id %||% 1L)
      }, integer(1)),
      stringsAsFactors = FALSE
    )
    fold_meta$fold[!is.finite(fold_meta$fold)] <- fold_meta$fold_seq[!is.finite(fold_meta$fold)]
    fold_idx <- match(fold_id, fold_meta$fold_seq)
    repeat_id <- fold_meta$repeat_id[fold_idx]
    fold_label <- fold_meta$fold[fold_idx]
  }

  if (is.null(coldata) && !is.null(fit@splits@info$coldata)) {
    coldata <- fit@splits@info$coldata
  }
  if (!is.null(coldata)) {
    rn <- rownames(coldata)
    aligned <- FALSE

    if (!is.null(rn) && all(ids_all_chr %in% rn)) {
      coldata <- coldata[ids_all_chr, , drop = FALSE]
      aligned <- TRUE
    } else if ("row_id" %in% names(coldata)) {
      rid_chr <- as.character(coldata[["row_id"]])
      if (!anyDuplicated(rid_chr) && all(ids_all_chr %in% rid_chr)) {
        coldata <- coldata[match(ids_all_chr, rid_chr), , drop = FALSE]
        aligned <- TRUE
      }
    }

    if (!aligned && is.numeric(ids_all) && max(ids_all, na.rm = TRUE) <= nrow(coldata)) {
      if (!is.null(rn) && !all(ids_all_chr %in% rn)) {
        warning("`coldata` row names do not match prediction ids; aligning by row order.")
      }
      coldata <- coldata[ids_all, , drop = FALSE]
      aligned <- TRUE
    }

    if (!aligned) {
      warning("`coldata` not aligned to predictions; batch association skipped.")
      coldata <- NULL
    }
  }

  coldata_pred <- NULL
  if (!is.null(coldata)) {
    idx_pred <- match(ids_pred_chr, ids_all_chr)
    if (anyNA(idx_pred)) {
      warning("`coldata` not aligned to predictions; batch association skipped.")
    } else {
      coldata_pred <- coldata[idx_pred, , drop = FALSE]
    }
  }

  batch_df <- data.frame()
  if (!skip_batch_assoc && !is.null(coldata_pred)) {
    if (is.null(batch_cols)) {
      batch_cols <- intersect(c("batch", "plate", "center", "site", "study"), colnames(coldata_pred))
    }
    if (length(batch_cols) > 0) {
      fold_valid <- !is.na(fold_label)
      repeat_valid <- !is.null(repeat_id) && !all(is.na(repeat_id))
      batch_results <- lapply(batch_cols, function(bc) {
        if (!bc %in% colnames(coldata_pred)) return(NULL)
        if (isTRUE(repeat_valid) && any(fold_valid, na.rm = TRUE)) {
          reps <- sort(unique(repeat_id[fold_valid]))
          reps <- reps[is.finite(reps)]
          rep_rows <- lapply(reps, function(r) {
            idx <- fold_valid & repeat_id == r
            if (!any(idx)) return(NULL)
            tab <- table(
              fold = fold_label[idx],
              batch = as.factor(coldata_pred[idx, bc])
            )
            df <- as.data.frame(.chisq_assoc(tab))
            df$repeat_id <- r
            df
          })
          rep_rows <- rep_rows[!vapply(rep_rows, is.null, logical(1))]
          if (!length(rep_rows)) return(NULL)
          do.call(rbind, rep_rows)
        } else {
          tab <- table(
            fold = fold_label[fold_valid],
            batch = as.factor(coldata_pred[fold_valid, bc])
          )
          df <- as.data.frame(.chisq_assoc(tab))
          df$repeat_id <- 1L
          df
        }
      })
      names(batch_results) <- batch_cols
      keep_idx <- !vapply(batch_results, is.null, logical(1))
      if (any(keep_idx)) {
        batch_results <- batch_results[keep_idx]
        batch_df <- do.call(rbind, Map(function(nm, df) { df$batch_col <- nm; df }, names(batch_results), batch_results))
        batch_df <- batch_df[, c("batch_col", "repeat_id", "stat", "df", "pval", "cramer_v")]
      }
    }
  }

  # --- Target leakage scan (feature-wise outcome association) ----------------
  target_df <- data.frame()
  if (isTRUE(target_scan)) {
    y_outcome <- NULL
    if (!is.null(coldata) && !is.null(outcome_col) && outcome_col %in% names(coldata)) {
      y_outcome <- coldata[[outcome_col]]
    }
    if (is.null(y_outcome) || !length(y_outcome)) {
      id_first <- !duplicated(as.character(all_pred$id))
      truth_map <- all_pred$truth[id_first]
      names(truth_map) <- as.character(all_pred$id[id_first])
      y_outcome <- truth_map[ids_all_chr]
    }
    if (is.null(X_ref) && !is.null(fit@info$X_ref)) X_ref <- fit@info$X_ref
    if (!is.null(X_ref) && !is.null(y_outcome)) {
      X_scan <- .align_by_ids(X_ref, ids_all_chr, sample_ids = sample_ids)
      if (!is.null(X_scan)) {
        target_df <- .target_assoc_scan(
          X_scan, y_outcome, task,
          positive_class = fit@info$positive_class,
          threshold = target_threshold
        )
      }
    }
  }

  # --- Multivariate/interaction target scan ----------------------------------
  target_multivariate <- data.frame()
  if (isTRUE(target_scan_multivariate)) {
    if (is.null(X_ref) && !is.null(fit@info$X_ref)) X_ref <- fit@info$X_ref
    max_idx <- suppressWarnings(max(vapply(folds, function(z) {
      idx <- c(z$train, z$test)
      if (length(idx)) max(idx) else 0L
    }, integer(1)), na.rm = TRUE))
    if (!is.finite(max_idx) || max_idx < 1L) max_idx <- NA_integer_
    sample_ids_scan <- resolve_sample_ids(fit, fallback_n = max_idx)
    if (!is.null(sample_ids_scan) && length(sample_ids_scan)) {
      sample_ids_scan <- as.character(sample_ids_scan)
    }
    coldata_scan <- fit@splits@info$coldata %||% NULL
    if (!is.null(coldata_scan) && !is.null(sample_ids_scan)) {
      coldata_scan <- .align_by_ids(coldata_scan, sample_ids_scan,
                                    sample_ids = sample_ids_scan, warn = FALSE)
    }
    y_outcome_scan <- NULL
    if (!is.null(coldata_scan) && !is.null(outcome_col) &&
        outcome_col %in% names(coldata_scan)) {
      y_outcome_scan <- coldata_scan[[outcome_col]]
    }
    if (is.null(y_outcome_scan) || !length(y_outcome_scan)) {
      id_first <- !duplicated(as.character(all_pred$id))
      truth_map <- all_pred$truth[id_first]
      names(truth_map) <- as.character(all_pred$id[id_first])
      if (!is.null(sample_ids_scan)) {
        y_outcome_scan <- truth_map[sample_ids_scan]
      }
    }
    if (!is.null(X_ref) && !is.null(y_outcome_scan) && length(y_outcome_scan) &&
        !is.null(sample_ids_scan) && length(sample_ids_scan)) {
      X_scan_mv <- .align_by_ids(X_ref, sample_ids_scan, sample_ids = sample_ids_scan)
      if (!is.null(X_scan_mv)) {
        permute_outcome <- NULL
        if (!is.null(coldata_scan) && !is.null(outcome_col) &&
            outcome_col %in% names(coldata_scan) &&
            perm_mode_use %in% c("subject_grouped", "batch_blocked", "study_loocv", "time_series")) {
          folds_perm <- list(list(test = seq_len(length(y_outcome_scan)), fold = 1L, repeat_id = 1L))
          perm_source <- .permute_labels_factory(
            cd = coldata_scan, outcome = outcome_col, mode = perm_mode_use,
            folds = folds_perm, perm_stratify = perm_stratify, time_block = time_block,
            block_len = block_len, seed = seed,
            group_col = fit@splits@info$group, batch_col = fit@splits@info$batch,
            study_col = fit@splits@info$study, time_col = fit@splits@info$time
          )
          permute_outcome <- function(b) perm_source(b)[[1]]
        }
        target_multivariate <- .target_scan_multivariate(
          X_scan_mv, y_outcome_scan, task, fit@splits, folds, coldata_scan, seed,
          B = target_scan_multivariate_B,
          max_components = target_scan_multivariate_components,
          interactions = target_scan_multivariate_interactions,
          positive_class = fit@info$positive_class,
          permute_outcome = permute_outcome
        )
      }
    }
  }

  # --- Duplicate / near-duplicate detection ---------------------------------
  dup_df <- data.frame()
  if (is.null(X_ref) && !is.null(fit@info$X_ref)) X_ref <- fit@info$X_ref
  if (!is.null(X_ref)) {
    n_samples <- NA_integer_
    if (!is.null(fit@info$sample_ids)) {
      n_samples <- length(fit@info$sample_ids)
    }
    if (!is.finite(n_samples) || n_samples < 1L) {
      if (!is.null(fit@splits@info$coldata)) {
        n_samples <- nrow(fit@splits@info$coldata)
      }
    }
    if (!is.finite(n_samples) || n_samples < 1L) {
      if (isTRUE(compact) && length(fold_assignments)) {
        n_samples <- max(vapply(fold_assignments, length, integer(1)), na.rm = TRUE)
      }
    }
    if (!is.finite(n_samples) || n_samples < 1L) {
      n_samples <- suppressWarnings(max(vapply(folds, function(z) {
        idx <- c(z$train, z$test)
        if (length(idx)) max(idx) else 0L
      }, integer(1)), na.rm = TRUE))
    }
    if (!is.finite(n_samples) || n_samples < 1L) n_samples <- NA_integer_

    sample_ids_all <- resolve_sample_ids(fit, fallback_n = n_samples)
    if (!is.null(sample_ids_all) && length(sample_ids_all)) {
      sample_ids_all <- as.character(sample_ids_all)
    }

    X_use <- NULL
    if (!is.null(sample_ids_all) && length(sample_ids_all)) {
      X_use <- .align_by_ids(
        X_ref, sample_ids_all, sample_ids = sample_ids_all,
        warn = identical(duplicate_scope, "train_test")
      )
    }
    if (is.null(X_use)) {
      if (identical(duplicate_scope, "train_test")) {
        warning("X_ref not aligned to splits; duplicate_scope='train_test' skipped.", call. = FALSE)
      } else {
        warning("X_ref not aligned to splits; proceeding with duplicate_scope='all'.", call. = FALSE)
        X_use <- X_ref
      }
    }

    if (is.null(X_use)) {
      # skip duplicate detection when train/test alignment is required
    } else {
      X <- as.matrix(X_use)
    # choose feature space
    if (feature_space == "rank") {
      X <- t(apply(X, 1, function(row) rank(row, ties.method = "average", na.last = "keep")))
      X[is.na(X)] <- 0
    }
    if (sim_method == "pearson") {
      X <- .row_center(X)
    }
    # standardize to enable similarity via Euclidean NN
    Xn <- .row_l2_normalize(X)

    n <- nrow(Xn)
    candidate_pairs <- NULL

    if (n > 3000 && requireNamespace("RANN", quietly = TRUE)) {
      # approximate k-NN search; cosine/pearson ~ 1 - 0.5*||u - v||^2 for unit vectors
      nn <- RANN::nn2(Xn, k = min(nn_k + 1, n))  # includes self
      idx <- rep(seq_len(n), times = ncol(nn$nn.idx))
      jdx <- as.vector(nn$nn.idx)
      mask <- (idx != jdx)
      idx <- idx[mask]
      jdx <- jdx[mask]
      # compute similarity for candidate pairs
      simv <- rowSums(Xn[idx, , drop = FALSE] * Xn[jdx, , drop = FALSE])
      keep <- which(simv >= sim_threshold)
      if (length(keep)) {
        candidate_pairs <- data.frame(i = idx[keep], j = jdx[keep], sim = simv[keep])
      }
    } else {
      # exact; cap output
      S <- .cosine_sim_block(Xn)
      S[upper.tri(S, TRUE)] <- 0
      which_dup <- which(S >= sim_threshold, arr.ind = TRUE)
      if (nrow(which_dup) > 0)
        candidate_pairs <- data.frame(i = which_dup[, 1], j = which_dup[, 2], sim = S[which_dup])
    }

    if (!is.null(candidate_pairs)) {
      compute_cross_fold <- function(pairs) {
        n_pairs <- nrow(pairs)
        if (!n_pairs) return(logical(0))
        if (!is.finite(n_samples) || n_samples < 1L) return(rep(NA, n_pairs))
        if (is.null(folds) || !length(folds)) return(rep(NA, n_pairs))
        split_mode <- fit@splits@mode %||% NA_character_

        if (!identical(split_mode, "time_series")) {
          fold_map_by_repeat <- NULL
          if (isTRUE(compact) && length(fold_assignments)) {
            ok_len <- all(vapply(fold_assignments, function(vec) {
              length(vec) == n_samples
            }, logical(1)))
            if (isTRUE(ok_len)) {
              fold_map_by_repeat <- lapply(fold_assignments, function(vec) as.integer(vec))
            }
          }
          if (is.null(fold_map_by_repeat)) {
            repeat_ids <- vapply(folds, function(f) {
              as.integer(f$repeat_id %||% 1L)
            }, integer(1))
            n_rep <- suppressWarnings(max(repeat_ids, na.rm = TRUE))
            if (!is.finite(n_rep) || n_rep < 1L) return(rep(NA, n_pairs))
            fold_map_by_repeat <- lapply(seq_len(n_rep), function(r) rep(NA_integer_, n_samples))
            for (i in seq_along(folds)) {
              r <- repeat_ids[[i]]
              test <- resolve_test_idx(folds[[i]])
              if (length(test)) {
                fold_id <- folds[[i]]$fold %||% i
                fold_map_by_repeat[[r]][test] <- fold_id
              }
            }
          }
          cross <- rep(FALSE, n_pairs)
          i_idx <- pairs$i
          j_idx <- pairs$j
          for (fm in fold_map_by_repeat) {
            if (length(fm) != n_samples) next
            fi <- fm[i_idx]
            fj <- fm[j_idx]
            valid <- !is.na(fi) & !is.na(fj)
            cross <- cross | (valid & fi != fj)
          }
          return(cross)
        }

        fold_map <- rep(NA_integer_, n_samples)
        time_vec <- NULL
        time_col <- fit@splits@info$time %||% NULL
        cd_time <- fit@splits@info$coldata %||% NULL
        if (!is.null(cd_time) && !is.null(time_col) && time_col %in% names(cd_time)) {
          time_vec <- cd_time[[time_col]]
        }
        fold_tmin <- NULL
        if (!is.null(time_vec) && length(time_vec) == n_samples) {
          fold_tmin <- rep(time_vec[1], length(folds))
          fold_tmin[] <- NA
        }
        for (i in seq_along(folds)) {
          test <- resolve_test_idx(folds[[i]])
          if (!length(test)) next
          fold_map[test] <- i
          if (!is.null(fold_tmin)) {
            fold_tmin[i] <- suppressWarnings(min(time_vec[test]))
          }
        }
        if (!is.null(time_vec) && length(time_vec) == n_samples && !is.null(fold_tmin)) {
          split_horizon <- fit@splits@info$horizon %||% 0
          cross <- rep(FALSE, n_pairs)
          i_idx <- pairs$i
          j_idx <- pairs$j
          fi <- fold_map[i_idx]
          fj <- fold_map[j_idx]
          idx_i <- which(!is.na(fi))
          if (length(idx_i)) {
            tmin_i <- fold_tmin[fi[idx_i]]
            valid_i <- which(!is.na(tmin_i))
            if (length(valid_i)) {
              idx_i <- idx_i[valid_i]
              tmin_i <- tmin_i[valid_i]
              if (split_horizon == 0) {
                cross[idx_i] <- time_vec[j_idx[idx_i]] < tmin_i
              } else {
                cutoff <- tmin_i - split_horizon
                cross[idx_i] <- time_vec[j_idx[idx_i]] <= cutoff
              }
            }
          }
          idx_j <- which(!is.na(fj))
          if (length(idx_j)) {
            tmin_j <- fold_tmin[fj[idx_j]]
            valid_j <- which(!is.na(tmin_j))
            if (length(valid_j)) {
              idx_j <- idx_j[valid_j]
              tmin_j <- tmin_j[valid_j]
              if (split_horizon == 0) {
                cross[idx_j] <- cross[idx_j] | (time_vec[i_idx[idx_j]] < tmin_j)
              } else {
                cutoff <- tmin_j - split_horizon
                cross[idx_j] <- cross[idx_j] | (time_vec[i_idx[idx_j]] <= cutoff)
              }
            }
          }
          return(cross)
        }

        has_train <- any(vapply(folds, function(f) !is.null(f$train), logical(1)))
        if (!has_train) return(rep(NA, n_pairs))
        cross <- rep(FALSE, n_pairs)
        for (i in seq_along(folds)) {
          fold <- folds[[i]]
          if (is.null(fold$train) || is.null(fold$test)) next
          train <- fold$train
          test <- fold$test
          if (!length(train) || !length(test)) next
          train_mask <- rep(FALSE, n_samples)
          test_mask <- rep(FALSE, n_samples)
          train_mask[train] <- TRUE
          test_mask[test] <- TRUE
          cross <- cross | (test_mask[pairs$i] & train_mask[pairs$j]) |
            (test_mask[pairs$j] & train_mask[pairs$i])
        }
        cross
      }

      cross_fold <- compute_cross_fold(candidate_pairs)
      if (length(cross_fold)) {
        candidate_pairs$cross_fold <- cross_fold
      }
      if (identical(duplicate_scope, "train_test")) {
        keep <- which(!is.na(cross_fold) & cross_fold)
        if (length(keep)) {
          candidate_pairs <- candidate_pairs[keep, , drop = FALSE]
        } else {
          candidate_pairs <- candidate_pairs[0, , drop = FALSE]
        }
      }

      if (nrow(candidate_pairs) > 0) {
        if (sim_method == "cosine") {
          candidate_pairs$cos_sim <- candidate_pairs$sim
        } else if (sim_method == "pearson") {
          candidate_pairs$pearson <- candidate_pairs$sim
        }
        # sort and cap
        ord <- order(candidate_pairs$sim, decreasing = TRUE)
        candidate_pairs <- candidate_pairs[ord, , drop = FALSE]
        if (nrow(candidate_pairs) > max_pairs)
          candidate_pairs <- candidate_pairs[seq_len(max_pairs), , drop = FALSE]
        dup_df <- candidate_pairs
      }
    }
    }
  }

  # --- Assemble S4 object ----------------------------------------------------
  perm_values <- if (isTRUE(return_perm)) perm_vals else numeric(0)

  new("LeakAudit",
      fit = fit,
      permutation_gap = perm_df,
      perm_values = perm_values,
      batch_assoc = batch_df,
      target_assoc = target_df,
      duplicates = dup_df,
      trail = trail,
      info = list(
        n_perm = B,
        perm_stratify = perm_stratify,
        perm_method = perm_method,
        perm_mode = perm_mode_use,
        perm_refit_mode = perm_refit_mode,
        perm_refit_reason = perm_refit_reason,
        perm_refit_auto_max = perm_refit_auto_max,
        parallel = parallel,
        target_scan = target_scan,
        target_scan_multivariate = target_scan_multivariate,
        target_scan_multivariate_B = target_scan_multivariate_B,
        target_scan_multivariate_components = target_scan_multivariate_components,
        target_scan_multivariate_interactions = target_scan_multivariate_interactions,
        target_multivariate = target_multivariate,
        target_threshold = target_threshold,
        sim_method = sim_method,
        feature_space = feature_space,
        duplicate_threshold = sim_threshold,
        duplicate_scope = duplicate_scope,
        nn_k = nn_k,
        max_pairs = max_pairs,
        batch_cols = batch_cols,
        permutation_se = p_se,
        ci_delta = seci$ci,
        se_delta = seci$se,
        ci_method = ci_method,
        include_z = include_z
      ))
}

#' Audit leakage per learner
#'
#' Runs [audit_leakage()] separately for each learner recorded in a [LeakFit]
#' and returns a named list of [LeakAudit] objects. Use this when a single fit
#' contains predictions for multiple models and you want model-specific audits.
#' If predictions do not include learner IDs, only a single audit can be run and
#' requesting multiple learners is an error.
#'
#' @param fit A [LeakFit] object produced by [fit_resample()]. It must contain
#'   predictions and split metadata. Learner IDs must be present in predictions
#'   to audit multiple models.
#' @param metric Character scalar. One of `"auc"`, `"pr_auc"`, `"accuracy"`,
#'   `"macro_f1"`, `"log_loss"`, `"rmse"`, or `"cindex"`. Controls which metric
#'   is audited for each learner.
#' @param learners Character vector or NULL. If NULL (default), audits all
#'   learners found in predictions. If provided, must match learner IDs stored in
#'   the predictions. Supplying more than one learner requires learner IDs.
#' @param parallel_learners Logical scalar. If TRUE, runs per-learner audits in
#'   parallel using `future.apply` (if installed). This changes runtime but not
#'   the audit results.
#' @param mc.cores Integer scalar or NULL. Number of workers used when
#'   `parallel_learners = TRUE`. Defaults to the minimum of available cores and
#'   the number of learners.
#' @param ... Additional named arguments forwarded to [audit_leakage()] for each
#'   learner. These control the audit itself. Common options include:
#'   `B` (integer permutations), `perm_stratify` (logical or `"auto"`),
#'   `perm_refit` (logical), `perm_refit_spec` (list),
#'   `time_block` (character), `block_len` (integer or NULL), `include_z`
#'   (logical), `ci_method` (character), `boot_B` (integer), `parallel`
#'   (logical), `seed` (integer), `return_perm` (logical), `batch_cols`
#'   (character vector), `coldata` (data.frame), `X_ref` (matrix/data.frame),
#'   `target_scan` (logical), `target_threshold` (numeric), `feature_space`
#'   (character), `sim_method` (character), `sim_threshold` (numeric),
#'   `nn_k` (integer), `max_pairs` (integer), and `duplicate_scope` (character).
#'   See [audit_leakage()] for full definitions; changing these values changes
#'   each learner's audit.
#' @return Named list of [LeakAudit] objects, keyed by learner ID.
#' @examples
#' \dontrun{
#' set.seed(1)
#' df <- data.frame(
#'   subject = rep(1:6, each = 2),
#'   outcome = factor(rep(c(0, 1), 6)),
#'   x1 = rnorm(12),
#'   x2 = rnorm(12)
#' )
#' splits <- make_split_plan(df, outcome = "outcome",
#'                       mode = "subject_grouped", group = "subject",
#'                       v = 3, progress = FALSE)
#' custom <- list(
#'   glm = list(
#'     fit = function(x, y, task, weights, ...) {
#'       stats::glm(y ~ ., data = data.frame(y = y, x),
#'                  family = stats::binomial(), weights = weights)
#'     },
#'     predict = function(object, newdata, task, ...) {
#'       as.numeric(stats::predict(object,
#'                                 newdata = as.data.frame(newdata),
#'                                 type = "response"))
#'     }
#'   )
#' )
#' custom$glm2 <- custom$glm
#' fit <- fit_resample(df, outcome = "outcome", splits = splits,
#'                     learner = c("glm", "glm2"), custom_learners = custom,
#'                     metrics = "auc", refit = FALSE, seed = 1)
#' audits <- audit_leakage_by_learner(fit, metric = "auc", B = 10,
#'                                    perm_stratify = FALSE)
#' names(audits)
#' }
#' @export
audit_leakage_by_learner <- function(fit,
                                     metric = c("auc", "pr_auc", "accuracy", "macro_f1", "log_loss", "rmse", "cindex"),
                                     learners = NULL,
                                     parallel_learners = FALSE,
                                     mc.cores = NULL,
                                     ...) {
  metric <- match.arg(metric)
  dots <- list(...)
  if ("learner" %in% names(dots)) {
    stop("Use `learners` to select models; do not pass `learner` to audit_leakage_by_learner().")
  }

  pred_list <- lapply(fit@predictions, function(df) data.frame(df, stringsAsFactors = FALSE))
  pred_df <- if (length(pred_list)) do.call(rbind, pred_list) else NULL
  if (is.null(pred_df) || !nrow(pred_df)) {
    stop("No predictions available in LeakFit object.")
  }

  has_learner <- "learner" %in% names(pred_df)
  if (is.null(learners)) {
    if (has_learner) {
      learners <- unique(as.character(pred_df$learner))
    } else if (!is.null(fit@metrics) && "learner" %in% names(fit@metrics)) {
      learners <- unique(as.character(fit@metrics$learner))
    } else {
      learners <- character(0)
    }
  }

  if (!length(learners)) {
    stop("No learner names found to audit.")
  }

  if (has_learner) {
    available <- unique(as.character(pred_df$learner))
    missing <- setdiff(learners, available)
    if (length(missing)) {
      stop(sprintf("Learner(s) not found in predictions: %s", paste(missing, collapse = ", ")))
    }
  } else if (length(learners) > 1L) {
    stop("Predictions do not include learner IDs; cannot audit per learner. Refit with updated bioLeak.")
  } else {
    warning("Predictions do not include learner IDs; returning a single audit without learner filtering.")
  }

  if (isTRUE(parallel_learners)) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      warning("parallel_learners=TRUE requires the 'future.apply' package; using sequential execution.")
      parallel_learners <- FALSE
    }
  }
  if (isTRUE(parallel_learners)) {
    if (is.null(mc.cores)) {
      detected <- parallel::detectCores()
      if (is.na(detected) || detected < 1L) detected <- 1L
      mc.cores <- min(length(learners), detected)
    }
    mc.cores <- max(1L, as.integer(mc.cores))
    if (requireNamespace("future", quietly = TRUE)) {
      old_plan <- future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan(future::multisession, workers = mc.cores)
    }
  }

  runner <- function(ln) {
    audit_leakage(fit, metric = metric, learner = if (has_learner) ln else NULL, ...)
  }
  audits <- if (isTRUE(parallel_learners)) {
    future.apply::future_lapply(learners, runner, future.seed = TRUE)
  } else {
    lapply(learners, runner)
  }
  names(audits) <- learners
  class(audits) <- c("LeakAuditList", class(audits))
  audits
}

#' @export
print.LeakAuditList <- function(x, digits = 3, ...) {
  n <- length(x)
  cat(sprintf("LeakAuditList with %d learner%s\n", n, ifelse(n == 1L, "", "s")))
  if (!n) return(invisible(x))

  learners <- names(x)
  if (is.null(learners) || any(!nzchar(learners))) {
    learners <- paste0("learner_", seq_len(n))
  }

  summary_rows <- lapply(seq_along(x), function(i) {
    aud <- x[[i]]
    pg <- if (inherits(aud, "LeakAudit")) aud@permutation_gap else NULL
    metric_obs <- NA_real_
    perm_mean <- NA_real_
    gap <- NA_real_
    p_value <- NA_real_
    z_val <- NA_real_
    n_perm <- NA_real_
    if (!is.null(pg) && nrow(pg) > 0) {
      metric_obs <- pg$metric_obs[1]
      perm_mean <- pg$perm_mean[1]
      gap <- pg$gap[1]
      p_value <- pg$p_value[1]
      if ("z" %in% names(pg)) z_val <- pg$z[1]
      if ("n_perm" %in% names(pg)) n_perm <- pg$n_perm[1]
    }

    batch_df <- if (inherits(aud, "LeakAudit")) aud@batch_assoc else NULL
    batch_p_min <- NA_real_
    batch_col_min_p <- NA_character_
    batch_v_max <- NA_real_
    batch_col_max_v <- NA_character_
    if (!is.null(batch_df) && nrow(batch_df) > 0) {
      pvals <- batch_df$pval
      if (any(is.finite(pvals))) {
        idx <- which.min(pvals)
        batch_p_min <- pvals[idx]
        batch_col_min_p <- as.character(batch_df$batch_col[idx])
      }
      vvals <- batch_df$cramer_v
      if (any(is.finite(vvals))) {
        idx <- which.max(vvals)
        batch_v_max <- vvals[idx]
        batch_col_max_v <- as.character(batch_df$batch_col[idx])
      }
    }

    dup_df <- if (inherits(aud, "LeakAudit")) aud@duplicates else NULL
    dup_pairs <- NA_real_
    dup_max_sim <- NA_real_
    if (!is.null(dup_df)) {
      if (nrow(dup_df) > 0) {
        dup_pairs <- nrow(dup_df)
        if ("sim" %in% names(dup_df)) {
          dup_max_sim <- max(dup_df$sim, na.rm = TRUE)
        } else if ("cos_sim" %in% names(dup_df)) {
          dup_max_sim <- max(dup_df$cos_sim, na.rm = TRUE)
        } else if ("pearson" %in% names(dup_df)) {
          dup_max_sim <- max(dup_df$pearson, na.rm = TRUE)
        }
      } else if (ncol(dup_df) > 0) {
        dup_pairs <- 0
      }
    }
    dup_threshold <- if (inherits(aud, "LeakAudit") &&
                         !is.null(aud@info$duplicate_threshold)) {
      aud@info$duplicate_threshold
    } else {
      NA_real_
    }

    data.frame(
      learner = learners[[i]],
      metric_obs = metric_obs,
      perm_mean = perm_mean,
      gap = gap,
      p_value = p_value,
      z = z_val,
      n_perm = n_perm,
      batch_p_min = batch_p_min,
      batch_col_min_p = batch_col_min_p,
      batch_v_max = batch_v_max,
      batch_col_max_v = batch_col_max_v,
      dup_pairs = dup_pairs,
      dup_threshold = dup_threshold,
      dup_max_sim = dup_max_sim,
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, summary_rows)
  if (nrow(summary_df)) {
    numeric_cols <- vapply(summary_df, is.numeric, logical(1))
    summary_df[numeric_cols] <- lapply(summary_df[numeric_cols], function(v) round(v, digits))
    print(summary_df, row.names = FALSE)
  }
  invisible(x)
}

.coerce_truth_like <- function(template, values) {
  if (inherits(template, "Surv")) {
    values
  } else if (is.factor(template)) {
    factor(values, levels = levels(template))
  } else if (is.logical(template)) {
    as.logical(values)
  } else if (is.numeric(template)) {
    as.numeric(values)
  } else {
    values
  }
}
