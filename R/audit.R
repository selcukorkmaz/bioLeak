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
#' (1) compares observed cross-validated performance to a permutation null to
#' detect label or split leakage, (2) tests whether fold assignments are
#' associated with batch or study metadata (confounding by design), (3) scans
#' features for unusually strong outcome proxies, and (4) flags duplicate or
#' near-duplicate samples in a reference feature matrix.
#'
#' The returned [LeakAudit] summarizes these diagnostics. It relies only on the
#' stored predictions, splits, and optional metadata; it does not refit models
#' or alter splits. Results are conditional on the chosen metric and supplied
#' metadata/features and should be interpreted as diagnostics, not proof of
#' leakage or its absence.
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
#'   distribution (default 1000). Larger values reduce Monte Carlo error but
#'   increase runtime.
#' @param perm_stratify Logical scalar or `"auto"`. If TRUE (default), permutations
#'   are stratified within each fold (factor levels; numeric outcomes are binned
#'   into quantiles when enough non-missing values are available). If FALSE, no
#'   stratification is used. Stratification only applies when `coldata` supplies
#'   the outcome; otherwise labels are shuffled within each fold.
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
#'   `fit@splits@info$coldata` when available. Misalignment triggers warnings
#'   and skips the affected checks.
#' @param X_ref Optional numeric matrix/data.frame (samples x features). Used for
#'   duplicate detection and the target leakage scan. If NULL, uses
#'   `fit@info$X_ref` when available. Rows must align to prediction ids via
#'   row names, sample ids, or row order; misalignment disables these checks.
#' @param target_scan Logical scalar. If TRUE (default), computes per-feature
#'   outcome associations on `X_ref` and flags proxy features; if FALSE, or if
#'   `X_ref`/outcomes are unavailable, `target_assoc` is empty. Not available
#'   for survival outcomes.
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
#'
#' `batch_assoc` contains chi-square tests between fold assignment and each
#' `batch_cols` variable (`stat`, `df`, `pval`, `cramer_v`). `target_assoc`
#' reports feature-wise outcome associations on `X_ref`; numeric features use
#' AUC (binomial), `eta_sq` (multiclass), or correlation (gaussian), while
#' categorical features use Cramer's V (binomial/multiclass) or `eta_sq` from a
#' one-way ANOVA (gaussian). The
#' `score` column is the scaled effect size used for flagging
#' (`flag = score >= target_threshold`).
#'
#' Duplicate detection compares rows of `X_ref` using the chosen `sim_method`
#' (cosine on L2-normalized rows, or Pearson via row-centering), optionally after
#' rank transformation (`feature_space = "rank"`). The `duplicates` slot returns
#' index pairs and similarity values for near-duplicate samples. Only duplicates
#' present in `X_ref` can be detected, and checks are skipped if inputs cannot
#' be aligned to prediction ids.
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   subject = rep(1:6, each = 2),
#'   outcome = rbinom(12, 1, 0.5),
#'   x1 = rnorm(12),
#'   x2 = rnorm(12)
#' )
#'
#' splits <- make_splits(df, outcome = "outcome",
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
#' @export
audit_leakage <- function(fit,
                          metric = c("auc", "pr_auc", "accuracy", "macro_f1", "log_loss", "rmse", "cindex"),
                          B = 1000,
                          perm_stratify = TRUE,
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
                          target_threshold = 0.9,
                          feature_space = c("raw", "rank"),
                          sim_method = c("cosine", "pearson"),
                          sim_threshold = 0.995,
                          nn_k = 50,
                          max_pairs = 5000,
                          learner = NULL) {

  metric <- match.arg(metric)
  feature_space <- match.arg(feature_space)
  sim_method    <- match.arg(sim_method)
  time_block <- match.arg(time_block)
  ci_method <- match.arg(ci_method)

  if (!is.numeric(target_threshold) || length(target_threshold) != 1L ||
      !is.finite(target_threshold) || target_threshold <= 0 || target_threshold >= 1) {
    stop("target_threshold must be a single numeric value in (0,1).")
  }

  set.seed(seed)

  # Trail / provenance
  trail <- list(
    indices_hash = .bio_hash_indices(fit@splits@indices),
    mode = fit@splits@mode,
    info = fit@splits@info,
    seed = seed
  )
  trail$metric <- metric
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
  perm_source <- NULL
  if (!is.null(coldata) && !is.null(outcome_col) && outcome_col %in% names(coldata)) {
    folds_perm <- folds
    if (isTRUE(compact)) {
      attr(folds_perm, "fold_assignments") <- fold_assignments
    }
    perm_source <- .permute_labels_factory(
      cd = coldata, outcome = outcome_col, mode = fit@splits@mode,
      folds = folds_perm, perm_stratify = perm_stratify, time_block = time_block,
      block_len = block_len, seed = seed,
      group_col = fit@splits@info$group, batch_col = fit@splits@info$batch,
      study_col = fit@splits@info$study
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
  ids_all <- unique(all_pred$id)
  ids_all_chr <- as.character(ids_all)
  fold_id <- rep(NA_integer_, length(ids_all_chr))
  names(fold_id) <- ids_all_chr

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
        fold_id[as.character(te_ids)] <- i
      }
    }
  } else if (is.numeric(ids_all)) {
    for (i in seq_along(fit@splits@indices)) {
      te_idx <- resolve_test_idx(fit@splits@indices[[i]])
      te_ids <- intersect(te_idx, ids_all)
      fold_id[as.character(te_ids)] <- i
    }
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

  batch_df <- data.frame()
  if (!is.null(coldata)) {
    if (is.null(batch_cols)) {
      batch_cols <- intersect(c("batch", "plate", "center", "site", "study"), colnames(coldata))
    }
    if (length(batch_cols) > 0) {
      fold_valid <- !is.na(fold_id)
      batch_results <- lapply(batch_cols, function(bc) {
        if (!bc %in% colnames(coldata)) return(NULL)
        tab <- table(
          fold = fold_id[fold_valid],
          batch = as.factor(coldata[fold_valid, bc])
        )
        as.data.frame(.chisq_assoc(tab))
      })
      names(batch_results) <- batch_cols
      keep_idx <- !vapply(batch_results, is.null, logical(1))
      if (any(keep_idx)) {
        batch_results <- batch_results[keep_idx]
        batch_df <- do.call(rbind, Map(function(nm, df) { df$batch_col <- nm; df }, names(batch_results), batch_results))
        batch_df <- batch_df[, c("batch_col", "stat", "df", "pval", "cramer_v")]
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

  # --- Duplicate / near-duplicate detection ---------------------------------
  dup_df <- data.frame()
  if (is.null(X_ref) && !is.null(fit@info$X_ref)) X_ref <- fit@info$X_ref
  if (!is.null(X_ref)) {
    X <- as.matrix(X_ref)
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
        parallel = parallel,
        target_scan = target_scan,
        target_threshold = target_threshold,
        sim_method = sim_method,
        feature_space = feature_space,
        duplicate_threshold = sim_threshold,
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
#'   `time_block` (character), `block_len` (integer or NULL), `include_z`
#'   (logical), `ci_method` (character), `boot_B` (integer), `parallel`
#'   (logical), `seed` (integer), `return_perm` (logical), `batch_cols`
#'   (character vector), `coldata` (data.frame), `X_ref` (matrix/data.frame),
#'   `target_scan` (logical), `target_threshold` (numeric), `feature_space`
#'   (character), `sim_method` (character), `sim_threshold` (numeric),
#'   `nn_k` (integer), and `max_pairs` (integer). See [audit_leakage()] for full
#'   definitions; changing these values changes each learner's audit.
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
#' splits <- make_splits(df, outcome = "outcome",
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
