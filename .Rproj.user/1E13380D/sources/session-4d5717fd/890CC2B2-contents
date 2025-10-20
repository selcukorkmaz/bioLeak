# audit_leakage(): permutation gap, batch association, duplicate detection, trail ----

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

# Safe χ² with Cramér's V
.chisq_assoc <- function(tab) {
  if (any(dim(tab) < 2) || any(tab == 0)) {
    return(list(stat = NA_real_, df = NA_integer_, pval = NA_real_, cramer_v = NA_real_))
  }
  cs <- suppressWarnings(stats::chisq.test(tab))
  n <- sum(tab)
  v <- sqrt(unname(cs$statistic) / (n * (min(dim(tab)) - 1)))
  list(stat = unname(cs$statistic), df = unname(cs$parameter),
       pval = unname(cs$p.value), cramer_v = as.numeric(v))
}

# Compute main metric from CV predictions
.main_metric <- function(task, truth, pred) {
  if (task == "binomial") {
    if (requireNamespace("pROC", quietly = TRUE)) {
      roc <- pROC::roc(truth, pred, quiet = TRUE)
      as.numeric(pROC::auc(roc))
    } else {
      yb <- if (is.factor(truth)) as.numeric(truth) - 1 else truth
      pos <- pred[yb == 1]; neg <- pred[yb == 0]
      if (length(pos) && length(neg)) mean(outer(pos, neg, ">")) else NA_real_
    }
  } else if (task == "gaussian") {
    sqrt(mean((as.numeric(truth) - pred)^2))  # RMSE (lower better)
  } else {
    NA_real_
  }
}

# Permute within folds (optionally stratified by outcome) and recompute global metric
.permute_once <- function(fit, task, stratified = TRUE) {
  preds <- fit@predictions
  permuted <- vector("list", length(preds))
  for (i in seq_along(preds)) {
    df <- preds[[i]]
    if (stratified && task == "binomial") {
      # permute labels within class inside fold (keeps class balance)
      if (is.factor(df$truth)) {
        levs <- levels(df$truth)
        for (lv in levs) {
          idx <- which(df$truth == lv)
          if (length(idx) > 1) df$truth[idx] <- sample(df$truth[idx])
        }
      } else {
        for (lv in c(0, 1)) {
          idx <- which(df$truth == lv)
          if (length(idx) > 1) df$truth[idx] <- sample(df$truth[idx])
        }
      }
    } else {
      df$truth <- sample(df$truth)
    }
    permuted[[i]] <- df
  }
  ap <- do.call(rbind, permuted)
  .main_metric(task, ap$truth, ap$pred)
}

#' Audit leakage and confounding
#'
#' @param fit LeakFit
#' @param n_perm integer, number of permutations for permutation test (default 200)
#' @param perm_stratify logical, stratify permutation by class within fold (binomial)
#' @param parallel logical, use future.apply for permutations
#' @param seed integer random seed
#' @param batch_cols character vector of metadata columns to test association with folds
#' @param coldata optional data.frame of sample-level metadata (rows align to original samples)
#' @param X_ref optional numeric matrix/data.frame (samples x features) used for duplicate detection
#' @param feature_space "raw" or "rank" (rank-normalize rows before similarity)
#' @param sim_method "cosine" or "pearson"
#' @param sim_threshold numeric in (0,1), similarity threshold for duplicates (default 0.995)
#' @param nn_k integer, if RANN available and n large, #nearest neighbors per row to check (default 50)
#' @param max_pairs cap on the number of duplicate pairs returned (default 5000)
#' @return LeakAudit
#' @export
audit_leakage <- function(fit,
                          n_perm = 200,
                          perm_stratify = TRUE,
                          parallel = FALSE,
                          seed = 1,
                          batch_cols = NULL,
                          coldata = NULL,
                          X_ref = NULL,
                          feature_space = c("raw", "rank"),
                          sim_method = c("cosine", "pearson"),
                          sim_threshold = 0.995,
                          nn_k = 50,
                          max_pairs = 5000) {

  feature_space <- match.arg(feature_space)
  sim_method    <- match.arg(sim_method)

  set.seed(seed)

  # Trail / provenance
  trail <- list(
    indices_hash = .bio_hash_indices(fit@splits@indices),
    mode = fit@splits@mode,
    info = fit@splits@info,
    seed = seed
  )

  # --- Reconstruct main metric from CV predictions --------------------------
  all_pred <- do.call(rbind, fit@predictions)
  task <- fit@task
  metric_obs <- .main_metric(task, all_pred$truth, all_pred$pred)

  # --- Permutation test (parallelizable) ------------------------------------
  perm_fun <- function(b) .permute_once(fit, task, stratified = perm_stratify)
  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    perm_vals <- future.apply::future_sapply(seq_len(n_perm), function(i) {
      set.seed(seed + i); perm_fun(i)
    }, future.seed = TRUE)
  } else {
    perm_vals <- sapply(seq_len(n_perm), function(i) { set.seed(seed + i); perm_fun(i) })
  }

  # Direction & p-value (exact, with +1 correction)
  if (task == "gaussian") {
    # lower is better
    pval <- (1 + sum(perm_vals <= metric_obs, na.rm = TRUE)) / (1 + sum(is.finite(perm_vals)))
    gap  <- mean(perm_vals, na.rm = TRUE) - metric_obs
    z    <- (mean(perm_vals, na.rm = TRUE) - metric_obs) / (sd(perm_vals, na.rm = TRUE) + 1e-12)
  } else {
    # higher is better
    pval <- (1 + sum(perm_vals >= metric_obs, na.rm = TRUE)) / (1 + sum(is.finite(perm_vals)))
    gap  <- metric_obs - mean(perm_vals, na.rm = TRUE)
    z    <- (metric_obs - mean(perm_vals, na.rm = TRUE)) / (sd(perm_vals, na.rm = TRUE) + 1e-12)
  }

  perm_df <- data.frame(
    metric_obs = metric_obs,
    perm_mean  = mean(perm_vals, na.rm = TRUE),
    perm_sd    = sd(perm_vals, na.rm = TRUE),
    gap        = gap,
    z          = z,
    p_value    = pval,
    n_perm     = length(perm_vals)
  )

  # --- Batch / study association with folds ---------------------------------
  # Need: a per-sample fold assignment (first time each sample appears in test)
  n_samples <- length(unique(do.call(c, lapply(fit@predictions, function(d) d$id))))
  fold_id <- rep(NA_integer_, n_samples)
  for (i in seq_along(fit@splits@indices)) {
    te <- fit@splits@indices[[i]]$test
    fold_id[te] <- ifelse(is.na(fold_id[te]), i, fold_id[te])
  }

  # Get metadata: prefer provided 'coldata', else try to retrieve from splits@info
  if (is.null(coldata) && !is.null(fit@splits@info$coldata)) {
    coldata <- fit@splits@info$coldata
  }

  batch_results <- NULL
  if (!is.null(coldata) && length(fold_id) == nrow(coldata)) {
    # If user didn't pass batch_cols, try plausible keys
    if (is.null(batch_cols)) {
      batch_cols <- intersect(c("batch", "plate", "center", "site", "study"), colnames(coldata))
    }
    if (length(batch_cols) > 0) {
      batch_results <- lapply(batch_cols, function(bc) {
        if (!bc %in% colnames(coldata)) return(NULL)
        tab <- table(fold = fold_id[!is.na(fold_id)],
                     batch = as.factor(coldata[[bc]][!is.na(fold_id)]))
        as.data.frame(.chisq_assoc(tab))
      })
      names(batch_results) <- batch_cols
      batch_df <- do.call(rbind, Map(function(nm, df) { df$batch_col <- nm; df }, names(batch_results), batch_results))
      batch_df <- batch_df[, c("batch_col", "stat", "df", "pval", "cramer_v")]
    } else {
      batch_df <- data.frame()
    }
  } else {
    batch_df <- data.frame()
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
    # standardize to enable cosine via Euclidean NN
    Xn <- .row_l2_normalize(X)

    n <- nrow(Xn)
    candidate_pairs <- NULL

    if (n > 3000 && requireNamespace("RANN", quietly = TRUE)) {
      # approximate k-NN search; cosine ~ 1 - 0.5*||u - v||^2 for unit vectors
      nn <- RANN::nn2(Xn, k = min(nn_k + 1, n))  # includes self
      idx <- rep(seq_len(n), times = ncol(nn$nn.idx))
      jdx <- as.vector(nn$nn.idx)
      idx <- idx[jdx != idx]   # drop self
      jdx <- jdx[jdx != idx]
      # compute cosine for candidate pairs
      cosv <- rowSums(Xn[idx, , drop = FALSE] * Xn[jdx, , drop = FALSE])
      keep <- which(cosv >= sim_threshold)
      if (length(keep)) {
        candidate_pairs <- data.frame(i = idx[keep], j = jdx[keep], cos_sim = cosv[keep])
      }
    } else {
      # exact; cap output
      S <- .cosine_sim_block(Xn)
      S[upper.tri(S, TRUE)] <- 0
      which_dup <- which(S >= sim_threshold, arr.ind = TRUE)
      if (nrow(which_dup) > 0)
        candidate_pairs <- data.frame(i = which_dup[, 1], j = which_dup[, 2], cos_sim = S[which_dup])
    }

    if (!is.null(candidate_pairs)) {
      # sort and cap
      ord <- order(candidate_pairs$cos_sim, decreasing = TRUE)
      candidate_pairs <- candidate_pairs[ord, , drop = FALSE]
      if (nrow(candidate_pairs) > max_pairs)
        candidate_pairs <- candidate_pairs[seq_len(max_pairs), , drop = FALSE]
      dup_df <- candidate_pairs
    }
  }

  # --- Assemble S4 object ----------------------------------------------------
  new("LeakAudit",
      fit = fit,
      permutation_gap = perm_df,
      perm_distribution = perm_vals,
      batch_assoc = batch_df,
      duplicates = dup_df,
      trail = trail,
      info = list(
        n_perm = n_perm,
        perm_stratify = perm_stratify,
        parallel = parallel,
        sim_method = sim_method,
        feature_space = feature_space,
        duplicate_threshold = sim_threshold,
        nn_k = nn_k,
        max_pairs = max_pairs,
        batch_cols = batch_cols
      ))
}
