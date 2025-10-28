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

# Compute metric from predictions
.metric_value <- function(metric, task, truth, pred) {
  if (metric == "rmse") {
    return(sqrt(mean((as.numeric(truth) - pred)^2)))
  }
  if (metric == "auc") {
    if (requireNamespace("pROC", quietly = TRUE)) {
      roc <- pROC::roc(truth, pred, quiet = TRUE)
      return(as.numeric(pROC::auc(roc)))
    }
    yb <- if (is.factor(truth)) as.numeric(truth) - 1 else truth
    pos <- pred[yb == 1]; neg <- pred[yb == 0]
    if (length(pos) && length(neg)) return(mean(outer(pos, neg, ">")))
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
  if (metric == "cindex" && requireNamespace("survival", quietly = TRUE)) {
    return(survival::concordance.index(pred, truth)$c.index)
  }
  NA_real_
}

#' Audit leakage and confounding
#'
#' @param fit LeakFit
#' @param metric performance metric ("auc","pr_auc","rmse","cindex")
#' @param B integer, number of permutations
#' @param perm_stratify logical, stratify permutation by class within fold (binomial)
#' @param time_block block resampling method for time-series modes
#' @param block_len integer block length (NULL for automatic)
#' @param include_z logical, whether to compute IF-based z-score for Δ
#' @param ci_method "if" or "bootstrap" for Δ confidence interval
#' @param boot_B bootstrap resamples when ci_method = "bootstrap"
#' @param parallel logical, use future.apply for permutations
#' @param seed integer random seed
#' @param return_perm logical, store permutation distribution
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
                          metric = c("auc", "pr_auc", "rmse", "cindex"),
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
                          feature_space = c("raw", "rank"),
                          sim_method = c("cosine", "pearson"),
                          sim_threshold = 0.995,
                          nn_k = 50,
                          max_pairs = 5000) {

  metric <- match.arg(metric)
  feature_space <- match.arg(feature_space)
  sim_method    <- match.arg(sim_method)
  time_block <- match.arg(time_block)
  ci_method <- match.arg(ci_method)

  set.seed(seed)

  # Trail / provenance
  trail <- list(
    indices_hash = .bio_hash_indices(fit@splits@indices),
    mode = fit@splits@mode,
    info = fit@splits@info,
    seed = seed
  )

  # --- Reconstruct main metric from CV predictions --------------------------
  task <- fit@task
  folds <- fit@splits@indices
  if (is.null(coldata) && !is.null(fit@splits@info$coldata)) {
    coldata <- fit@splits@info$coldata
  }
  outcome_col <- fit@splits@info$outcome
  if (is.null(outcome_col)) outcome_col <- fit@outcome

  pred_list <- lapply(fit@predictions, function(df) data.frame(df, stringsAsFactors = FALSE))
  all_pred <- do.call(rbind, pred_list)

  metric_obs <- .metric_value(metric, task, all_pred$truth, all_pred$pred)

  # Weights per fold for IF aggregation
  weights <- vapply(folds, function(f) length(f$test), integer(1))
  weights <- weights / sum(weights)

  se_obs <- NA_real_
  if (metric == "auc") {
    obs_truth <- lapply(pred_list, `[[`, "truth")
    obs_pred <- lapply(pred_list, `[[`, "pred")
    if_stat <- .cvauc_if(pred = obs_pred, truth = obs_truth, weights = weights)
    se_obs <- if_stat$se_auc
  }

  higher_better <- metric %in% c("auc", "pr_auc", "cindex")

  # --- Permutations ----------------------------------------------------------
  perm_source <- NULL
  if (!is.null(coldata) && !is.null(outcome_col) && outcome_col %in% names(coldata)) {
    perm_source <- .permute_labels_factory(
      cd = coldata, outcome = outcome_col, mode = fit@splits@mode,
      folds = folds, perm_stratify = perm_stratify, time_block = time_block,
      block_len = block_len, seed = seed,
      group_col = fit@splits@info$group, batch_col = fit@splits@info$batch,
      study_col = fit@splits@info$study
    )
  }
  if (is.null(perm_source)) {
    perm_source <- function(b) {
      lapply(pred_list, function(df) sample(df$truth))
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
    new_preds <- Map(function(df, tr) {
      df$truth <- .coerce_truth_like(df$truth, tr)
      df
    }, pred_list, truths)
    agg <- do.call(rbind, new_preds)
    .metric_value(metric, task, agg$truth, agg$pred)
  }

  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    perm_vals <- future.apply::future_sapply(seq_len(B), function(i) {
      set.seed(seed + i); perm_eval(i)
    }, future.seed = TRUE)
  } else {
    perm_vals <- sapply(seq_len(B), function(i) { set.seed(seed + i); perm_eval(i) })
  }

  perm_mean <- mean(perm_vals, na.rm = TRUE)
  delta <- if (higher_better) metric_obs - perm_mean else perm_mean - metric_obs
  pval <- if (higher_better) {
    (1 + sum(perm_vals >= metric_obs, na.rm = TRUE)) / (1 + sum(is.finite(perm_vals)))
  } else {
    (1 + sum(perm_vals <= metric_obs, na.rm = TRUE)) / (1 + sum(is.finite(perm_vals)))
  }
  perm_sd <- stats::sd(perm_vals, na.rm = TRUE)
  p_se <- sqrt(pval * (1 - pval) / (B + 1))

  seci <- if (ci_method == "if") {
    .se_ci_delta(delta, se_obs, perm_vals)
  } else {
    boot_vals <- replicate(boot_B, {
      sample_vals <- sample(perm_vals, replace = TRUE)
      if (higher_better) metric_obs - mean(sample_vals) else mean(sample_vals) - metric_obs
    })
    alpha <- (1 - 0.95) / 2
    list(
      se = stats::sd(boot_vals),
      ci = stats::quantile(boot_vals, probs = c(alpha, 1 - alpha)),
      z = ifelse(stats::sd(boot_vals) > 0, delta / stats::sd(boot_vals), NA_real_)
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

  # --- Batch / study association with folds ---------------------------------
  # Need: a per-sample fold assignment (first time each sample appears in test)
  n_samples <- length(unique(do.call(c, lapply(fit@predictions, function(d) d$id))))
  # fold_id <- rep(NA_integer_, n_samples)
  # for (i in seq_along(fit@splits@indices)) {
  #   te <- fit@splits@indices[[i]]$test
  #   fold_id[te] <- ifelse(is.na(fold_id[te]), i, fold_id[te])
  # }

  # Extract all unique sample IDs appearing in predictions
  ids_all <- sort(unique(do.call(c, lapply(fit@predictions, `[[`, "id"))))
  fold_id <- rep(NA_integer_, length(ids_all))
  names(fold_id) <- ids_all

  # Assign fold numbers to each test ID
  for (i in seq_along(fit@splits@indices)) {
    te_ids <- fit@splits@indices[[i]]$test
    te_ids <- intersect(te_ids, ids_all)
    fold_id[as.character(te_ids)] <- i
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
  perm_distribution <- if (isTRUE(return_perm)) perm_vals else numeric(0)

  new("LeakAudit",
      fit = fit,
      permutation_gap = perm_df,
      perm_distribution = perm_distribution,
      batch_assoc = batch_df,
      duplicates = dup_df,
      trail = trail,
      info = list(
        n_perm = B,
        perm_stratify = perm_stratify,
        parallel = parallel,
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
.coerce_truth_like <- function(template, values) {
  if (is.factor(template)) {
    factor(values, levels = levels(template))
  } else if (is.logical(template)) {
    as.logical(values)
  } else if (is.numeric(template)) {
    as.numeric(values)
  } else {
    values
  }
}

