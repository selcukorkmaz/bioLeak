# Delta Leakage Sensitivity Index (ΔLSI) ────────────────────────────────────
# Compares a naive (potentially leaky) CV pipeline against a guarded pipeline
# and quantifies leakage-induced performance inflation.

# ── Internal helpers ──────────────────────────────────────────────────────────

# Resolve which learner to use from a LeakFit's @metrics.
.dlsi_resolve_learner <- function(fit, learner) {
  if (!is.null(learner)) return(as.character(learner[1L]))
  if ("learner" %in% names(fit@metrics)) {
    lrns <- unique(as.character(fit@metrics$learner))
    lrns <- lrns[nzchar(lrns)]
    if (!length(lrns)) return(NULL)
    if (length(lrns) > 1L)
      warning(sprintf(
        "[delta_lsi] Multiple learners in fit@metrics (%s); using '%s'.",
        paste(lrns, collapse = ", "), lrns[1L]
      ))
    return(lrns[1L])
  }
  NULL
}

# Extract per-fold metric values from fit@metrics.
# Returns data.frame(fold, metric).
.dlsi_fold_metrics <- function(fit, metric, learner) {
  mdf <- fit@metrics

  if (!is.null(learner) && "learner" %in% names(mdf))
    mdf <- mdf[as.character(mdf$learner) == learner, , drop = FALSE]

  if (metric %in% names(mdf) && "fold" %in% names(mdf)) {
    out <- mdf[, c("fold", metric), drop = FALSE]
    names(out)[2L] <- "metric"
    return(out)
  }

  # Fallback: recompute fold metrics from @predictions
  .dlsi_fold_metrics_from_preds(fit, metric, learner)
}

# Recompute fold metrics from @predictions when not in @metrics.
.dlsi_fold_metrics_from_preds <- function(fit, metric, learner) {
  task    <- if (length(fit@task)) fit@task else "binomial"
  pred_df <- .dlsi_combine_preds(fit, learner)
  if (is.null(pred_df) || !"fold" %in% names(pred_df))
    stop(sprintf(
      "[delta_lsi] Metric '%s' not in fit@metrics and no usable predictions.",
      metric
    ))
  folds <- sort(unique(pred_df$fold))
  vals  <- vapply(folds, function(f) {
    sub <- pred_df[pred_df$fold == f, , drop = FALSE]
    tryCatch(
      .metric_value(metric, task, sub$truth, sub$pred, pred_df = sub),
      error = function(e) NA_real_
    )
  }, numeric(1L))
  data.frame(fold = folds, metric = vals, stringsAsFactors = FALSE)
}

# Combine all per-fold prediction data frames into one.
.dlsi_combine_preds <- function(fit, learner = NULL) {
  preds <- fit@predictions
  if (!length(preds)) return(NULL)
  pred_df <- tryCatch(
    do.call(rbind, Filter(is.data.frame, preds)),
    error = function(e) NULL
  )
  if (is.null(pred_df) || !nrow(pred_df)) return(NULL)
  if (!is.null(learner) && "learner" %in% names(pred_df)) {
    pred_df <- pred_df[as.character(pred_df$learner) == learner, , drop = FALSE]
  } else if ("learner" %in% names(pred_df)) {
    first_lrn <- unique(as.character(pred_df$learner))[1L]
    pred_df   <- pred_df[as.character(pred_df$learner) == first_lrn, , drop = FALSE]
  }
  pred_df
}

# Get repeat_id for each sequential fold position (1..N).
.dlsi_repeat_ids <- function(fit) {
  vapply(fit@splits@indices, function(z) {
    as.integer(z$repeat_id %||% 1L)
  }, integer(1L))
}

# Get test-set size for each sequential fold position.
# Falls back to predictions-based count when splits use compact form.
.dlsi_fold_sizes <- function(fit, pred_df = NULL) {
  sizes <- vapply(fit@splits@indices, function(z) {
    if (!is.null(z$test)) as.integer(length(z$test)) else NA_integer_
  }, integer(1L))

  if (any(is.na(sizes)) && !is.null(pred_df) && "fold" %in% names(pred_df)) {
    tab      <- table(pred_df$fold)
    fold_seq <- seq_along(fit@splits@indices)
    for (i in which(is.na(sizes))) {
      k <- as.character(fold_seq[i])
      if (k %in% names(tab)) sizes[i] <- as.integer(tab[[k]])
    }
  }
  sizes[is.na(sizes)] <- 1L
  sizes
}

# Check whether two fits share exactly the same fold-level test indices.
# Returns TRUE only when every fold's test set is identical (same observations,
# same order after sorting).  This guards against pairing fits that happen to
# have the same repeat count but were run on different data splits.
.dlsi_splits_match <- function(fit_naive, fit_guarded) {
  idx_n <- fit_naive@splits@indices
  idx_g <- fit_guarded@splits@indices
  if (length(idx_n) != length(idx_g)) return(FALSE)
  all(vapply(seq_along(idx_n), function(i) {
    te_n <- sort(as.integer(idx_n[[i]]$test %||% integer(0)))
    te_g <- sort(as.integer(idx_g[[i]]$test %||% integer(0)))
    identical(te_n, te_g)
  }, logical(1L)))
}

# Huber M-estimator of location (iterative; Huber k = 1.345).
.dlsi_huber <- function(x, k = 1.345, tol = 1e-8, max_iter = 100L) {
  x <- as.numeric(x[is.finite(x)])
  n <- length(x)
  if (!n)    return(NA_real_)
  if (n == 1L) return(x[1L])
  mu    <- stats::median(x)
  sigma <- stats::mad(x, constant = 1.4826)
  if (!is.finite(sigma) || sigma < tol) return(mu)
  for (iter in seq_len(max_iter)) {
    r     <- (x - mu) / sigma
    w     <- ifelse(abs(r) <= k, 1.0, k / abs(r))
    mu_new <- sum(w * x) / sum(w)
    if (abs(mu_new - mu) < tol * (1.0 + abs(mu))) break
    mu <- mu_new
  }
  mu
}

# BCa bootstrap confidence interval (Efron, 1987).
.dlsi_bca_ci <- function(x, statfn, M = 2000L, conf = 0.95, seed = 42L) {
  x <- as.numeric(x[is.finite(x)])
  n <- length(x)
  if (n < 2L) return(c(NA_real_, NA_real_))

  set.seed(seed)
  theta_hat <- statfn(x)
  if (!is.finite(theta_hat)) return(c(NA_real_, NA_real_))

  boot_stats <- replicate(M, statfn(sample(x, n, replace = TRUE)))
  boot_stats <- as.numeric(boot_stats[is.finite(boot_stats)])
  M_eff      <- length(boot_stats)
  if (M_eff < 10L) return(c(NA_real_, NA_real_))

  # Bias correction z0
  prop <- mean(boot_stats < theta_hat, na.rm = TRUE)
  prop <- max(1.0 / M_eff, min(1.0 - 1.0 / M_eff, prop))
  z0   <- stats::qnorm(prop)

  # Acceleration a_hat via jackknife
  jack   <- vapply(seq_len(n),
                   function(i) tryCatch(statfn(x[-i]), error = function(e) NA_real_),
                   numeric(1L))
  jack   <- as.numeric(jack[is.finite(jack)])
  a_hat  <- if (length(jack) >= 2L) {
    j_bar <- mean(jack)
    diff  <- j_bar - jack
    denom <- 6.0 * sum(diff^2)^1.5
    if (is.finite(denom) && denom > 1e-12) sum(diff^3) / denom else 0.0
  } else 0.0

  # BCa adjusted quantiles
  alpha    <- 1.0 - conf
  adj_p    <- function(z_a) {
    denom <- 1.0 - a_hat * (z0 + z_a)
    if (abs(denom) < 1e-10) return(alpha / 2.0)
    stats::pnorm(z0 + (z0 + z_a) / denom)
  }
  p_lo <- adj_p(stats::qnorm(alpha / 2.0))
  p_hi <- adj_p(stats::qnorm(1.0 - alpha / 2.0))
  eps  <- 1.0 / M_eff
  p_lo <- max(eps, min(1.0 - eps, p_lo))
  p_hi <- max(eps, min(1.0 - eps, p_hi))

  unname(stats::quantile(boot_stats, c(p_lo, p_hi), na.rm = TRUE))
}

# Sign-flip randomization test (Phipson & Smyth, 2010).
# H0: mean(delta_r) = 0 (no leakage inflation).
# Test statistic: arithmetic mean (well-defined under sign-flip null).
#
# Exact enumeration (R <= 15): all 2^R sign combinations are evaluated.
#   p = count(|null| >= |obs|) / 2^R  — no continuity correction.
#
# Monte Carlo (R > 15): Phipson & Smyth (2010) correction:
#   p = (count + 1) / (M + 1)  — guards against p = 0 from sampling.
.dlsi_sign_flip <- function(delta_r, M_flip = 10000L, seed = 42L) {
  delta_r  <- as.numeric(delta_r[is.finite(delta_r)])
  R        <- length(delta_r)
  if (R < 2L) return(NA_real_)

  obs_stat <- mean(delta_r)
  if (!is.finite(obs_stat)) return(NA_real_)

  set.seed(seed)

  if (R <= 15L) {
    # ── Exact enumeration ──────────────────────────────────────────────────
    n_combos   <- 2L^R
    null_stats <- vapply(seq_len(n_combos), function(j) {
      bits  <- as.integer(intToBits(j - 1L)[seq_len(R)])
      signs <- ifelse(bits == 1L, 1.0, -1.0)
      mean(signs * delta_r)
    }, numeric(1L))
    null_stats <- as.numeric(null_stats[is.finite(null_stats)])
    if (!length(null_stats)) return(NA_real_)
    # Exact p-value: proportion of null stats at least as extreme as observed.
    # No continuity correction — the null distribution is fully enumerated.
    sum(abs(null_stats) >= abs(obs_stat) - .Machine$double.eps * 100) /
      length(null_stats)
  } else {
    # ── Monte Carlo sign-flip ──────────────────────────────────────────────
    null_stats <- replicate(M_flip, {
      signs <- sample(c(-1.0, 1.0), R, replace = TRUE)
      mean(signs * delta_r)
    })
    null_stats <- as.numeric(null_stats[is.finite(null_stats)])
    if (!length(null_stats)) return(NA_real_)
    # Phipson & Smyth (2010) +1 correction prevents p = 0 from Monte Carlo.
    (sum(abs(null_stats) >= abs(obs_stat) - .Machine$double.eps * 100) + 1L) /
      (length(null_stats) + 1L)
  }
}

# Auto-estimate block size from AR(1) autocorrelation of Δ_r for blocked_time
# sign-flip.  With small R (5-20 values) the AR(1) estimate is noisy; the result
# is capped at floor(R/3) to guarantee at least three independent blocks.
.dlsi_resolve_block_size <- function(delta_r, block_size = NULL) {
  R <- length(delta_r)
  if (!is.null(block_size)) return(max(1L, as.integer(block_size)))
  if (R < 4L) return(1L)
  rho1 <- tryCatch(
    stats::cor(delta_r[-R], delta_r[-1L]),
    error = function(e) 0.0
  )
  rho1 <- if (is.finite(rho1)) max(0.0, rho1) else 0.0
  # Effective block size from AR(1) rate: b ≈ 1/(1-rho1)
  bs <- if (rho1 >= 0.99) R else max(1L, ceiling(1.0 / (1.0 - rho1)))
  # Cap so there are at least 3 blocks; floor(R/3) is the upper limit
  as.integer(min(bs, max(1L, floor(R / 3L))))
}

# Block sign-flip randomization test for temporally structured Δ_r.
# Flips contiguous blocks of repeats together to preserve serial autocorrelation
# under the null.  Returns a list(p_value, n_blocks, block_size_used).
#
# Exact enumeration (n_blocks <= 15): all 2^n_blocks block-sign vectors.
#   p = count(|null| >= |obs|) / 2^n_blocks  — no continuity correction.
#
# Monte Carlo (n_blocks > 15): Phipson & Smyth (2010) correction:
#   p = (count + 1) / (M + 1).
.dlsi_sign_flip_blocked <- function(delta_r, block_size, M_flip = 10000L,
                                    seed = 42L) {
  delta_r    <- as.numeric(delta_r[is.finite(delta_r)])
  R          <- length(delta_r)
  block_size <- max(1L, as.integer(block_size))

  if (R < 2L) return(list(p_value = NA_real_, n_blocks = 0L,
                           block_size_used = block_size))

  block_id <- ceiling(seq_len(R) / block_size)
  n_blocks <- max(block_id)

  # Single block: every sign vector flips all Δ_r simultaneously; p is always 1.
  if (n_blocks < 2L)
    return(list(p_value = 1.0, n_blocks = n_blocks,
                block_size_used = block_size))

  obs_stat <- mean(delta_r)
  if (!is.finite(obs_stat))
    return(list(p_value = NA_real_, n_blocks = n_blocks,
                block_size_used = block_size))

  set.seed(seed)

  if (n_blocks <= 15L) {
    # ── Exact enumeration over 2^n_blocks block-sign vectors ─────────────────
    combos     <- as.matrix(expand.grid(rep(list(c(-1.0, 1.0)), n_blocks)))
    null_stats <- apply(combos, 1L, function(bsigns) {
      mean(bsigns[block_id] * delta_r)
    })
    null_stats <- as.numeric(null_stats[is.finite(null_stats)])
    if (!length(null_stats))
      return(list(p_value = NA_real_, n_blocks = n_blocks,
                  block_size_used = block_size))
    p_val <- sum(abs(null_stats) >= abs(obs_stat) - .Machine$double.eps * 100) /
      length(null_stats)
  } else {
    # ── Monte Carlo with Phipson & Smyth (2010) correction ───────────────────
    null_stats <- replicate(M_flip, {
      bsigns <- sample(c(-1.0, 1.0), n_blocks, replace = TRUE)
      mean(bsigns[block_id] * delta_r)
    })
    null_stats <- as.numeric(null_stats[is.finite(null_stats)])
    if (!length(null_stats))
      return(list(p_value = NA_real_, n_blocks = n_blocks,
                  block_size_used = block_size))
    p_val <- (sum(abs(null_stats) >= abs(obs_stat) - .Machine$double.eps * 100) +
                1L) / (length(null_stats) + 1L)
  }

  list(p_value = p_val, n_blocks = n_blocks, block_size_used = block_size)
}

# ── Main function ─────────────────────────────────────────────────────────────

#' Delta Leakage Sensitivity Index (Delta LSI)
#'
#' Compares a naive (potentially leaky) cross-validation pipeline against a
#' guarded (leakage-protected) pipeline and quantifies leakage-induced
#' performance inflation using the Leakage Sensitivity Index (LSI).
#'
#' @details
#' \subsection{Method}{
#' For each fit, per-fold metric values are extracted from \code{fit@@metrics}
#' (or recomputed from \code{fit@@predictions} if necessary). Fold test-set
#' sizes are used as weights to aggregate fold metrics into per-repeat
#' estimates \eqn{\mu_r}. The repeat-level delta
#' \eqn{\Delta_r = s \cdot (\mu_r^{\text{naive}} - \mu_r^{\text{guarded}})}
#' captures leakage-induced performance inflation for each CV repeat, where
#' \eqn{s = +1} for higher-is-better metrics (e.g., AUC) and \eqn{s = -1}
#' for lower-is-better metrics (e.g., RMSE), so that \eqn{\Delta_r > 0}
#' always indicates the naive pipeline is more optimistic than the guarded one.
#'
#' The \strong{delta_lsi} point estimate is the Huber M-estimator (k = 1.345)
#' applied to \eqn{\{\Delta_r\}}, which is robust to occasional outlier
#' repeats. \strong{delta_metric} is the arithmetic mean of \eqn{\{\Delta_r\}}
#' for easy interpretation in the original metric's units.
#'
#' Pairing requires that \code{fit_leaky} and \code{fit_guarded} share
#' \emph{identical fold structures} (same test-set membership per fold) in
#' addition to the same number of repeats.  When repeat counts match but fold
#' structures differ, a warning is issued and the fits are treated as unpaired.
#'
#' When \eqn{R_{\text{eff}} \geq 5} (equal, paired repeats), a sign-flip
#' randomization test (Phipson & Smyth, 2010) is performed: under
#' \eqn{H_0} (no leakage) the sign of each \eqn{\Delta_r} is exchangeable.
#' All \eqn{2^R} sign combinations are enumerated exactly for
#' \eqn{R \leq 15} (no continuity correction); Monte Carlo sampling is used
#' for larger \eqn{R} with the Phipson & Smyth (2010) correction.
#'
#' BCa bootstrap confidence intervals (Efron, 1987) require
#' \eqn{R_{\text{eff}} \geq 10}.
#' }
#'
#' \subsection{Inference tiers}{
#' \describe{
#'   \item{\code{"A_full_inference"}}{R_eff >= 20: point + BCa CI + sign-flip p-value; \code{inference_ok = TRUE}}
#'   \item{\code{"B_signflip_ci"}}{10 <= R_eff < 20: point + sign-flip p-value + BCa CI}
#'   \item{\code{"C_signflip"}}{5 <= R_eff < 10: point + sign-flip p-value (no CI)}
#'   \item{\code{"D_insufficient"}}{R_eff < 5 or unpaired: point estimate only}
#' }}
#'
#' @param fit_leaky A \code{\linkS4class{LeakFit}} object from the leaky
#'   (unprotected) evaluation pipeline.
#' @param fit_guarded A \code{\linkS4class{LeakFit}} object from the guarded
#'   (leakage-protected) evaluation pipeline.
#' @param metric Character. Performance metric to compare. Must appear in
#'   \code{fit@@metrics} of both fits (e.g., \code{"auc"}, \code{"rmse"}).
#' @param exchangeability Character. Exchangeability assumption for the
#'   sign-flip test. One of \code{"iid"} (default), \code{"by_group"},
#'   \code{"within_batch"}, \code{"blocked_time"}.
#'   \code{"blocked_time"} activates a block sign-flip procedure that flips
#'   contiguous groups of repeats together, preserving serial autocorrelation
#'   under the null; see \code{block_size}.  \code{"by_group"} and
#'   \code{"within_batch"} are stored and reported but inference still uses the
#'   iid sign-flip procedure (a warning is issued; contributions welcome).
#'   \code{"iid"} (default) applies the standard independent sign-flip test.
#' @param block_size Integer or \code{NULL}. Block length for the block
#'   sign-flip test, used only when \code{exchangeability = "blocked_time"}.
#'   When \code{NULL} (default), the block size is auto-estimated from the
#'   first-order autocorrelation of \eqn{\{\Delta_r\}} and capped at
#'   \code{floor(R/3)} to ensure at least three independent blocks.  A warning
#'   is issued when the estimate is used with \code{R_eff < 20} because the
#'   AR(1) estimate is noisy at small sample sizes.  Provide an explicit integer
#'   to override auto-estimation.
#' @param learner Optional character. Learner name to select from multi-learner
#'   fits. If \code{NULL}, the first learner found in \code{fit@@metrics} is used.
#' @param higher_is_better Logical or \code{NULL}. Whether a higher value of
#'   \code{metric} indicates better performance. When \code{NULL} (default),
#'   auto-detected from the metric name: \code{"rmse"}, \code{"mse"},
#'   \code{"mae"}, \code{"log_loss"}, \code{"brier"}, \code{"error"},
#'   \code{"loss"}, and \code{"deviance"} are treated as lower-is-better;
#'   all others default to higher-is-better. Setting this correctly ensures
#'   that a positive \code{delta_lsi} always indicates leakage inflation
#'   (the naive pipeline is artificially more optimistic than the guarded one).
#' @param M_boot Integer. Number of bootstrap samples for BCa CI (default 2000).
#' @param M_flip Integer. Maximum Monte Carlo samples for sign-flip test when
#'   R_eff > 15 (default 10000).
#' @param strict Logical. If \code{TRUE}, error on insufficient R_eff instead
#'   of a warning.
#' @param return_details Logical. If \code{TRUE}, include the per-repeat
#'   \eqn{\Delta_r} vector and the original fit objects in the \code{info} slot.
#' @param seed Integer. Random seed for bootstrap and sign-flip test.
#' @param ... Unused. Reserved for deprecated aliases such as
#'   \code{fit_naive}.
#'
#' @return A \code{\linkS4class{LeakDeltaLSI}} object.
#'
#' @seealso \code{\link{audit_leakage}}, \code{\link{fit_resample}}
#' @export
delta_lsi <- function(
    fit_leaky,
    fit_guarded,
    metric          = "auc",
    exchangeability = c("iid", "by_group", "within_batch", "blocked_time"),
    learner         = NULL,
    higher_is_better = NULL,
    block_size      = NULL,
    M_boot          = 2000L,
    M_flip          = 10000L,
    strict          = FALSE,
    return_details  = FALSE,
    seed            = 42L,
    ...
) {
  extras <- list(...)
  if ("fit_naive" %in% names(extras)) {
    if (!missing(fit_leaky)) {
      stop("Use either fit_leaky or the deprecated fit_naive argument, not both.")
    }
    fit_leaky <- extras$fit_naive
    extras$fit_naive <- NULL
    warning("delta_lsi(): fit_naive is deprecated; use fit_leaky instead.")
  }
  if (length(extras)) {
    stop(
      sprintf(
        "Unused argument(s): %s",
        paste(names(extras), collapse = ", ")
      )
    )
  }

  exchangeability <- match.arg(exchangeability)

  if (!inherits(fit_leaky,   "LeakFit")) stop("fit_leaky must be a LeakFit object.")
  if (!inherits(fit_guarded, "LeakFit")) stop("fit_guarded must be a LeakFit object.")

  # ── Resolve metric directionality ────────────────────────────────────────────
  if (is.null(higher_is_better)) {
    lower_is_better <- grepl(
      "rmse|mse|mae|log_loss|logloss|brier|error|loss|deviance",
      tolower(metric)
    )
    higher_is_better <- !lower_is_better
  }
  # sign_factor > 0 means "naive - guarded > 0 is leakage"; flip for lower-is-better
  sign_factor <- if (isTRUE(higher_is_better)) 1.0 else -1.0

  # Resolve learners (one per fit; may differ)
  lrn_n <- .dlsi_resolve_learner(fit_leaky,   learner)
  lrn_g <- .dlsi_resolve_learner(fit_guarded, learner)

  # Per-fold metric data frames: columns fold, metric
  fm_n <- .dlsi_fold_metrics(fit_leaky,   metric, lrn_n)
  fm_g <- .dlsi_fold_metrics(fit_guarded, metric, lrn_g)

  # Sequential fold positions (1..N) for each fit
  fold_seq_n <- seq_along(fit_leaky@splits@indices)
  fold_seq_g <- seq_along(fit_guarded@splits@indices)

  # Repeat IDs and fold sizes
  rep_ids_n <- .dlsi_repeat_ids(fit_leaky)
  rep_ids_g <- .dlsi_repeat_ids(fit_guarded)

  pred_n <- .dlsi_combine_preds(fit_leaky,   lrn_n)
  pred_g <- .dlsi_combine_preds(fit_guarded, lrn_g)

  sizes_n <- .dlsi_fold_sizes(fit_leaky,   pred_n)
  sizes_g <- .dlsi_fold_sizes(fit_guarded, pred_g)

  # Build fold-level data frames
  .build_folds <- function(fm, rep_ids, sizes, fold_seq) {
    idx <- match(fm$fold, fold_seq)
    data.frame(
      fold      = fm$fold,
      metric    = fm$metric,
      repeat_id = rep_ids[idx],
      n         = sizes[idx],
      stringsAsFactors = FALSE
    )
  }

  folds_n <- .build_folds(fm_n, rep_ids_n, sizes_n, fold_seq_n)
  folds_g <- .build_folds(fm_g, rep_ids_g, sizes_g, fold_seq_g)
  folds_n$n[is.na(folds_n$n)] <- 1L
  folds_g$n[is.na(folds_g$n)] <- 1L

  # Aggregate folds → per-repeat size-weighted mean metric
  .agg_repeats <- function(folds_df) {
    reps <- sort(unique(folds_df$repeat_id[is.finite(folds_df$metric)]))
    rows <- lapply(reps, function(r) {
      sub <- folds_df[folds_df$repeat_id == r & is.finite(folds_df$metric), ]
      if (!nrow(sub)) return(NULL)
      w  <- pmax(as.numeric(sub$n), 1.0)
      mu <- sum(w * sub$metric) / sum(w)
      data.frame(
        repeat_id = r, metric = mu,
        n_folds = nrow(sub), total_n = sum(sub$n),
        stringsAsFactors = FALSE
      )
    })
    out <- do.call(rbind, Filter(Negate(is.null), rows))
    if (is.null(out)) {
      out <- data.frame(repeat_id = integer(0), metric = numeric(0),
                        n_folds = integer(0), total_n = integer(0),
                        stringsAsFactors = FALSE)
    }
    out
  }

  reps_n <- .agg_repeats(folds_n)
  reps_g <- .agg_repeats(folds_g)

  R_naive   <- nrow(reps_n)
  R_guarded <- nrow(reps_g)

  # ── Pairing: equal count AND matching fold structures ─────────────────────
  same_count <- (R_naive > 0L && R_guarded > 0L && R_naive == R_guarded)
  if (same_count) {
    splits_ok <- .dlsi_splits_match(fit_leaky, fit_guarded)
    if (!splits_ok) {
      warning(sprintf(
        paste0("[delta_lsi] fit_leaky and fit_guarded have R=%d repeats each ",
               "but different fold structures; treating as unpaired."),
        R_naive
      ))
    }
    paired <- splits_ok
    if (paired) {
      # Normalize guarded repeat_ids to match naive's positional mapping.
      # .dlsi_splits_match() verified that sequential position i has identical
      # test members in both fits, so naive's repeat_id labels are canonical.
      rid_lookup <- setNames(folds_n$repeat_id, folds_n$fold)
      folds_g$repeat_id <- unname(rid_lookup[as.character(folds_g$fold)])
      reps_g <- .agg_repeats(folds_g)
    }
  } else {
    paired <- FALSE
  }
  R_eff <- if (paired) as.integer(R_naive) else 0L

  # Per-repeat Δ values — sign_factor applied so positive = leakage inflation.
  # Align by repeat_id (both are sorted by .agg_repeats; guard against rare
  # case where different repeats yield all-NA metrics in each fit).
  delta_r <- if (paired) {
    if (identical(reps_n$repeat_id, reps_g$repeat_id)) {
      sign_factor * (reps_n$metric - reps_g$metric)
    } else {
      common <- intersect(reps_n$repeat_id, reps_g$repeat_id)
      rn <- reps_n[match(common, reps_n$repeat_id), ]
      rg <- reps_g[match(common, reps_g$repeat_id), ]
      sign_factor * (rn$metric - rg$metric)
    }
  } else NULL

  # Update R_eff if intersection reduced pairing
  if (!is.null(delta_r) && length(delta_r) < R_eff) {
    warning(sprintf(
      "[delta_lsi] %d of %d repeats dropped (all-NA metrics); R_eff reduced from %d to %d.",
      R_eff - length(delta_r), R_eff, R_eff, length(delta_r)
    ))
    R_eff <- length(delta_r)
  }

  # Inference tier (computed after R_eff adjustment)
  tier <- if (R_eff >= 20L) "A_full_inference" else
          if (R_eff >= 10L) "B_signflip_ci"    else
          if (R_eff >=  5L) "C_signflip"       else
                            "D_insufficient"

  if (tier == "D_insufficient") {
    msg <- sprintf(
      paste0("[delta_lsi] R_eff = %d (R_leaky=%d, R_guarded=%d, paired=%s). ",
             "Need paired R_eff >= 5 for sign-flip inference; reporting point estimate only."),
      R_eff, R_naive, R_guarded, paired
    )
    if (isTRUE(strict)) stop(msg) else warning(msg)
  }

  # Point estimates
  delta_lsi_est    <- if (!is.null(delta_r)) {
    .dlsi_huber(delta_r)
  } else {
    sign_factor * (.dlsi_huber(reps_n$metric) - .dlsi_huber(reps_g$metric))
  }
  delta_metric_est <- if (!is.null(delta_r)) {
    mean(delta_r, na.rm = TRUE)
  } else {
    sign_factor * (mean(reps_n$metric, na.rm = TRUE) - mean(reps_g$metric, na.rm = TRUE))
  }

  # BCa CI (requires R_eff >= 10 and paired)
  ci_lsi <- c(NA_real_, NA_real_)
  ci_dm  <- c(NA_real_, NA_real_)
  if (!is.null(delta_r) && R_eff >= 10L) {
    ci_lsi <- tryCatch(
      .dlsi_bca_ci(delta_r, .dlsi_huber,
                   M = M_boot, seed = seed),
      error = function(e) c(NA_real_, NA_real_)
    )
    ci_dm  <- tryCatch(
      .dlsi_bca_ci(delta_r, function(x) mean(x, na.rm = TRUE),
                   M = M_boot, seed = seed + 1L),
      error = function(e) c(NA_real_, NA_real_)
    )
  }

  # Sign-flip test (requires R_eff >= 5 and paired).
  # Routing depends on exchangeability:
  #   "blocked_time"       -> block sign-flip (contiguous blocks of Δ_r)
  #   "by_group" / "within_batch" -> iid sign-flip with an explanatory warning
  #   "iid"                -> standard independent sign-flip
  p_val            <- NA_real_
  block_size_used  <- NA_integer_
  n_blocks_used    <- NA_integer_

  if (!is.null(delta_r) && R_eff >= 5L) {
    if (exchangeability == "blocked_time") {
      bs     <- .dlsi_resolve_block_size(delta_r, block_size)
      if (is.null(block_size) && R_eff < 20L)
        warning(sprintf(
          paste0("[delta_lsi] exchangeability='blocked_time': block_size ",
                 "auto-estimated as %d from AR(1) of \u0394r (R_eff=%d; ",
                 "estimate is noisy at small R). Supply block_size explicitly ",
                 "to override."),
          bs, R_eff
        ))
      sf_out         <- .dlsi_sign_flip_blocked(delta_r, block_size = bs,
                                                M_flip = M_flip,
                                                seed = seed + 2L)
      p_val          <- sf_out$p_value
      block_size_used <- as.integer(sf_out$block_size_used)
      n_blocks_used   <- as.integer(sf_out$n_blocks)
      # If too few independent blocks for p < 0.05, set p to NA and warn.
      if (!is.na(n_blocks_used) && n_blocks_used < 5L) {
        min_p <- if (n_blocks_used >= 1L) 1.0 / 2.0^n_blocks_used else 1.0
        warning(sprintf(
          paste0("[delta_lsi] blocked_time sign-flip: n_blocks = %d < 5; ",
                 "minimum achievable p = %.4f > 0.05. Setting p_value = NA. ",
                 "Increase R_eff or reduce block_size."),
          n_blocks_used, min_p
        ))
        p_val <- NA_real_
      }
    } else {
      if (exchangeability %in% c("by_group", "within_batch"))
        warning(sprintf(
          paste0("[delta_lsi] exchangeability = '%s' is stored but the ",
                 "sign-flip test still assumes iid repeat structure. ",
                 "Use 'blocked_time' for temporal autocorrelation or 'iid' ",
                 "to suppress this warning."),
          exchangeability
        ))
      p_val <- tryCatch(
        .dlsi_sign_flip(delta_r, M_flip = M_flip, seed = seed + 2L),
        error = function(e) NA_real_
      )
    }
  }

  inference_ok <- (R_eff >= 20L && is.finite(p_val) && all(is.finite(ci_lsi)))

  # Metadata
  info_list <- list(
    paired           = paired,
    R_naive          = R_naive,
    R_guarded        = R_guarded,
    higher_is_better = higher_is_better,
    metric_naive     = .dlsi_huber(reps_n$metric),
    metric_guarded   = .dlsi_huber(reps_g$metric),
    block_size_used  = block_size_used,
    n_blocks         = n_blocks_used
  )
  if (return_details) {
    info_list$delta_r     <- delta_r
    info_list$fit_leaky   <- fit_leaky
    info_list$fit_naive   <- fit_leaky
    info_list$fit_guarded <- fit_guarded
  }

  methods::new("LeakDeltaLSI",
    metric          = metric,
    exchangeability = exchangeability,
    tier            = tier,
    strict          = strict,
    R_eff           = R_eff,
    delta_lsi       = delta_lsi_est,
    delta_lsi_ci    = as.numeric(ci_lsi),
    delta_metric    = delta_metric_est,
    delta_metric_ci = as.numeric(ci_dm),
    p_value         = p_val,
    inference_ok    = inference_ok,
    folds_naive     = folds_n,
    folds_guarded   = folds_g,
    repeats_naive   = reps_n,
    repeats_guarded = reps_g,
    info            = info_list
  )
}

# ── S4 show method ─────────────────────────────────────────────────────────────

#' @rdname LeakClasses
#' @param object A \code{\linkS4class{LeakDeltaLSI}} object.
setMethod("show", "LeakDeltaLSI", function(object) {
  hib <- object@info[["higher_is_better"]]
  hib_str <- if (is.null(hib)) "" else
             if (isTRUE(hib)) " [higher=better]" else " [lower=better]"
  cat("A LeakDeltaLSI object\n")
  cat(sprintf("  Metric:         %s%s\n", object@metric, hib_str))
  cat(sprintf("  Tier:           %s  (R_eff = %d)\n", object@tier, object@R_eff))
  fmt <- function(x) if (is.finite(x)) sprintf("%.4f", x) else "NA"
  cat(sprintf("  delta_metric:   %s", fmt(object@delta_metric)))
  if (all(is.finite(object@delta_metric_ci)))
    cat(sprintf("  [%.4f, %.4f] 95%% CI", object@delta_metric_ci[1L], object@delta_metric_ci[2L]))
  cat("\n")
  cat(sprintf("  delta_lsi:      %s", fmt(object@delta_lsi)))
  if (all(is.finite(object@delta_lsi_ci)))
    cat(sprintf("  [%.4f, %.4f] 95%% CI", object@delta_lsi_ci[1L], object@delta_lsi_ci[2L]))
  cat("\n")
  if (is.finite(object@p_value))
    cat(sprintf("  p-value:        %.4f  (sign-flip; tests mean \u0394r)\n",
                object@p_value))
  invisible(object)
})

# ── summary method ─────────────────────────────────────────────────────────────

#' Summarize a LeakDeltaLSI object
#'
#' Prints a human-readable summary of the Delta LSI analysis comparing leaky vs
#' guarded evaluation pipelines.
#'
#' @param object A \code{LeakDeltaLSI} object from \code{\link{delta_lsi}}.
#' @param digits Integer. Number of decimal places to show (default 3).
#' @param ... Unused.
#' @return Invisibly returns \code{object}.
#' @export
summary.LeakDeltaLSI <- function(object, digits = 3L, ...) {
  fmt <- function(x) if (is.finite(x)) formatC(x, digits = digits, format = "f") else "NA"
  hib <- object@info[["higher_is_better"]] %||% TRUE

  cat("\n=====================================\n")
  cat(" bioLeak Delta LSI Summary\n")
  cat("=====================================\n\n")
  cat(sprintf("Metric:            %s\n", object@metric))
  exch <- object@exchangeability
  cat(sprintf("Exchangeability:   %s\n", exch))
  if (identical(exch, "blocked_time")) {
    bs   <- object@info$block_size_used
    nblk <- object@info$n_blocks
    if (!is.null(bs) && !is.na(bs))
      cat(sprintf("  Block sign-flip: block_size = %d, n_blocks = %d\n", bs, nblk))
    else
      cat("  Block sign-flip: not computed (R_eff < 5)\n")
  } else if (exch %in% c("by_group", "within_batch")) {
    cat(sprintf(
      "  [Note: '%s' is stored; inference uses iid sign-flip.]\n", exch
    ))
  }
  cat(sprintf("Inference tier:    %s\n", object@tier))
  cat(sprintf("R_eff:             %d  (R_leaky=%s, R_guarded=%s, paired=%s)\n",
              object@R_eff,
              object@info[["R_naive"]]   %||% "?",
              object@info[["R_guarded"]] %||% "?",
              object@info[["paired"]]    %||% "?"))

  cat("\nLeaky pipeline:    ", fmt(object@info[["metric_naive"]]   %||% NA_real_), "\n", sep = "")
  cat("Guarded pipeline:  ", fmt(object@info[["metric_guarded"]] %||% NA_real_), "\n", sep = "")

  direction <- if (isTRUE(hib)) "leaky - guarded" else "-(leaky - guarded)"
  cat(sprintf("\nPoint estimates (%s; positive = leakage inflation):\n", direction))
  cat(sprintf("  delta_metric:  %s  (raw metric difference)\n", fmt(object@delta_metric)))
  ci_dm <- object@delta_metric_ci
  if (all(is.finite(ci_dm)))
    cat(sprintf("  95%% BCa CI:    [%s, %s]\n", fmt(ci_dm[1L]), fmt(ci_dm[2L])))
  cat(sprintf("  delta_lsi:     %s  (Huber-robust)\n", fmt(object@delta_lsi)))
  ci_lsi <- object@delta_lsi_ci
  if (all(is.finite(ci_lsi)))
    cat(sprintf("  95%% BCa CI:    [%s, %s]\n", fmt(ci_lsi[1L]), fmt(ci_lsi[2L])))

  cat("\nHypothesis test  [H0: no systematic inflation; paired repeat signs exchangeable]:\n")
  if (is.finite(object@p_value)) {
    p   <- object@p_value
    sig <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
    cat(sprintf("  Sign-flip p:   %.4f %s\n", p, sig))
    # Diagnostic: warn when p-value and BCa CI lead to different conclusions.
    # This happens when the arithmetic mean (tested by p) and the Huber estimate
    # (covered by the CI) diverge — a sign of outlier repeats skewing the mean.
    if (all(is.finite(object@delta_lsi_ci))) {
      pval_sig    <- p < 0.05
      ci_exc_zero <- (object@delta_lsi_ci[1L] > 0 || object@delta_lsi_ci[2L] < 0)
      if (pval_sig != ci_exc_zero)
        cat(paste0("  [Note: p-value and BCa CI disagree. The p-value tests delta_metric\n",
                   "   (mean); the CI covers delta_lsi (Huber robust). Outlier repeats\n",
                   "   may be pulling the mean. Prefer delta_lsi for point estimates.]\n"))
    }
  } else {
    reason <- if (identical(object@exchangeability, "blocked_time") &&
                  !is.null(object@info$n_blocks) &&
                  !is.na(object@info$n_blocks) &&
                  object@info$n_blocks < 5L)
      sprintf("blocked_time: n_blocks = %d < 5", object@info$n_blocks)
    else
      "requires paired R_eff >= 5"
    cat(sprintf("  Sign-flip p:   NA  (%s)\n", reason))
  }

  cat(sprintf("\nInference valid:   %s\n",
              if (object@inference_ok) "YES (tier A)" else
                sprintf("NO (%s; need paired R_eff >= 20)", object@tier)))
  cat("\n")
  invisible(object)
}
