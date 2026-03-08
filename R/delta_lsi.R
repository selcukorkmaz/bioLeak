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

  if (metric %in% names(mdf)) {
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

# Two-pass τ stabilisation: MAD of per-fold predictions.
# Returns a named numeric vector  fold_seq_char -> tau.
.dlsi_tau <- function(pred_df, fold_seq) {
  if (is.null(pred_df) || !"fold" %in% names(pred_df) || !"pred" %in% names(pred_df)) {
    return(setNames(rep(1e-3, length(fold_seq)), as.character(fold_seq)))
  }

  # Pass 1: MAD per fold
  s <- vapply(fold_seq, function(f) {
    p <- as.numeric(pred_df$pred[pred_df$fold == f])
    p <- p[is.finite(p)]
    if (length(p) < 2L) return(NA_real_)
    stats::mad(p, constant = 1.0, na.rm = TRUE)
  }, numeric(1L))

  s_clean   <- s[is.finite(s) & s > 0]
  tau_floor <- if (length(s_clean)) max(1e-3, stats::median(s_clean) / 10) else 1e-3

  # Pass 2: apply floor
  tau                <- pmax(s, tau_floor)
  tau[!is.finite(tau)] <- tau_floor
  setNames(tau, as.character(fold_seq))
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
# Uses exact enumeration for R <= 15, Monte Carlo otherwise.
.dlsi_sign_flip <- function(delta_r, M_flip = 10000L, seed = 42L) {
  delta_r  <- as.numeric(delta_r[is.finite(delta_r)])
  R        <- length(delta_r)
  if (R < 2L) return(NA_real_)

  obs_stat <- mean(delta_r)
  if (!is.finite(obs_stat)) return(NA_real_)

  set.seed(seed)

  null_stats <- if (R <= 15L) {
    # Exact enumeration over all 2^R sign combinations
    n_combos <- 2L^R
    vapply(seq_len(n_combos), function(j) {
      bits  <- as.integer(intToBits(j - 1L)[seq_len(R)])
      signs <- ifelse(bits == 1L, 1.0, -1.0)
      mean(signs * delta_r)
    }, numeric(1L))
  } else {
    # Monte Carlo sign-flip
    replicate(M_flip, {
      signs <- sample(c(-1.0, 1.0), R, replace = TRUE)
      mean(signs * delta_r)
    })
  }

  null_stats <- as.numeric(null_stats[is.finite(null_stats)])
  if (!length(null_stats)) return(NA_real_)

  # Two-sided p-value with +1 continuity correction (Phipson & Smyth)
  (sum(abs(null_stats) >= abs(obs_stat) - .Machine$double.eps * 100) + 1L) /
    (length(null_stats) + 1L)
}

# ── Main function ─────────────────────────────────────────────────────────────

#' Delta Leakage Sensitivity Index (ΔLSI)
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
#' \eqn{\Delta_r = \mu_r^{\text{naive}} - \mu_r^{\text{guarded}}} captures
#' leakage-induced performance inflation for each CV repeat.
#'
#' The \strong{delta_lsi} point estimate is the Huber M-estimator (k = 1.345)
#' applied to \eqn{\{\Delta_r\}}, which is robust to occasional outlier
#' repeats. \strong{delta_metric} is the arithmetic mean of \eqn{\Delta_r}
#' for easy interpretation in the original metric's units.
#'
#' When \eqn{R_{\text{eff}} \geq 5} (equal, paired repeats), a sign-flip
#' randomization test (Phipson & Smyth, 2010) is performed: under
#' \eqn{H_0} (no leakage) the sign of each \eqn{\Delta_r} is exchangeable.
#' All \eqn{2^R} sign combinations are enumerated exactly for
#' \eqn{R \leq 15}; Monte Carlo sampling is used for larger \eqn{R}.
#'
#' BCa bootstrap confidence intervals (Efron, 1987) require
#' \eqn{R_{\text{eff}} \geq 10}.
#' }
#'
#' \subsection{Inference tiers}{
#' \describe{
#'   \item{\code{"A_full_inference"}}{R_eff >= 20: point + BCa CI + p-value}
#'   \item{\code{"B_ci_only"}}{10 <= R_eff < 20: point + BCa CI}
#'   \item{\code{"C_point_only"}}{5 <= R_eff < 10: point + coarse sign-flip p}
#'   \item{\code{"D_insufficient"}}{R_eff < 5 or unpaired: point estimate only}
#' }}
#'
#' @param fit_naive A \code{\link{LeakFit}} object from the naive (unprotected)
#'   evaluation pipeline.
#' @param fit_guarded A \code{\link{LeakFit}} object from the guarded
#'   (leakage-protected) evaluation pipeline.
#' @param metric Character. Performance metric to compare. Must appear in
#'   \code{fit@@metrics} of both fits (e.g., \code{"auc"}, \code{"rmse"}).
#' @param exchangeability Character. Exchangeability assumption used to justify
#'   the sign-flip test. Informational; does not change the computation. One of
#'   \code{"iid"} (default), \code{"by_group"}, \code{"within_batch"},
#'   \code{"blocked_time"}.
#' @param learner Optional character. Learner name to select from multi-learner
#'   fits. If \code{NULL}, the first learner found in \code{fit@@metrics} is used.
#' @param M_boot Integer. Number of bootstrap samples for BCa CI (default 2000).
#' @param M_flip Integer. Maximum Monte Carlo samples for sign-flip test when
#'   R_eff > 15 (default 10000).
#' @param strict Logical. If \code{TRUE}, error on insufficient R_eff instead
#'   of a warning.
#' @param return_details Logical. If \code{TRUE}, include the per-repeat
#'   \eqn{\Delta_r} vector and the original fit objects in the \code{info} slot.
#' @param seed Integer. Random seed for bootstrap and sign-flip test.
#'
#' @return A \code{\link{LeakDeltaLSI}} object.
#'
#' @seealso \code{\link{audit_leakage}}, \code{\link{fit_resample}}
#' @export
delta_lsi <- function(
    fit_naive,
    fit_guarded,
    metric          = "auc",
    exchangeability = c("iid", "by_group", "within_batch", "blocked_time"),
    learner         = NULL,
    M_boot          = 2000L,
    M_flip          = 10000L,
    strict          = FALSE,
    return_details  = FALSE,
    seed            = 42L
) {
  exchangeability <- match.arg(exchangeability)

  if (!inherits(fit_naive,   "LeakFit")) stop("fit_naive must be a LeakFit object.")
  if (!inherits(fit_guarded, "LeakFit")) stop("fit_guarded must be a LeakFit object.")

  # Resolve learners (one per fit; may differ)
  lrn_n <- .dlsi_resolve_learner(fit_naive,   learner)
  lrn_g <- .dlsi_resolve_learner(fit_guarded, learner)

  # Per-fold metric data frames: columns fold, metric
  fm_n <- .dlsi_fold_metrics(fit_naive,   metric, lrn_n)
  fm_g <- .dlsi_fold_metrics(fit_guarded, metric, lrn_g)

  # Sequential fold positions (1..N) for each fit
  fold_seq_n <- seq_along(fit_naive@splits@indices)
  fold_seq_g <- seq_along(fit_guarded@splits@indices)

  # Repeat IDs and fold sizes
  rep_ids_n <- .dlsi_repeat_ids(fit_naive)
  rep_ids_g <- .dlsi_repeat_ids(fit_guarded)

  pred_n <- .dlsi_combine_preds(fit_naive,   lrn_n)
  pred_g <- .dlsi_combine_preds(fit_guarded, lrn_g)

  sizes_n <- .dlsi_fold_sizes(fit_naive,   pred_n)
  sizes_g <- .dlsi_fold_sizes(fit_guarded, pred_g)

  # τ per fold (two-pass MAD stabilisation)
  tau_n <- .dlsi_tau(pred_n, fold_seq_n)
  tau_g <- .dlsi_tau(pred_g, fold_seq_g)

  # Build fold-level data frames
  .build_folds <- function(fm, rep_ids, sizes, tau_vec, fold_seq) {
    idx <- match(fm$fold, fold_seq)
    data.frame(
      fold      = fm$fold,
      metric    = fm$metric,
      repeat_id = rep_ids[idx],
      n         = sizes[idx],
      tau       = tau_vec[as.character(fm$fold)],
      stringsAsFactors = FALSE
    )
  }

  folds_n <- .build_folds(fm_n, rep_ids_n, sizes_n, tau_n, fold_seq_n)
  folds_g <- .build_folds(fm_g, rep_ids_g, sizes_g, tau_g, fold_seq_g)
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
    do.call(rbind, Filter(Negate(is.null), rows))
  }

  reps_n <- .agg_repeats(folds_n)
  reps_g <- .agg_repeats(folds_g)

  R_naive   <- if (!is.null(reps_n)) nrow(reps_n) else 0L
  R_guarded <- if (!is.null(reps_g)) nrow(reps_g) else 0L

  # Paired design requires equal number of repeats
  paired <- (R_naive > 0L && R_guarded > 0L && R_naive == R_guarded)
  R_eff  <- if (paired) as.integer(R_naive) else 0L

  # Inference tier
  tier <- if (R_eff >= 20L) "A_full_inference" else
          if (R_eff >= 10L) "B_ci_only"        else
          if (R_eff >=  5L) "C_point_only"     else
                            "D_insufficient"

  if (tier == "D_insufficient") {
    msg <- sprintf(
      paste0("[delta_lsi] R_eff = %d (R_naive=%d, R_guarded=%d, paired=%s). ",
             "Need paired R_eff >= 5 for sign-flip inference; reporting point estimate only."),
      R_eff, R_naive, R_guarded, paired
    )
    if (isTRUE(strict)) stop(msg) else warning(msg)
  }

  # Per-repeat ΔLSI values (paired: element-wise)
  delta_r <- if (paired) reps_n$metric - reps_g$metric else NULL

  # Point estimates
  delta_lsi_est    <- if (!is.null(delta_r)) {
    .dlsi_huber(delta_r)
  } else {
    .dlsi_huber(reps_n$metric) - .dlsi_huber(reps_g$metric)
  }
  delta_metric_est <- if (!is.null(delta_r)) {
    mean(delta_r, na.rm = TRUE)
  } else {
    mean(reps_n$metric, na.rm = TRUE) - mean(reps_g$metric, na.rm = TRUE)
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

  # Sign-flip test (requires R_eff >= 5 and paired)
  p_val <- NA_real_
  if (!is.null(delta_r) && R_eff >= 5L) {
    p_val <- tryCatch(
      .dlsi_sign_flip(delta_r, M_flip = M_flip, seed = seed + 2L),
      error = function(e) NA_real_
    )
  }

  inference_ok <- (R_eff >= 20L && is.finite(p_val) && all(is.finite(ci_lsi)))

  # Metadata
  info_list <- list(
    paired         = paired,
    R_naive        = R_naive,
    R_guarded      = R_guarded,
    metric_naive   = .dlsi_huber(reps_n$metric),
    metric_guarded = .dlsi_huber(reps_g$metric)
  )
  if (return_details) {
    info_list$delta_r     <- delta_r
    info_list$fit_naive   <- fit_naive
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
setMethod("show", "LeakDeltaLSI", function(object) {
  cat("A LeakDeltaLSI object\n")
  cat(sprintf("  Metric:         %s\n", object@metric))
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
    cat(sprintf("  p-value:        %.4f  (sign-flip)\n", object@p_value))
  invisible(object)
})

# ── summary method ─────────────────────────────────────────────────────────────

#' Summarize a LeakDeltaLSI object
#'
#' Prints a human-readable summary of the ΔLSI analysis comparing naive vs
#' guarded evaluation pipelines.
#'
#' @param object A \code{LeakDeltaLSI} object from \code{\link{delta_lsi}}.
#' @param digits Integer. Number of decimal places to show (default 3).
#' @param ... Unused.
#' @return Invisibly returns \code{object}.
#' @export
summary.LeakDeltaLSI <- function(object, digits = 3L, ...) {
  fmt <- function(x) if (is.finite(x)) formatC(x, digits = digits, format = "f") else "NA"

  cat("\n=====================================\n")
  cat(" bioLeak Delta LSI Summary\n")
  cat("=====================================\n\n")
  cat(sprintf("Metric:            %s\n", object@metric))
  cat(sprintf("Exchangeability:   %s\n", object@exchangeability))
  cat(sprintf("Inference tier:    %s\n", object@tier))
  cat(sprintf("R_eff:             %d  (R_naive=%s, R_guarded=%s, paired=%s)\n",
              object@R_eff,
              object@info[["R_naive"]]   %||% "?",
              object@info[["R_guarded"]] %||% "?",
              object@info[["paired"]]    %||% "?"))

  cat("\nNaive pipeline:    ", fmt(object@info[["metric_naive"]]   %||% NA_real_), "\n", sep = "")
  cat("Guarded pipeline:  ", fmt(object@info[["metric_guarded"]] %||% NA_real_), "\n", sep = "")

  cat("\nPoint estimates (naive - guarded):\n")
  cat(sprintf("  delta_metric:  %s  (raw metric difference)\n", fmt(object@delta_metric)))
  ci_dm <- object@delta_metric_ci
  if (all(is.finite(ci_dm)))
    cat(sprintf("  95%% BCa CI:    [%s, %s]\n", fmt(ci_dm[1L]), fmt(ci_dm[2L])))
  cat(sprintf("  delta_lsi:     %s  (Huber-robust)\n", fmt(object@delta_lsi)))
  ci_lsi <- object@delta_lsi_ci
  if (all(is.finite(ci_lsi)))
    cat(sprintf("  95%% BCa CI:    [%s, %s]\n", fmt(ci_lsi[1L]), fmt(ci_lsi[2L])))

  cat("\nHypothesis test (H0: no leakage inflation):\n")
  if (is.finite(object@p_value)) {
    p <- object@p_value
    sig <- if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
    cat(sprintf("  Sign-flip p:   %.4f %s\n", p, sig))
  } else {
    cat("  Sign-flip p:   NA  (requires paired R_eff >= 5)\n")
  }

  cat(sprintf("\nInference valid:   %s\n",
              if (object@inference_ok) "YES (tier A)" else
                sprintf("NO (%s; need paired R_eff >= 20)", object@tier)))
  cat("\n")
  invisible(object)
}
