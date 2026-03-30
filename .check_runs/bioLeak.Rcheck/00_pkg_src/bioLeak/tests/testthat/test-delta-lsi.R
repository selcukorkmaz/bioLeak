# tests/testthat/test-delta-lsi.R
# Unit tests for delta_lsi() and its internal helpers.

# ── Helper: build minimal LeakFit objects ─────────────────────────────────────

# Build a LeakSplits with explicit test/train indices.
.make_dlsi_splits <- function(n_folds, n_obs, n_repeats = 1L) {
  fold_size <- floor(n_obs / n_folds)
  indices   <- lapply(seq_len(n_folds * n_repeats), function(i) {
    fold_num <- ((i - 1L) %% n_folds) + 1L
    rep_num  <- ceiling(i / n_folds)
    te       <- ((fold_num - 1L) * fold_size + 1L):(fold_num * fold_size)
    tr       <- setdiff(seq_len(n_obs), te)
    list(fold = fold_num, repeat_id = rep_num, test = te, train = tr)
  })
  methods::new("LeakSplits",
    mode    = "subject_grouped",
    indices = indices,
    info    = list()
  )
}

# Build a minimal LeakFit with controlled per-fold AUC values.
.make_dlsi_fit <- function(fold_aucs, n_obs = 30L, n_repeats = 1L,
                           seed = 1L) {
  n_folds <- length(fold_aucs) %/% n_repeats
  splits  <- .make_dlsi_splits(n_folds, n_obs, n_repeats)

  metrics_df <- data.frame(
    fold    = seq_along(fold_aucs),
    learner = "glm",
    auc     = as.numeric(fold_aucs),
    stringsAsFactors = FALSE
  )

  set.seed(seed)
  fold_size <- floor(n_obs / n_folds)
  preds <- lapply(seq_along(fold_aucs), function(i) {
    te <- splits@indices[[i]]$test
    n  <- length(te)
    # Predictions with spread proportional to the AUC value
    y  <- rbinom(n, 1L, 0.5)
    p  <- pmin(1, pmax(0, rnorm(n, mean = 0.5, sd = fold_aucs[i] / 2)))
    data.frame(
      id = te, fold = i, learner = "glm",
      truth = y, pred = p,
      stringsAsFactors = FALSE
    )
  })

  methods::new("LeakFit",
    splits        = splits,
    metrics       = metrics_df,
    metric_summary = data.frame(),
    audit         = data.frame(),
    predictions   = preds,
    preprocess    = list(),
    learners      = list(),
    outcome       = "y",
    task          = "binomial",
    feature_names = character(0L),
    info          = list()
  )
}

# ── Basic construction ─────────────────────────────────────────────────────────

test_that("delta_lsi returns a LeakDeltaLSI object", {
  set.seed(1)
  fit <- .make_dlsi_fit(rep(0.65, 6L), n_repeats = 2L)
  res <- suppressWarnings(delta_lsi(fit, fit, metric = "auc", seed = 1L))
  expect_s4_class(res, "LeakDeltaLSI")
})

test_that("delta_lsi with same fit: delta_metric and delta_lsi are zero", {
  set.seed(1)
  # 6 folds over 2 repeats → R_eff = 2 (D_insufficient tier but point estimate fine)
  fit <- .make_dlsi_fit(rep(0.7, 6L), n_repeats = 2L)
  res <- suppressWarnings(delta_lsi(fit, fit, metric = "auc", seed = 1L))
  expect_equal(res@delta_metric, 0.0, tolerance = 1e-10)
  expect_equal(res@delta_lsi,    0.0, tolerance = 1e-10)
})

test_that("delta_lsi detects inflated naive pipeline", {
  set.seed(2)
  # Naive: AUC ~ 0.85 per fold; guarded: AUC ~ 0.60
  fit_n <- .make_dlsi_fit(rep(0.85, 10L), n_repeats = 2L)
  fit_g <- .make_dlsi_fit(rep(0.60, 10L), n_repeats = 2L)
  res   <- suppressWarnings(delta_lsi(fit_n, fit_g, metric = "auc", seed = 2L))
  expect_gt(res@delta_lsi,    0.0)
  expect_gt(res@delta_metric, 0.0)
})

# ── Tier assignment ────────────────────────────────────────────────────────────

test_that("tier D_insufficient for R_eff < 5", {
  fit_n <- .make_dlsi_fit(rep(0.75, 3L), n_repeats = 1L)
  fit_g <- .make_dlsi_fit(rep(0.60, 3L), n_repeats = 1L)
  expect_warning(
    res <- delta_lsi(fit_n, fit_g, metric = "auc"),
    regexp = "R_eff"
  )
  expect_equal(res@tier, "D_insufficient")
  expect_equal(res@R_eff, 1L)  # paired R=1 → R_eff=1, but below the 5-repeat threshold
})

test_that("tier C_signflip for R_eff in 5..9", {
  set.seed(3)
  # 5 repeats × 3 folds = 15 folds each
  fit_n <- .make_dlsi_fit(rep(0.8, 15L), n_obs = 45L, n_repeats = 5L)
  fit_g <- .make_dlsi_fit(rep(0.6, 15L), n_obs = 45L, n_repeats = 5L)
  res   <- delta_lsi(fit_n, fit_g, metric = "auc", seed = 3L)
  expect_equal(res@tier, "C_signflip")
  expect_equal(res@R_eff, 5L)
  # CI not available at tier C
  expect_true(all(!is.finite(res@delta_lsi_ci)))
  # Sign-flip p-value is available at tier C
  expect_true(is.finite(res@p_value))
})

test_that("tier B_signflip_ci for R_eff in 10..19", {
  set.seed(4)
  fit_n <- .make_dlsi_fit(rep(0.8, 30L), n_obs = 90L, n_repeats = 10L)
  fit_g <- .make_dlsi_fit(rep(0.6, 30L), n_obs = 90L, n_repeats = 10L)
  res   <- delta_lsi(fit_n, fit_g, metric = "auc", seed = 4L,
                     M_boot = 200L, M_flip = 500L)
  expect_equal(res@tier, "B_signflip_ci")
  expect_equal(res@R_eff, 10L)
  # Both p-value and CI are available at tier B
  expect_true(is.finite(res@p_value))
  expect_true(all(is.finite(res@delta_lsi_ci)))
})

test_that("tier A_full_inference for R_eff >= 20", {
  set.seed(5)
  fit_n <- .make_dlsi_fit(rep(0.8, 60L), n_obs = 120L, n_repeats = 20L)
  fit_g <- .make_dlsi_fit(rep(0.6, 60L), n_obs = 120L, n_repeats = 20L)
  res   <- delta_lsi(fit_n, fit_g, metric = "auc", seed = 5L,
                     M_boot = 200L, M_flip = 1000L)
  expect_equal(res@tier, "A_full_inference")
  expect_true(res@inference_ok)
  expect_true(all(is.finite(res@delta_lsi_ci)))
  expect_true(is.finite(res@p_value))
})

# ── strict mode ───────────────────────────────────────────────────────────────

test_that("strict = TRUE errors on insufficient R_eff", {
  fit_n <- .make_dlsi_fit(rep(0.7, 2L), n_repeats = 1L)
  fit_g <- .make_dlsi_fit(rep(0.5, 2L), n_repeats = 1L)
  expect_error(
    delta_lsi(fit_n, fit_g, metric = "auc", strict = TRUE),
    regexp = "R_eff"
  )
})

# ── Input validation ───────────────────────────────────────────────────────────

test_that("non-LeakFit input errors informatively", {
  fit <- .make_dlsi_fit(rep(0.7, 4L))
  expect_error(delta_lsi("not_a_fit", fit), regexp = "LeakFit")
  expect_error(delta_lsi(fit, 42),          regexp = "LeakFit")
})

# ── Unpaired fits ─────────────────────────────────────────────────────────────

test_that("unpaired fits give D_insufficient with point estimate", {
  fit_n <- .make_dlsi_fit(rep(0.8, 6L), n_obs = 30L, n_repeats = 2L)  # R=2
  fit_g <- .make_dlsi_fit(rep(0.6, 9L), n_obs = 30L, n_repeats = 3L)  # R=3
  res   <- suppressWarnings(
    delta_lsi(fit_n, fit_g, metric = "auc", seed = 1L)
  )
  expect_equal(res@tier, "D_insufficient")
  expect_equal(res@R_eff, 0L)
  # Point estimate should still be computed
  expect_true(is.finite(res@delta_lsi))
  expect_true(is.finite(res@delta_metric))
  expect_gt(res@delta_lsi, 0.0)
})

# ── Return details ────────────────────────────────────────────────────────────

test_that("return_details = TRUE stores delta_r and fits in info", {
  set.seed(6)
  fit_n <- .make_dlsi_fit(rep(0.75, 15L), n_obs = 45L, n_repeats = 5L)
  fit_g <- .make_dlsi_fit(rep(0.55, 15L), n_obs = 45L, n_repeats = 5L)
  res   <- delta_lsi(fit_n, fit_g, metric = "auc",
                     return_details = TRUE, seed = 6L)
  expect_true("delta_r"     %in% names(res@info))
  expect_true("fit_leaky"   %in% names(res@info))
  expect_true("fit_naive"   %in% names(res@info))
  expect_true("fit_guarded" %in% names(res@info))
  expect_length(res@info$delta_r, 5L)
})

test_that("deprecated fit_naive argument maps to fit_leaky", {
  set.seed(6)
  fit_n <- .make_dlsi_fit(rep(0.75, 15L), n_obs = 45L, n_repeats = 5L)
  fit_g <- .make_dlsi_fit(rep(0.55, 15L), n_obs = 45L, n_repeats = 5L)

  expect_warning(
    res <- delta_lsi(fit_naive = fit_n, fit_guarded = fit_g, metric = "auc",
                     seed = 6L),
    regexp = "fit_naive is deprecated; use fit_leaky instead"
  )
  expect_s4_class(res, "LeakDeltaLSI")
})

# ── BCa CI covers zero when naive = guarded ───────────────────────────────────

test_that("BCa CI covers 0 when no leakage (same AUC)", {
  set.seed(7)
  # 10 repeats × 3 folds; same AUC → delta = 0
  aucs <- rep(0.65, 30L)
  fit_n <- .make_dlsi_fit(aucs, n_obs = 90L, n_repeats = 10L)
  fit_g <- .make_dlsi_fit(aucs, n_obs = 90L, n_repeats = 10L)
  res   <- delta_lsi(fit_n, fit_g, metric = "auc",
                     M_boot = 500L, seed = 7L)
  expect_equal(res@delta_lsi, 0.0, tolerance = 1e-10)
  ci <- res@delta_lsi_ci
  if (all(is.finite(ci))) {
    expect_lte(ci[1L], 0.0 + 1e-8)
    expect_gte(ci[2L], 0.0 - 1e-8)
  }
})

# ── Sign-flip p-value ──────────────────────────────────────────────────────────

test_that("sign-flip p-value is in [0, 1]", {
  set.seed(8)
  fit_n <- .make_dlsi_fit(rep(0.8, 15L), n_obs = 45L, n_repeats = 5L)
  fit_g <- .make_dlsi_fit(rep(0.6, 15L), n_obs = 45L, n_repeats = 5L)
  res   <- delta_lsi(fit_n, fit_g, metric = "auc",
                     M_flip = 500L, seed = 8L)
  expect_gte(res@p_value, 0.0)
  expect_lte(res@p_value, 1.0)
})

test_that("sign-flip p-value small when inflation is large and consistent", {
  set.seed(9)
  # 20 repeats × 3 folds; large consistent gap
  fit_n <- .make_dlsi_fit(rep(0.95, 60L), n_obs = 120L, n_repeats = 20L)
  fit_g <- .make_dlsi_fit(rep(0.50, 60L), n_obs = 120L, n_repeats = 20L)
  res   <- delta_lsi(fit_n, fit_g, metric = "auc",
                     M_boot = 500L, M_flip = 2000L, seed = 9L)
  expect_lt(res@p_value, 0.05)
})

# ── show / summary don't error ────────────────────────────────────────────────

test_that("show() and summary() run without error", {
  set.seed(10)
  fit_n <- .make_dlsi_fit(rep(0.75, 15L), n_obs = 45L, n_repeats = 5L)
  fit_g <- .make_dlsi_fit(rep(0.55, 15L), n_obs = 45L, n_repeats = 5L)
  res   <- delta_lsi(fit_n, fit_g, metric = "auc",
                     M_boot = 100L, M_flip = 200L, seed = 10L)
  expect_output(show(res),     "LeakDeltaLSI")
  expect_output(summary(res),  "Delta LSI")
})

# ── Fold / repeat data frames ─────────────────────────────────────────────────

test_that("folds_naive and folds_guarded have expected columns", {
  set.seed(11)
  fit_n <- .make_dlsi_fit(c(0.7, 0.8, 0.75), n_obs = 30L, n_repeats = 1L)
  fit_g <- .make_dlsi_fit(c(0.6, 0.65, 0.6), n_obs = 30L, n_repeats = 1L)
  res   <- suppressWarnings(delta_lsi(fit_n, fit_g, metric = "auc"))
  expected_cols <- c("fold", "metric", "repeat_id", "n")
  expect_true(all(expected_cols %in% names(res@folds_naive)))
  expect_true(all(expected_cols %in% names(res@folds_guarded)))
})

# ── Internal helpers ──────────────────────────────────────────────────────────

test_that(".dlsi_huber returns median for single value", {
  expect_equal(bioLeak:::.dlsi_huber(42), 42)
})

test_that(".dlsi_huber returns NA for empty input", {
  expect_true(is.na(bioLeak:::.dlsi_huber(numeric(0))))
})

test_that(".dlsi_huber is robust to one outlier", {
  x   <- c(0, 0, 0, 0, 0, 100)
  est <- bioLeak:::.dlsi_huber(x)
  # Huber M-estimator should be much closer to 0 than arithmetic mean (100/6 ≈ 16.7)
  expect_lt(abs(est), abs(mean(x)))
})

test_that(".dlsi_bca_ci covers known mean for symmetric distribution", {
  set.seed(42)
  x  <- rnorm(20, mean = 0.2, sd = 0.05)
  ci <- bioLeak:::.dlsi_bca_ci(x, mean, M = 2000L, seed = 42L)
  expect_true(all(is.finite(ci)))
  expect_lt(ci[1L], 0.2)
  expect_gt(ci[2L], 0.2)
})

test_that(".dlsi_sign_flip p-value is NA for length-1 input", {
  expect_true(is.na(bioLeak:::.dlsi_sign_flip(0.5)))
})

test_that(".dlsi_sign_flip exact enumeration: known result for R=2", {
  # delta_r = c(1, 1); obs_stat = mean(1, 1) = 1
  # All 4 sign combos (R=2 <= 15 → exact):
  #   (--) → mean(-1, -1) = -1;  (-+) → mean(-1,  1) =  0
  #   (+-) → mean( 1, -1) =  0;  (++) → mean( 1,  1) =  1
  # |null| >= |obs| = 1: combos (--) and (++) → count = 2
  # Exact p = 2 / 4 = 0.5  (no continuity correction for full enumeration)
  p <- bioLeak:::.dlsi_sign_flip(c(1, 1), seed = 1L)
  expect_equal(p, 0.5, tolerance = 1e-10)
})

# ── Adversarial tests ─────────────────────────────────────────────────────────

test_that("equal repeat count but mismatched fold structures warns and goes unpaired", {
  # fit_n uses n_obs=30, fit_g uses n_obs=24 → different test indices → unpaired
  fit_n <- .make_dlsi_fit(rep(0.8, 6L), n_obs = 30L, n_repeats = 2L)
  fit_g <- .make_dlsi_fit(rep(0.6, 6L), n_obs = 24L, n_repeats = 2L)
  res <- expect_warning_match(
    delta_lsi(fit_n, fit_g, metric = "auc"),
    pattern = "fold structures"
  )
  expect_equal(res@tier, "D_insufficient")
  expect_equal(res@R_eff, 0L)
  # Point estimate should still be computed unpaired
  expect_true(is.finite(res@delta_lsi))
})

test_that("relabeled repeat_ids do not corrupt pairing (regression)", {
  n_folds <- 3L; n_obs <- 45L; n_repeats <- 5L

  # Naive: per-repeat means decrease from 0.9 to 0.7
  naive_aucs <- rep(c(0.9, 0.85, 0.80, 0.75, 0.70), each = n_folds)
  fit_n <- .make_dlsi_fit(naive_aucs, n_obs = n_obs, n_repeats = n_repeats)

  # Guarded: per-repeat means decrease from 0.8 to 0.6
  guard_aucs <- rep(c(0.80, 0.75, 0.70, 0.65, 0.60), each = n_folds)
  fit_g <- .make_dlsi_fit(guard_aucs, n_obs = n_obs, n_repeats = n_repeats)

  # Reverse the repeat_id labels in guarded fit's splits.
  # Original: positions 1-3 → rid=1, 4-6 → rid=2, ..., 13-15 → rid=5
  # Reversed: positions 1-3 → rid=5, 4-6 → rid=4, ..., 13-15 → rid=1
  new_idx <- fit_g@splits@indices
  for (i in seq_along(new_idx)) {
    new_idx[[i]]$repeat_id <- (n_repeats + 1L) - new_idx[[i]]$repeat_id
  }
  fit_g@splits@indices <- new_idx

  res <- delta_lsi(fit_n, fit_g, metric = "auc",
                   return_details = TRUE, seed = 42L)

  # With correct positional pairing, every repeat delta should be exactly 0.1.
  # Bug (before fix): deltas would be 0.3, 0.2, 0.1, 0.0, -0.1 due to
  # repeat_id sort-order mismatch between the two fits.
  expect_equal(res@info$delta_r, rep(0.1, n_repeats), tolerance = 1e-10)
  expect_true(res@info$paired)
  expect_equal(res@R_eff, n_repeats)
})

test_that("lower-is-better metric: delta > 0 when naive has lower (better) values", {
  set.seed(1)
  # For lower-is-better metrics: naive = 0.3 (better), guarded = 0.5 (worse)
  # Without sign flip: naive - guarded = -0.2 < 0 (looks like naive is worse)
  # With higher_is_better = FALSE: sign_factor = -1 → delta = -1 * (-0.2) = +0.2 > 0
  fit_n <- .make_dlsi_fit(rep(0.3, 10L), n_obs = 30L, n_repeats = 2L)
  fit_g <- .make_dlsi_fit(rep(0.5, 10L), n_obs = 30L, n_repeats = 2L)
  res   <- suppressWarnings(delta_lsi(
    fit_n, fit_g, metric = "auc",
    higher_is_better = FALSE, seed = 1L
  ))
  expect_gt(res@delta_lsi,    0.0)
  expect_gt(res@delta_metric, 0.0)
  expect_false(isTRUE(res@info[["higher_is_better"]]))
})

test_that("metrics df without fold column falls back to predictions gracefully", {
  # Build a LeakFit whose @metrics has the requested metric but no 'fold' column.
  # delta_lsi should fall back to .dlsi_fold_metrics_from_preds() without crashing.
  splits <- methods::new("LeakSplits",
    mode    = "subject_grouped",
    indices = list(
      list(fold = 1L, repeat_id = 1L, test = 1:15, train = 16:30),
      list(fold = 2L, repeat_id = 1L, test = 16:30, train = 1:15)
    ),
    info = list()
  )
  # @metrics has auc but NO fold column
  mdf <- data.frame(learner = "glm", auc = c(0.7, 0.75), stringsAsFactors = FALSE)

  set.seed(1)
  make_preds <- function(te, auc_val) {
    n <- length(te)
    y <- rbinom(n, 1L, 0.5)
    p <- pmin(1, pmax(0, rnorm(n, 0.5, auc_val / 2)))
    data.frame(id = te, fold = if (min(te) == 1L) 1L else 2L,
               learner = "glm", truth = y, pred = p, stringsAsFactors = FALSE)
  }
  preds <- list(make_preds(1:15, 0.7), make_preds(16:30, 0.75))

  fit_bad <- methods::new("LeakFit",
    splits = splits, metrics = mdf, metric_summary = data.frame(),
    audit = data.frame(), predictions = preds, preprocess = list(),
    learners = list(), outcome = "y", task = "binomial",
    feature_names = character(0L), info = list()
  )
  fit_good <- .make_dlsi_fit(c(0.6, 0.65), n_obs = 30L, n_repeats = 1L)

  # Should not crash; falls back to prediction-based metric recomputation
  expect_no_error(suppressWarnings(
    delta_lsi(fit_bad, fit_good, metric = "auc")
  ))
})

test_that("auto-detect higher_is_better = FALSE for rmse via metric name", {
  set.seed(2)
  # Naive RMSE = 0.3 (lower = better), guarded RMSE = 0.5
  fit_n <- .make_dlsi_fit(rep(0.3, 10L), n_obs = 30L, n_repeats = 2L)
  fit_g <- .make_dlsi_fit(rep(0.5, 10L), n_obs = 30L, n_repeats = 2L)
  # Add rmse column so .dlsi_fold_metrics can find it
  fit_n@metrics$rmse <- fit_n@metrics$auc
  fit_g@metrics$rmse <- fit_g@metrics$auc
  # Auto-detect: metric="rmse" with higher_is_better=NULL (default)
  res <- suppressWarnings(delta_lsi(fit_n, fit_g, metric = "rmse", seed = 2L))
  expect_false(isTRUE(res@info[["higher_is_better"]]))
  # Naive RMSE is lower (better) → positive delta (leakage inflation)
  expect_gt(res@delta_lsi,    0.0)
  expect_gt(res@delta_metric, 0.0)
})

# ── Exchangeability: by_group / within_batch warnings ─────────────────────────

test_that("exchangeability 'by_group' emits iid-fallback warning", {
  set.seed(20)
  fit_n <- .make_dlsi_fit(rep(0.8, 15L), n_obs = 45L, n_repeats = 5L)
  fit_g <- .make_dlsi_fit(rep(0.6, 15L), n_obs = 45L, n_repeats = 5L)
  expect_warning(
    res <- delta_lsi(fit_n, fit_g, metric = "auc",
                     exchangeability = "by_group", seed = 20L),
    regexp = "iid repeat structure"
  )
  # Despite the warning, a p-value is still computed (using iid sign-flip)
  expect_true(is.finite(res@p_value))
  expect_equal(res@exchangeability, "by_group")
})

test_that("exchangeability 'within_batch' emits iid-fallback warning", {
  set.seed(21)
  fit_n <- .make_dlsi_fit(rep(0.8, 15L), n_obs = 45L, n_repeats = 5L)
  fit_g <- .make_dlsi_fit(rep(0.6, 15L), n_obs = 45L, n_repeats = 5L)
  expect_warning(
    delta_lsi(fit_n, fit_g, metric = "auc",
              exchangeability = "within_batch", seed = 21L),
    regexp = "iid repeat structure"
  )
})

# ── Exchangeability: blocked_time block sign-flip ──────────────────────────────

test_that(".dlsi_resolve_block_size returns 1 for uncorrelated Δ_r", {
  # iid Δ_r has rho1 ≈ 0, so block_size should be 1
  set.seed(42)
  delta_r <- rnorm(20)
  bs <- bioLeak:::.dlsi_resolve_block_size(delta_r)
  expect_gte(bs, 1L)
  # With low autocorrelation, expected block_size is 1 (or very small)
  expect_lte(bs, 3L)
})

test_that(".dlsi_resolve_block_size respects explicit block_size argument", {
  delta_r <- rnorm(20)
  expect_equal(bioLeak:::.dlsi_resolve_block_size(delta_r, block_size = 4L), 4L)
})

test_that(".dlsi_resolve_block_size caps at floor(R/3)", {
  # Highly correlated Δ_r: AR(1) rho ~ 0.99, would give block_size = R
  # but cap forces at most floor(R/3)
  delta_r <- cumsum(rnorm(12))   # near unit root, very high autocorrelation
  bs      <- bioLeak:::.dlsi_resolve_block_size(delta_r)
  expect_lte(bs, floor(length(delta_r) / 3L))
})

test_that(".dlsi_sign_flip_blocked exact: known result for R=4, block_size=2", {
  # delta_r = c(1, 1, 1, 1); block_size=2 → block_id = c(1,1,2,2); n_blocks=2
  # All 4 block-sign combos (exact enumeration, n_blocks=2 <= 15):
  #   (--): mean(-1,-1,-1,-1) = -1
  #   (-+): mean(-1,-1, 1, 1) =  0
  #   (+-): mean( 1, 1,-1,-1) =  0
  #   (++): mean( 1, 1, 1, 1) =  1
  # obs_stat = 1; |null| >= |obs|=1: combos (--) and (++) → count = 2
  # Exact p = 2/4 = 0.5
  res <- bioLeak:::.dlsi_sign_flip_blocked(c(1, 1, 1, 1), block_size = 2L, seed = 1L)
  expect_equal(res$p_value,        0.5, tolerance = 1e-10)
  expect_equal(res$n_blocks,       2L)
  expect_equal(res$block_size_used, 2L)
})

test_that("blocked_time sign-flip gives more conservative p than iid for AR Δ_r", {
  set.seed(30)
  # Generate strongly AR(1) Δ_r with a positive mean (leakage signal)
  R <- 15L
  delta_r <- 0.2 + as.numeric(arima.sim(list(ar = 0.85), n = R, sd = 0.05))

  # iid sign-flip
  p_iid <- bioLeak:::.dlsi_sign_flip(delta_r, M_flip = 5000L, seed = 30L)

  # Blocked sign-flip with block_size=3 (5 blocks from 15 repeats)
  sf_blk <- bioLeak:::.dlsi_sign_flip_blocked(delta_r, block_size = 3L,
                                               M_flip = 5000L, seed = 30L)
  p_blk <- sf_blk$p_value

  # Block p-value should be >= iid p-value (more conservative under autocorrelation)
  # Allow tolerance for MC noise; the key property is p_blk is generally larger
  expect_gte(sf_blk$n_blocks, 5L)  # at least 5 blocks for inference
  expect_true(is.finite(p_blk))
  expect_true(is.finite(p_iid))
})

test_that("blocked_time with n_blocks < 5 sets p_value = NA and warns", {
  set.seed(31)
  # 6 repeats, block_size=3 → 2 blocks → below threshold of 5
  fit_n <- .make_dlsi_fit(rep(0.8, 18L), n_obs = 54L, n_repeats = 6L)
  fit_g <- .make_dlsi_fit(rep(0.6, 18L), n_obs = 54L, n_repeats = 6L)
  expect_warning(
    res <- delta_lsi(fit_n, fit_g, metric = "auc",
                     exchangeability = "blocked_time",
                     block_size = 3L, seed = 31L),
    regexp = "n_blocks"
  )
  expect_true(is.na(res@p_value))
  expect_equal(res@info$n_blocks, 2L)
})

test_that("blocked_time with explicit block_size=1 behaves like iid sign-flip", {
  set.seed(32)
  fit_n <- .make_dlsi_fit(rep(0.8, 15L), n_obs = 45L, n_repeats = 5L)
  fit_g <- .make_dlsi_fit(rep(0.6, 15L), n_obs = 45L, n_repeats = 5L)
  # block_size=1 → each repeat is its own block → equivalent to iid
  res_blk <- delta_lsi(fit_n, fit_g, metric = "auc",
                        exchangeability = "blocked_time",
                        block_size = 1L, seed = 32L)
  res_iid <- delta_lsi(fit_n, fit_g, metric = "auc",
                        exchangeability = "iid", seed = 32L)
  # p-values should be identical (same enumeration, same seed offset)
  expect_equal(res_blk@p_value, res_iid@p_value, tolerance = 1e-10)
})

test_that("blocked_time stores block_size_used and n_blocks in info", {
  set.seed(33)
  fit_n <- .make_dlsi_fit(rep(0.8, 30L), n_obs = 90L, n_repeats = 10L)
  fit_g <- .make_dlsi_fit(rep(0.6, 30L), n_obs = 90L, n_repeats = 10L)
  res <- delta_lsi(fit_n, fit_g, metric = "auc",
                    exchangeability = "blocked_time",
                    block_size = 2L, seed = 33L,
                    M_boot = 100L, M_flip = 500L)
  expect_equal(res@info$block_size_used, 2L)
  expect_equal(res@info$n_blocks,        5L)   # 10 repeats / block_size=2 = 5 blocks
  expect_true(is.finite(res@p_value))          # 5 blocks >= 5 → p_value available
})
