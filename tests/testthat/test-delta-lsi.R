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
  res <- delta_lsi(fit, fit, metric = "auc", seed = 1L)
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
  res   <- delta_lsi(fit_n, fit_g, metric = "auc", seed = 2L)
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

test_that("tier C_point_only for R_eff in 5..9", {
  set.seed(3)
  # 5 repeats × 3 folds = 15 folds each
  fit_n <- .make_dlsi_fit(rep(0.8, 15L), n_obs = 45L, n_repeats = 5L)
  fit_g <- .make_dlsi_fit(rep(0.6, 15L), n_obs = 45L, n_repeats = 5L)
  res   <- delta_lsi(fit_n, fit_g, metric = "auc", seed = 3L)
  expect_equal(res@tier, "C_point_only")
  expect_equal(res@R_eff, 5L)
  # CI should be NA (not yet)
  expect_true(all(!is.finite(res@delta_lsi_ci)))
  # But p-value should be finite (sign-flip with R=5)
  expect_true(is.finite(res@p_value))
})

test_that("tier B_ci_only for R_eff in 10..19", {
  set.seed(4)
  fit_n <- .make_dlsi_fit(rep(0.8, 30L), n_obs = 90L, n_repeats = 10L)
  fit_g <- .make_dlsi_fit(rep(0.6, 30L), n_obs = 90L, n_repeats = 10L)
  res   <- delta_lsi(fit_n, fit_g, metric = "auc", seed = 4L,
                     M_boot = 200L, M_flip = 500L)
  expect_equal(res@tier, "B_ci_only")
  expect_equal(res@R_eff, 10L)
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
  expect_true("fit_naive"   %in% names(res@info))
  expect_true("fit_guarded" %in% names(res@info))
  expect_length(res@info$delta_r, 5L)
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
  expected_cols <- c("fold", "metric", "repeat_id", "n", "tau")
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
  # delta_r = c(1, 1); all four sign combos give means: 1, 0, 0, -1
  # obs_stat = mean(1,1) = 1; |null| >= 1 for (++): count=1 and for (--): count=1
  # p = (2 + 1) / (4 + 1) = 0.6
  p <- bioLeak:::.dlsi_sign_flip(c(1, 1), seed = 1L)
  expect_equal(p, 3/5, tolerance = 1e-10)
})
