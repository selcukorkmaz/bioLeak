# tests/testthat/test-delta-lsi.R
# Unit tests for delta_lsi() and its internal helpers.

# ── Helper: build minimal LeakFit objects ─────────────────────────────────────

# Build a LeakSplits with explicit test/train indices.
# distinct = TRUE gives each repeat its own partition (as repeated random CV
# does); distinct = FALSE re-uses one partition in every repeat (as
# leave-one-group-out does), which delta_lsi() collapses to one unit.
.make_dlsi_splits <- function(n_folds, n_obs, n_repeats = 1L, distinct = TRUE) {
  fold_size <- floor(n_obs / n_folds)
  perms <- lapply(seq_len(n_repeats), function(r) {
    if (!distinct || r == 1L) return(seq_len(n_obs))
    old <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv())
    set.seed(7919L * r)
    p <- sample.int(n_obs)
    if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
    p
  })
  indices   <- lapply(seq_len(n_folds * n_repeats), function(i) {
    fold_num <- ((i - 1L) %% n_folds) + 1L
    rep_num  <- ceiling(i / n_folds)
    te       <- sort(perms[[rep_num]][((fold_num - 1L) * fold_size + 1L):(fold_num * fold_size)])
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
                           seed = 1L, distinct = TRUE) {
  n_folds <- length(fold_aucs) %/% n_repeats
  splits  <- .make_dlsi_splits(n_folds, n_obs, n_repeats, distinct = distinct)

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

# ── Tier assignment and the corrected repeated-CV t-test ─────────────────────
# R_eff is the number of folds per repeat; repeats add no independent units.

# Paired fits with k folds per repeat, R repeats, noisy fold metrics and a
# leaky-minus-guarded gap.
.noisy_pair <- function(k, R, gap, sd = 0.02, n_obs = 120L, seed = 1L) {
  set.seed(seed)
  base  <- runif(k * R, 0.6, 0.7)
  aux_n <- base + gap + rnorm(k * R, 0, sd)
  aux_g <- base + rnorm(k * R, 0, sd)
  list(n = .make_dlsi_fit(aux_n, n_obs = n_obs, n_repeats = R),
       g = .make_dlsi_fit(aux_g, n_obs = n_obs, n_repeats = R))
}

test_that("tier D_insufficient with fewer than 3 folds per repeat", {
  fit_n <- .make_dlsi_fit(c(0.75, 0.70), n_repeats = 1L)
  fit_g <- .make_dlsi_fit(c(0.60, 0.62), n_repeats = 1L)
  expect_warning(res <- delta_lsi(fit_n, fit_g, metric = "auc"), regexp = "R_eff")
  expect_equal(res@tier, "D_insufficient")
  expect_equal(res@R_eff, 2L)
  expect_true(is.na(res@p_value))
  expect_false(res@inference_ok)
})

test_that("tier C for 3-9 folds per repeat, with p-value and CI", {
  f <- .noisy_pair(k = 5L, R = 4L, gap = 0.05)
  res <- delta_lsi(f$n, f$g, metric = "auc")
  expect_equal(res@tier, "C_signflip")
  expect_equal(res@R_eff, 5L)
  expect_true(is.finite(res@p_value))
  expect_true(all(is.finite(res@delta_metric_ci)))
  expect_true(all(is.finite(res@delta_lsi_ci)))
  expect_true(res@inference_ok)
})

test_that("tier B for 10-19 folds per repeat", {
  f <- .noisy_pair(k = 10L, R = 2L, gap = 0.05)
  res <- delta_lsi(f$n, f$g, metric = "auc")
  expect_equal(res@tier, "B_signflip_ci")
  expect_equal(res@R_eff, 10L)
})

test_that("tier A for >= 20 folds per repeat", {
  f <- .noisy_pair(k = 20L, R = 1L, gap = 0.05)
  res <- delta_lsi(f$n, f$g, metric = "auc")
  expect_equal(res@tier, "A_full_inference")
  expect_true(res@inference_ok)
})

test_that("repeats do not raise R_eff or the tier", {
  f <- .noisy_pair(k = 5L, R = 20L, gap = 0.01)
  res <- delta_lsi(f$n, f$g, metric = "auc")
  expect_equal(res@R_eff, 5L)
  expect_equal(res@tier, "C_signflip")
  expect_equal(res@info$df, 4L)
  expect_equal(res@info$n_fold_deltas, 100L)
})

test_that("p-value and CI equal the corrected repeated-CV t computed by hand", {
  f <- .noisy_pair(k = 5L, R = 3L, gap = 0.03, n_obs = 100L)
  res <- delta_lsi(f$n, f$g, metric = "auc", return_details = TRUE)
  fd  <- res@info$fold_deltas
  J   <- nrow(fd)
  expect_equal(res@info$n_test_sets, J)   # random repeats: no repeated test sets
  rho <- mean(fd$n_test) / mean(fd$n_train)
  expect_equal(rho, 20 / 80)
  se  <- sqrt((1 / J + rho) * var(fd$d))
  est <- mean(res@info$delta_r)
  expect_equal(res@delta_metric, est)
  expect_equal(res@p_value, 2 * pt(-abs(est / se), df = 4))
  expect_equal(res@delta_metric_ci, est + c(-1, 1) * qt(0.975, 4) * se)
  expect_equal(res@delta_lsi_ci, res@delta_lsi + c(-1, 1) * qt(0.975, 4) * se)
})

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
  expect_true(is.finite(res@delta_lsi))
  expect_true(is.finite(res@delta_metric))
  expect_gt(res@delta_lsi, 0.0)
})

# ── Return details ────────────────────────────────────────────────────────────

test_that("return_details = TRUE stores delta_r, fold deltas and fits in info", {
  set.seed(6)
  fit_n <- .make_dlsi_fit(rep(0.75, 15L), n_obs = 45L, n_repeats = 5L)
  fit_g <- .make_dlsi_fit(rep(0.55, 15L), n_obs = 45L, n_repeats = 5L)
  res   <- delta_lsi(fit_n, fit_g, metric = "auc",
                     return_details = TRUE, seed = 6L)
  expect_true(all(c("delta_r", "fold_deltas", "fit_leaky", "fit_naive",
                    "fit_guarded") %in% names(res@info)))
  expect_length(res@info$delta_r, 5L)
  expect_equal(nrow(res@info$fold_deltas), 15L)
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

# ── CI and p-value behaviour ──────────────────────────────────────────────────

test_that("CI covers 0 and p is large when there is no leakage", {
  f <- .noisy_pair(k = 5L, R = 10L, gap = 0)
  res <- delta_lsi(f$n, f$g, metric = "auc")
  expect_lte(res@delta_metric_ci[1L], 0)
  expect_gte(res@delta_metric_ci[2L], 0)
  expect_gt(res@p_value, 0.05)
})

test_that("identical fits give zero estimate, p = 1 and a zero-width CI", {
  aucs  <- rep(0.65, 30L)
  fit_n <- .make_dlsi_fit(aucs, n_obs = 90L, n_repeats = 10L)
  fit_g <- .make_dlsi_fit(aucs, n_obs = 90L, n_repeats = 10L)
  res   <- delta_lsi(fit_n, fit_g, metric = "auc")
  expect_equal(res@delta_lsi, 0.0, tolerance = 1e-10)
  expect_equal(res@p_value, 1)
  expect_equal(res@delta_metric_ci, c(0, 0))
})

test_that("p-value is in [0, 1]", {
  f <- .noisy_pair(k = 3L, R = 5L, gap = 0.02)
  res <- delta_lsi(f$n, f$g, metric = "auc")
  expect_gte(res@p_value, 0.0)
  expect_lte(res@p_value, 1.0)
})

test_that("p-value small when inflation is large and consistent", {
  f <- .noisy_pair(k = 10L, R = 5L, gap = 0.2, sd = 0.01)
  res <- delta_lsi(f$n, f$g, metric = "auc")
  expect_lt(res@p_value, 0.01)
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
  expect_equal(res@R_eff, n_folds)
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


# ── Exchangeability settings are reported, not separate tests ─────────────────

test_that("all exchangeability settings give the same corrected test", {
  f <- .noisy_pair(k = 5L, R = 4L, gap = 0.03)
  ps <- vapply(c("iid", "by_group", "within_batch", "blocked_time"), function(e)
    delta_lsi(f$n, f$g, metric = "auc", exchangeability = e)@p_value, numeric(1))
  expect_true(all(abs(ps - ps[1]) < 1e-12))
})

test_that("by_group / within_batch report the number of group / batch levels", {
  cd <- data.frame(study = rep(paste0("S", 1:6), each = 5),
                   batch = rep(paste0("B", 1:3), each = 10), stringsAsFactors = FALSE)
  f <- .noisy_pair(k = 6L, R = 1L, gap = 0.03, n_obs = 30L)
  for (nm in c("n", "g")) {
    f[[nm]]@splits@mode <- "study_loocv"
    f[[nm]]@splits@info <- list(study = "study", batch = "batch", coldata = cd)
  }
  rg <- delta_lsi(f$n, f$g, metric = "auc", exchangeability = "by_group")
  rb <- delta_lsi(f$n, f$g, metric = "auc", exchangeability = "within_batch")
  expect_equal(rg@info$n_groups_used, 6L)
  expect_equal(rg@R_eff, 6L)
  expect_equal(rb@info$n_batches_used, 3L)
  expect_equal(rg@p_value, rb@p_value)
})

test_that("identical repeated test sets count once in the variance", {
  cd <- data.frame(b = rep(paste0("B", 1:6), each = 5), stringsAsFactors = FALSE)
  set.seed(4)
  m1 <- runif(6, 0.6, 0.8); m2 <- runif(6, 0.6, 0.8)
  one <- .make_dlsi_fit(m1, n_obs = 30L); one_g <- .make_dlsi_fit(m2, n_obs = 30L)
  ten <- .make_dlsi_fit(rep(m1, 10), n_obs = 30L, n_repeats = 10L, distinct = FALSE)
  ten_g <- .make_dlsi_fit(rep(m2, 10), n_obs = 30L, n_repeats = 10L, distinct = FALSE)
  r1  <- delta_lsi(one, one_g, metric = "auc")
  r10 <- delta_lsi(ten, ten_g, metric = "auc")
  expect_equal(r10@info$n_test_sets, 6L)
  expect_equal(r10@info$n_fold_deltas, 60L)
  expect_equal(r10@p_value, r1@p_value)
  expect_equal(r10@delta_metric_ci, r1@delta_metric_ci)
})

test_that("block_size is ignored with a warning", {
  f <- .noisy_pair(k = 5L, R = 2L, gap = 0.03)
  expect_warning(delta_lsi(f$n, f$g, metric = "auc", exchangeability = "blocked_time",
                           block_size = 2L), regexp = "block_size")
})

# ── Missing metrics (BL-10, BL-12) ────────────────────────────────────────────

test_that("an all-NA repeat in one fit keeps the comparison paired", {
  f <- .noisy_pair(k = 5L, R = 6L, gap = 0.03)
  f$n@metrics$auc[1:5] <- NA          # repeat 1 missing in the leaky fit only
  w <- testthat::capture_warnings(res <- delta_lsi(f$n, f$g, metric = "auc"))
  expect_true(any(grepl("dropped from both fits", w)))
  expect_true(res@info$paired)
  expect_equal(res@R_eff, 5L)
  expect_equal(res@info$n_fold_deltas, 25L)
  expect_true(is.finite(res@p_value))
})

test_that("a fold missing in one fit is dropped from both before pairing", {
  aux_n <- c(0.80, 0.80, 0.80, 0.80)
  aux_g <- c(0.60, 0.60, 0.60, 0.20)   # fold 4 differs only in guarded
  fit_n <- .make_dlsi_fit(aux_n, n_obs = 40L)
  fit_g <- .make_dlsi_fit(aux_g, n_obs = 40L)
  fit_n@metrics$auc[4] <- NA
  res <- suppressWarnings(delta_lsi(fit_n, fit_g, metric = "auc", return_details = TRUE))
  expect_equal(res@info$delta_r, 0.2)   # fold 4 left out of both, not just one
  expect_equal(res@info$n_fold_deltas, 3L)
})

# ── Compact splits (BL-11) ────────────────────────────────────────────────────

test_that("compact and explicit splits give identical results", {
  set.seed(3)
  df <- data.frame(g = rep(sprintf("G%02d", 1:20), each = 3), y = rbinom(60, 1, 0.5), x = rnorm(60))
  sp_e <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "g",
                          v = 5, repeats = 3, seed = 9, progress = FALSE)
  sp_c <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "g",
                          v = 5, repeats = 3, seed = 9, compact = TRUE, progress = FALSE)
  k <- length(sp_e@indices)
  m1 <- runif(k, 0.6, 0.8); m2 <- runif(k, 0.6, 0.8)
  mk <- function(sp, m) methods::new("LeakFit", splits = sp,
    metrics = data.frame(fold = seq_len(k), learner = "glm", auc = m),
    outcome = "y", task = "binomial")
  re <- delta_lsi(mk(sp_e, m1), mk(sp_e, m2), metric = "auc")
  rc <- delta_lsi(mk(sp_c, m1), mk(sp_c, m2), metric = "auc")
  expect_equal(rc@delta_metric, re@delta_metric)
  expect_equal(rc@p_value, re@p_value)
  expect_equal(rc@info$rho, re@info$rho)
})

test_that("delta_lsi() leaves the caller's RNG state unchanged", {
  f <- .noisy_pair(k = 5L, R = 3L, gap = 0.03)
  set.seed(10); a <- runif(1)
  set.seed(10); invisible(delta_lsi(f$n, f$g, metric = "auc")); b <- runif(1)
  expect_identical(a, b)
})

# ── Review round 4 (corrected-t engine) ───────────────────────────────────────

test_that("compact time_series and combined plans use the real training sizes", {
  set.seed(3)
  n <- 100
  dat <- data.frame(t = seq_len(n), y = rbinom(n, 1, 0.5),
                    s = rep(paste0("S", 1:5), each = 20),   # subjects nested in studies
                    p = rep(paste0("P", 1:25), each = 4))
  mk <- function(sp, m) methods::new("LeakFit", splits = sp,
    metrics = data.frame(fold = seq_along(sp@indices), learner = "glm", auc = m),
    outcome = "y", task = "binomial")
  plans <- list(
    ts = list(mode = "time_series", time = "t", v = 5, horizon = 10, purge = 5),
    cb = list(mode = "combined", v = 3,
              constraints = list(list(type = "subject", col = "p"), list(type = "study", col = "s")))
  )
  for (pl in plans) {
    sp_e <- do.call(make_split_plan, c(list(dat, outcome = "y", progress = FALSE, seed = 2), pl))
    sp_c <- do.call(make_split_plan, c(list(dat, outcome = "y", progress = FALSE, seed = 2, compact = TRUE), pl))
    k  <- length(sp_e@indices)
    m1 <- runif(k, 0.5, 0.8); m2 <- runif(k, 0.5, 0.8)
    re <- delta_lsi(mk(sp_e, m1), mk(sp_e, m2), metric = "auc")
    rc <- delta_lsi(mk(sp_c, m1), mk(sp_c, m2), metric = "auc")
    expect_equal(bioLeak:::.dlsi_fold_train_sizes(sp_c),
                 vapply(sp_e@indices, function(z) as.numeric(length(z$train)), numeric(1)))
    expect_equal(rc@info$rho, re@info$rho)
    expect_equal(rc@p_value, re@p_value)
  }
})

test_that("metric direction recognises lower-is-better names and warns on unknown ones", {
  dir <- bioLeak:::.dlsi_metric_direction
  for (m in c("rmse", "mape", "smape", "mase", "mn_log_loss", "huber_loss",
              "classification_cost", "brier"))
    expect_false(dir(m), info = m)
  for (m in c("auc", "roc_auc", "pr_auc", "accuracy", "macro_f1", "kappa", "rsq", "cindex"))
    expect_true(suppressWarnings(dir(m)), info = m)
  expect_warning(dir("my_score"), "Cannot tell")
  for (m in c("ece", "fpr", "nll", "cross_entropy", "misclassification_rate"))
    expect_false(dir(m), info = m)
  for (m in c("brier_skill_score", "deviance_explained", "kap", "ccc", "roc_aunp", "gain_capture"))
    expect_true(expect_silent(dir(m)), info = m)
})

test_that("a fold missing from one fit's metrics is counted as dropped", {
  f <- .noisy_pair(k = 5L, R = 3L, gap = 0.03)
  f$n@metrics <- f$n@metrics[-2, ]            # fold 2 absent from the leaky fit
  w <- testthat::capture_warnings(res <- delta_lsi(f$n, f$g, metric = "auc"))
  expect_true(any(grepl("1 fold\\(s\\)", w)))
  expect_equal(res@info$n_fold_deltas, 14L)
  expect_true(is.na(res@folds_guarded$metric[res@folds_guarded$fold == 2]))
})

test_that("repeat counts in info follow the paired folds", {
  f <- .noisy_pair(k = 5L, R = 3L, gap = 0.03)
  f$n@metrics$auc[1:5] <- NA                   # whole repeat 1 fails in the leaky fit
  res <- suppressWarnings(delta_lsi(f$n, f$g, metric = "auc"))
  expect_equal(res@info$R_naive, 2L)
  expect_equal(res@info$R_guarded, 2L)
  expect_equal(nrow(res@repeats_guarded), 2L)
})

test_that("per-pipeline summaries use the same paired folds as the delta", {
  fit_n <- .make_dlsi_fit(c(0.6, 0.6, 0.6, 0.9), n_obs = 40L)
  fit_g <- .make_dlsi_fit(c(0.6, 0.6, 0.6, 0.6), n_obs = 40L)
  fit_g@metrics$auc[4] <- NA
  res <- suppressWarnings(delta_lsi(fit_n, fit_g, metric = "auc"))
  expect_equal(res@delta_metric, 0)
  expect_equal(res@info$metric_naive, res@info$metric_guarded)
})
