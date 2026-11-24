# Refit permutation nulls for stratified, grouped designs.
#
# make_split_plan(stratify = TRUE) balances folds on the TRUE labels. Refitting
# on group-permuted labels while keeping those folds makes the permuted class
# balance vary across folds; each fold model's baseline tracks its training
# balance, which is anti-correlated with its test balance, so pooled AUC under
# the null is biased below 0.5 (Parker, Gunter & Bedo 2007). perm_folds =
# "redraw" re-draws the plan on each permuted outcome and perm_summary =
# "fold_mean" averages per-fold AUCs; together the null centres at 0.5.

# Pure noise: 40 families of 6 sharing a family-level outcome, plus 24
# singletons; features are independent of the outcome.
make_noise_family_df <- function(seed = 1) {
  set.seed(seed)
  fy <- stats::rbinom(40, 1, 0.3)
  y <- c(rep(fy, each = 6), stats::rbinom(24, 1, 0.3))
  n <- length(y)
  X <- matrix(stats::rnorm(n * 5), n, dimnames = list(NULL, paste0("x", 1:5)))
  data.frame(grp = c(rep(sprintf("F%02d", 1:40), each = 6), sprintf("S%02d", 1:24)),
             outcome = factor(y, levels = 0:1), X, stringsAsFactors = FALSE)
}

fit_noise_family <- function(df, stratify = TRUE, group = "grp") {
  sp <- make_split_plan(df, outcome = "outcome", mode = "subject_grouped",
                        group = group, v = 5, stratify = stratify, seed = 1,
                        progress = FALSE)
  out <- NULL
  utils::capture.output(
    out <- fit_resample(df, outcome = "outcome", splits = sp, learner = "glm",
                        custom_learners = make_custom_learners(),
                        metrics = "auc", refit = FALSE, seed = 1)
  )
  out
}

audit_quiet <- function(...) {
  out <- NULL
  utils::capture.output(
    out <- audit_leakage(..., metric = "auc", target_scan = FALSE,
                         target_scan_multivariate = FALSE)
  )
  out
}

test_that("re-drawn fold-mean null centres at 0.5; fixed-fold pooled null is biased low", {
  skip_on_cran()
  B <- 50
  fixed_pooled <- numeric(0)
  redrawn_fold_mean <- numeric(0)
  for (s in 1:2) {
    fit <- fit_noise_family(make_noise_family_df(s))
    a_fixed <- audit_quiet(fit, B = B, perm_refit = TRUE, perm_folds = "fixed")
    a_redraw <- audit_quiet(fit, B = B, perm_refit = TRUE, perm_folds = "redraw",
                            perm_summary = "fold_mean")
    expect_identical(audit_info(a_fixed)$perm_scheme, "group_restricted")
    expect_identical(audit_info(a_fixed)$perm_null, "refit_fixed_folds")
    expect_identical(audit_info(a_redraw)$perm_null, "refit_redrawn_folds")
    fixed_pooled <- c(fixed_pooled, a_fixed@perm_values)
    redrawn_fold_mean <- c(redrawn_fold_mean, a_redraw@perm_values)
  }
  expect_true(all(is.finite(fixed_pooled)))
  expect_true(all(is.finite(redrawn_fold_mean)))

  mc_se <- function(v) stats::sd(v) / sqrt(length(v))
  # Re-drawn folds + per-fold-mean AUC: centred at 0.5 within Monte Carlo error.
  expect_lt(abs(mean(redrawn_fold_mean) - 0.5), 4 * mc_se(redrawn_fold_mean))
  # Fixed folds + pooled AUC: the stratification bias pulls the null below 0.5.
  expect_lt(mean(fixed_pooled), 0.5 - 4 * mc_se(fixed_pooled))
  expect_gt(mean(fixed_pooled < 0.5), 0.6)
  expect_gt(mean(redrawn_fold_mean) - mean(fixed_pooled), 0.03)
})

test_that("perm_folds = 'auto' re-draws stratified plans and keeps unstratified ones", {
  df <- make_noise_family_df(1)
  a_strat <- audit_quiet(fit_noise_family(df, stratify = TRUE), B = 3, perm_refit = TRUE)
  info <- audit_info(a_strat)
  expect_identical(info$perm_folds, "redrawn")
  expect_identical(info$perm_null, "refit_redrawn_folds")
  expect_match(info$perm_folds_reason, "^auto: stratified")

  a_plain <- audit_quiet(fit_noise_family(df, stratify = FALSE), B = 3, perm_refit = TRUE)
  info <- audit_info(a_plain)
  expect_identical(info$perm_folds, "fixed")
  expect_identical(info$perm_null, "refit_fixed_folds")
  expect_match(info$perm_folds_reason, "not stratified")
})

test_that("re-drawn permutations differ from fixed-fold permutations", {
  fit <- fit_noise_family(make_noise_family_df(1))
  a_fixed <- audit_quiet(fit, B = 5, perm_refit = TRUE, perm_folds = "fixed", seed = 3)
  a_redraw <- audit_quiet(fit, B = 5, perm_refit = TRUE, perm_folds = "redraw", seed = 3)
  a_redraw2 <- audit_quiet(fit, B = 5, perm_refit = TRUE, perm_folds = "redraw", seed = 3)
  expect_false(isTRUE(all.equal(a_fixed@perm_values, a_redraw@perm_values)))
  # Re-drawing is reproducible for a given seed.
  expect_identical(a_redraw@perm_values, a_redraw2@perm_values)
  # The observed statistic does not depend on the null.
  expect_identical(audit_perm_gap(a_fixed)$metric_obs, audit_perm_gap(a_redraw)$metric_obs)
})

test_that("perm_folds = 'redraw' errors when the plan cannot be re-drawn", {
  fit <- fit_noise_family(make_noise_family_df(1))
  fit@splits@info$coldata <- NULL
  expect_error(
    audit_quiet(fit, B = 2, perm_refit = TRUE, perm_folds = "redraw"),
    "perm_folds = \"redraw\" is not possible"
  )
  # auto falls back to the observed folds and says so.
  expect_message(
    a <- audit_quiet(fit, B = 2, perm_refit = TRUE),
    "cannot be re-drawn"
  )
  expect_identical(audit_info(a)$perm_folds, "fixed")
})

test_that("perm_summary = 'fold_mean' uses the mean of per-fold metrics", {
  fit <- fit_noise_family(make_noise_family_df(1))
  a <- audit_quiet(fit, B = 4, perm_refit = TRUE, perm_summary = "fold_mean")
  per_fold <- vapply(fit@predictions, function(d) {
    bioLeak:::.auc_binary(d$truth, d$pred)
  }, numeric(1))
  expect_equal(audit_perm_gap(a)$metric_obs, signif(mean(per_fold), 6))
  expect_identical(audit_info(a)$perm_summary, "fold_mean")

  # Refit audits report both summaries whichever is chosen.
  gs <- audit_info(a)$perm_gap_summaries
  expect_s3_class(gs, "data.frame")
  expect_identical(gs$summary, c("pooled", "fold_mean"))
  expect_equal(gs$metric_obs[gs$summary == "fold_mean"], mean(per_fold))
  expect_equal(gs$perm_mean[gs$summary == "fold_mean"],
               mean(a@perm_values), tolerance = 1e-6)

  # The fold mean is not calibrated under fixed predictions: pooled is used.
  expect_warning(
    a_fixed <- audit_quiet(fit, B = 20, perm_refit = FALSE, perm_summary = "fold_mean"),
    "not calibrated for fixed-prediction", class = "bioLeak_input_warning")
  expect_identical(audit_info(a_fixed)$perm_summary, "pooled")
  expect_true(all(is.finite(a_fixed@perm_values)))
  expect_null(audit_info(a_fixed)$perm_gap_summaries)
  expect_true(is.na(audit_info(a_fixed)$perm_folds))
})

test_that("an unrestricted refit null on a grouped design raises a warning", {
  df <- make_noise_family_df(1)
  fit <- fit_noise_family(df)
  spec <- fit@info$perm_refit_spec
  # Refit data without the grouping column: the permutation cannot respect groups.
  spec$x <- df[, c("outcome", paste0("x", 1:5))]
  expect_warning(
    a <- audit_quiet(fit, B = 2, perm_refit = TRUE, perm_refit_spec = spec),
    class = "bioLeak_permutation_warning"
  )
  expect_identical(audit_info(a)$perm_scheme, "unrestricted")

  # Supplying coldata with the outcome and the group column restores the
  # group-restricted null and silences the warning.
  spec$coldata <- df[, c("outcome", "grp")]
  expect_no_warning(
    a <- audit_quiet(fit, B = 2, perm_refit = TRUE, perm_refit_spec = spec)
  )
  expect_identical(audit_info(a)$perm_scheme, "group_restricted")
})

test_that("no permutation warning when the grouping is one sample per group", {
  df <- make_noise_family_df(1)
  df$grp <- NULL
  fit <- fit_noise_family(df, group = "row_id")
  expect_no_warning(a <- audit_quiet(fit, B = 2, perm_refit = TRUE))
  expect_identical(audit_info(a)$perm_scheme, "unrestricted")
})

test_that("re-drawn plans keep as_leaksplits()' sample-id column out of the predictors", {
  set.seed(2); n <- 120
  sd <- data.frame(sample_id = paste0("S", 1:n), group_id = paste0("G", rep(1:30, each = 4)),
                   stringsAsFactors = FALSE)
  spec <- structure(list(sample_data = sd, group_var = "group_id", constraint_mode = "subject",
                         ordering_required = FALSE), class = "split_spec")
  dat <- data.frame(sample_id = sd$sample_id, y = factor(rbinom(n, 1, .4)),
                    x1 = rnorm(n), x2 = rnorm(n))
  sp <- as_leaksplits(spec, dat, "y", v = 3, stratify = TRUE, progress = FALSE)
  capture.output(fit <- fit_resample(dat, outcome = "y", splits = sp, learner = "glm",
                                     custom_learners = make_custom_learners(), metrics = "auc",
                                     refit = FALSE, seed = 1))
  ws <- character()
  capture.output(a <- withCallingHandlers(
    audit_leakage(fit, metric = "auc", B = 5, perm_refit = TRUE, perm_folds = "auto",
                  target_scan = FALSE, target_scan_multivariate = FALSE),
    warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") }))
  expect_identical(audit_info(a)$perm_folds, "redrawn")
  expect_true(all(is.finite(a@perm_values)))
  expect_false(any(grepl("sample_id", ws)))
})

test_that("refit permutations of as_leaksplits() plans keep the spec's design", {
  set.seed(4); n <- 120
  sd <- data.frame(sample_id = paste0("S", 1:n),
                   group_id = paste0("P", rep(1:30, each = 4)),
                   batch_group = paste0("b", rep(1:6, each = 20)), stringsAsFactors = FALSE)
  dat <- data.frame(sample_id = sd$sample_id, y = factor(rbinom(n, 1, .4)),
                    x1 = rnorm(n), x2 = rnorm(n))
  for (cm in c("batch", "subject")) {
    spec <- structure(list(sample_data = sd,
                           group_var = if (cm == "batch") "batch_group" else "group_id",
                           constraint_mode = cm, ordering_required = FALSE), class = "split_spec")
    sp <- suppressWarnings(as_leaksplits(spec, dat, "y", v = 3, progress = FALSE))
    capture.output(fit <- suppressWarnings(fit_resample(dat, "y", sp, learner = "glm",
      custom_learners = make_custom_learners(), metrics = "auc", seed = 1)))
    ws <- character()
    capture.output(a <- withCallingHandlers(
      audit_leakage(fit, metric = "auc", B = 5, perm_refit = TRUE,
                    target_scan = FALSE, target_scan_multivariate = FALSE),
      warning = function(w) { ws <<- c(ws, conditionMessage(w)); invokeRestart("muffleWarning") }))
    expect_true(all(is.finite(a@perm_values)), label = cm)
    expect_false(identical(audit_info(a)$perm_scheme, "unrestricted"), label = cm)
    expect_false(any(grepl("unrestricted label shuffle", ws)), label = cm)
  }
})

test_that("supplied coldata without the spec's design columns and perm_stratify on batch plans", {
  set.seed(4); n <- 120
  sd <- data.frame(sample_id = paste0("S", 1:n), batch_group = paste0("b", rep(1:6, each = 20)),
                   stringsAsFactors = FALSE)
  dat <- data.frame(sample_id = sd$sample_id, y = factor(rbinom(n, 1, .4)),
                    x1 = rnorm(n), x2 = rnorm(n))
  spec <- structure(list(sample_data = sd, group_var = "batch_group", constraint_mode = "batch",
                         ordering_required = FALSE), class = "split_spec")
  sp <- suppressWarnings(as_leaksplits(spec, dat, "y", v = 3, progress = FALSE))
  capture.output(fit <- suppressWarnings(fit_resample(dat, "y", sp, learner = "glm",
    custom_learners = make_custom_learners(), metrics = "auc", seed = 1)))
  capture.output(a <- suppressWarnings(audit_leakage(fit, metric = "auc", B = 6, perm_refit = TRUE,
    perm_refit_spec = list(x = dat, learner = "glm", custom_learners = make_custom_learners(),
                           coldata = dat[, c("sample_id", "y")]),
    target_scan = FALSE, target_scan_multivariate = FALSE)))
  expect_identical(audit_info(a)$perm_scheme, "within_batch")
  expect_true(all(is.finite(a@perm_values)))
  capture.output(expect_warning(
    b <- audit_leakage(fit, metric = "auc", B = 6, perm_refit = TRUE, perm_stratify = TRUE,
                       perm_folds = "fixed", target_scan = FALSE, target_scan_multivariate = FALSE),
    "perm_stratify is ignored", class = "bioLeak_input_warning"))
  expect_gt(length(unique(round(b@perm_values, 8))), 1L)
})
