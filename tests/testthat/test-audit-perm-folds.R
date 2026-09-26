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

  # The fixed-prediction null supports the fold mean too.
  a_fixed <- audit_quiet(fit, B = 20, perm_refit = FALSE, perm_summary = "fold_mean")
  expect_equal(audit_perm_gap(a_fixed)$metric_obs, signif(mean(per_fold), 6))
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
