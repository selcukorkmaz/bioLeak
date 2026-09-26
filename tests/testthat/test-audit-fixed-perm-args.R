# Regression tests: with perm_refit = FALSE the permutation-gap null is a
# global label shuffle. Before 0.3.8 a restricted permutation source was
# built and silently discarded, so perm_stratify / time_block / block_len had
# no effect without telling the user.

make_fixed_fit <- function(seed = 1) {
  set.seed(seed)
  df <- make_class_df(40)
  sp <- make_split_plan(df, outcome = "outcome", mode = "subject_grouped",
                        group = "subject", v = 4, progress = FALSE)
  custom <- make_custom_learners()
  fit_resample(df, outcome = "outcome", splits = sp, learner = "glm",
               custom_learners = custom, metrics = "auc", refit = FALSE, seed = 1)
}

test_that("default fixed-prediction audit raises no restricted-permutation warning", {
  fit <- make_fixed_fit()
  expect_no_warning(
    aud <- audit_leakage(fit, metric = "auc", B = 20, perm_refit = FALSE,
                         target_scan = FALSE, target_scan_multivariate = FALSE)
  )
  expect_identical(aud@info$perm_null, "global_shuffle")
})

test_that("perm_stratify with perm_refit = FALSE warns that it has no effect", {
  fit <- make_fixed_fit()
  expect_warning(
    audit_leakage(fit, metric = "auc", B = 20, perm_refit = FALSE,
                  perm_stratify = TRUE,
                  target_scan = FALSE, target_scan_multivariate = FALSE),
    "perm_stratify has no effect on the fixed-prediction permutation null"
  )
  expect_warning(
    audit_leakage(fit, metric = "auc", B = 20, perm_refit = FALSE,
                  perm_stratify = "auto",
                  target_scan = FALSE, target_scan_multivariate = FALSE),
    "perm_stratify"
  )
})

test_that("time_block and block_len with perm_refit = FALSE warn", {
  fit <- make_fixed_fit()
  expect_warning(
    audit_leakage(fit, metric = "auc", B = 20, perm_refit = FALSE,
                  time_block = "stationary", block_len = 5,
                  target_scan = FALSE, target_scan_multivariate = FALSE),
    "time_block, block_len have no effect"
  )
})

test_that("the fixed-prediction null does not depend on perm_stratify", {
  fit <- make_fixed_fit()
  a1 <- audit_leakage(fit, metric = "auc", B = 30, perm_refit = FALSE, seed = 7,
                      target_scan = FALSE, target_scan_multivariate = FALSE)
  a2 <- suppressWarnings(
    audit_leakage(fit, metric = "auc", B = 30, perm_refit = FALSE, seed = 7,
                  perm_stratify = TRUE,
                  target_scan = FALSE, target_scan_multivariate = FALSE)
  )
  expect_identical(a1@perm_values, a2@perm_values)
})
