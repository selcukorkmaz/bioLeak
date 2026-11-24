# K1 (JSS review): perm_refit = "auto" says what it costs and when it falls back.

k1_setup <- function() {
  set.seed(11)
  df <- make_class_df(40)
  sp <- make_split_plan(df, outcome = "outcome", mode = "subject_grouped",
                        group = "subject", v = 3, seed = 1, progress = FALSE)
  custom <- make_custom_learners()
  fit <- NULL
  utils::capture.output(fit <- fit_resample(df, outcome = "outcome", splits = sp,
    learner = "glm", custom_learners = custom, metrics = "auc", refit = FALSE, seed = 1))
  list(fit = fit, spec = list(x = df, outcome = "outcome", learner = "glm",
                              custom_learners = custom))
}

k1_messages <- function(...) {
  msgs <- character()
  withCallingHandlers(
    suppressWarnings(utils::capture.output(audit_leakage(..., metric = "auc",
      target_scan = FALSE, target_scan_multivariate = FALSE))),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
  )
  msgs
}

test_that("auto refit reports the number of model fits", {
  s <- k1_setup()
  m <- k1_messages(s$fit, B = 4, perm_refit = "auto", perm_refit_spec = s$spec)
  expect_true(any(grepl("4 permutations x 3 folds = 12 model fits", m, fixed = TRUE)))
})

test_that("auto falls back with a message when B exceeds perm_refit_auto_max", {
  s <- k1_setup()
  m <- k1_messages(s$fit, B = 6, perm_refit = "auto", perm_refit_auto_max = 5,
                   perm_refit_spec = s$spec)
  expect_true(any(grepl("exceeds perm_refit_auto_max", m)))
  expect_false(any(grepl("model fits", m)))
})

test_that("the default refits with the spec stored in the fit and reports the cost", {
  s <- k1_setup()
  m <- k1_messages(s$fit, B = 4)
  expect_true(any(grepl("12 model fits", m, fixed = TRUE)))
})

test_that("auto without a refit spec says so, without a calibration claim", {
  s <- k1_setup()
  s$fit@info$perm_refit_spec <- NULL
  m <- k1_messages(s$fit, B = 4, perm_refit = "auto")
  expect_true(any(grepl("perm_refit_spec not supplied", m)))
  expect_false(any(grepl("liberal", m)))
  # explicit choices stay quiet about the cost
  m2 <- k1_messages(s$fit, B = 4, perm_refit = FALSE)
  expect_false(any(grepl("model fits|perm_refit_spec not supplied", m2)))
})
