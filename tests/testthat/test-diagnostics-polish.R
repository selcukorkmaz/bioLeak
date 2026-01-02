test_that("calibration_summary computes metrics for binomial fits", {
  df <- make_class_df(30)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = custom,
                            metrics = "accuracy", refit = FALSE, seed = 1)
  cal <- calibration_summary(fit, bins = 5)
  expect_true(is.data.frame(cal$curve))
  expect_true(is.data.frame(cal$metrics))
  expect_true(all(c("ece", "mce", "brier") %in% names(cal$metrics)))
})

test_that("confounder_sensitivity summarizes metrics by confounder", {
  df <- make_class_df(30)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = custom,
                            metrics = "accuracy", refit = FALSE, seed = 1)
  cs <- confounder_sensitivity(fit, confounders = "batch", metric = "accuracy", min_n = 2)
  expect_true(nrow(cs) > 0)
  expect_true(all(c("confounder", "level", "metric", "value", "n") %in% names(cs)))
})

test_that("diagnostics require learner selection when multiple models are present", {
  df <- make_class_df(40)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 4, seed = 2)
  custom <- make_custom_learners()
  custom$glm2 <- custom$glm
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = c("glm", "glm2"), custom_learners = custom,
                            metrics = "accuracy", refit = FALSE, seed = 2)

  expect_error(calibration_summary(fit, bins = 5), "Multiple learners found")
  expect_error(confounder_sensitivity(fit, confounders = "batch", metric = "accuracy", min_n = 2),
               "Multiple learners found")

  cal <- calibration_summary(fit, bins = 5, learner = "glm")
  expect_true(is.data.frame(cal$curve))
  cs <- confounder_sensitivity(fit, confounders = "batch", metric = "accuracy",
                               min_n = 2, learner = "glm")
  expect_true(nrow(cs) > 0)
})
