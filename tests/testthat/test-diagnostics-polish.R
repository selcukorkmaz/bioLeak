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
