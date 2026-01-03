test_that("tune_resample runs nested tuning with workflows", {
  skip_if_not_installed("tune")
  skip_if_not_installed("dials")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("recipes")
  skip_if_not_installed("rsample")
  skip_if_not_installed("yardstick")
  skip_if_not_installed("workflows")

  df <- make_class_df(24)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, nested = TRUE, stratify = TRUE, seed = 1)

  spec <- parsnip::logistic_reg(penalty = tune::tune(), mixture = 1) |>
    parsnip::set_engine("glmnet")
  rec <- recipes::recipe(outcome ~ x1 + x2, data = df)

  tuned <- tune_resample(df, outcome = "outcome", splits = splits,
                         learner = spec, preprocess = rec, grid = 3,
                         metrics = c("auc", "accuracy"), seed = 1)

  expect_s3_class(tuned, "LeakTune")
  expect_true(nrow(tuned$metrics) > 0)
  expect_true(nrow(tuned$best_params) > 0)
})
