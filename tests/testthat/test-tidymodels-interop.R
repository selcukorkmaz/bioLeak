test_that("fit_resample accepts rsample splits", {
  skip_if_not_installed("rsample")
  df <- make_class_df(12)
  rs <- rsample::vfold_cv(df, v = 2)
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = rs,
                            learner = "glm", custom_learners = make_custom_learners(),
                            metrics = "accuracy", refit = FALSE)
  expect_true(nrow(fit@metrics) > 0)
})

test_that("as_rsample converts LeakSplits", {
  skip_if_not_installed("rsample")
  df <- make_class_df(12)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject", v = 3, seed = 1)
  rs <- as_rsample(splits, data = df)
  expect_true(inherits(rs, "rset"))
  expect_equal(nrow(rs), length(splits@indices))
})

test_that("fit_resample accepts recipe preprocessing", {
  skip_if_not_installed("recipes")
  df <- make_class_df(12)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject", v = 3, seed = 1)
  rec <- recipes::recipe(outcome ~ ., data = df) |>
    recipes::step_normalize(recipes::all_numeric_predictors())
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            preprocess = rec,
                            learner = "glm", custom_learners = make_custom_learners(),
                            metrics = "accuracy", refit = FALSE)
  expect_true(nrow(fit@metrics) > 0)
})

test_that("fit_resample accepts workflow learners", {
  skip_if_not_installed("workflows")
  skip_if_not_installed("parsnip")
  df <- make_class_df(12)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject", v = 3, seed = 1)
  spec <- parsnip::logistic_reg(mode = "classification") |>
    parsnip::set_engine("glm")
  wf <- workflows::workflow() |>
    workflows::add_model(spec) |>
    workflows::add_formula(outcome ~ .)
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = wf, metrics = "accuracy", refit = FALSE)
  expect_true(nrow(fit@metrics) > 0)
})

test_that("fit_resample accepts yardstick metric sets", {
  skip_if_not_installed("yardstick")
  df <- make_class_df(12)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject", v = 3, seed = 1)
  ys <- yardstick::metric_set(yardstick::accuracy, yardstick::roc_auc)
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = make_custom_learners(),
                            metrics = ys, refit = FALSE)
  expect_true(all(c("accuracy", "roc_auc") %in% colnames(fit@metrics)))
})
