test_that("fit_resample validates custom learner specifications", {
  df <- make_class_df(10)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 2, seed = 1)

  expect_error(fit_resample_quiet(df, outcome = "outcome", splits = splits,
                                  learner = "glm", custom_learners = "bad"),
               "custom_learners must be a named list")

  expect_error(fit_resample_quiet(df, outcome = "outcome", splits = splits,
                                  learner = "glm", custom_learners = list(list())),
               "custom_learners must be a named list")

  bad_custom <- list(glm = list(fit = function(...) NULL))
  expect_error(fit_resample_quiet(df, outcome = "outcome", splits = splits,
                                  learner = "glm", custom_learners = bad_custom),
               "fit.*predict")
})

test_that("fit_resample validates class weights and positive class", {
  df <- make_class_df(12)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()

  expect_error(fit_resample_quiet(df, outcome = "outcome", splits = splits,
                                  learner = "glm", custom_learners = custom,
                                  class_weights = c(a = 1),
                                  metrics = "auc", refit = FALSE),
               "missing levels")

  expect_error(fit_resample_quiet(df, outcome = "outcome", splits = splits,
                                  learner = "glm", custom_learners = custom,
                                  positive_class = c(0, 1),
                                  metrics = "auc", refit = FALSE),
               "single value")
})

test_that("fit_resample supports class weights with parsnip learners", {
  skip_if_not_installed("parsnip")
  df <- make_class_df(12)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  spec <- parsnip::logistic_reg() |> parsnip::set_engine("glm")
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = spec,
                            class_weights = c("0" = 1, "1" = 1),
                            metrics = "auc", refit = FALSE)
  expect_true(nrow(fit@metrics) > 0)
})

test_that("fit_resample drops invalid metrics with warnings", {
  df <- make_class_df(10)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 2, seed = 1)
  custom <- make_custom_learners()

  expect_warning({
    fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                              learner = "glm", custom_learners = custom,
                              metrics = c("auc", "rmse"), refit = FALSE)
  }, "Dropping metrics")
  expect_true(nrow(fit@metrics) > 0)
})

test_that("fit_resample summarizes metrics with all-NA columns", {
  df <- make_class_df(10)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 2, seed = 1)
  custom <- make_custom_learners()
  na_metric <- function(y, pred) NA_real_

  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = custom,
                            metrics = list(auc = "auc", na_metric = na_metric),
                            refit = FALSE)
  expect_true(nrow(fit@metric_summary) > 0)
  expect_true(any(grepl("^na_metric\\.", colnames(fit@metric_summary))))
})

test_that("fit_resample handles gaussian tasks and ignores classification options", {
  df <- make_regression_df(12)
  splits <- make_splits_quiet(df, outcome = "y",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()

  expect_warning({
    fit <- fit_resample_quiet(df, outcome = "y", splits = splits,
                              learner = "glm", custom_learners = custom,
                              class_weights = c(a = 1, b = 2),
                              positive_class = "a",
                              metrics = c("rmse", "cindex"), refit = FALSE)
  }, "ignored", all = TRUE)
  expect_equal(fit@task, "gaussian")
})

test_that("fit_resample errors on unsupported outcomes", {
  df <- make_class_df(10)
  df$outcome <- factor(c("a", "b", "c", "a", "b", "c", "a", "b", "c", "a"))
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 2, seed = 1)
  custom <- make_custom_learners()
  expect_error(fit_resample_quiet(df, outcome = "outcome", splits = splits,
                                  learner = "glm", custom_learners = custom),
               "require binomial")
})

test_that("fit_resample respects positive_class releveling", {
  df <- make_class_df(12)
  df$outcome <- factor(ifelse(df$outcome == 1, "yes", "no"), levels = c("no", "yes"))
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = custom,
                            positive_class = "yes", metrics = "auc",
                            refit = FALSE)
  expect_equal(fit@info$positive_class, "yes")
})

test_that("fit_resample surfaces custom learner prediction length errors", {
  df <- make_class_df(10)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 2, seed = 1)
  custom <- list(
    bad = list(
      fit = function(x, y, task, weights, ...) stats::glm(y ~ ., data = data.frame(y = y, x),
                                                         family = stats::binomial()),
      predict = function(object, newdata, task, ...) rep(0.5, nrow(newdata) + 1)
    )
  )
  expect_warning({
    expect_error(fit_resample_quiet(df, outcome = "outcome", splits = splits,
                                    learner = "bad", custom_learners = custom,
                                    metrics = "auc", refit = FALSE),
                 "No successful folds were completed")
  }, "Custom learner 'bad' returned")
})

test_that("fit_resample warns when a fold lacks both classes", {
  df <- make_class_df(8)
  df$outcome <- c(rep(0, 4), rep(1, 4))
  indices <- list(
    list(train = 1:4, test = 5:6, fold = 1, repeat_id = 1),
    list(train = c(1:3, 5:8), test = 4, fold = 2, repeat_id = 1)
  )
  splits <- LeakSplits(mode = "custom", indices = indices,
                       info = list(outcome = "outcome", coldata = df))
  custom <- make_custom_learners()
  expect_warning({
    fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                              learner = "glm", custom_learners = custom,
                              metrics = "auc", refit = FALSE)
  }, "only one class")
  expect_true(nrow(fit@metrics) > 0)
})

test_that("fit_resample errors for compact time_series without time metadata", {
  df <- make_class_df(10)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "time_series", time = "time",
                              v = 2, seed = 1, compact = TRUE)
  splits@info$time <- NULL
  custom <- make_custom_learners()
  expect_error(fit_resample_quiet(df, outcome = "outcome", splits = splits,
                                  learner = "glm", custom_learners = custom,
                                  metrics = "auc", refit = FALSE),
               "time_series compact splits")
})

test_that("fit_resample supports parsnip learners when available", {
  skip_if_not_installed("parsnip")
  df <- make_class_df(12)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  spec <- parsnip::logistic_reg() |>
    parsnip::set_engine("glm")
  expect_warning({
    fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                              learner = spec, custom_learners = make_custom_learners(),
                              learner_args = list(alpha = 1),
                              metrics = "auc", refit = FALSE, seed = 1)
  }, "ignored")
  expect_true(nrow(fit@metrics) > 0)
  expect_true(any(nzchar(fit@metrics$learner)))
})
