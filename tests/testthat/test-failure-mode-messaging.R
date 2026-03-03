test_that("fold_status contains elapsed_sec column", {
  set.seed(1)
  df <- data.frame(
    subject = rep(1:10, each = 2),
    outcome = rbinom(20, 1, 0.5),
    x1 = rnorm(20),
    x2 = rnorm(20)
  )
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 3, progress = FALSE)
  fit <- fit_resample(df, outcome = "outcome", splits = splits,
                      learner = "glmnet",
                      preprocess = list(
                        impute = list(method = "median"),
                        normalize = list(method = "zscore")
                      ))
  fs <- fit@info$fold_status
  expect_true("elapsed_sec" %in% names(fs))
})

test_that("successful folds have non-negative elapsed_sec", {
  set.seed(1)
  df <- data.frame(
    subject = rep(1:10, each = 2),
    outcome = rbinom(20, 1, 0.5),
    x1 = rnorm(20),
    x2 = rnorm(20)
  )
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 3, progress = FALSE)
  fit <- fit_resample(df, outcome = "outcome", splits = splits,
                      learner = "glmnet",
                      preprocess = list(
                        impute = list(method = "median"),
                        normalize = list(method = "zscore")
                      ))
  fs <- fit@info$fold_status
  success_rows <- fs[fs$status == "success", ]
  expect_true(all(!is.na(success_rows$elapsed_sec)))
  expect_true(all(success_rows$elapsed_sec >= 0))
})

test_that("schema consistency between successful and failed folds", {
  set.seed(1)
  df <- data.frame(
    subject = rep(1:10, each = 2),
    outcome = rbinom(20, 1, 0.5),
    x1 = rnorm(20),
    x2 = rnorm(20)
  )
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 3, progress = FALSE)
  fit <- fit_resample(df, outcome = "outcome", splits = splits,
                      learner = "glmnet",
                      preprocess = list(
                        impute = list(method = "median"),
                        normalize = list(method = "zscore")
                      ))
  fs <- fit@info$fold_status
  expected_cols <- c("fold", "stage", "status", "reason", "notes", "elapsed_sec")
  expect_true(all(expected_cols %in% names(fs)))
})
