skip_if_not_installed("recipes")

test_that("guard_to_recipe returns a recipe object", {
  df <- data.frame(outcome = c(0, 1, 0, 1), x1 = rnorm(4), x2 = rnorm(4))
  steps <- list(impute = list(method = "median"),
                normalize = list(method = "zscore"))
  rec <- guard_to_recipe(steps, outcome ~ ., df)
  expect_s3_class(rec, "recipe")
})

test_that("median + zscore adds correct step classes", {
  df <- data.frame(outcome = c(0, 1, 0, 1), x1 = rnorm(4), x2 = rnorm(4))
  steps <- list(impute = list(method = "median"),
                normalize = list(method = "zscore"))
  rec <- guard_to_recipe(steps, outcome ~ ., df)
  step_classes <- vapply(rec$steps, function(s) class(s)[1], character(1))
  expect_true("step_impute_median" %in% step_classes)
  expect_true("step_normalize" %in% step_classes)
})

test_that("kNN imputation adds step_impute_knn", {
  df <- data.frame(outcome = c(0, 1, 0, 1), x1 = rnorm(4), x2 = rnorm(4))
  steps <- list(impute = list(method = "knn", k = 3))
  rec <- guard_to_recipe(steps, outcome ~ ., df)
  step_classes <- vapply(rec$steps, function(s) class(s)[1], character(1))
  expect_true("step_impute_knn" %in% step_classes)
})

test_that("warns for unsupported FS methods (ttest, lasso)", {
  df <- data.frame(outcome = c(0, 1, 0, 1), x1 = rnorm(4), x2 = rnorm(4))
  expect_warning(
    guard_to_recipe(list(fs = list(method = "ttest")), outcome ~ ., df),
    "ttest"
  )
  expect_warning(
    guard_to_recipe(list(fs = list(method = "lasso")), outcome ~ ., df),
    "lasso"
  )
})

test_that("warns for missForest imputation fallback", {
  df <- data.frame(outcome = c(0, 1, 0, 1), x1 = rnorm(4), x2 = rnorm(4))
  expect_warning(
    guard_to_recipe(list(impute = list(method = "missForest")), outcome ~ ., df),
    "missForest"
  )
})

test_that("output recipe can be prep'd and bake'd without error", {
  set.seed(1)
  df <- data.frame(outcome = c(0, 1, 0, 1, 0, 1),
                   x1 = c(1, NA, 3, 4, 5, 6),
                   x2 = c(10, 20, 30, 40, 50, 60))
  steps <- list(impute = list(method = "median"),
                normalize = list(method = "zscore"))
  rec <- guard_to_recipe(steps, outcome ~ ., df)
  prepped <- recipes::prep(rec, training = df)
  baked <- recipes::bake(prepped, new_data = df)
  expect_true(is.data.frame(baked))
  expect_false(any(is.na(baked$x1)))
})
