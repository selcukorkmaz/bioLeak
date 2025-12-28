test_that("knn imputation does not depend on test data values", {
  skip_if_not_installed("VIM")

  train <- data.frame(
    a = c(1, 2, NA, 4, 5),
    b = c(2, NA, 2, 2, 2)
  )

  test1 <- data.frame(
    a = c(NA, 1000, 3),
    b = c(NA, 2, NA)
  )
  test2 <- data.frame(
    a = c(NA, -999, 3),
    b = c(NA, 2, NA)
  )

  res1 <- impute_guarded(train, test1, method = "knn", k = 3,
                         winsor = FALSE, seed = 1)
  res2 <- impute_guarded(train, test2, method = "knn", k = 3,
                         winsor = FALSE, seed = 1)

  expect_identical(res1$train, res2$train)
  na_idx <- is.na(test1)
  expect_true(any(na_idx))
  expect_identical(res1$test[na_idx], res2$test[na_idx])
  expect_identical(res1$model$impute_values, res2$model$impute_values)
})

test_that("missForest imputation does not depend on test data values", {
  skip_if_not_installed("missForest")

  train <- data.frame(
    a = c(1, 2, NA, 4, 5),
    b = c(2, NA, 2, 2, 2)
  )

  test1 <- data.frame(
    a = c(NA, 1000, 3),
    b = c(NA, 2, NA)
  )
  test2 <- data.frame(
    a = c(NA, -999, 3),
    b = c(NA, 2, NA)
  )

  res1 <- impute_guarded(train, test1, method = "missForest",
                         winsor = FALSE, seed = 1, parallel = FALSE)
  res2 <- impute_guarded(train, test2, method = "missForest",
                         winsor = FALSE, seed = 1, parallel = FALSE)

  expect_identical(res1$train, res2$train)
  na_idx <- is.na(test1)
  expect_true(any(na_idx))
  expect_identical(res1$test[na_idx], res2$test[na_idx])
  expect_identical(res1$model$impute_values, res2$model$impute_values)
})
