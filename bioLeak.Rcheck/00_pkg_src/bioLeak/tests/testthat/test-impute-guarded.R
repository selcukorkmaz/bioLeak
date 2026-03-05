test_that("impute_guarded matches guarded preprocessing output", {
  train <- data.frame(
    a = c(1, 2, NA, 4),
    b = c(NA, 1, 2, 3)
  )
  test <- data.frame(
    a = c(NA, 5),
    b = c(1, NA)
  )

  res <- impute_guarded(train, test, method = "median", winsor = FALSE, seed = 1)
  guard <- bioLeak:::.guard_fit(
    train,
    y = NULL,
    steps = list(
      impute = list(method = "median", winsor = FALSE, winsor_k = 3),
      normalize = list(method = "none"),
      filter = list(var_thresh = 0, iqr_thresh = 0),
      fs = list(method = "none")
    ),
    task = "gaussian"
  )

  expect_equal(res$train, guard$transform(train))
  expect_equal(res$test, guard$transform(test))
})

test_that("impute_guarded validates inputs and emits warnings", {
  train <- data.frame(a = c(1, NA), b = c(2, 3))
  test <- data.frame(a = c(NA, 4), b = c(1, NA))

  expect_error(impute_guarded(1, test), "data frames")
  expect_error(impute_guarded(train, 1), "data frames")

  res <- expect_warning_match(
    impute_guarded(train, test, method = "median",
                   parallel = TRUE, return_outliers = TRUE,
                   constant_value = 1),
    "unused|not supported",
    all = TRUE
  )
  expect_equal(class(res), "LeakImpute")
})

test_that("impute_guarded respects vars selection", {
  train <- data.frame(a = c(1, NA), b = c(2, 3), c = c(4, 5))
  test <- data.frame(a = c(NA, 4), b = c(1, NA), c = c(6, 7))
  res <- impute_guarded(train, test, method = "median", vars = c("a", "b"))
  expect_equal(names(res$train), c("a", "b"))
  expect_equal(names(res$test), c("a", "b"))
})
