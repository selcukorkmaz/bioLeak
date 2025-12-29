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
