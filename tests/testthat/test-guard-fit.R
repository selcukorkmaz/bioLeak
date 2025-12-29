test_that("guard_fit produces a GuardFit pipeline with median imputation", {
  X <- data.frame(a = c(1, 2, NA, 4), b = c(NA, 1, 2, 3))
  fit <- bioLeak:::.guard_fit(
    X, y = c(0, 1, 0, 1),
    steps = list(impute = list(method = "median"),
                 normalize = list(method = "zscore"),
                 filter = list(var_thresh = 0, iqr_thresh = 0),
                 fs = list(method = "none")),
    task = "binomial"
  )

  expect_true(inherits(fit, "GuardFit"))
  Xtr <- fit$transform(X)
  expect_true(is.data.frame(Xtr))
  expect_false(anyNA(Xtr))
})

test_that("guard_fit warns on impute none with missing values", {
  X <- data.frame(a = c(1, NA, 3), b = c(4, 5, NA))
  expect_warning({
    fit <- bioLeak:::.guard_fit(
      X, y = c(1, 0, 1),
      steps = list(impute = list(method = "none"),
                   normalize = list(method = "none"),
                   filter = list(var_thresh = 0, iqr_thresh = 0),
                   fs = list(method = "none")),
      task = "binomial"
    )
  }, "impute\\$method='none'")
  out <- fit$transform(X)
  expect_true(any(grepl("__missing", names(out))))
})

test_that("guard_fit supports knn imputation and fallback warnings", {
  X <- data.frame(a = c(1, NA, 3), b = c(4, 5, NA))
  fit <- bioLeak:::.guard_fit(
    X, y = c(1, 0, 1),
    steps = list(impute = list(method = "knn", k = 1),
                 normalize = list(method = "none"),
                 filter = list(var_thresh = 0, iqr_thresh = 0),
                 fs = list(method = "none")),
    task = "binomial"
  )

  fit$state$impute$knn_ref <- NULL
  expect_warning({
    out <- fit$transform(X)
  }, "KNN imputation state missing")
})

test_that("guard_fit validates imputation method values", {
  X <- data.frame(a = c(1, NA), b = c(2, 3))
  expect_error(bioLeak:::.guard_fit(X, y = c(1, 0),
                                    steps = list(impute = list(method = "mice")),
                                    task = "binomial"),
               "mice")
  expect_error(bioLeak:::.guard_fit(X, y = c(1, 0),
                                    steps = list(impute = list(method = "bad")),
                                    task = "binomial"),
               "Unknown impute method")
})

test_that("guard_fit validates normalization and feature selection modes", {
  X <- data.frame(a = 1:4, b = 2:5)
  expect_error(bioLeak:::.guard_fit(X, y = 1:4,
                                    steps = list(normalize = list(method = "bad")),
                                    task = "gaussian"),
               "Unknown normalize method")

  expect_error(bioLeak:::.guard_fit(X, y = c(1, 1, 1, 1),
                                    steps = list(fs = list(method = "ttest")),
                                    task = "binomial"),
               "requires a binary outcome")
})

test_that("guard_fit filter respects min_keep", {
  X <- data.frame(a = rep(1, 5), b = 1:5, c = 2:6)
  fit <- bioLeak:::.guard_fit(
    X, y = c(0, 1, 0, 1, 0),
    steps = list(impute = list(method = "median"),
                 normalize = list(method = "none"),
                 filter = list(var_thresh = 10, iqr_thresh = 10, min_keep = 1),
                 fs = list(method = "none")),
    task = "binomial"
  )
  expect_true(fit$state$diagnostics$p_out >= 1)
})

test_that("guard_fit supports PCA feature selection", {
  X <- data.frame(a = 1:10, b = 2:11, c = 3:12)
  fit <- bioLeak:::.guard_fit(
    X, y = c(rep(0, 5), rep(1, 5)),
    steps = list(impute = list(method = "median"),
                 normalize = list(method = "none"),
                 filter = list(var_thresh = 0, iqr_thresh = 0),
                 fs = list(method = "pca", pca_comp = 2)),
    task = "binomial"
  )
  out <- fit$transform(X)
  expect_equal(ncol(out), 2)
})

test_that("guard_fit handles missForest when available", {
  X <- data.frame(a = c(1, NA, 3), b = c(4, 5, NA))
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    expect_error(bioLeak:::.guard_fit(
      X, y = c(1, 0, 1),
      steps = list(impute = list(method = "missForest")),
      task = "binomial"
    ), "randomForest")
  } else {
    fit <- bioLeak:::.guard_fit(
      X, y = c(1, 0, 1),
      steps = list(impute = list(method = "missForest", maxiter = 1, ntree = 5),
                   normalize = list(method = "none"),
                   filter = list(var_thresh = 0, iqr_thresh = 0),
                   fs = list(method = "none")),
      task = "binomial"
    )
    expect_true(inherits(fit, "GuardFit"))
  }
})

test_that("guard_fit handles lasso feature selection when available", {
  X <- data.frame(a = rnorm(20), b = rnorm(20), c = rnorm(20))
  y <- factor(sample(c(0, 1), 20, replace = TRUE))
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    expect_error(bioLeak:::.guard_fit(
      X, y = y,
      steps = list(fs = list(method = "lasso")),
      task = "binomial"
    ), "glmnet")
  } else {
    fit <- bioLeak:::.guard_fit(
      X, y = y,
      steps = list(impute = list(method = "median"),
                   normalize = list(method = "none"),
                   filter = list(var_thresh = 0, iqr_thresh = 0),
                   fs = list(method = "lasso")),
      task = "binomial"
    )
    expect_true(fit$state$fs$method == "lasso")
  }
})

test_that("guard_fit coerces character columns to numeric or factor", {
  X <- data.frame(a = c("1", "2", "3"), b = c("x", "y", "x"),
                  stringsAsFactors = FALSE)
  fit <- bioLeak:::.guard_fit(
    X, y = c(0, 1, 0),
    steps = list(impute = list(method = "median"),
                 normalize = list(method = "none"),
                 filter = list(var_thresh = 0, iqr_thresh = 0),
                 fs = list(method = "none")),
    task = "binomial"
  )
  out <- fit$transform(X)
  expect_true(ncol(out) >= 2)
})

test_that("predict_guard requires GuardFit objects", {
  X <- data.frame(a = 1:3, b = 2:4)
  fit <- bioLeak:::.guard_fit(X, y = c(1, 0, 1),
                              steps = list(impute = list(method = "median")),
                              task = "binomial")
  out <- predict_guard(fit, X)
  expect_true(is.data.frame(out))
  expect_error(predict_guard(list(), X))
})
