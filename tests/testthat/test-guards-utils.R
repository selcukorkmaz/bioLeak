test_that("basic guard utility helpers behave as expected", {
  expect_equal(bioLeak:::`%||%`(1, 2), 1)
  expect_equal(bioLeak:::`%||%`(NULL, 2), 2)

  expect_true(bioLeak:::.guard_is_binary(factor(c("a", "b"))))
  expect_true(bioLeak:::.guard_is_binary(c(0, 1, 1)))
  expect_false(bioLeak:::.guard_is_binary(c(0, 1, 2)))

  expect_equal(bioLeak:::.guard_majority(c(1, 1, 2, NA)), 1)
  expect_true(is.na(bioLeak:::.guard_majority(c(NA, NA))))
})

test_that("winsorization guards numeric inputs and ignores non-numeric", {
  x <- c(-100, 0, 1, 2, 100)
  out <- bioLeak:::.guard_mad_winsorize(x, k = 1)
  expect_true(min(out) >= 0)
  expect_true(max(out) <= 2)

  y <- c("a", "b")
  expect_equal(bioLeak:::.guard_mad_winsorize(y), y)
})

test_that("column alignment and NA filling behave consistently", {
  X <- data.frame(a = c(1, NA), b = c(2, 3))
  aligned <- bioLeak:::.guard_align_cols(X, c("b", "a", "c"))
  expect_equal(names(aligned), c("b", "a", "c"))
  expect_true(all(is.na(aligned$c)))

  fill <- list(a = 5, b = 6, c = 7)
  filled <- bioLeak:::.guard_fill_na(aligned, fill)
  expect_equal(filled$a[2], 5)
  expect_equal(filled$c[1], 7)
})

test_that("kNN search and imputation functions return valid shapes", {
  train <- data.frame(a = c(1, 2, 3), b = c(4, 5, 6))
  query <- data.frame(a = c(2, 3), b = c(5, 6))
  idx <- bioLeak:::.guard_knn_search(train, query, k = 2)
  expect_equal(dim(idx), c(2, 2))
  expect_true(all(idx >= 1 & idx <= nrow(train)))

  empty_idx <- bioLeak:::.guard_knn_search(train[0, ], query, k = 2)
  expect_equal(dim(empty_idx), c(2, 0))

  X <- data.frame(a = c(1, NA, 3), b = c(4, 5, NA))
  knn_fit <- bioLeak:::.guard_knn_impute_train(X, k = 1)
  expect_equal(dim(knn_fit$imp), dim(X))
  expect_true(all(is.finite(as.matrix(knn_fit$imp))))

  Xnew <- data.frame(a = c(NA, 2), b = c(5, NA))
  imputed <- bioLeak:::.guard_knn_impute_new(Xnew, knn_fit$imp, k = 1,
                                             med = knn_fit$med, scale = knn_fit$scale)
  expect_equal(dim(imputed), dim(Xnew))
  expect_true(all(is.finite(as.matrix(imputed))))
})

test_that("random-forest imputation guard handles missing dependency", {
  X <- data.frame(a = c(1, NA, 3), b = c(4, 5, NA))
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    expect_error(bioLeak:::.guard_rf_impute_train(X, maxiter = 1, ntree = 10),
                 "randomForest")
  } else {
    rf_fit <- suppressWarnings(bioLeak:::.guard_rf_impute_train(X, maxiter = 1,
                                                                ntree = 10, seed = 1))
    expect_equal(dim(rf_fit$imp), dim(X))
    out <- suppressWarnings(bioLeak:::.guard_rf_impute_new(X, rf_fit$models, rf_fit$med))
    expect_equal(dim(out), dim(X))
  }
})

test_that("guard hashing and level normalization behave predictably", {
  h <- bioLeak:::.guard_hash(list(a = 1))
  expect_true(is.character(h))
  expect_true(nchar(h) > 0)

  df <- data.frame(flag = c(TRUE, FALSE), group = c("a", "a"))
  out <- bioLeak:::.guard_ensure_levels(df)
  expect_true(is.factor(out$data$flag))
  expect_true(is.factor(out$data$group))
  expect_true(length(out$levels$group) >= 2)

  out2 <- bioLeak:::.guard_ensure_levels(data.frame(flag = c(TRUE, FALSE)),
                                         levels_map = out$levels)
  expect_true("group" %in% names(out2$data))
  expect_true(is.factor(out2$data$group))
})
