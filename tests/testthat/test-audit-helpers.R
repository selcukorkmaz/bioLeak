test_that("audit helper utilities behave as expected", {
  A <- matrix(c(1, 0, 0, 1), nrow = 2, byrow = TRUE)
  S <- bioLeak:::.cosine_sim_block(A)
  expect_equal(dim(S), c(2, 2))
  expect_equal(diag(S), c(1, 1))

  B <- matrix(c(3, 4, 0, 0), nrow = 2, byrow = TRUE)
  Bn <- bioLeak:::.row_l2_normalize(B)
  expect_equal(round(rowSums(Bn * Bn), 6), c(1, 0))

  C <- matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE)
  Cc <- bioLeak:::.row_center(C)
  expect_equal(rowMeans(Cc), c(0, 0))
})

test_that("chi-squared association returns NA for small or sparse tables", {
  tab <- matrix(c(1, 0), nrow = 1)
  res <- bioLeak:::.chisq_assoc(tab)
  expect_true(all(is.na(unlist(res))))

  tab2 <- matrix(c(5, 2, 1, 3), nrow = 2)
  res2 <- bioLeak:::.chisq_assoc(tab2)
  expect_true(is.finite(res2$stat))
  expect_true(is.finite(res2$pval))
})

test_that("metric_value computes supported metrics", {
  truth <- factor(c(0, 1, 0, 1))
  pred <- c(0.1, 0.9, 0.2, 0.8)
  rmse <- bioLeak:::.metric_value("rmse", "gaussian", c(1, 2, 3), c(1, 2, 4))
  expect_equal(rmse, sqrt(1 / 3))

  auc <- bioLeak:::.metric_value("auc", "binomial", truth, pred)
  expect_true(is.finite(auc))

  pr <- bioLeak:::.metric_value("pr_auc", "binomial", truth, pred)
  if (requireNamespace("PRROC", quietly = TRUE)) {
    expect_true(is.finite(pr))
  } else {
    expect_true(is.na(pr))
  }

  cidx <- bioLeak:::.metric_value("cindex", "gaussian", c(1, 2, 3), c(0.2, 0.3, 0.1))
  expect_true(is.finite(cidx))
  expect_true(is.na(bioLeak:::.metric_value("unknown", "gaussian", 1:3, 1:3)))
})

test_that("coerce_truth_like preserves types", {
  expect_true(is.factor(bioLeak:::.coerce_truth_like(factor(c("a", "b")), c("b", "a"))))
  expect_true(is.logical(bioLeak:::.coerce_truth_like(c(TRUE, FALSE), c(1, 0))))
  expect_true(is.numeric(bioLeak:::.coerce_truth_like(c(1, 2), c("3", "4"))))
})

