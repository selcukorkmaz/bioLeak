test_that("cindex_pairwise computes concordance correctly", {
  pred <- c(0.1, 0.2, 0.3)
  truth <- c(1, 2, 3)
  expect_equal(bioLeak:::.cindex_pairwise(pred, truth), 1)

  pred2 <- c(0.3, 0.2, 0.1)
  expect_equal(bioLeak:::.cindex_pairwise(pred2, truth), 0)

  expect_true(is.na(bioLeak:::.cindex_pairwise(1, 1)))
  expect_true(is.na(bioLeak:::.cindex_pairwise(c(0.1, 0.2), c(1, 1))))
})
