test_that("time-series block permutation utilities return valid indices", {
  set.seed(1)
  idx <- 1:10
  out1 <- bioLeak:::.circular_block_permute(idx, block_len = 3)
  expect_equal(length(out1), length(idx))
  expect_true(all(out1 %in% idx))

  set.seed(2)
  out2 <- bioLeak:::.stationary_bootstrap(idx, mean_block = 3)
  expect_equal(length(out2), length(idx))
  expect_true(all(out2 %in% idx))
})

test_that("time-series block permutation utilities validate inputs", {
  expect_error(bioLeak:::.circular_block_permute(integer(0), block_len = 1))
  expect_error(bioLeak:::.circular_block_permute(1:5, block_len = 0))
  expect_error(bioLeak:::.stationary_bootstrap(integer(0), mean_block = 1))
  expect_error(bioLeak:::.stationary_bootstrap(1:5, mean_block = 0))
})
