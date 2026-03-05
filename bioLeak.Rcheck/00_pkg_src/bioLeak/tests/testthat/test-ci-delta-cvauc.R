test_that("se_ci_delta handles empty and zero-variance permutations", {
  out_empty <- bioLeak:::.se_ci_delta(0.1, 0.2, numeric(0))
  expect_true(all(is.na(unlist(out_empty))))

  perm <- rep(0.2, 10)
  out_zero <- bioLeak:::.se_ci_delta(0.1, 0.05, perm)
  expect_equal(out_zero$se, 0.05)
  expect_equal(length(out_zero$ci), 2)
  expect_true(is.finite(out_zero$z))
})

test_that("se_ci_delta handles finite variance permutations", {
  perm <- seq(0.1, 0.5, length.out = 25)
  out <- bioLeak:::.se_ci_delta(0.2, 0.1, perm, level = 0.9)
  expect_true(is.finite(out$se))
  expect_equal(length(out$ci), 2)
})

test_that("psi_auc and cvauc_if compute fold statistics", {
  pos <- c(0.9, 0.8)
  neg <- c(0.1, 0.2)
  psi <- bioLeak:::.psi_auc(pos, neg)
  expect_equal(dim(psi), c(2, 2))
  expect_true(all(psi %in% c(0, 0.5, 1)))

  pred <- list(c(0.9, 0.1, 0.8, 0.2), c(0.7, 0.3))
  truth <- list(factor(c(1, 0, 1, 0)), factor(c(1, 0)))
  res <- bioLeak:::.cvauc_if(pred, truth, weights = c(0.7, 0.3))
  expect_true(is.finite(res$mean_auc))
  expect_true(is.finite(res$se_auc))
  expect_equal(nrow(res$fold_stats), 2)
})

test_that("cvauc_if returns NA when no valid folds", {
  pred <- list(c(0.1, 0.2), c(0.3, 0.4))
  truth <- list(factor(c(1, 1)), factor(c(0, 0)))
  res <- bioLeak:::.cvauc_if(pred, truth)
  expect_true(all(is.na(c(res$mean_auc, res$var_auc, res$se_auc))))
  expect_equal(nrow(res$fold_stats), 0)
})
