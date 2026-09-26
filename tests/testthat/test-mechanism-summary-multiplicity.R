# Regression tests: the confounding_alignment rule must correct for the
# number of repeats. Before 0.3.8 it took the raw minimum p over repeats, so
# one chance p = 0.038 among 10 repeats flagged an unrelated batch.

batch_rows <- function(pvals, v = 0.3) {
  data.frame(batch_col = "batch", repeat_id = seq_along(pvals),
             stat = 1, df = 4, pval = pvals, cramer_v = v)
}

confounding_row <- function(batch_df) {
  ms <- bioLeak:::.bio_mechanism_summary(
    perm_df = data.frame(), batch_df = batch_df,
    target_df = data.frame(), dup_df = data.frame()
  )
  ms[ms$mechanism_class == "confounding_alignment", ]
}

test_that("a single nominal p < 0.05 across 10 repeats is not flagged", {
  p <- c(0.61, 0.44, 0.83, 0.29, 0.51, 0.72, 0.038, 0.35, 0.92, 0.18)
  row <- confounding_row(batch_rows(p))
  expect_false(row$flagged)
  expect_equal(row$p_value, min(stats::p.adjust(p, "holm")))
  expect_equal(row$p_value, 0.38)
})

test_that("a strong association survives the Holm correction", {
  p <- c(0.001, rep(0.5, 9))
  expect_true(confounding_row(batch_rows(p))$flagged)
})

test_that("the effect-size floor still applies after correction", {
  p <- c(0.0001, rep(0.5, 9))
  expect_false(confounding_row(batch_rows(p, v = 0.05))$flagged)
})

test_that("a single repeat is unaffected by the correction", {
  row <- confounding_row(batch_rows(0.03))
  expect_true(row$flagged)
  expect_equal(row$p_value, 0.03)
})

test_that("NA p-values do not enlarge the family", {
  expect_equal(bioLeak:::.holm_adjust(c(0.01, NA, 0.04)), c(0.02, NA, 0.04))
})

test_that("audit_leakage reports Holm-adjusted batch p-values", {
  set.seed(5)
  df <- make_class_df(60)
  sp <- make_split_plan(df, outcome = "outcome", mode = "subject_grouped",
                        group = "subject", v = 3, repeats = 4, progress = FALSE)
  fit <- fit_resample(df, outcome = "outcome", splits = sp, learner = "glm",
                      custom_learners = make_custom_learners(), metrics = "auc",
                      refit = FALSE, seed = 1)
  aud <- suppressMessages(audit_leakage(fit, metric = "auc", B = 5,
                                        coldata = df, batch_cols = "batch",
                                        target_scan = FALSE,
                                        target_scan_multivariate = FALSE))
  ba <- aud@batch_assoc
  expect_true("pval_adj" %in% names(ba))
  expect_equal(nrow(ba), 4L)
  expect_equal(ba$pval_adj, stats::p.adjust(ba$pval, "holm"))
  ms <- aud@info$mechanism_summary
  expect_equal(ms$p_value[ms$mechanism_class == "confounding_alignment"],
               min(ba$pval_adj))
})
