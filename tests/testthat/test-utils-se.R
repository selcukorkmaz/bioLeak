test_that("SummarizedExperiment helpers behave for data.frame and matrix inputs", {
  df <- make_class_df(6)

  expect_false(bioLeak:::.bio_is_se(df))
  expect_equal(bioLeak:::.bio_get_x(df)$x1, df$x1)
  expect_equal(bioLeak:::.bio_get_y(df, "outcome"), df$outcome)

  expect_error(bioLeak:::.bio_get_y(df, "missing"), "Outcome not found")
  expect_error(bioLeak:::.bio_get_y(df, NULL), "Provide outcome")

  mat <- as.matrix(df[, c("x1", "x2")])
  mat_out <- bioLeak:::.bio_get_x(mat)
  expect_true(is.data.frame(mat_out))
  expect_true(all(vapply(mat_out, is.numeric, logical(1))))
})

test_that("SummarizedExperiment helpers behave for SE inputs", {
  skip_if_not_installed("SummarizedExperiment")
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:9, nrow = 3)),
    colData = data.frame(outcome = c(0, 1, 0))
  )
  expect_true(bioLeak:::.bio_is_se(se))
  x <- bioLeak:::.bio_get_x(se)
  expect_equal(dim(x), c(3, 3))
  expect_equal(bioLeak:::.bio_get_y(se, "outcome"), c(0, 1, 0))
  expect_error(bioLeak:::.bio_get_y(se, "missing"), "Outcome column not in colData")
})

test_that("metadata and hash helpers are stable", {
  df <- make_class_df(6)
  meta <- bioLeak:::.bio_get_meta(df, c("batch", "study"))
  expect_equal(names(meta), c("batch", "study"))
  expect_true(all(vapply(meta, is.null, logical(1))))

  idx <- list(list(train = 1:3, test = 4:6))
  h <- bioLeak:::.bio_hash_indices(idx)
  expect_true(is.character(h))
  expect_true(nchar(h) >= 2)
})

test_that("classification helper predicates behave as expected", {
  expect_true(bioLeak:::.bio_is_binomial(c(0, 1, 0)))
  expect_false(bioLeak:::.bio_is_binomial(c(0, 1, 2)))
  expect_true(bioLeak:::.bio_is_classification(factor(c("a", "b"))))
  expect_true(bioLeak:::.bio_is_regression(c(1.1, 2.2, 3.3)))
  expect_false(bioLeak:::.bio_is_regression(factor(c("a", "b"))))
})
