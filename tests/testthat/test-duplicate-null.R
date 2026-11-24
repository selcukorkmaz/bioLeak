# L3.4 (JSS review): expected number of chance pairs and enrichment ratio.

.l34_fit <- function(n = 200, p = 8, seed = 1) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, dimnames = list(NULL, paste0("f", seq_len(p))))
  df <- data.frame(subject = paste0("s", seq_len(n)),
                   y = factor(sample(c("a", "b"), n, TRUE)), X)
  sp <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject",
                        v = 5, seed = 1, progress = FALSE)
  fit <- NULL
  utils::capture.output(fit <- suppressWarnings(fit_resample(
    df[, c("subject", "y", "f1", "f2")], outcome = "y", splits = sp, learner = "glm",
    custom_learners = make_custom_learners(), metrics = "auc", refit = FALSE, seed = 1)))
  list(fit = fit, X = X)
}

.l34_audit <- function(fit, X, ...) {
  a <- NULL
  utils::capture.output(a <- suppressWarnings(audit_leakage(
    fit, metric = "auc", B = 5, X_ref = X, target_scan = FALSE,
    target_scan_multivariate = FALSE, ...)))
  a
}

test_that("duplicate_null reports observed, expected and enrichment", {
  ff <- .l34_fit()
  a <- .l34_audit(ff$fit, ff$X)
  dn <- a@info$duplicate_null
  expect_s3_class(dn, "data.frame")
  expect_named(dn, c("scope", "observed", "expected", "enrichment", "null_rows", "reps"))
  expect_identical(dn$scope, "train_test")
  expect_equal(dn$observed, nrow(a@duplicates))
  expect_true(is.finite(dn$expected) && dn$expected >= 0)
})

test_that("copied rows are enriched over the chance count", {
  ff <- .l34_fit(seed = 2)
  X <- ff$X
  X[1:15, ] <- X[101:115, ] + 1e-4
  a <- .l34_audit(ff$fit, X, duplicate_null_reps = 5)
  dn <- a@info$duplicate_null
  expect_gte(dn$observed, 8)
  expect_lt(dn$expected, 1)
  expect_true(is.na(dn$enrichment) || dn$enrichment > 8)
  out <- capture.output(summary(a))
  expect_true(any(grepl("Expected by chance", out)))
  expect_true(any(grepl("Flagged leakage mechanisms: duplicate_overlap", out)))
})

test_that("duplicate_null_reps = 0 skips the estimate; bad values error", {
  ff <- .l34_fit(seed = 3)
  a <- .l34_audit(ff$fit, ff$X, duplicate_null_reps = 0)
  expect_null(a@info$duplicate_null)
  expect_error(.l34_audit(ff$fit, ff$X, duplicate_null_reps = -1), "duplicate_null_reps")
})

test_that("the null draws leave the user's RNG state unchanged", {
  ff <- .l34_fit(seed = 4)
  set.seed(123); before <- .Random.seed
  a1 <- .l34_audit(ff$fit, ff$X, duplicate_null_reps = 3)
  set.seed(123)
  a0 <- .l34_audit(ff$fit, ff$X, duplicate_null_reps = 0)
  after_rng <- stats::runif(1)
  set.seed(123)
  a1b <- .l34_audit(ff$fit, ff$X, duplicate_null_reps = 3)
  expect_identical(stats::runif(1), after_rng)
  expect_equal(a1@info$duplicate_null, a1b@info$duplicate_null)
})
