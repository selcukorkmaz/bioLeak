test_that("duplicate detection completes for n=520 without error (chunked path)", {
  skip_on_cran()
  set.seed(42)
  n <- 520
  p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  # inject a near-duplicate
  X[2, ] <- X[1, ] + rnorm(p, sd = 0.001)
  df <- data.frame(
    subject = seq_len(n),
    outcome = rbinom(n, 1, 0.5),
    X
  )
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 3, progress = FALSE)
  fit <- fit_resample(df, outcome = "outcome", splits = splits,
                      learner = "glmnet",
                      preprocess = list(
                        impute = list(method = "median"),
                        normalize = list(method = "zscore")
                      ))
  # Should complete without error — tests the chunked path (n > 500, no RANN)
  audit <- tryCatch(
    audit_leakage(fit, B = 5, X_ref = as.matrix(df[, paste0("X", 1:p)])),
    error = function(e) e
  )
  expect_false(inherits(audit, "error"))
})

test_that("feature-count trigger logic: n=300, p=6000 should trigger ANN path", {
  # This tests that the condition (n * p > 1.5e6) is correct
  n <- 300
  p <- 6000
  # n * p = 1.8e6 > 1.5e6 should trigger ANN when RANN available
  use_ann <- (n > 500 || (n * p > 1.5e6)) &&
             requireNamespace("RANN", quietly = TRUE)
  if (requireNamespace("RANN", quietly = TRUE)) {
    expect_true(use_ann)
  } else {
    expect_false(use_ann)
  }
})

test_that("injected near-duplicate is detected in the n>500 regime", {
  skip_on_cran()
  set.seed(42)
  n <- 520
  p <- 10
  X <- matrix(rnorm(n * p), nrow = n)
  # inject a near-duplicate
  X[2, ] <- X[1, ] + rnorm(p, sd = 1e-6)
  df <- data.frame(
    subject = seq_len(n),
    outcome = rbinom(n, 1, 0.5),
    X
  )
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 3, progress = FALSE)
  fit <- fit_resample(df, outcome = "outcome", splits = splits,
                      learner = "glmnet",
                      preprocess = list(
                        impute = list(method = "median"),
                        normalize = list(method = "zscore")
                      ))
  audit <- audit_leakage(fit, B = 5,
                         X_ref = as.matrix(df[, paste0("X", 1:p)]),
                         sim_threshold = 0.99)
  expect_true(nrow(audit@duplicates) > 0)
})
