test_that("benchmark_leakage_suite returns scenario matrix and summary", {
  learner <- NULL
  if (requireNamespace("glmnet", quietly = TRUE)) {
    learner <- "glmnet"
  } else if (requireNamespace("ranger", quietly = TRUE)) {
    learner <- "ranger"
  } else {
    skip("benchmark_leakage_suite requires glmnet or ranger.")
  }

  out <- benchmark_leakage_suite(
    modalities = "omics",
    leakages = c("none", "subject_overlap"),
    modes = "subject_grouped",
    learner = learner,
    seeds = 1,
    B = 2,
    alpha = 0.1,
    parallel = FALSE
  )

  expect_true(is.data.frame(out))
  expect_true(nrow(out) >= 1)
  expect_true(all(c("modality", "leakage", "mode", "seed", "metric_obs", "gap", "p_value", "detected") %in% names(out)))
  summ <- attr(out, "summary")
  expect_true(is.data.frame(summ))
  expect_true(all(c("modality", "leakage", "mode", "detection_rate") %in% names(summ)))
})

test_that("benchmark_leakage_suite validates modalities", {
  expect_error(
    benchmark_leakage_suite(modalities = "unknown_modality", seeds = 1, B = 1, parallel = FALSE),
    "Unsupported modalities"
  )
})
