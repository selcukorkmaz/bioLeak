test_that("LeakSplits validity and show method work", {
  expect_error(LeakSplits(mode = "x", indices = list(), info = list()),
               "indices list cannot be empty")

  df <- make_class_df(8)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 2, seed = 1)
  out <- capture.output(show(splits))
  expect_true(any(grepl("LeakSplits object", out)))
})

test_that("LeakFit and LeakAudit validity checks enforce class constraints", {
  splits <- make_splits_quiet(make_class_df(8), outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 2, seed = 1)
  expect_error(LeakFit(splits = list(), metrics = data.frame(),
                       metric_summary = data.frame(), audit = data.frame(),
                       predictions = list(), preprocess = list(), learners = list(),
                       outcome = "outcome", task = "binomial",
                       feature_names = character(), info = list()),
               "splits must be a LeakSplits")

  fit <- LeakFit(splits = splits, metrics = data.frame(),
                 metric_summary = data.frame(), audit = data.frame(),
                 predictions = list(), preprocess = list(), learners = list(),
                 outcome = "outcome", task = "binomial",
                 feature_names = character(), info = list())
  expect_error(LeakAudit(fit = list(), permutation_gap = data.frame(),
                         perm_values = numeric(), batch_assoc = data.frame(),
                         target_assoc = data.frame(), duplicates = data.frame(),
                         trail = list(), info = list()),
               "fit must be a LeakFit")
})

test_that("summary methods emit expected headers", {
  df <- make_class_df(10)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 2, seed = 1)
  custom <- make_custom_learners()
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = custom,
                            metrics = "auc", refit = FALSE, seed = 1)
  audit <- audit_leakage(fit, metric = "auc", B = 3)

  fit_out <- capture.output(summary(fit))
  expect_true(any(grepl("Model Fit Summary", fit_out)))

  aud_out <- capture.output(summary(audit))
  expect_true(any(grepl("Leakage Audit Summary", aud_out)))
})

test_that("GuardFit print and summary methods work", {
  X <- data.frame(a = c(1, 2, NA), b = c(3, 4, 5))
  fit <- bioLeak:::.guard_fit(X, y = c(1, 2, 3),
                              steps = list(impute = list(method = "median")),
                              task = "gaussian")
  out <- capture.output(print(fit))
  expect_true(any(grepl("Guarded preprocessing pipeline", out)))

  sm <- summary(fit)
  expect_true(inherits(sm, "summary.GuardFit"))
  out2 <- capture.output(print(sm))
  expect_true(any(grepl("GuardFit summary", out2)))
})
