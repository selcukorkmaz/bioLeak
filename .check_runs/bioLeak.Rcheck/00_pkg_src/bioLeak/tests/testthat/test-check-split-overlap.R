test_that("check_split_overlap passes for correct subject_grouped splits", {
  set.seed(1)
  df <- data.frame(
    subject = rep(1:10, each = 2),
    outcome = rbinom(20, 1, 0.5),
    x1 = rnorm(20)
  )
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 3, progress = FALSE)
  result <- check_split_overlap(splits)
  expect_true(all(result$pass))
})

test_that("check_split_overlap detects injected overlap", {
  # Manually construct a splits object with deliberate subject overlap
  indices <- list(
    list(train = 1:8, test = 7:10, fold = 1L, repeat_id = 1L)
  )
  cd <- data.frame(subject = c(1,1,2,2,3,3,4,4,5,5))
  splits_bad <- bioLeak:::LeakSplits(
    mode = "subject_grouped",
    indices = indices,
    info = list(v = 1L, repeats = 1L, seed = 1L, mode = "subject_grouped",
                group = "subject", batch = NULL, study = NULL, time = NULL,
                primary_axis = NULL, secondary_axis = NULL,
                stratify = FALSE, nested = FALSE,
                horizon = 0, purge = 0, embargo = 0,
                summary = data.frame(fold=1L, repeat_id=1L, train_n=8L, test_n=4L),
                hash = NA_character_, inner = NULL,
                compact = FALSE, fold_assignments = NULL,
                coldata = cd)
  )
  expect_error(
    check_split_overlap(splits_bad, coldata = cd, cols = "subject"),
    "Overlap detected"
  )
})

test_that("check_split_overlap stop_on_fail=FALSE returns warning not error", {
  indices <- list(
    list(train = 1:8, test = 7:10, fold = 1L, repeat_id = 1L)
  )
  cd <- data.frame(subject = c(1,1,2,2,3,3,4,4,5,5))
  splits_bad <- bioLeak:::LeakSplits(
    mode = "subject_grouped",
    indices = indices,
    info = list(v = 1L, repeats = 1L, seed = 1L, mode = "subject_grouped",
                group = "subject", batch = NULL, study = NULL, time = NULL,
                primary_axis = NULL, secondary_axis = NULL,
                stratify = FALSE, nested = FALSE,
                horizon = 0, purge = 0, embargo = 0,
                summary = data.frame(fold=1L, repeat_id=1L, train_n=8L, test_n=4L),
                hash = NA_character_, inner = NULL,
                compact = FALSE, fold_assignments = NULL,
                coldata = cd)
  )
  expect_warning(
    result <- check_split_overlap(splits_bad, coldata = cd, cols = "subject",
                                  stop_on_fail = FALSE),
    "Overlap detected"
  )
  expect_false(all(result$pass))
})

test_that("check_split_overlap works for combined mode and checks both axes", {
  set.seed(42)
  df <- data.frame(
    subject = rep(1:20, each = 3),
    batch = rep(rep(c("A","B","C","D"), each=5), 3),
    outcome = rbinom(60, 1, 0.5),
    x1 = rnorm(60)
  )
  splits <- suppressWarnings(make_split_plan(
    df, outcome = "outcome",
    mode = "combined",
    primary_axis = list(type = "subject", col = "subject"),
    secondary_axis = list(type = "batch", col = "batch"),
    v = 3, progress = FALSE
  ))
  result <- check_split_overlap(splits)
  # Both subject and batch columns should show no overlap
  expect_true(all(result$pass))
  expect_true(all(c("subject", "batch") %in% result$col))
})

test_that("check_split_overlap returns a data.frame with expected columns", {
  set.seed(1)
  df <- data.frame(
    subject = rep(1:10, each = 2),
    outcome = rbinom(20, 1, 0.5),
    x1 = rnorm(20)
  )
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 3, progress = FALSE)
  result <- check_split_overlap(splits)
  expect_true(is.data.frame(result))
  expect_true(all(c("fold", "repeat_id", "col", "n_overlap", "pass") %in% names(result)))
})

test_that("LeakAudit info slot carries provenance", {
  set.seed(42)
  df <- data.frame(
    subject = rep(1:30, each = 3),
    outcome = rep(c(0, 1), length.out = 90),
    x1 = rnorm(90),
    x2 = rnorm(90)
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
  audit <- audit_leakage(fit, B = 5)
  expect_true("provenance" %in% names(audit@info))
  expect_true(is.list(audit@info$provenance))
})
