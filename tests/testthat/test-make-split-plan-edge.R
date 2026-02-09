test_that("make_split_plan validates required columns and time types", {
  df <- make_class_df(6)
  expect_error(make_split_plan_quiet(df, outcome = "outcome",
                                 mode = "subject_grouped", group = NULL,
                                 v = 2, seed = 1),
               "requires a non-NULL 'group'")
  expect_error(make_split_plan_quiet(df, outcome = "outcome",
                                 mode = "batch_blocked", batch = "missing"),
               "Column 'missing' not found")
  expect_error(make_split_plan_quiet(df, outcome = "outcome",
                                 mode = "study_loocv", study = "missing"),
               "Column 'missing' not found")
  expect_error(make_split_plan_quiet(df, outcome = "outcome",
                                 mode = "time_series", time = "missing"),
               "Column 'missing' not found")

  df$time <- as.character(df$time)
  expect_error(make_split_plan_quiet(df, outcome = "outcome",
                                 mode = "time_series", time = "time"),
               "time' column must be numeric")
})

test_that("make_split_plan warns on stratification issues and compact/nested conflicts", {
  df <- make_class_df(10)
  expect_warning_match(
    make_split_plan_quiet(df, outcome = "missing",
                      mode = "subject_grouped", group = "subject",
                      v = 2, stratify = TRUE),
    "Stratification requested"
  )

  df$outcome <- 1
  expect_warning_match(
    make_split_plan_quiet(df, outcome = "outcome",
                      mode = "subject_grouped", group = "subject",
                      v = 2, stratify = TRUE),
    "Stratification ignored"
  )

  splits <- expect_warning_match(
    make_split_plan_quiet(df, outcome = "outcome",
                      mode = "subject_grouped", group = "subject",
                      v = 2, nested = TRUE, compact = TRUE),
    "compact=TRUE is not supported"
  )
  expect_false(splits@info$compact)
})

test_that("make_split_plan supports compact assignments and nested splits", {
  df <- make_class_df(12)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, compact = TRUE, seed = 1)
  expect_true(isTRUE(splits@info$compact))
  expect_true(length(splits@info$fold_assignments) > 0)
  expect_true(all(vapply(splits@indices, function(x) !is.null(x$fold), logical(1))))

  nested <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, nested = TRUE, seed = 1)
  expect_true(is.list(nested@info$inner))
  expect_true(length(nested@info$inner) == length(nested@indices))
})

test_that("make_split_plan errors when no valid folds are produced", {
  df <- make_class_df(4)
  expect_error(make_split_plan_quiet(df, outcome = "outcome",
                                 mode = "time_series", time = "time",
                                 v = 2, horizon = 100, seed = 1),
               "No valid folds generated")
})

test_that("make_split_plan time_series includes tail samples when n not divisible by v", {
  df <- make_class_df(10)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "time_series", time = "time",
                              v = 3, horizon = 0, seed = 1)
  test_idx <- sort(unique(unlist(lapply(splits@indices, `[[`, "test"))))
  expect_true(max(test_idx) == nrow(df))
})

test_that("make_split_plan time_series applies purge and embargo gaps", {
  df <- make_class_df(12)

  base <- make_split_plan_quiet(df, outcome = "outcome",
                           mode = "time_series", time = "time",
                           v = 3, horizon = 0, seed = 1)
  with_purge <- make_split_plan_quiet(df, outcome = "outcome",
                                 mode = "time_series", time = "time",
                                 v = 3, horizon = 0, purge = 2, seed = 1)
  with_embargo <- make_split_plan_quiet(df, outcome = "outcome",
                                   mode = "time_series", time = "time",
                                   v = 3, horizon = 0, embargo = 5, seed = 1)

  fold2_base <- base@indices[[which(vapply(base@indices, function(z) z$fold == 2L, logical(1)))]]
  fold2_purge <- with_purge@indices[[which(vapply(with_purge@indices, function(z) z$fold == 2L, logical(1)))]]
  fold2_embargo <- with_embargo@indices[[which(vapply(with_embargo@indices, function(z) z$fold == 2L, logical(1)))]]

  expect_equal(max(df$time[fold2_base$train]), 4)
  expect_equal(max(df$time[fold2_purge$train]), 3)
  expect_equal(max(df$time[fold2_embargo$train]), 3)
  expect_equal(with_purge@info$purge, 2)
  expect_equal(with_purge@info$embargo, 0)
  expect_equal(with_embargo@info$embargo, 5)
})

test_that("make_split_plan supports matrix and SummarizedExperiment inputs", {
  mat <- matrix(rnorm(20), nrow = 10)
  splits <- make_split_plan(mat, outcome = NULL,
                        mode = "subject_grouped", group = "row_id",
                        v = 2, seed = 1)
  expect_true(inherits(splits, "LeakSplits"))
  expect_true("row_id" %in% names(splits@info$coldata))

  skip_if_not_installed("SummarizedExperiment")
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(expr = matrix(rnorm(40), nrow = 4)),
    colData = data.frame(outcome = rep(0:1, each = 5),
                         subject = rep(1:5, 2))
  )
  splits_se <- make_split_plan(se, outcome = "outcome",
                           mode = "subject_grouped", group = "subject",
                           v = 2, seed = 1, progress = FALSE)
  expect_true(inherits(splits_se, "LeakSplits"))
})
