test_that(".bio_stop raises error inheriting from bioLeak_error", {
  expect_error(
    bioLeak:::.bio_stop("test error", "bioLeak_input_error"),
    "test error"
  )
  err <- tryCatch(
    bioLeak:::.bio_stop("test error", "bioLeak_input_error"),
    error = function(e) e
  )
  expect_s3_class(err, "bioLeak_error")
  expect_s3_class(err, "bioLeak_input_error")
  expect_s3_class(err, "error")
  expect_s3_class(err, "condition")
})

test_that(".bio_warn raises warning inheriting from bioLeak_warning", {
  expect_warning(
    bioLeak:::.bio_warn("test warning", "bioLeak_fold_warning"),
    "test warning"
  )
  w <- tryCatch(
    bioLeak:::.bio_warn("test warning", "bioLeak_fold_warning"),
    warning = function(w) w
  )
  expect_s3_class(w, "bioLeak_warning")
  expect_s3_class(w, "bioLeak_fold_warning")
  expect_s3_class(w, "warning")
  expect_s3_class(w, "condition")
})

test_that("bioLeak_column_error is catchable from make_split_plan", {
  df <- data.frame(x1 = rnorm(10), outcome = rbinom(10, 1, 0.5))
  caught <- tryCatch(
    make_split_plan(df, outcome = "outcome",
                    mode = "subject_grouped", group = "nonexistent",
                    v = 3, progress = FALSE),
    bioLeak_column_error = function(e) "caught_column"
  )
  expect_equal(caught, "caught_column")
})

test_that("bioLeak_input_error is catchable from make_split_plan", {
  df <- data.frame(x1 = rnorm(10), outcome = rbinom(10, 1, 0.5))
  caught <- tryCatch(
    make_split_plan(df, outcome = "outcome",
                    mode = "subject_grouped", group = NULL,
                    v = 3, progress = FALSE),
    bioLeak_input_error = function(e) "caught_input"
  )
  expect_equal(caught, "caught_input")
})

test_that("check_split_overlap raises bioLeak_overlap_error", {
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
  caught <- tryCatch(
    check_split_overlap(splits_bad, coldata = cd, cols = "subject"),
    bioLeak_overlap_error = function(e) "caught_overlap"
  )
  expect_equal(caught, "caught_overlap")
})

test_that("fit_resample with bad input raises bioLeak_input_error", {
  caught <- tryCatch(
    fit_resample(data.frame(x = 1:5), outcome = "x",
                 splits = "not_a_split",
                 learner = "glmnet"),
    bioLeak_input_error = function(e) "caught_input"
  )
  expect_equal(caught, "caught_input")
})

test_that("missing dependency raises bioLeak_dependency_error", {
  caught <- tryCatch(
    bioLeak:::.bio_stop("Install 'foo' package", "bioLeak_dependency_error"),
    bioLeak_dependency_error = function(e) "caught_dep"
  )
  expect_equal(caught, "caught_dep")
})

test_that("guard_to_recipe fallback raises bioLeak_fallback_warning", {
  skip_if_not_installed("recipes")
  df <- data.frame(outcome = factor(c(0,1,0,1)), x1 = rnorm(4))
  w <- tryCatch(
    guard_to_recipe(
      steps = list(impute = list(method = "missForest")),
      formula = outcome ~ .,
      training_data = df
    ),
    bioLeak_fallback_warning = function(w) w
  )
  expect_s3_class(w, "bioLeak_fallback_warning")
})
