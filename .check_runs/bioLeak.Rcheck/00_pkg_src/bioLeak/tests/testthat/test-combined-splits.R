test_that("combined mode produces valid splits with no subject overlap AND no batch overlap", {
  set.seed(42)
  df <- data.frame(
    subject = rep(1:20, each = 3),
    batch = rep(rep(c("A", "B", "C", "D"), each = 5), 3),
    outcome = rep(c(0, 1), length.out = 60),
    x1 = rnorm(60),
    x2 = rnorm(60)
  )
  splits <- suppressWarnings(make_split_plan(
    df, outcome = "outcome",
    mode = "combined",
    primary_axis = list(type = "subject", col = "subject"),
    secondary_axis = list(type = "batch", col = "batch"),
    v = 3, progress = FALSE
  ))
  expect_s4_class(splits, "LeakSplits")
  expect_equal(splits@mode, "combined")

  # Check no subject overlap between train and test
  for (fold in splits@indices) {
    if (is.null(fold$train) || is.null(fold$test)) next
    train_subjects <- unique(df$subject[fold$train])
    test_subjects <- unique(df$subject[fold$test])
    expect_length(intersect(train_subjects, test_subjects), 0)

    # Check no batch overlap between train and test
    train_batches <- unique(df$batch[fold$train])
    test_batches <- unique(df$batch[fold$test])
    expect_length(intersect(train_batches, test_batches), 0)
  }
})

test_that("error when primary_axis or secondary_axis missing", {
  df <- data.frame(subject = 1:10, batch = rep("A", 10),
                   outcome = rbinom(10, 1, 0.5), x1 = rnorm(10))
  expect_error(
    make_split_plan(df, outcome = "outcome", mode = "combined",
                    primary_axis = list(type = "subject", col = "subject"),
                    v = 3, progress = FALSE),
    "secondary_axis"
  )
  expect_error(
    make_split_plan(df, outcome = "outcome", mode = "combined",
                    secondary_axis = list(type = "batch", col = "batch"),
                    v = 3, progress = FALSE),
    "primary_axis"
  )
})

test_that("error when axis list is malformed", {
  df <- data.frame(subject = 1:10, batch = rep("A", 10),
                   outcome = rbinom(10, 1, 0.5), x1 = rnorm(10))
  # not a list
  expect_error(
    make_split_plan(df, outcome = "outcome", mode = "combined",
                    primary_axis = "subject",
                    secondary_axis = list(type = "batch", col = "batch"),
                    v = 3, progress = FALSE),
    "primary_axis"
  )
  # missing type
  expect_error(
    make_split_plan(df, outcome = "outcome", mode = "combined",
                    primary_axis = list(col = "subject"),
                    secondary_axis = list(type = "batch", col = "batch"),
                    v = 3, progress = FALSE),
    "primary_axis"
  )
})

test_that("compact = TRUE produces fold_assignments", {
  set.seed(42)
  df <- data.frame(
    subject = rep(1:20, each = 3),
    batch = rep(rep(c("A", "B", "C", "D"), each = 5), 3),
    outcome = rbinom(60, 1, 0.5),
    x1 = rnorm(60),
    x2 = rnorm(60)
  )
  splits <- suppressWarnings(make_split_plan(
    df, outcome = "outcome",
    mode = "combined",
    primary_axis = list(type = "subject", col = "subject"),
    secondary_axis = list(type = "batch", col = "batch"),
    v = 3, compact = TRUE, progress = FALSE
  ))
  expect_true(isTRUE(splits@info$compact))
  expect_true(length(splits@info$fold_assignments) > 0)
})

test_that("combined splits pass through to fit_resample successfully", {
  set.seed(42)
  df <- data.frame(
    subject = rep(1:20, each = 3),
    batch = rep(rep(c("A", "B", "C", "D"), each = 5), 3),
    outcome = rbinom(60, 1, 0.5),
    x1 = rnorm(60),
    x2 = rnorm(60)
  )
  splits <- suppressWarnings(make_split_plan(
    df, outcome = "outcome",
    mode = "combined",
    primary_axis = list(type = "subject", col = "subject"),
    secondary_axis = list(type = "batch", col = "batch"),
    v = 3, progress = FALSE
  ))
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glmnet",
                            preprocess = list(
                              impute = list(method = "median"),
                              normalize = list(method = "zscore")
                            ))
  expect_s4_class(fit, "LeakFit")
  expect_gt(nrow(fit@metrics), 0)
})

test_that("stratification works with combined mode", {
  set.seed(42)
  df <- data.frame(
    subject = rep(1:20, each = 3),
    batch = rep(rep(c("A", "B", "C", "D"), each = 5), 3),
    outcome = rep(c(0, 1), 30),
    x1 = rnorm(60),
    x2 = rnorm(60)
  )
  splits <- suppressWarnings(make_split_plan(
    df, outcome = "outcome",
    mode = "combined",
    primary_axis = list(type = "subject", col = "subject"),
    secondary_axis = list(type = "batch", col = "batch"),
    v = 3, stratify = TRUE, progress = FALSE
  ))
  expect_s4_class(splits, "LeakSplits")
  expect_gt(length(splits@indices), 0)
})
