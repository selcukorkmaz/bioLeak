# tests/testthat/test-strict-mode.R
# Tests for options(bioLeak.strict = TRUE) strict leakage mode

test_that("strict mode forces bioLeak:::.bio_validation_mode() to return 'error'", {
  withr::local_options(list(bioLeak.strict = TRUE))
  expect_equal(bioLeak:::.bio_validation_mode(), "error")
  expect_equal(bioLeak:::.bio_validation_mode("warn"), "error")
  expect_equal(bioLeak:::.bio_validation_mode("off"), "error")
})

test_that("strict mode warns when seed is NULL", {
  withr::local_options(list(bioLeak.strict = TRUE))
  expect_warning(
    bioLeak:::.bio_strict_checks(context = "test", seed = NULL),
    "\\[strict\\]:.*seed"
  )
  expect_warning(
    bioLeak:::.bio_strict_checks(context = "test", seed = NA),
    "\\[strict\\]:.*seed"
  )
})

test_that("strict mode warns on nested + combined", {
  withr::local_options(list(bioLeak.strict = TRUE))
  expect_warning(
    bioLeak:::.bio_strict_checks(context = "test", seed = 1, nested = TRUE, mode = "combined"),
    "\\[strict\\]:.*nested.*combined"
  )
})

test_that("strict mode auto-runs check_split_overlap (correct splits succeed)", {
  withr::local_options(list(bioLeak.strict = TRUE))
  set.seed(1)
  df <- data.frame(
    subject = rep(1:10, each = 2),
    outcome = rbinom(20, 1, 0.5),
    x1 = rnorm(20)
  )
  # Should succeed without error — splits are correct
  expect_no_error(
    make_split_plan(df, outcome = "outcome",
                    mode = "subject_grouped", group = "subject",
                    v = 5, seed = 1, progress = FALSE)
  )
})

test_that("strict mode rejects trained recipes in fit_resample", {
  skip_if_not_installed("recipes")
  withr::local_options(list(bioLeak.strict = TRUE))

  set.seed(1)
  df <- data.frame(
    subject = rep(1:6, each = 2),
    outcome = factor(rep(c(0, 1), each = 6)),
    x1 = rnorm(12),
    x2 = rnorm(12)
  )
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 3, seed = 1, progress = FALSE)

  rec <- recipes::recipe(outcome ~ x1 + x2, data = df)
  trained_rec <- recipes::prep(rec, training = df)

  custom <- list(
    glm = list(
      fit = function(x, y, task, weights, ...) {
        stats::glm(y ~ ., data = as.data.frame(x),
                   family = stats::binomial(), weights = weights)
      },
      predict = function(object, newdata, task, ...) {
        as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
                                  type = "response"))
      }
    )
  )

  # Strict mode forces validation_mode to "error", so trained recipe should error

  expect_error(
    fit_resample(df, outcome = "outcome", splits = splits,
                 learner = "glm", custom_learners = custom,
                 preprocess = trained_rec, metrics = "auc",
                 refit = FALSE, seed = 1),
    "trained recipe"
  )
})
