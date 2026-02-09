test_that("compact splits are resolved during fitting", {
  set.seed(1)
  df <- data.frame(
    outcome = rep(c(0, 1), times = 10),
    subject = rep(1:10, each = 2),
    x1 = rnorm(20),
    x2 = rnorm(20)
  )
  rownames(df) <- paste0("S", seq_len(nrow(df)))

  splits <- make_split_plan(df,
                        outcome = "outcome",
                        mode = "subject_grouped",
                        group = "subject",
                        v = 5,
                        compact = TRUE,
                        progress = FALSE)

  expect_true(isTRUE(splits@info$compact))
  expect_equal(length(splits@info$fold_assignments), 1)
  expect_true(all(!is.na(splits@info$fold_assignments[[1]])))

  custom <- list(
    glm = list(
      fit = function(x, y, task, weights, ...) {
        suppressWarnings(stats::glm(y ~ ., data = as.data.frame(x),
                                    family = stats::binomial(), weights = weights))
      },
      predict = function(object, newdata, task, ...) {
        as.numeric(suppressWarnings(stats::predict(object,
                                                   newdata = as.data.frame(newdata),
                                                   type = "response")))
      }
    )
  )

  fit <- fit_resample(
    df,
    outcome = "outcome",
    splits = splits,
    learner = "glm",
    custom_learners = custom,
    metrics = "accuracy",
    refit = FALSE,
    seed = 1
  )

  pred_ids <- unique(do.call(rbind, fit@predictions)$id)
  expect_true(all(pred_ids %in% rownames(df)))
})

test_that("compact time_series splits keep purge/embargo behavior in as_rsample", {
  skip_if_not_installed("rsample")
  df <- make_class_df(12)
  splits <- make_split_plan_quiet(
    df,
    outcome = "outcome",
    mode = "time_series",
    time = "time",
    v = 3,
    compact = TRUE,
    purge = 2,
    embargo = 5,
    seed = 1
  )
  rs <- as_rsample(splits, data = df)

  expect_equal(nrow(rs), 2)
  expect_equal(max(rsample::analysis(rs$splits[[1]])$time), 3)
  expect_equal(max(rsample::analysis(rs$splits[[2]])$time), 7)
})
