test_that("compact splits are resolved during fitting", {
  set.seed(1)
  df <- data.frame(
    outcome = rep(c(0, 1), times = 10),
    subject = rep(1:10, each = 2),
    x1 = rnorm(20),
    x2 = rnorm(20)
  )
  rownames(df) <- paste0("S", seq_len(nrow(df)))

  splits <- make_splits(df,
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
        stats::glm(y ~ ., data = as.data.frame(x),
                   family = stats::binomial(), weights = weights)
      },
      predict = function(object, newdata, task, ...) {
        as.numeric(stats::predict(object, newdata = as.data.frame(newdata), type = "response"))
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
