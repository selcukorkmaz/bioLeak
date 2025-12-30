test_that("fit_resample supports custom learners", {
  set.seed(1)
  df <- data.frame(
    subject = rep(1:10, each = 2),
    outcome = rep(c(0, 1), times = 10),
    x1 = rnorm(20),
    x2 = rnorm(20)
  )

  splits <- make_splits(
    df,
    outcome = "outcome",
    mode = "subject_grouped",
    group = "subject",
    v = 5,
    stratify = FALSE
  )

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

  expect_s4_class(fit, "LeakFit")
  expect_true(nrow(fit@metrics) > 0)
})
