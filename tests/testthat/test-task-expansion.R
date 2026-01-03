test_that("fit_resample supports multiclass custom learners", {
  df <- make_multiclass_df(15, k = 3)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- list(
    freq = list(
      fit = function(x, y, task, weights, ...) {
        probs <- prop.table(table(y))
        list(probs = probs, levels = names(probs))
      },
      predict = function(object, newdata, task, ...) {
        n <- nrow(newdata)
        prob_mat <- matrix(rep(object$probs, each = n), nrow = n, byrow = FALSE)
        colnames(prob_mat) <- object$levels
        prob_mat
      }
    )
  )
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "freq", custom_learners = custom,
                            metrics = c("accuracy", "macro_f1", "log_loss"),
                            refit = FALSE)
  expect_true(all(c("accuracy", "macro_f1", "log_loss") %in% colnames(fit@metrics)))
})

test_that("fit_resample supports survival tasks with custom learners", {
  skip_if_not_installed("survival")
  n <- 20
  df <- data.frame(
    subject = rep(seq_len(10), each = 2),
    time = rexp(n, rate = 0.1),
    status = rbinom(n, 1, 0.7),
    x1 = rnorm(n),
    x2 = rnorm(n),
    stringsAsFactors = FALSE
  )
  df$surv <- survival::Surv(df$time, df$status)
  splits <- make_split_plan_quiet(df, outcome = "surv",
                              mode = "subject_grouped", group = "subject",
                              v = 2, seed = 1, stratify = FALSE)
  custom <- list(
    cox = list(
      fit = function(x, y, task, weights, ...) {
        df_fit <- data.frame(y = y, x, check.names = FALSE)
        survival::coxph(y ~ ., data = df_fit, weights = weights)
      },
      predict = function(object, newdata, task, ...) {
        as.numeric(stats::predict(object, newdata = as.data.frame(newdata),
                                  type = "lp"))
      }
    )
  )
  fit <- fit_resample_quiet(df, outcome = "surv", splits = splits,
                            learner = "cox", custom_learners = custom,
                            metrics = "cindex", refit = FALSE)
  expect_true(nrow(fit@metrics) > 0)
})
