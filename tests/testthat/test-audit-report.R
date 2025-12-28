test_that("audit_report renders an HTML report", {
  skip_if_not_installed("rmarkdown")
  if (!rmarkdown::pandoc_available()) {
    skip("pandoc not available")
  }

  set.seed(3)
  df <- data.frame(
    subject = rep(1:6, each = 2),
    outcome = rep(c(0, 1), 6),
    x1 = rnorm(12),
    x2 = rnorm(12)
  )
  splits <- make_splits(df, outcome = "outcome",
                        mode = "subject_grouped", group = "subject",
                        v = 3, seed = 1, progress = FALSE)

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
    metrics = "auc",
    refit = FALSE,
    seed = 1
  )
  audit <- audit_leakage(fit, metric = "auc", B = 5, perm_stratify = FALSE)

  out_file <- audit_report(audit, output_dir = tempdir(),
                           output_file = "bioLeak_audit_report_test.html",
                           quiet = TRUE, open = FALSE)
  expect_true(file.exists(out_file))
})
