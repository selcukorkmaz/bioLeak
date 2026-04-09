test_that("cv_ci returns correct column names", {
  df <- data.frame(
    fold = 1:5,
    learner = rep("glmnet", 5),
    auc = runif(5, 0.6, 0.9),
    accuracy = runif(5, 0.7, 0.95)
  )
  res <- cv_ci(df)
  expect_true("learner" %in% names(res))
  expect_true(all(c("auc_mean", "auc_sd", "auc_ci_lo", "auc_ci_hi") %in% names(res)))
  expect_true(all(c("accuracy_mean", "accuracy_sd", "accuracy_ci_lo", "accuracy_ci_hi") %in% names(res)))
})

test_that("CI bounds satisfy ci_lo < mean < ci_hi", {
  set.seed(1)
  df <- data.frame(
    fold = 1:10,
    learner = rep("glmnet", 10),
    auc = runif(10, 0.6, 0.9)
  )
  res <- cv_ci(df)
  expect_true(res$auc_ci_lo < res$auc_mean)
  expect_true(res$auc_mean < res$auc_ci_hi)
})

test_that("Nadeau-Bengio CIs are wider than normal CIs", {
  set.seed(1)
  df <- data.frame(
    fold = 1:10,
    learner = rep("glmnet", 10),
    auc = runif(10, 0.6, 0.9)
  )
  res_nb <- cv_ci(df, method = "nadeau_bengio", n_train = 90, n_test = 10)
  res_normal <- cv_ci(df, method = "normal")
  nb_width <- res_nb$auc_ci_hi - res_nb$auc_ci_lo
  normal_width <- res_normal$auc_ci_hi - res_normal$auc_ci_lo
  expect_gt(nb_width, normal_width)
})

test_that("single-fold input returns NA for CIs", {
  df <- data.frame(fold = 1, learner = "glmnet", auc = 0.8)
  res <- cv_ci(df)
  expect_true(is.na(res$auc_ci_lo))
  expect_true(is.na(res$auc_ci_hi))
})

test_that("multiple learners produce separate rows", {
  df <- data.frame(
    fold = rep(1:5, 2),
    learner = rep(c("glmnet", "ranger"), each = 5),
    auc = runif(10, 0.6, 0.9)
  )
  res <- cv_ci(df)
  expect_equal(nrow(res), 2)
  expect_true(all(c("glmnet", "ranger") %in% res$learner))
})

test_that("fit_resample metric_summary gains _ci_lo/_ci_hi columns", {
  set.seed(1)
  df <- data.frame(
    subject = rep(1:30, each = 2),
    outcome = rep(c(0, 1), length.out = 60),
    x1 = rnorm(60),
    x2 = rnorm(60)
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
  ms <- fit@metric_summary
  ci_cols <- grep("_ci_lo$|_ci_hi$", names(ms), value = TRUE)
  expect_gt(length(ci_cols), 0)
})

test_that(".nb_corrected_var returns NA for single fold", {
  expect_true(is.na(bioLeak:::.nb_corrected_var(0.8)))
})

test_that(".nb_corrected_var without n_train falls back to var/K", {
  vals <- c(0.7, 0.8, 0.9, 0.75, 0.85)
  K <- length(vals)
  res <- bioLeak:::.nb_corrected_var(vals, n_train = NULL, n_test = NULL)
  expected <- var(vals) / K
  expect_equal(res, expected)
})

test_that("cv_ci handles NA fold metrics correctly", {
  df <- data.frame(
    fold = 1:5,
    learner = rep("glm", 5),
    auc = c(0.9, NA, NA, NA, 0.7)
  )

  # Normal method: CI should be based on K_eff = 2, not K = 5
  res_n <- cv_ci(df, method = "normal")
  expect_equal(res_n$auc_mean, 0.8, tolerance = 1e-10)
  # With only 2 informative folds, df = 1, CI should be wide
  expect_true(all(is.finite(c(res_n$auc_ci_lo, res_n$auc_ci_hi))))

  # Nadeau-Bengio method: should also produce finite CIs
  res_nb <- cv_ci(df, method = "nadeau_bengio", n_train = 80, n_test = 20)
  expect_true(all(is.finite(c(res_nb$auc_ci_lo, res_nb$auc_ci_hi))))
  expect_equal(res_nb$auc_mean, 0.8, tolerance = 1e-10)
})

test_that("cv_ci returns NA CIs when only one finite fold", {
  df <- data.frame(
    fold = 1:3,
    learner = rep("glm", 3),
    auc = c(0.8, NA, NA)
  )
  res <- cv_ci(df, method = "normal")
  expect_true(is.na(res$auc_ci_lo))
  expect_true(is.na(res$auc_ci_hi))
})
