test_that("cv_ci returns correct column names", {
  df <- data.frame(
    fold = 1:5,
    learner = rep("glmnet", 5),
    auc = runif(5, 0.6, 0.9),
    accuracy = runif(5, 0.7, 0.95)
  )
  res <- cv_ci(df, n_train = 80, n_test = 20)
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
  res <- cv_ci(df, n_train = 80, n_test = 20)
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
  res <- cv_ci(df, n_train = 80, n_test = 20)
  expect_true(is.na(res$auc_ci_lo))
  expect_true(is.na(res$auc_ci_hi))
})

test_that("multiple learners produce separate rows", {
  df <- data.frame(
    fold = rep(1:5, 2),
    learner = rep(c("glmnet", "ranger"), each = 5),
    auc = runif(10, 0.6, 0.9)
  )
  res <- cv_ci(df, n_train = 80, n_test = 20)
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

test_that("cv_ci nadeau_bengio uses folds_per_repeat - 1 degrees of freedom (BL-14)", {
  set.seed(1)
  df <- data.frame(fold = 1:50, learner = "glm", auc = rnorm(50, 0.7, 0.05))
  r5 <- cv_ci(df, method = "nadeau_bengio", n_train = 80, n_test = 20, folds_per_repeat = 5)
  m <- mean(df$auc); se <- sqrt((1 / 50 + 20 / 80) * var(df$auc))
  expect_equal(r5$auc_ci_hi, m + qt(0.975, 4) * se)
  # without folds_per_repeat: one repeat of 50 folds, still a t reference (no qnorm switch)
  r1 <- cv_ci(df, method = "nadeau_bengio", n_train = 80, n_test = 20)
  expect_equal(r1$auc_ci_hi, m + qt(0.975, 49) * se)
  expect_error(cv_ci(df, method = "nadeau_bengio", n_train = 80, n_test = 20, folds_per_repeat = 1),
               "folds_per_repeat")
})

test_that("cv_ci warns when nadeau_bengio lacks n_train / n_test (BL-15)", {
  df <- data.frame(fold = 1:5, learner = "glm", auc = c(0.7, 0.72, 0.68, 0.71, 0.69))
  expect_warning(cv_ci(df, method = "nadeau_bengio"), class = "bioLeak_cv_ci_warning")
  expect_silent(cv_ci(df, method = "normal"))
})

test_that("cv_ci validates sizes and folds_per_repeat", {
  df <- data.frame(fold = 1:10, learner = "glm", auc = seq(0.6, 0.8, length.out = 10))
  expect_error(cv_ci(df, method = "nadeau_bengio", n_train = 80, n_test = -5), "n_test")
  expect_error(cv_ci(df, method = "nadeau_bengio", n_train = c(80, 90), n_test = 20), "n_train")
  expect_error(cv_ci(df, method = "nadeau_bengio", n_train = 80, n_test = 20, folds_per_repeat = 4.5),
               "whole number")
  expect_silent(cv_ci(df, method = "normal", folds_per_repeat = 1))   # unused by "normal"
})

test_that("cv_ci defaults to nadeau_bengio and takes sizes from a LeakFit (BL-17)", {
  set.seed(1)
  df <- data.frame(id = 1:60, y = factor(rep(c("a", "b"), 30)), x1 = rnorm(60), x2 = rnorm(60))
  sp <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "id",
                        v = 5, repeats = 2, progress = FALSE)
  fit <- suppressWarnings(fit_resample(df, outcome = "y", splits = sp, learner = "glm",
                                       metrics = "accuracy"))
  expect_silent(res_fit <- cv_ci(fit))
  manual <- cv_ci(fit@metrics, method = "nadeau_bengio",
                  n_train = mean(fit@audit$n_train), n_test = mean(fit@audit$n_test),
                  folds_per_repeat = 5)
  expect_equal(res_fit, manual)
  # a bare metrics data.frame without sizes warns under the default
  expect_warning(cv_ci(fit@metrics), class = "bioLeak_cv_ci_warning")
  expect_identical(formals(cv_ci)$method[[2]], "nadeau_bengio")
})

test_that("cv_ci counts repeated identical test sets once (verification finding)", {
  set.seed(5)
  df <- data.frame(b = rep(paste0("B", 1:5), each = 12), y = factor(rep(c("a", "b"), 30)),
                   x1 = rnorm(60), x2 = rnorm(60))
  sp1 <- make_split_plan(df, outcome = "y", mode = "batch_blocked", batch = "b", v = 5, progress = FALSE)
  sp4 <- make_split_plan(df, outcome = "y", mode = "batch_blocked", batch = "b", v = 5, repeats = 4,
                         progress = FALSE)
  m1 <- data.frame(fold = 1:5, learner = "glm", auc = c(0.6, 0.7, 0.8, 0.65, 0.75))
  m4 <- data.frame(fold = 1:20, learner = "glm", auc = rep(m1$auc, 4))
  # same metrics per identical test set: identical intervals
  key_order <- vapply(sp4@indices, function(z) paste(sort(z$test), collapse = ","), "")
  key1 <- vapply(sp1@indices, function(z) paste(sort(z$test), collapse = ","), "")
  m4$auc <- m1$auc[match(key_order, key1)]
  r1 <- cv_ci(methods::new("LeakFit", splits = sp1, metrics = m1))
  r4 <- cv_ci(methods::new("LeakFit", splits = sp4, metrics = m4))
  expect_equal(r4$auc_ci_lo, r1$auc_ci_lo)
  expect_equal(r4$auc_ci_hi, r1$auc_ci_hi)
})

test_that("cv_ci accepts a LeakTune-shaped object", {
  df <- data.frame(id = 1:40, y = factor(rep(c("a", "b"), 20)), x = rnorm(40))
  sp <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "id", v = 4, progress = FALSE)
  lt <- structure(list(splits = sp, metrics = data.frame(fold = 1:4, learner = "glm",
                                                         auc = c(0.6, 0.7, 0.65, 0.72))),
                  class = "LeakTune")
  expect_silent(r <- cv_ci(lt))
  expect_true(all(is.finite(c(r$auc_ci_lo, r$auc_ci_hi))))
})

test_that("metric_summary CI columns belong to the right learner", {
  skip_if_not_installed("ranger")
  skip_if_not_installed("glmnet")
  set.seed(1)
  d <- data.frame(id = 1:120, x1 = rnorm(120), x2 = rnorm(120))
  d$y <- factor(rbinom(120, 1, plogis(1.5 * d$x1)))
  sp <- make_split_plan(d, outcome = "y", mode = "subject_grouped", group = "id", v = 5,
                        progress = FALSE)
  fit <- suppressWarnings(fit_resample(d, outcome = "y", splits = sp, learner = c("ranger", "glmnet"),
                                       metrics = "auc", refit = FALSE))
  ms <- fit@metric_summary
  ci <- cv_ci(fit)
  ci <- ci[match(ms$learner, ci$learner), ]
  expect_equal(ms$auc_ci_lo, ci$auc_ci_lo)
  expect_equal(ms$auc_ci_hi, ci$auc_ci_hi)
  expect_true(all(ms$auc_ci_lo < ms$auc_mean & ms$auc_mean < ms$auc_ci_hi))
})
