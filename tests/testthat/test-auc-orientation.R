# Regression tests: AUC must be oriented (higher prediction = positive class).
# Before 0.3.8, pROC::roc() was called with direction = "auto", so an
# anti-correlated predictor scored AUC ~ 1 and permutation nulls centred
# above 0.5.

make_anti_fixture <- function(n = 120, seed = 1) {
  set.seed(seed)
  df <- data.frame(subject = rep(seq_len(n / 3), each = 3),
                   x1 = rnorm(n), x2 = rnorm(n))
  df$y <- factor(ifelse(df$x1 + rnorm(n, sd = 0.5) > 0, "case", "ctrl"),
                 levels = c("ctrl", "case"))
  df
}

score_learner <- function(sign = 1) {
  list(score = list(
    fit = function(x, y, task, weights, ...) NULL,
    predict = function(object, newdata, task, ...) sign * as.numeric(newdata[, "x1"])
  ))
}

test_that(".auc_binary is oriented and matches the Mann-Whitney statistic", {
  truth <- factor(c("a", "a", "b", "b"), levels = c("a", "b"))
  expect_equal(bioLeak:::.auc_binary(truth, c(0.1, 0.2, 0.8, 0.9)), 1)
  expect_equal(bioLeak:::.auc_binary(truth, c(0.9, 0.8, 0.2, 0.1)), 0)
  expect_equal(bioLeak:::.auc_binary(truth, c(0.5, 0.5, 0.5, 0.5)), 0.5)
  # positive_class overrides the level order
  expect_equal(bioLeak:::.auc_binary(truth, c(0.1, 0.2, 0.8, 0.9), positive_class = "a"), 0)
  # numeric 0/1 and logical outcomes
  expect_equal(bioLeak:::.auc_binary(c(0, 0, 1, 1), c(0.9, 0.8, 0.2, 0.1)), 0)
  expect_equal(bioLeak:::.auc_binary(c(FALSE, FALSE, TRUE, TRUE), c(0.1, 0.2, 0.8, 0.9)), 1)
  # a single class present gives NA, not an error
  expect_true(is.na(bioLeak:::.auc_binary(factor(c("b", "b"), levels = c("a", "b")), c(0.1, 0.2))))

  set.seed(3)
  y <- rbinom(60, 1, 0.4)
  p <- rnorm(60) - y
  mw <- mean(outer(p[y == 1], p[y == 0], function(a, b) (a > b) + 0.5 * (a == b)))
  expect_equal(bioLeak:::.auc_binary(y, p), mw)
  expect_lt(bioLeak:::.auc_binary(y, p), 0.5)
})

test_that("fit_resample reports AUC < 0.5 for anti-correlated predictions", {
  df <- make_anti_fixture()
  sp <- make_split_plan(df, outcome = "y", mode = "subject_grouped",
                        group = "subject", v = 4, progress = FALSE)
  pre <- list(normalize = list(method = "none"))
  fit_pos <- fit_resample(df, outcome = "y", splits = sp, learner = "score",
                          custom_learners = score_learner(1), metrics = "auc",
                          refit = FALSE, preprocess = pre)
  fit_neg <- fit_resample(df, outcome = "y", splits = sp, learner = "score",
                          custom_learners = score_learner(-1), metrics = "auc",
                          refit = FALSE, preprocess = pre)
  expect_gt(fit_pos@metric_summary$auc_mean, 0.8)
  expect_lt(fit_neg@metric_summary$auc_mean, 0.2)
  expect_equal(fit_pos@metrics$auc + fit_neg@metrics$auc,
               rep(1, nrow(fit_pos@metrics)), tolerance = 1e-12)
})

test_that("audit_leakage AUC is oriented and its fixed-prediction null centres at 0.5", {
  df <- make_anti_fixture()
  sp <- make_split_plan(df, outcome = "y", mode = "subject_grouped",
                        group = "subject", v = 4, progress = FALSE)
  fit <- fit_resample(df, outcome = "y", splits = sp, learner = "score",
                      custom_learners = score_learner(-1), metrics = "auc",
                      refit = FALSE, preprocess = list(normalize = list(method = "none")))
  aud <- suppressMessages(audit_leakage(fit, metric = "auc", B = 300,
                                        target_scan_multivariate = FALSE))
  pg <- aud@permutation_gap
  expect_lt(pg$metric_obs, 0.2)
  expect_lt(pg$gap, 0)
  # Under direction = "auto" every permuted AUC was >= 0.5 and the mean was
  # well above 0.5; oriented AUCs are symmetric around 0.5.
  expect_equal(mean(aud@perm_values), 0.5, tolerance = 0.02)
  expect_true(any(aud@perm_values < 0.5))
})

test_that("positive_class flips the reported AUC", {
  df <- make_anti_fixture()
  sp <- make_split_plan(df, outcome = "y", mode = "subject_grouped",
                        group = "subject", v = 4, progress = FALSE)
  pre <- list(normalize = list(method = "none"))
  fit_case <- fit_resample(df, outcome = "y", splits = sp, learner = "score",
                           custom_learners = score_learner(1), metrics = "auc",
                           refit = FALSE, preprocess = pre, positive_class = "case")
  fit_ctrl <- fit_resample(df, outcome = "y", splits = sp, learner = "score",
                           custom_learners = score_learner(1), metrics = "auc",
                           refit = FALSE, preprocess = pre, positive_class = "ctrl")
  expect_gt(fit_case@metric_summary$auc_mean, 0.8)
  expect_lt(fit_ctrl@metric_summary$auc_mean, 0.2)
})

test_that("target scans report oriented AUC values", {
  set.seed(4)
  y <- factor(rep(c("neg", "pos"), each = 30), levels = c("neg", "pos"))
  X <- data.frame(up = as.numeric(y == "pos") + rnorm(60, sd = 0.3),
                  down = -as.numeric(y == "pos") + rnorm(60, sd = 0.3))
  res <- bioLeak:::.target_assoc_scan(X, y, task = "binomial", positive_class = "pos")
  expect_gt(res$value[res$feature == "up"], 0.9)
  expect_lt(res$value[res$feature == "down"], 0.1)
  # the flagging score stays symmetric: a strong inverse proxy is still a proxy
  expect_equal(res$score[res$feature == "up"], res$score[res$feature == "down"],
               tolerance = 0.1)
})
