# L1.2-L1.3 (JSS review): confounding_alignment flags only when batch is
# associated with the outcome AND batches are shared by training and test rows.

.l1_data <- function(seed, confounded = TRUE, n_batch = 6, per_batch = 30,
                     gaussian = FALSE) {
  set.seed(seed)
  batch <- rep(paste0("b", seq_len(n_batch)), each = per_batch)
  p <- if (confounded) rep(c(0.15, 0.85), length.out = n_batch) else rep(0.5, n_batch)
  pb <- p[as.integer(factor(batch))]
  df <- data.frame(subject = paste0("s", seq_along(batch)), batch = batch,
                   x1 = stats::rnorm(length(batch)), x2 = stats::rnorm(length(batch)))
  if (gaussian) {
    df$y <- stats::rnorm(length(batch), mean = if (confounded) 2 * pb else 0)
  } else {
    df$y <- factor(ifelse(stats::runif(length(batch)) < pb, "case", "ctrl"),
                   levels = c("ctrl", "case"))
  }
  df
}

.l1_audit <- function(df, mode = "subject_grouped", compact = FALSE) {
  sp <- switch(mode,
    subject_grouped = make_split_plan(df, outcome = "y", mode = mode, group = "subject",
                                      v = 5, seed = 1, compact = compact, progress = FALSE),
    batch_blocked = make_split_plan(df, outcome = "y", mode = mode, batch = "batch",
                                    v = 3, seed = 1, progress = FALSE),
    study_loocv = make_split_plan(transform(df, study = batch), outcome = "y", mode = mode,
                                  study = "study", seed = 1, progress = FALSE))
  keep <- c(switch(mode, subject_grouped = "subject", batch_blocked = "batch",
                   study_loocv = "study"), "y", "x1", "x2")
  if (identical(mode, "study_loocv")) df$study <- df$batch
  fit <- NULL
  utils::capture.output(
    fit <- fit_resample(df[, keep], outcome = "y", splits = sp,
                        learner = "glm", custom_learners = make_custom_learners(),
                        metrics = if (is.factor(df$y)) "auc" else "rmse", refit = FALSE, seed = 1)
  )
  aud <- NULL
  utils::capture.output(suppressWarnings(
    aud <- audit_leakage(fit, metric = if (is.factor(df$y)) "auc" else "rmse", B = 20,
                         coldata = df, batch_cols = "batch",
                         target_scan = FALSE, target_scan_multivariate = FALSE)
  ))
  aud
}

.l1_flag <- function(aud) {
  ms <- aud@info$mechanism_summary
  ms$flagged[ms$mechanism_class == "confounding_alignment"]
}

test_that("batch associated with the outcome and shared by train/test is flagged", {
  aud <- .l1_audit(.l1_data(1))
  ba <- aud@batch_assoc
  expect_identical(ba$outcome_test[1], "chisq")
  expect_equal(ba$outcome_df[1], 5)                 # 6 batches x 2 classes
  expect_lt(ba$outcome_pval_adj[1], 0.05)
  expect_gt(ba$span_frac[1], 0.9)
  expect_true(.l1_flag(aud))
})

test_that("blocked designs are not flagged even when batch predicts the outcome", {
  skip_on_cran()
  aud_b <- .l1_audit(.l1_data(1), mode = "batch_blocked")
  expect_equal(aud_b@batch_assoc$span_frac[1], 0)
  expect_lt(aud_b@batch_assoc$outcome_pval_adj[1], 0.05)
  expect_false(.l1_flag(aud_b))
  aud_s <- .l1_audit(.l1_data(1), mode = "study_loocv")
  expect_equal(aud_s@batch_assoc$span_frac[1], 0)
  expect_false(.l1_flag(aud_s))
})

test_that("batch unrelated to the outcome is not flagged", {
  skip_on_cran()
  aud <- .l1_audit(.l1_data(2, confounded = FALSE))
  expect_gt(aud@batch_assoc$span_frac[1], 0.9)
  expect_false(.l1_flag(aud))
})

test_that("regression outcomes use a Kruskal-Wallis batch-outcome test", {
  skip_on_cran()
  aud <- .l1_audit(.l1_data(3, gaussian = TRUE))
  expect_identical(aud@batch_assoc$outcome_test[1], "kruskal")
  expect_true(.l1_flag(aud))
  aud0 <- .l1_audit(.l1_data(3, confounded = FALSE, gaussian = TRUE))
  expect_false(.l1_flag(aud0))
})

test_that("compact plans give the same span and flag", {
  skip_on_cran()
  df <- .l1_data(4)
  a1 <- .l1_audit(df)
  a2 <- .l1_audit(df, compact = TRUE)
  expect_equal(a1@batch_assoc$span_frac, a2@batch_assoc$span_frac)
  expect_identical(.l1_flag(a1), .l1_flag(a2))
})

test_that(".batch_outcome_assoc handles one batch, NA and survival", {
  expect_true(is.na(.batch_outcome_assoc(rep("a", 10), factor(rep(c("x", "y"), 5)), "binomial")$outcome_pval))
  r <- .batch_outcome_assoc(c(rep("a", 10), rep("b", 10), NA),
                            factor(c(rep("x", 9), "y", rep("y", 9), "x", "x")), "binomial")
  expect_lt(r$outcome_pval, 0.05)
  expect_true(is.na(.batch_outcome_assoc(rep(c("a", "b"), 5), 1:10, "survival")$outcome_pval))
})

test_that("repeated rows of a subject count once in the batch-outcome test", {
  skip_on_cran()
  set.seed(12)
  sb <- sample(paste0("b", 1:4), 40, TRUE); sy <- sample(c("ctrl", "case"), 40, TRUE)
  s <- rep(1:40, each = 5)
  df <- data.frame(subject = paste0("s", s), batch = sb[s],
                   y = factor(sy[s], levels = c("ctrl", "case")),
                   x1 = stats::rnorm(200), x2 = stats::rnorm(200))
  aud <- .l1_audit(df)
  ba <- aud@batch_assoc
  expect_equal(ba$outcome_n[1], 40)
  ref <- .chisq_assoc(table(sb, sy))
  expect_equal(ba$outcome_pval[1], ref$pval)
})

test_that("within-subject confounding is reported, not tested (L1.2 contract)", {
  skip_on_cran()
  set.seed(13)
  s <- rep(1:30, each = 3)
  df <- data.frame(patient = paste0("p", s), y = factor(rep(c("ctrl", "case", "case"), 30),
                                                        levels = c("ctrl", "case")),
                   x1 = stats::rnorm(90), x2 = stats::rnorm(90))
  df$batch <- ifelse(df$y == "ctrl", "N", "T")
  sp <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "patient",
                        v = 5, seed = 1, progress = FALSE)
  fit <- NULL
  utils::capture.output(fit <- fit_resample(df[, c("patient", "y", "x1", "x2")], outcome = "y",
    splits = sp, learner = "glm", custom_learners = make_custom_learners(),
    metrics = "auc", refit = FALSE, seed = 1))
  aud <- NULL
  utils::capture.output(suppressWarnings(aud <- audit_leakage(fit, metric = "auc", B = 5,
    coldata = df, batch_cols = "batch", target_scan = FALSE, target_scan_multivariate = FALSE,
    duplicate_null_reps = 0)))
  ba <- aud@batch_assoc
  expect_identical(ba$outcome_method[1], "units")
  expect_gt(ba$within_unit_effect[1], 0.9)
  expect_false(.l1_flag(aud))
  out <- capture.output(summary(aud))
  expect_true(any(grepl("within-subject association not tested", out)))
})

test_that("rsample grouped plans test subjects, not rows", {
  skip_on_cran()
  skip_if_not_installed("rsample")
  set.seed(14)
  sb <- sample(c("A", "B", "C"), 40, TRUE); sy <- sample(c("n", "p"), 40, TRUE)
  u <- rep(1:40, sample(c(1, 3, 8), 40, TRUE))
  df <- data.frame(subject = paste0("s", u), batch = sb[u], y = factor(sy[u]),
                   x1 = stats::rnorm(length(u)), x2 = stats::rnorm(length(u)))
  rs <- rsample::group_vfold_cv(df, group = "subject", v = 5)
  fit <- NULL
  utils::capture.output(fit <- suppressWarnings(fit_resample(df, outcome = "y", splits = rs,
    learner = "glm", custom_learners = make_custom_learners(), metrics = "auc",
    split_cols = list(group = "subject", batch = "batch"))))
  aud <- NULL
  utils::capture.output(suppressWarnings(aud <- audit_leakage(fit, metric = "auc", B = 5,
    perm_refit = FALSE, target_scan = FALSE, target_scan_multivariate = FALSE,
    duplicate_null_reps = 0)))
  expect_identical(aud@batch_assoc$outcome_method[1], "units")
  expect_equal(aud@batch_assoc$outcome_n[1], 40)
})

test_that("a subject constraint axis of a combined plan is the unit", {
  skip_on_cran()
  set.seed(18)
  n_subj <- 60
  s <- rep(seq_len(n_subj), each = 4)
  df <- data.frame(subject = paste0("s", s),
                   proc = sample(paste0("p", 1:6), length(s), TRUE),
                   site = sample(paste0("site", 1:3), n_subj, TRUE)[s],
                   y = factor(sample(c("ctrl", "case"), n_subj, TRUE)[s], levels = c("ctrl", "case")),
                   x1 = stats::rnorm(length(s)), x2 = stats::rnorm(length(s)))
  sp <- suppressWarnings(make_split_plan(df, outcome = "y", mode = "combined",
    constraints = list(list(type = "batch", col = "proc"), list(type = "subject", col = "subject")),
    v = 3, seed = 1, progress = FALSE))
  fit <- NULL
  utils::capture.output(fit <- suppressWarnings(fit_resample(df[, c("proc", "subject", "y", "x1", "x2")],
    outcome = "y", splits = sp, learner = "glm", custom_learners = make_custom_learners(),
    metrics = "auc", refit = FALSE, seed = 1)))
  aud <- NULL
  utils::capture.output(suppressWarnings(aud <- audit_leakage(fit, metric = "auc", B = 5,
    coldata = df, batch_cols = "site", target_scan = FALSE, target_scan_multivariate = FALSE,
    duplicate_null_reps = 0)))
  expect_identical(aud@batch_assoc$outcome_method[1], "units")
  expect_equal(aud@batch_assoc$outcome_n[1], n_subj)
})

test_that("audit_leakage_by_learner prints batch-outcome evidence, not the fold-batch table", {
  skip_on_cran()
  df <- .l1_data(1)
  sp <- make_split_plan(df, outcome = "y", mode = "batch_blocked", batch = "batch", v = 3,
                        seed = 1, progress = FALSE)
  custom <- make_custom_learners()
  fit <- NULL
  utils::capture.output(fit <- suppressWarnings(fit_resample(df[, c("batch", "y", "x1", "x2")],
    outcome = "y", splits = sp, learner = c("glm", "glm2"),
    custom_learners = c(custom, list(glm2 = custom$glm)), metrics = "auc", refit = FALSE, seed = 1)))
  al <- NULL
  utils::capture.output(suppressWarnings(al <- audit_leakage_by_learner(fit, metric = "auc", B = 5,
    coldata = df, batch_cols = "batch", target_scan = FALSE, target_scan_multivariate = FALSE,
    duplicate_null_reps = 0)))
  out <- capture.output(tab <- print(al))
  expect_false(any(grepl("batch_v_max|batch_p_min", out)))
  expect_true(any(grepl("batch_outcome_p", out)))
  expect_true(any(grepl("batch_span_max", out)))
})

test_that("a recorded group is the unit in study_loocv plans", {
  skip_on_cran()
  set.seed(16)
  s <- rep(1:40, each = 5)
  df <- data.frame(subject = paste0("s", s), study = rep(paste0("st", 1:4), each = 50),
                   batch = sample(paste0("b", 1:4), 40, TRUE)[s],
                   y = factor(sample(c("ctrl", "case"), 40, TRUE)[s], levels = c("ctrl", "case")),
                   x1 = stats::rnorm(200), x2 = stats::rnorm(200))
  sp <- make_split_plan(df, outcome = "y", mode = "study_loocv", study = "study",
                        group = "subject", seed = 1, progress = FALSE)
  fit <- NULL
  utils::capture.output(fit <- suppressWarnings(fit_resample(df[, c("study", "y", "x1", "x2")],
    outcome = "y", splits = sp, learner = "glm", custom_learners = make_custom_learners(),
    metrics = "auc", refit = FALSE, seed = 1)))
  aud <- NULL
  utils::capture.output(suppressWarnings(aud <- audit_leakage(fit, metric = "auc", B = 5,
    coldata = df, batch_cols = "batch", target_scan = FALSE,
    target_scan_multivariate = FALSE, duplicate_null_reps = 0)))
  expect_identical(aud@batch_assoc$outcome_method[1], "units")
  expect_equal(aud@batch_assoc$outcome_n[1], 40)
})

test_that("time-series plans without a group do not test dependent rows", {
  skip_on_cran()
  set.seed(17)
  n <- 150
  df <- data.frame(t = seq_len(n), batch = rep(paste0("r", 1:6), each = 25),
                   y = as.numeric(stats::arima.sim(list(ar = 0.9), n)),
                   x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  sp <- make_split_plan(df, outcome = "y", mode = "time_series", time = "t", v = 4,
                        progress = FALSE)
  fit <- NULL
  utils::capture.output(fit <- suppressWarnings(fit_resample(df[, c("t", "y", "x1", "x2")],
    outcome = "y", splits = sp, learner = "glm", custom_learners = make_custom_learners(),
    metrics = "rmse", refit = FALSE, seed = 1)))
  aud <- NULL
  utils::capture.output(suppressWarnings(aud <- audit_leakage(fit, metric = "rmse", B = 5,
    coldata = df, batch_cols = "batch", target_scan = FALSE,
    target_scan_multivariate = FALSE, duplicate_null_reps = 0)))
  expect_identical(aud@batch_assoc$outcome_method[1], "not_tested")
  expect_true(is.na(aud@batch_assoc$outcome_pval[1]))
  expect_false(.l1_flag(aud))
  out <- capture.output(summary(aud))
  expect_true(any(grepl("not tested", out)))
})

test_that("ties in the per-unit summary are broken at random, not by name", {
  set.seed(19)
  rej <- vapply(1:150, function(b) {
    set.seed(b)
    u <- rep(1:300, sample(1:2, 300, TRUE))
    r <- .batch_outcome_assoc(sample(paste0("b", 1:4), length(u), TRUE),
                              factor(sample(c("case", "ctrl"), length(u), TRUE)),
                              "binomial", unit = u, seed = b)
    r$outcome_pval <= 0.05
  }, logical(1))
  expect_lt(mean(rej), 0.05 + 3 * sqrt(0.05 * 0.95 / 150))
  set.seed(5); before <- stats::runif(1)
  set.seed(5)
  invisible(.batch_outcome_assoc(c("a", "b", "a", "b"), factor(c("x", "y", "y", "x")),
                                 "binomial", unit = c(1, 1, 2, 2)))
  expect_identical(stats::runif(1), before)
})
