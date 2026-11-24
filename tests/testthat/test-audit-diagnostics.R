# The class fixture gives every subject the labels (0, 1) in the same order,
# so subject-level permutations cannot move a label and audit_leakage() warns
# that the null is degenerate. These tests check other behaviour with B = 2.
audit_leakage_nd <- function(...) {
  withCallingHandlers(audit_leakage(...), bioLeak_input_warning = function(w) {
    if (grepl("every permutation gave the same", conditionMessage(w), fixed = TRUE)) {
      invokeRestart("muffleWarning")
    }
  })
}

test_that("audit_leakage reports batch association and duplicates", {
  set.seed(2)
  X <- matrix(rnorm(24), nrow = 12, ncol = 2)
  X[7, ] <- X[1, ]
  X[8, ] <- X[2, ]

  df <- data.frame(
    outcome = rep(c(0, 1), 6),
    batch = rep(c("A", "B"), each = 6),
    x1 = X[, 1],
    x2 = X[, 2]
  )

  fold1 <- c(1L, 2L, 3L, 7L, 8L, 9L)
  fold2 <- c(4L, 5L, 6L, 10L, 11L, 12L)
  indices <- list(
    list(train = setdiff(seq_len(12), fold1), test = fold1, fold = 1, repeat_id = 1),
    list(train = setdiff(seq_len(12), fold2), test = fold2, fold = 2, repeat_id = 1)
  )
  splits <- bioLeak:::LeakSplits(mode = "custom", indices = indices,
                                 info = list(outcome = "outcome", coldata = df))

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
    metrics = "auc",
    refit = FALSE,
    seed = 1
  )

  audit <- audit_leakage(
    fit,
    metric = "auc",
    B = 10,
    perm_stratify = FALSE,
    batch_cols = "batch",
    X_ref = df[, c("x1", "x2")],
    sim_threshold = 0.999,
    duplicate_scope = "all",
    # ADD THIS LINE TO FIX WARNINGS:
    target_scan_multivariate = FALSE
  )

  expect_true(nrow(audit@batch_assoc) >= 1)
  expect_true(is.finite(audit@batch_assoc$pval[1]))
  expect_true(nrow(audit@duplicates) >= 1)
})

test_that("audit_leakage batch association respects repeated CV folds", {
  set.seed(3)
  df <- make_class_df(8)
  df$batch <- rep(c("A", "B"), each = 4)
  df <- df[, c("batch", "outcome", "x1", "x2")]

  fold1 <- c(1L, 2L, 5L, 6L)
  fold2 <- c(3L, 4L, 7L, 8L)
  fold3 <- c(1L, 2L, 7L, 8L)
  fold4 <- c(3L, 4L, 5L, 6L)
  indices <- list(
    list(train = setdiff(seq_len(8), fold1), test = fold1, fold = 1, repeat_id = 1),
    list(train = setdiff(seq_len(8), fold2), test = fold2, fold = 2, repeat_id = 1),
    list(train = setdiff(seq_len(8), fold3), test = fold3, fold = 1, repeat_id = 2),
    list(train = setdiff(seq_len(8), fold4), test = fold4, fold = 2, repeat_id = 2)
  )
  splits <- bioLeak:::LeakSplits(mode = "custom", indices = indices,
                                 info = list(outcome = "outcome", coldata = df, batch = "batch"))

  custom <- make_custom_learners()
  fit <- fit_resample_quiet(
    df,
    outcome = "outcome",
    splits = splits,
    learner = "glm",
    custom_learners = custom,
    metrics = "auc",
    refit = FALSE,
    seed = 1
  )

  audit <- audit_leakage(
    fit,
    metric = "auc",
    B = 2,
    perm_stratify = FALSE,
    batch_cols = "batch",
    return_perm = FALSE,
    target_scan = FALSE,
    target_scan_multivariate = FALSE
  )

  expect_equal(sort(unique(audit@batch_assoc$repeat_id)), c(1, 2))
  expect_true(all(audit@batch_assoc$df == 1))
})

test_that("audit_leakage supports refit-based permutations", {
  df <- make_class_df(20)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()

  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = custom,
                            metrics = "auc", refit = FALSE, seed = 1)

  audit <- audit_leakage_nd(
    fit,
    metric = "auc",
    B = 2,
    perm_refit = TRUE,
    perm_refit_spec = list(
      x = df,
      outcome = "outcome",
      learner = "glm",
      custom_learners = custom
    ),
    perm_stratify = FALSE
  )

  expect_true(identical(audit@info$perm_method, "refit"))
  expect_equal(audit@permutation_gap$n_perm[1], 2)
})

test_that("audit_leakage supports perm_refit auto mode", {
  df <- make_class_df(20)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = custom,
                            metrics = "auc", refit = FALSE, seed = 1)

  audit_fixed <- audit_leakage_nd(
    fit,
    metric = "auc",
    B = 5,
    perm_refit = "auto",
    perm_refit_auto_max = 2,
    perm_stratify = FALSE
  )
  expect_true(identical(audit_fixed@info$perm_method, "fixed"))
  expect_true(grepl("^auto", audit_fixed@info$perm_refit_mode))

  audit_refit <- audit_leakage_nd(
    fit,
    metric = "auc",
    B = 2,
    perm_refit = "auto",
    perm_refit_auto_max = 5,
    perm_refit_spec = list(
      x = df,
      outcome = "outcome",
      learner = "glm",
      custom_learners = custom
    ),
    perm_stratify = FALSE
  )
  expect_true(identical(audit_refit@info$perm_method, "refit"))
  expect_true(grepl("^auto", audit_refit@info$perm_refit_mode))
})

test_that("audit_leakage multivariate target scan runs", {
  df <- make_class_df(24)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = custom,
                            metrics = "auc", refit = FALSE, seed = 1)
  aud <- audit_leakage_nd(
    fit,
    metric = "auc",
    B = 3,
    perm_stratify = FALSE,
    X_ref = df[, c("x1", "x2")],
    target_scan_multivariate = TRUE,
    target_scan_multivariate_B = 2,
    target_scan_multivariate_components = 2
  )
  mv <- aud@info$target_multivariate
  expect_true(is.data.frame(mv))
  expect_true(nrow(mv) > 0)
  expect_true(all(c("metric", "score", "p_value") %in% names(mv)))
})

test_that("audit_leakage adds mechanism taxonomy and FDR-adjusted target flags", {
  df <- make_class_df(24)
  splits <- make_split_plan_quiet(
    df,
    outcome = "outcome",
    mode = "subject_grouped",
    group = "subject",
    v = 3,
    seed = 7
  )
  custom <- make_custom_learners()
  fit <- fit_resample_quiet(
    df,
    outcome = "outcome",
    splits = splits,
    learner = "glm",
    custom_learners = custom,
    metrics = "auc",
    refit = FALSE,
    seed = 7
  )

  aud <- audit_leakage_nd(
    fit,
    metric = "auc",
    B = 3,
    perm_stratify = FALSE,
    X_ref = df[, c("x1", "x2")],
    target_scan_multivariate = FALSE,
    target_p_adjust = "BH",
    target_alpha = 0.2
  )

  expect_true("mechanism_class" %in% names(aud@permutation_gap))
  if (nrow(aud@batch_assoc)) {
    expect_true("mechanism_class" %in% names(aud@batch_assoc))
  }
  if (nrow(aud@target_assoc)) {
    expect_true("mechanism_class" %in% names(aud@target_assoc))
    expect_true("p_value_adj" %in% names(aud@target_assoc))
    expect_true("flag_fdr" %in% names(aud@target_assoc))
  }
  if (nrow(aud@duplicates)) {
    expect_true("mechanism_class" %in% names(aud@duplicates))
  }

  expect_true(is.data.frame(aud@info$taxonomy))
  expect_true(all(c("mechanism_class", "evidence_slot") %in% names(aud@info$taxonomy)))
  expect_true(is.data.frame(aud@info$mechanism_summary))
  expect_true(all(c("mechanism_class", "flagged", "evidence") %in% names(aud@info$mechanism_summary)))
})

test_that("duplicate detection ignores column scale and still finds real duplicates", {
  set.seed(5); n <- 120
  cl <- data.frame(id = 1:n, age = round(rnorm(n, 55, 10)), bmi = round(rnorm(n, 27, 4), 1),
                   sbp = round(rnorm(n, 130, 15)), glucose = round(rnorm(n, 100, 15)),
                   chol = round(rnorm(n, 200, 30)))
  cl$y <- factor(rbinom(n, 1, 0.4))
  X <- cl[, c("age", "bmi", "sbp", "glucose", "chol")]
  sp <- make_split_plan(cl, outcome = "y", mode = "subject_grouped", group = "id", v = 5,
                        progress = FALSE)
  capture.output(fit <- suppressWarnings(fit_resample(cl, outcome = "y", splits = sp, learner = "glm",
    custom_learners = make_custom_learners(), metrics = "auc", seed = 1)))
  capture.output(a <- suppressWarnings(audit_leakage(fit, metric = "auc", B = 10, X_ref = X,
                                                     target_scan = FALSE)))
  expect_equal(nrow(a@duplicates), 0L)  # unscaled tabular data, no duplicates
  X2 <- X; X2[3, ] <- X2[110, ]
  capture.output(a2 <- suppressWarnings(audit_leakage(fit, metric = "auc", B = 10, X_ref = X2,
                                                      target_scan = FALSE)))
  expect_gte(nrow(a2@duplicates), 1L)
})

test_that("batch alignment the split creates on purpose is not flagged", {
  set.seed(6); n <- 120
  d <- data.frame(batch = rep(paste0("B", 1:6), each = 20), x1 = rnorm(n), x2 = rnorm(n))
  d$y <- factor(rbinom(n, 1, .5)); d$subj <- rep(1:6, each = 20)
  sp <- make_split_plan(d, outcome = "y", mode = "batch_blocked", batch = "batch", v = 6,
                        progress = FALSE)
  capture.output(f <- suppressWarnings(fit_resample(d[, c("y", "x1", "x2", "batch")], outcome = "y",
    splits = sp, learner = "glm", custom_learners = make_custom_learners(), metrics = "auc", seed = 1)))
  capture.output(a <- suppressWarnings(audit_leakage(f, metric = "auc", B = 10, coldata = d,
                                                     target_scan = FALSE)))
  ms <- audit_info(a)$mechanism_summary
  expect_true(all(a@batch_assoc$by_design))
  expect_false(ms$flagged[ms$mechanism_class == "confounding_alignment"])
  # a subject split whose batches each span two subjects (batch not nested in
  # subject): batches are shared by train and test, but the outcome is random,
  # so fold-batch alignment alone is not confounding (JSS review, L1.3)
  set.seed(7)
  d2 <- data.frame(subj = rep(1:30, each = 4), batch = rep(paste0("b", 1:15), each = 8),
                   x1 = rnorm(120), x2 = rnorm(120))
  d2$y <- factor(rbinom(120, 1, .5))
  sp2 <- make_split_plan(d2, outcome = "y", mode = "subject_grouped", group = "subj", v = 5,
                         progress = FALSE)
  capture.output(f2 <- suppressWarnings(fit_resample(d2[, c("y", "x1", "x2", "subj")], outcome = "y",
    splits = sp2, learner = "glm", custom_learners = make_custom_learners(), metrics = "auc", seed = 1)))
  capture.output(a2 <- suppressWarnings(audit_leakage(f2, metric = "auc", B = 10, coldata = d2,
                                                      target_scan = FALSE)))
  ms2 <- audit_info(a2)$mechanism_summary
  expect_false(any(a2@batch_assoc$by_design))
  expect_true(all(a2@batch_assoc$span_frac > 0))
  expect_false(ms2$flagged[ms2$mechanism_class == "confounding_alignment"])
  # a batch nested in subjects is kept apart by the subject split: by design
  d2$batch <- paste0("b", d2$subj)
  capture.output(a3 <- suppressWarnings(audit_leakage(f2, metric = "auc", B = 10, coldata = d2,
                                                      target_scan = FALSE)))
  expect_true(all(a3@batch_assoc$by_design))
  out <- capture.output(summary(a))
  expect_true(any(grepl("permutation p =", out)))
  expect_true(any(grepl("blocked by the split", out)))
  expect_false(any(grepl("Strong non-random signal", out)))
})

test_that("design permutations move whole same-size subjects and keep class balance", {
  key <- c("1", "1", "0", "0", "1", "0", "1", "1", "0", "0", "0")
  grp <- c("a", "a", "b", "b", "c", "c", "d", "e", "e", "e", "e")
  set.seed(1)
  for (i in 1:30) {
    d <- bioLeak:::.design_perm_donor(key, grp)
    # each subject takes the whole label vector of one subject of its size
    for (g in unique(grp)) {
      donors <- unique(grp[d[grp == g]])
      expect_length(donors, 1L)
      expect_equal(sum(grp == donors), sum(grp == g))
    }
    # row-level class balance is unchanged
    expect_equal(sort(key[d]), sort(key))
    # a subject without a same-size partner stays within itself
    expect_setequal(d[grp == "d"], which(grp == "d"))
  }
  # paired design (one sample of each class per subject): labels flip within pairs
  pk <- rep(c("0", "1"), 10); pg <- rep(1:10, each = 2)
  flips <- vapply(1:50, function(i) {
    d <- bioLeak:::.design_perm_donor(pk, pg)
    expect_equal(sort(pk[d]), sort(pk))
    !identical(pk[d], pk)
  }, logical(1))
  expect_true(any(flips))
  if (FALSE) {
  }
  expect_equal(bioLeak:::.design_perm_coverage(grp), 6 / 11)
})

test_that("duplicate detection with few features reports only identical rows", {
  set.seed(9); n <- 150
  d <- data.frame(id = 1:n, x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  d$y <- factor(rbinom(n, 1, .5))
  sp <- make_split_plan(d, outcome = "y", mode = "subject_grouped", group = "id", v = 5,
                        progress = FALSE)
  capture.output(fit <- suppressWarnings(fit_resample(d, outcome = "y", splits = sp, learner = "glm",
    custom_learners = make_custom_learners(), metrics = "auc", seed = 1)))
  X <- d[, c("x1", "x2", "x3")]
  capture.output(a <- suppressWarnings(audit_leakage(fit, metric = "auc", B = 10, X_ref = X,
                                                     target_scan = FALSE, duplicate_scope = "all")))
  expect_equal(nrow(a@duplicates), 0L)
  X2 <- X; X2[40, ] <- X2[7, ]
  capture.output(a2 <- suppressWarnings(audit_leakage(fit, metric = "auc", B = 10, X_ref = X2,
                                                      target_scan = FALSE, duplicate_scope = "all")))
  expect_equal(nrow(a2@duplicates), 1L)
  # the summary reports repeats from the plan
  sp5 <- make_split_plan(d, outcome = "y", mode = "subject_grouped", group = "id", v = 5,
                         repeats = 3, progress = FALSE)
  capture.output(f5 <- suppressWarnings(fit_resample(d, outcome = "y", splits = sp5, learner = "glm",
    custom_learners = make_custom_learners(), metrics = "auc", seed = 1)))
  capture.output(a5 <- suppressWarnings(audit_leakage(f5, metric = "auc", B = 5, target_scan = FALSE,
                                                      target_scan_multivariate = FALSE)))
  expect_true(any(grepl("Repeats: 3", capture.output(summary(a5)), fixed = TRUE)))
})
