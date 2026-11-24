test_that("make_split_plan respects grouping constraints", {
  set.seed(1)
  df <- data.frame(
    outcome = rep(c(0, 1), each = 10),
    subject = rep(1:10, each = 2),
    batch = rep(letters[1:4], length.out = 20),
    study = rep(LETTERS[1:5], length.out = 20),
    time = seq_len(20),
    x1 = rnorm(20),
    x2 = rnorm(20)
  )

  splits_subject <- make_split_plan(df, outcome = "outcome",
                                mode = "subject_grouped", group = "subject",
                                v = 5, repeats = 1, stratify = TRUE,
                                seed = 1, progress = FALSE)
  for (fold in splits_subject@indices) {
    tr <- df$subject[fold$train]
    te <- df$subject[fold$test]
    expect_equal(length(intersect(unique(tr), unique(te))), 0)
  }

  splits_batch <- make_split_plan(df, outcome = "outcome",
                              mode = "batch_blocked", batch = "batch",
                              v = 4, repeats = 1, stratify = FALSE,
                              seed = 1, progress = FALSE)
  for (fold in splits_batch@indices) {
    tr <- df$batch[fold$train]
    te <- df$batch[fold$test]
    expect_equal(length(intersect(unique(tr), unique(te))), 0)
  }

  splits_study <- make_split_plan(df, outcome = "outcome",
                              mode = "study_loocv", study = "study",
                              seed = 1, progress = FALSE)
  for (fold in splits_study@indices) {
    te_study <- unique(df$study[fold$test])
    expect_equal(length(te_study), 1)
    expect_false(te_study %in% unique(df$study[fold$train]))
  }

  splits_time <- make_split_plan(df, outcome = "outcome",
                             mode = "time_series", time = "time",
                             v = 4, horizon = 1, seed = 1, progress = FALSE)
  expect_true(length(splits_time@indices) > 0)
  for (fold in splits_time@indices) {
    tmin <- min(df$time[fold$test])
    expect_true(all(df$time[fold$train] <= (tmin - 1)))
  }

  splits_grouped_time <- make_split_plan(df, outcome = "outcome",
                                         mode = "time_series", time = "time",
                                         group = "subject",
                                         v = 4, horizon = 1, seed = 1, progress = FALSE)
  expect_true(length(splits_grouped_time@indices) > 0)
  for (fold in splits_grouped_time@indices) {
    tmin <- min(df$time[fold$test])
    expect_true(all(df$time[fold$train] <= (tmin - 1)))
    expect_equal(length(intersect(df$subject[fold$train], df$subject[fold$test])), 0L)
  }
})

test_that("grouped time-series panel fallback tests every group and reports its training share (BL-5)", {
  grid <- expand.grid(subj = 1:30, visit = 1:4)
  df <- data.frame(subj = paste0("P", grid$subj), t = grid$visit,
                   y = factor(rep(c("a", "b"), 60)), x = rnorm(120))
  expect_warning(
    sp <- make_split_plan(df, outcome = "y", mode = "time_series", time = "t", group = "subj",
                          v = 4, progress = FALSE),
    "panel fallback.*% of the available rows", class = "bioLeak_fold_warning")
  tested <- unique(unlist(lapply(sp@indices, function(f) f$test)))
  expect_setequal(tested, which(df$t > 1))                # every row after the first time
  expect_setequal(unique(df$subj[tested]), unique(df$subj))
  for (f in sp@indices) {
    expect_length(intersect(df$subj[f$train], df$subj[f$test]), 0L)
    expect_lt(max(df$t[f$train]), min(df$t[f$test]))
  }
})

test_that("grouped time series reports dropped later blocks and test coverage (BL-18)", {
  # 4 subjects at t = 1..10, 4 more at t = 1..6 only: later blocks hold few groups
  df <- rbind(
    data.frame(subj = rep(paste0("L", 1:4), each = 10), t = rep(1:10, 4)),
    data.frame(subj = rep(paste0("S", 1:4), each = 6),  t = rep(1:6, 4))
  )
  df$y <- factor(rep(c("a", "b"), length.out = nrow(df))); df$x <- rnorm(nrow(df))
  w <- testthat::capture_warnings(
    sp <- make_split_plan(df, outcome = "y", mode = "time_series", time = "t", group = "subj",
                          v = 4, progress = FALSE))
  tested <- length(unique(unlist(lapply(sp@indices, function(f) f$test))))
  # Regular folds lose blocks here; whichever plan is kept (scored by
  # sum(test_n * train_n)), exactly one warning reports it with its coverage.
  plan_w <- grep("panel fallback|could not be used", w, value = TRUE)
  expect_length(plan_w, 1L)
  if (grepl("could not be used", plan_w)) {
    expect_match(plan_w, sprintf("%d of %d rows are tested", tested, nrow(df)))
  }
  for (f in sp@indices) {
    expect_length(intersect(df$subj[f$train], df$subj[f$test]), 0L)
    expect_lt(max(df$t[f$train]), min(df$t[f$test]))
  }
  # ungrouped time series on the same data: only the first block is untested, no warning
  expect_no_warning(make_split_plan(df, outcome = "y", mode = "time_series", time = "t",
                                    v = 4, progress = FALSE))
})

test_that("grouped time series prefers the panel fallback over a degenerate regular split", {
  grid <- expand.grid(subj = 1:36, visit = 1:4)
  df <- data.frame(subj = paste0("P", grid$subj), t = grid$visit,
                   y = factor(rep(c("a", "b"), 72)), x = rnorm(144))
  for (v in 3:5) {
    sp <- suppressWarnings(make_split_plan(df, outcome = "y", mode = "time_series", time = "t",
                                           group = "subj", v = v, progress = FALSE))
    tr <- vapply(sp@indices, function(f) length(f$train), integer(1))
    expect_gte(min(tr), 20L, label = paste("v =", v))
    expect_equal(length(unique(unlist(lapply(sp@indices, function(f) f$test)))), 108L)
  }
})

test_that("tied times are named as the reason for a dropped block", {
  df <- data.frame(t = rep(1:4, each = 36), y = factor(rep(c("a", "b"), 72)), x = rnorm(144))
  expect_warning(make_split_plan(df, outcome = "y", mode = "time_series", time = "t", v = 5,
                                 progress = FALSE), "tied times")
  w <- testthat::capture_warnings(make_split_plan(df, outcome = "y", mode = "time_series",
                                                  time = "t", v = 5, progress = FALSE))
  expect_false(any(grepl("With group = ", w)))
})

test_that("a missing stratify_by column is named in the warning", {
  df <- data.frame(id = 1:20, y = factor(rep(c("a", "b"), 10)), x = rnorm(20))
  expect_warning(make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "id",
                                 stratify = TRUE, stratify_by = "nope", v = 4, progress = FALSE),
                 "stratify_by column 'nope'")
})

test_that("nested = TRUE on a short panel does not fail the whole plan", {
  grid <- expand.grid(subj = 1:30, visit = 1:4)
  df <- data.frame(subj = paste0("P", grid$subj), t = grid$visit,
                   y = factor(rep(c("a", "b"), 60)), x = rnorm(120))
  w <- testthat::capture_warnings(
    sp <- make_split_plan(df, outcome = "y", mode = "time_series", time = "t", group = "subj",
                          v = 4, nested = TRUE, progress = FALSE))
  expect_s4_class(sp, "LeakSplits")
  expect_length(sp@info$inner, length(sp@indices))
  expect_true(sum(grepl("^nested = TRUE", w)) <= 2L)   # inner issues reported once each
})

test_that("the panel fallback does not replace a well-trained regular plan", {
  # 15 groups entering one per timepoint, each followed for 8 timepoints
  P <- do.call(rbind, lapply(1:15, function(i) data.frame(g = paste0("S", i), t = i + seq_len(8) - 1)))
  P$y <- rnorm(nrow(P)); P$x <- rnorm(nrow(P))
  w <- testthat::capture_warnings(
    sp <- make_split_plan(P, outcome = "y", mode = "time_series", time = "t", group = "g",
                          v = 6, progress = FALSE))
  tr <- vapply(sp@indices, function(f) length(f$train), integer(1))
  expect_false(any(grepl("panel fallback", w)))
  expect_gte(min(tr), 10L)
})

test_that("panel fallback folds keep the 3-test-row minimum and the warning names the plan", {
  df <- expand.grid(g = 1:3, t = 1:2)
  df$y <- rnorm(nrow(df)); df$x <- rnorm(nrow(df))
  w <- testthat::capture_warnings(
    sp <- try(make_split_plan(df, outcome = "y", mode = "time_series", time = "t", group = "g",
                              v = 3, progress = FALSE), silent = TRUE))
  if (!inherits(sp, "try-error")) {
    expect_true(all(vapply(sp@indices, function(f) length(f$test), integer(1)) >= 3L))
  }
  # balanced panel: fallback kept, warning reports both plans' coverage
  P <- expand.grid(g = paste0("S", 1:36), t = 1:4)
  P$y <- rnorm(nrow(P)); P$x <- rnorm(nrow(P))
  w <- testthat::capture_warnings(
    make_split_plan(P, outcome = "y", mode = "time_series", time = "t", group = "g",
                    v = 4, progress = FALSE))
  fb <- grep("panel fallback", w, value = TRUE)
  expect_length(fb, 1L)
  expect_match(fb, "scores higher on the sum of test rows x training rows; it tests 108 of 144 rows")
  expect_false(grepl("tests more rows with comparably", fb))
})

test_that("an empty grouped time-series plan explains why", {
  df <- expand.grid(g = 1:3, t = 1:2)
  df$y <- rnorm(nrow(df)); df$x <- rnorm(nrow(df))
  expect_error(suppressWarnings(make_split_plan(df, outcome = "y", mode = "time_series", time = "t",
                                                group = "g", v = 3, progress = FALSE)),
               "no regular fold was usable.*fewer than 3 test rows.*smaller v")
})

test_that("info$stratify is TRUE only when stratification was applied", {
  d <- data.frame(t = 1:100, y = rep(0:1, c(80, 20)), x = rnorm(100))
  expect_warning(sp <- make_split_plan(d, outcome = "y", mode = "time_series", time = "t", v = 5,
                                       stratify = TRUE, progress = FALSE),
                 "not applied to time_series folds", class = "bioLeak_fold_warning")
  expect_false(sp@info$stratify)
  d$g <- rep(1:20, 5)
  d$y <- rep(0:1, 50)  # half the groups are class 1
  expect_warning(sp2 <- make_split_plan(d, outcome = "y", mode = "subject_grouped", group = "g",
                                        v = 5, stratify = TRUE, stratify_by = "nope",
                                        progress = FALSE),
                 "was not found", class = "bioLeak_fold_warning")
  expect_false(sp2@info$stratify)
  expect_identical(sp2@info$stratify_by, "nope")
  sp3 <- make_split_plan(d, outcome = "y", mode = "subject_grouped", group = "g", v = 5,
                         stratify = TRUE, progress = FALSE)
  expect_true(sp3@info$stratify)
  # grouped panel fallback does stratify its groups
  P <- expand.grid(g = paste0("S", 1:36), t = 1:4)
  P$y <- rep(rbinom(36, 1, 0.4), 4); P$x <- rnorm(nrow(P))
  w <- testthat::capture_warnings(
    sp4 <- make_split_plan(P, outcome = "y", mode = "time_series", time = "t", group = "g",
                           v = 4, stratify = TRUE, progress = FALSE))
  expect_true(any(grepl("panel fallback", w)))
  expect_true(sp4@info$stratify)
  expect_false(any(grepl("not applied to time_series", w)))
})

test_that("a stratify_by column stays out of the predictors when stratification is not applied", {
  skip_if_not_installed("glmnet")
  set.seed(7); n <- 200
  df <- data.frame(t = 1:n, g = rep(1:40, 5), x1 = rnorm(n), x2 = rnorm(n))
  df$y <- factor(ifelse(0.3 * df$x1 + rnorm(n) > 0, "b", "a"))
  df$subtype <- as.numeric(df$y == "b")  # derived from the outcome
  plans <- list(
    time_series = suppressWarnings(make_split_plan(df, outcome = "y", mode = "time_series",
      time = "t", v = 5, stratify = TRUE, stratify_by = "subtype", progress = FALSE)),
    not_requested = make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "g",
      v = 5, stratify = FALSE, stratify_by = "subtype", progress = FALSE))
  for (nm in names(plans)) {
    expect_identical(plans[[nm]]@info$stratify_by, "subtype", label = nm)
    fit <- suppressWarnings(fit_resample(df[, c("y", "x1", "x2", "subtype", "t", "g")],
                                         outcome = "y", splits = plans[[nm]], learner = "glmnet",
                                         metrics = "auc", refit = FALSE, seed = 1))
    expect_lt(fit@metric_summary$auc_mean, 0.8, label = nm)
  }
})

test_that("stratified plans keep every requested fold when classes have few groups", {
  set.seed(1)
  d <- data.frame(g = rep(1:24, each = 3))
  d$y <- factor(rep(rep(c("a", "b", "c", "d"), each = 6), each = 3)); d$x <- rnorm(nrow(d))
  sp <- make_split_plan(d, outcome = "y", mode = "subject_grouped", group = "g", v = 10,
                        stratify = TRUE, progress = FALSE)
  sizes <- vapply(sp@indices, function(f) length(f$test), integer(1))
  expect_length(sp@indices, 10L)
  expect_lte(diff(range(sizes)), 3L)  # 24 groups of 3 rows over 10 folds: 2 or 3 groups each
  d2 <- data.frame(g = rep(1:12, each = 2))
  d2$y <- factor(rep(c(rep("a", 6), rep("b", 3), rep("c", 3)), each = 2)); d2$x <- rnorm(24)
  sp2 <- make_split_plan(d2, outcome = "y", mode = "subject_grouped", group = "g", v = 5,
                         stratify = TRUE, progress = FALSE, seed = 3)
  sizes2 <- vapply(sp2@indices, function(f) length(f$test), integer(1))
  expect_lte(diff(range(sizes2)), 2L)
  # every fold has a test set without groups shared with training
  for (f in sp@indices) expect_length(intersect(d$g[f$train], d$g[f$test]), 0L)
})
