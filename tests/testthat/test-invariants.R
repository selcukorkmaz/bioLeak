# tests/testthat/test-invariants.R
# Invariant tests: seed reproducibility, nested-CV independence,
# transform locality, split completeness, multi-axis edge cases

# --- Helper data ---
.make_invariant_df <- function(n_subjects = 10, per_subject = 2) {
  n <- n_subjects * per_subject
  set.seed(999)
  data.frame(
    subject = rep(seq_len(n_subjects), each = per_subject),
    batch   = rep(letters[seq_len(min(n_subjects, 5))], length.out = n),
    outcome = rbinom(n, 1, 0.5),
    x1      = rnorm(n),
    x2      = rnorm(n),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Seed reproducibility
# ============================================================

test_that("same seed -> identical splits (subject_grouped)", {
  df <- .make_invariant_df()
  s1 <- make_split_plan(df, outcome = "outcome",
                        mode = "subject_grouped", group = "subject",
                        v = 5, seed = 42, progress = FALSE)
  s2 <- make_split_plan(df, outcome = "outcome",
                        mode = "subject_grouped", group = "subject",
                        v = 5, seed = 42, progress = FALSE)
  expect_identical(s1@indices, s2@indices)
})

test_that("different seeds -> different splits", {
  df <- .make_invariant_df()
  s1 <- make_split_plan(df, outcome = "outcome",
                        mode = "subject_grouped", group = "subject",
                        v = 5, seed = 1, progress = FALSE)
  s2 <- make_split_plan(df, outcome = "outcome",
                        mode = "subject_grouped", group = "subject",
                        v = 5, seed = 99, progress = FALSE)
  # At least one fold should differ
  any_diff <- !identical(s1@indices, s2@indices)
  expect_true(any_diff)
})

test_that("same seed -> identical splits (batch_blocked)", {
  df <- .make_invariant_df()
  s1 <- make_split_plan(df, outcome = "outcome",
                        mode = "batch_blocked", batch = "batch",
                        v = 3, seed = 7, progress = FALSE)
  s2 <- make_split_plan(df, outcome = "outcome",
                        mode = "batch_blocked", batch = "batch",
                        v = 3, seed = 7, progress = FALSE)
  expect_identical(s1@indices, s2@indices)
})

test_that("same seed -> identical splits (combined mode)", {
  df <- .make_invariant_df()
  s1 <- suppressWarnings(make_split_plan(
    df, outcome = "outcome",
    mode = "combined",
    constraints = list(
      list(type = "subject", col = "subject"),
      list(type = "batch", col = "batch")
    ),
    v = 3, seed = 5, progress = FALSE
  ))
  s2 <- suppressWarnings(make_split_plan(
    df, outcome = "outcome",
    mode = "combined",
    constraints = list(
      list(type = "subject", col = "subject"),
      list(type = "batch", col = "batch")
    ),
    v = 3, seed = 5, progress = FALSE
  ))
  expect_identical(s1@indices, s2@indices)
})

test_that("same seed -> identical fit_resample metrics (custom GLM)", {
  df <- .make_invariant_df()
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 3, seed = 1, progress = FALSE)
  custom <- list(
    glm = list(
      fit = function(x, y, task, weights, ...) {
        suppressWarnings(
          stats::glm(y ~ ., data = as.data.frame(x),
                     family = stats::binomial(), weights = weights)
        )
      },
      predict = function(object, newdata, task, ...) {
        as.numeric(suppressWarnings(
          stats::predict(object, newdata = as.data.frame(newdata),
                         type = "response")
        ))
      }
    )
  )
  f1 <- fit_resample(df, outcome = "outcome", splits = splits,
                     learner = "glm", custom_learners = custom,
                     metrics = "auc", refit = FALSE, seed = 10)
  f2 <- fit_resample(df, outcome = "outcome", splits = splits,
                     learner = "glm", custom_learners = custom,
                     metrics = "auc", refit = FALSE, seed = 10)
  expect_identical(f1@metrics, f2@metrics)
})

# ============================================================
# Nested CV independence
# ============================================================

test_that("inner train indices do not intersect outer test indices", {
  df <- .make_invariant_df(n_subjects = 12, per_subject = 2)
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 4, nested = TRUE, seed = 1, progress = FALSE)
  inner <- splits@info$inner
  expect_true(!is.null(inner))

  for (i in seq_along(splits@indices)) {
    outer_test <- splits@indices[[i]]$test
    outer_train <- splits@indices[[i]]$train
    for (j in seq_along(inner[[i]])) {
      # Inner indices are local to outer_train; map to global
      inner_train_global <- outer_train[inner[[i]][[j]]$train]
      inner_test_global  <- outer_train[inner[[i]][[j]]$test]
      expect_length(intersect(inner_train_global, outer_test), 0)
    }
  }
})

test_that("inner test indices do not intersect outer test indices", {
  df <- .make_invariant_df(n_subjects = 12, per_subject = 2)
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 4, nested = TRUE, seed = 1, progress = FALSE)
  inner <- splits@info$inner

  for (i in seq_along(splits@indices)) {
    outer_test <- splits@indices[[i]]$test
    outer_train <- splits@indices[[i]]$train
    for (j in seq_along(inner[[i]])) {
      inner_test_global <- outer_train[inner[[i]][[j]]$test]
      expect_length(intersect(inner_test_global, outer_test), 0)
    }
  }
})

# ============================================================
# Transform locality
# ============================================================

test_that("guard_fit zscore uses train stats, not test stats", {
  set.seed(1)
  n_train <- 50
  n_test <- 50
  # Use train mean=0, test mean=5 (within winsorization bounds)
  train_x <- matrix(rnorm(n_train * 2, mean = 0, sd = 1), ncol = 2,
                     dimnames = list(NULL, c("x1", "x2")))
  test_x  <- matrix(rnorm(n_test * 2, mean = 5, sd = 1), ncol = 2,
                     dimnames = list(NULL, c("x1", "x2")))
  train_y <- rbinom(n_train, 1, 0.5)

  guard_obj <- guard_fit(
    X = train_x, y = train_y,
    steps = list(
      impute = list(winsor = FALSE),
      normalize = list(method = "zscore")
    )
  )

  pred <- predict_guard(guard_obj, newdata = test_x)
  # After zscore normalization using train stats (mean~0, sd~1),
  # test values (mean~5) should be transformed to ~(5-0)/1 = ~5
  # If test stats were used, test mean would be ~0 after transform
  expect_true(mean(pred[, 1]) > 2)
  expect_true(mean(pred[, 2]) > 2)
})

# ============================================================
# Split completeness
# ============================================================

test_that("every sample appears in at least one test fold (subject_grouped)", {
  df <- .make_invariant_df()
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 5, seed = 1, progress = FALSE)
  all_test <- sort(unique(unlist(lapply(splits@indices, `[[`, "test"))))
  expect_equal(all_test, seq_len(nrow(df)))
})

test_that("every sample appears in at least one test fold (batch_blocked)", {
  df <- .make_invariant_df()
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "batch_blocked", batch = "batch",
                            v = 5, seed = 1, progress = FALSE)
  all_test <- sort(unique(unlist(lapply(splits@indices, `[[`, "test"))))
  expect_equal(all_test, seq_len(nrow(df)))
})

# ============================================================
# Multi-axis edge cases
# ============================================================

test_that("combined mode with 3 constraint axes works and has zero overlap on all 3", {
  set.seed(1)
  # Need enough unique levels per axis so that removing test-axis levels

  # still leaves training samples
  n_subj <- 30
  df <- data.frame(
    subject = rep(seq_len(n_subj), each = 3),
    batch   = rep(paste0("B", seq_len(n_subj)), each = 3),
    site    = rep(paste0("S", seq_len(n_subj)), each = 3),
    outcome = rbinom(n_subj * 3, 1, 0.5),
    x1      = rnorm(n_subj * 3),
    stringsAsFactors = FALSE
  )
  constraints <- list(
    list(type = "subject", col = "subject"),
    list(type = "batch", col = "batch"),
    list(type = "batch", col = "site")
  )
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "combined", constraints = constraints,
                            v = 3, seed = 1, progress = FALSE)
  expect_true(length(splits@indices) > 0)

  for (fold in splits@indices) {
    for (col in c("subject", "batch", "site")) {
      tr_vals <- unique(df[[col]][fold$train])
      te_vals <- unique(df[[col]][fold$test])
      overlap <- intersect(tr_vals, te_vals)
      expect_equal(length(overlap), 0,
                   info = sprintf("overlap on '%s'", col))
    }
  }
})

test_that("perfectly confounded axes (1:1 subject-batch) still produce valid folds", {
  set.seed(1)
  n_subj <- 10
  df <- data.frame(
    subject = rep(1:n_subj, each = 2),
    batch   = rep(paste0("B", 1:n_subj), each = 2),
    outcome = rbinom(n_subj * 2, 1, 0.5),
    x1      = rnorm(n_subj * 2),
    stringsAsFactors = FALSE
  )
  constraints <- list(
    list(type = "subject", col = "subject"),
    list(type = "batch", col = "batch")
  )
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "combined", constraints = constraints,
                            v = 5, seed = 1, progress = FALSE)
  expect_true(length(splits@indices) > 0)
  for (fold in splits@indices) {
    expect_length(intersect(df$subject[fold$train], df$subject[fold$test]), 0)
    expect_length(intersect(df$batch[fold$train], df$batch[fold$test]), 0)
  }
})

# ============================================================
# Mechanism summary in audit output
# ============================================================

test_that("summary(audit) output contains Mechanism Risk Assessment", {
  set.seed(1)
  df <- data.frame(
    subject = rep(1:6, each = 2),
    outcome = rbinom(12, 1, 0.5),
    x1 = rnorm(12),
    x2 = rnorm(12)
  )
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 3, seed = 1, progress = FALSE)
  custom <- list(
    glm = list(
      fit = function(x, y, task, weights, ...) {
        suppressWarnings(
          stats::glm(y ~ ., data = as.data.frame(x),
                     family = stats::binomial(), weights = weights)
        )
      },
      predict = function(object, newdata, task, ...) {
        as.numeric(suppressWarnings(
          stats::predict(object, newdata = as.data.frame(newdata),
                         type = "response")
        ))
      }
    )
  )
  fit <- fit_resample(df, outcome = "outcome", splits = splits,
                      learner = "glm", custom_learners = custom,
                      metrics = "auc", refit = FALSE, seed = 1)
  audit <- audit_leakage(fit, metric = "auc", B = 5,
                         X_ref = df[, c("x1", "x2")], seed = 1)

  out <- capture.output(summary(audit))
  expect_true(any(grepl("Mechanism Risk Assessment", out)))
  expect_true(any(grepl("non_random_signal", out)))
})
