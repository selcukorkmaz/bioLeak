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
    duplicate_scope = "all"
  )

  expect_true(nrow(audit@batch_assoc) >= 1)
  expect_true(is.finite(audit@batch_assoc$pval[1]))
  expect_true(nrow(audit@duplicates) >= 1)
})

test_that("audit_leakage supports refit-based permutations", {
  df <- make_class_df(20)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()

  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = custom,
                            metrics = "auc", refit = FALSE, seed = 1)

  audit <- audit_leakage(
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
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = custom,
                            metrics = "auc", refit = FALSE, seed = 1)

  audit_fixed <- audit_leakage(
    fit,
    metric = "auc",
    B = 5,
    perm_refit = "auto",
    perm_refit_auto_max = 2,
    perm_stratify = FALSE
  )
  expect_true(identical(audit_fixed@info$perm_method, "fixed"))
  expect_true(grepl("^auto", audit_fixed@info$perm_refit_mode))

  audit_refit <- audit_leakage(
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
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = custom,
                            metrics = "auc", refit = FALSE, seed = 1)
  aud <- audit_leakage(
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
