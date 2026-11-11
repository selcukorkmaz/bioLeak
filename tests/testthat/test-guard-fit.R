test_that("guard_fit produces a GuardFit pipeline with median imputation", {
  X <- data.frame(a = c(1, 2, NA, 4), b = c(NA, 1, 2, 3))
  fit <- guard_fit(
    X, y = c(0, 1, 0, 1),
    steps = list(impute = list(method = "median"),
                 normalize = list(method = "zscore"),
                 filter = list(var_thresh = 0, iqr_thresh = 0),
                 fs = list(method = "none")),
    task = "binomial"
  )

  expect_true(inherits(fit, "GuardFit"))
  Xtr <- fit$transform(X)
  expect_true(is.data.frame(Xtr))
  expect_false(anyNA(Xtr))
})

test_that("guard_fit warns on impute none with missing values", {
  X <- data.frame(a = c(1, NA, 3), b = c(4, 5, NA))
  fit <- expect_warning_match(
    guard_fit(
      X, y = c(1, 0, 1),
      steps = list(impute = list(method = "none"),
                   normalize = list(method = "none"),
                   filter = list(var_thresh = 0, iqr_thresh = 0),
                   fs = list(method = "none")),
      task = "binomial"
    ),
    "impute\\$method='none'"
  )
  out <- fit$transform(X)
  expect_true(any(grepl("__missing", names(out))))
})

test_that("guard_fit supports knn imputation and fallback warnings", {
  X <- data.frame(a = c(1, NA, 3), b = c(4, 5, NA))
  fit <- guard_fit(
    X, y = c(1, 0, 1),
    steps = list(impute = list(method = "knn", k = 1),
                 normalize = list(method = "none"),
                 filter = list(var_thresh = 0, iqr_thresh = 0),
                 fs = list(method = "none")),
    task = "binomial"
  )

  env <- environment(fit$transform)
  env$state$impute$knn_ref <- NULL
  out <- expect_warning_match(fit$transform(X), "KNN imputation state missing")
})

test_that("guard_fit validates imputation method values", {
  X <- data.frame(a = c(1, NA), b = c(2, 3))
  expect_error(guard_fit(X, y = c(1, 0),
                                    steps = list(impute = list(method = "mice")),
                                    task = "binomial"),
               "mice")
  expect_error(guard_fit(X, y = c(1, 0),
                                    steps = list(impute = list(method = "bad")),
                                    task = "binomial"),
               "Unknown impute method")
})

test_that("guard_fit validates normalization and feature selection modes", {
  X <- data.frame(a = 1:4, b = 2:5)
  expect_error(guard_fit(X, y = 1:4,
                                    steps = list(normalize = list(method = "bad")),
                                    task = "gaussian"),
               "Unknown normalize method")

  expect_error(guard_fit(X, y = c(1, 1, 1, 1),
                                    steps = list(fs = list(method = "ttest")),
                                    task = "binomial"),
               "requires a binary outcome")
})

test_that("guard_fit filter respects min_keep", {
  X <- data.frame(a = rep(1, 5), b = 1:5, c = 2:6)
  fit <- guard_fit(
    X, y = c(0, 1, 0, 1, 0),
    steps = list(impute = list(method = "median"),
                 normalize = list(method = "none"),
                 filter = list(var_thresh = 10, iqr_thresh = 10, min_keep = 1),
                 fs = list(method = "none")),
    task = "binomial"
  )
  expect_true(fit$state$diagnostics$p_out >= 1)
})

test_that("guard_fit supports PCA feature selection", {
  X <- data.frame(a = 1:10, b = 2:11, c = 3:12)
  fit <- guard_fit(
    X, y = c(rep(0, 5), rep(1, 5)),
    steps = list(impute = list(method = "median"),
                 normalize = list(method = "none"),
                 filter = list(var_thresh = 0, iqr_thresh = 0),
                 fs = list(method = "pca", pca_comp = 2)),
    task = "binomial"
  )
  out <- fit$transform(X)
  expect_equal(ncol(out), 2)
})

test_that("guard_fit handles missForest when available", {
  X <- data.frame(a = c(1, NA, 3), b = c(4, 5, NA))
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    expect_error(guard_fit(
      X, y = c(1, 0, 1),
      steps = list(impute = list(method = "missForest")),
      task = "binomial"
    ), "randomForest")
  } else {
    fit <- suppressWarnings(guard_fit(
      X, y = c(1, 0, 1),
      steps = list(impute = list(method = "missForest", maxiter = 1, ntree = 5),
                   normalize = list(method = "none"),
                   filter = list(var_thresh = 0, iqr_thresh = 0),
                   fs = list(method = "none")),
      task = "binomial"
    ))
    expect_true(inherits(fit, "GuardFit"))
  }
})

test_that("guard_fit handles lasso feature selection when available", {
  X <- data.frame(a = rnorm(20), b = rnorm(20), c = rnorm(20))
  y <- factor(rep(c(0, 1), each = 10))
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    expect_error(guard_fit(
      X, y = y,
      steps = list(fs = list(method = "lasso")),
      task = "binomial"
    ), "glmnet")
  } else {
    fit <- suppressWarnings(guard_fit(
      X, y = y,
      steps = list(impute = list(method = "median"),
                   normalize = list(method = "none"),
                   filter = list(var_thresh = 0, iqr_thresh = 0),
                   fs = list(method = "lasso")),
      task = "binomial"
    ))
    expect_true(fit$state$fs$method == "lasso")
  }
})

test_that("guard_fit coerces character columns to numeric or factor", {
  X <- data.frame(a = c("1", "2", "3"), b = c("x", "y", "x"),
                  stringsAsFactors = FALSE)
  fit <- guard_fit(
    X, y = c(0, 1, 0),
    steps = list(impute = list(method = "median"),
                 normalize = list(method = "none"),
                 filter = list(var_thresh = 0, iqr_thresh = 0),
                 fs = list(method = "none")),
    task = "binomial"
  )
  out <- fit$transform(X)
  expect_true(ncol(out) >= 2)
})

test_that("predict_guard requires GuardFit objects", {
  X <- data.frame(a = 1:3, b = 2:4)
  fit <- guard_fit(X, y = c(1, 0, 1),
                              steps = list(impute = list(method = "median")),
                              task = "binomial")
  out <- predict_guard(fit, X)
  expect_true(is.data.frame(out))
  expect_error(predict_guard(list(), X))
})

test_that("variance/IQR filter uses the pre-normalization scale", {
  # Before 0.3.8 the filter ran on z-scored data, where every non-constant
  # column has variance 1, so var_thresh could not separate these features.
  set.seed(8)
  X <- data.frame(low = rnorm(50, sd = 0.1), mid = rnorm(50, sd = 1),
                  high = rnorm(50, sd = 5))
  for (norm in c("zscore", "robust", "none")) {
    fit <- guard_fit(
      X, y = rep(0:1, 25),
      steps = list(impute = list(method = "median"),
                   normalize = list(method = norm),
                   filter = list(var_thresh = 0.5),
                   fs = list(method = "none")),
      task = "binomial"
    )
    expect_identical(names(which(fit$state$filter$keep)), c("mid", "high"), info = norm)
    out <- fit$transform(X)
    expect_identical(colnames(out), c("mid", "high"), info = norm)
    if (norm == "zscore") {
      expect_equal(unname(vapply(out, stats::sd, numeric(1))), c(1, 1))
    }
  }

  fit_iqr <- guard_fit(
    X, y = rep(0:1, 25),
    steps = list(normalize = list(method = "zscore"),
                 filter = list(iqr_thresh = 3)),
    task = "binomial"
  )
  expect_identical(names(which(fit_iqr$state$filter$keep)), "high")

  # min_keep ranks by the original-scale spread
  fit_mk <- guard_fit(
    X, y = rep(0:1, 25),
    steps = list(normalize = list(method = "zscore"),
                 filter = list(var_thresh = 1e6, min_keep = 1)),
    task = "binomial"
  )
  expect_identical(names(which(fit_mk$state$filter$keep)), "high")
})

test_that("Phase 2 (2.3): guard_fit supports frozen quantile, log_cpm, clr, and combat_frozen batch adjustment", {
  set.seed(101)
  X_train <- data.frame(
    g1 = c(10, 20, 30, 100, 110, 120),
    g2 = c(20, 30, 40, 200, 210, 220),
    g3 = c(40, 50, 60, 300, 310, 320),
    batch = c("B1", "B1", "B1", "B2", "B2", "B2"),
    stringsAsFactors = FALSE
  )
  X_test <- data.frame(
    g1 = c(15, 115, 50),
    g2 = c(25, 215, 60),
    g3 = c(45, 315, 70),
    batch = c("B1", "B2", "B_unseen"),
    stringsAsFactors = FALSE
  )

  # 1) Frozen quantile normalization maps new rows solely against training q_ref
  fit_q <- guard_fit(
    X_train[, c("g1", "g2", "g3")],
    steps = list(impute = list(winsor = FALSE), normalize = list(method = "quantile")),
    task = "gaussian"
  )
  out_q1 <- fit_q$transform(X_test[1, c("g1", "g2", "g3"), drop = FALSE])
  out_q_all <- fit_q$transform(X_test[, c("g1", "g2", "g3"), drop = FALSE])
  expect_equal(as.numeric(out_q1[1, ]), as.numeric(out_q_all[1, ]))

  # 2) log_cpm and clr transforms
  fit_cpm <- guard_fit(
    X_train[, c("g1", "g2", "g3")],
    steps = list(impute = list(winsor = FALSE), normalize = list(method = "log_cpm")),
    task = "gaussian"
  )
  out_cpm <- fit_cpm$transform(X_test[, c("g1", "g2", "g3")])
  expect_false(anyNA(out_cpm))

  fit_clr <- guard_fit(
    X_train[, c("g1", "g2", "g3")],
    steps = list(impute = list(winsor = FALSE), normalize = list(method = "clr")),
    task = "gaussian"
  )
  out_clr <- fit_clr$transform(X_test[, c("g1", "g2", "g3")])
  expect_equal(unname(rowMeans(as.matrix(out_clr))), rep(0, nrow(X_test)), tolerance = 1e-8)

  # 3) Frozen ComBat / EB batch adjustment strips batch_col and shrinks batch shift
  fit_eb <- guard_fit(
    X_train,
    steps = list(
      impute = list(winsor = FALSE),
      normalize = list(method = "none"),
      batch_adjust = list(method = "combat_frozen", batch_col = "batch")
    ),
    task = "gaussian"
  )
  expect_false("batch" %in% colnames(fit_eb$transform(X_train)))
  adj_train <- fit_eb$transform(X_train)
  raw_gap <- abs(mean(X_train$g1[4:6]) - mean(X_train$g1[1:3]))
  adj_gap <- abs(mean(adj_train$g1[4:6]) - mean(adj_train$g1[1:3]))
  expect_lt(adj_gap, raw_gap * 0.25)

  # Held-out transform works with known and unseen batches without leaking test rows
  adj_test_single <- fit_eb$transform(X_test[1, , drop = FALSE])
  adj_test_full <- fit_eb$transform(X_test)
  expect_equal(as.numeric(adj_test_single[1, ]), as.numeric(adj_test_full[1, ]))

  # Single-column quantile normalization (p = 1), predict_guard(..., batch=), hash, and fit_resample integration
  fit_q1 <- guard_fit(
    data.frame(g1 = c(10, 20, 30, 40)),
    steps = list(impute = list(winsor = FALSE), normalize = list(method = "quantile")),
    task = "gaussian"
  )
  out_q1_single <- predict(fit_q1, data.frame(g1 = c(15, 35)))
  expect_equal(nrow(out_q1_single), 2L)
  expect_false(anyNA(out_q1_single))

  fit_no_eb <- guard_fit(
    X_train[, c("g1", "g2", "g3")],
    steps = list(impute = list(winsor = FALSE), normalize = list(method = "none")),
    task = "gaussian"
  )
  expect_false(identical(fit_eb$hash, fit_no_eb$hash))
  expect_equal(
    as.numeric(predict_guard(fit_eb, X_test[, c("g1", "g2", "g3")], batch = X_test$batch)[1, ]),
    as.numeric(adj_test_full[1, ])
  )

  # End-to-end fit_resample with batch_adjust when batch is also a split blocking column
  set.seed(202)
  df_fr <- data.frame(
    batch = rep(c("B1", "B2", "B3"), length.out = 48),
    y = factor(rep(c("ctrl", "case"), 24)),
    g1 = rnorm(48) + ifelse(rep(c("B1", "B2", "B3"), length.out = 48) == "B2", 5, 0),
    g2 = rnorm(48),
    g3 = rnorm(48),
    stringsAsFactors = FALSE
  )
  custom_glm <- list(
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
  sp_fr <- make_split_plan(df_fr, outcome = "y", mode = "batch_blocked", batch = "batch", v = 3, progress = FALSE)
  fit_fr <- fit_resample(
    df_fr, outcome = "y", splits = sp_fr,
    learner = "glm", custom_learners = custom_glm, metrics = "auc",
    preprocess = list(
      impute = list(method = "median"),
      normalize = list(method = "quantile"),
      batch_adjust = list(method = "combat_frozen", batch_col = "batch")
    )
  )
  expect_s4_class(fit_fr, "LeakFit")
  expect_identical(fit_fr@feature_names, c("g1", "g2", "g3"))
})
