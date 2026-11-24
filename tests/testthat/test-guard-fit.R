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
  expect_warning(adj_test_full <- fit_eb$transform(X_test),
                 class = "bioLeak_unseen_batch_warning")   # B_unseen
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
    as.numeric(suppressWarnings(
      predict_guard(fit_eb, X_test[, c("g1", "g2", "g3")], batch = X_test$batch))[1, ]),
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
  expect_warning(fit_fr <- fit_resample(
    df_fr, outcome = "y", splits = sp_fr,
    learner = "glm", custom_learners = custom_glm, metrics = "auc",
    preprocess = list(
      impute = list(method = "median"),
      normalize = list(method = "quantile"),
      batch_adjust = list(method = "combat_frozen", batch_col = "batch")
    )
  ), class = "bioLeak_unseen_batch_warning")   # every held-out batch is unseen
  expect_s4_class(fit_fr, "LeakFit")
  expect_identical(fit_fr@feature_names, c("g1", "g2", "g3"))
})

# ── Batch adjustment: unseen batches (BL-3) ───────────────────────────────────

.ba_data <- function(seed = 2) {
  set.seed(seed)
  list(X = data.frame(a = c(rnorm(20, 0), rnorm(20, 3)), b = rnorm(40)),
       batch = rep(c("A", "B"), each = 20),
       y = factor(rep(0:1, 20)))
}
.ba_fit <- function(d, ...) {
  guard_fit(d$X, d$y,
            steps = list(batch_adjust = list(method = "combat_frozen", batch = d$batch, ...),
                         impute = list(winsor = FALSE),
                         normalize = list(method = "none")),
            task = "binomial")
}

test_that("unseen batches at transform warn by default", {
  d <- .ba_data(); gf <- .ba_fit(d)
  Xnew <- data.frame(a = rnorm(5, 10), b = rnorm(5))
  expect_warning(gf$transform(Xnew, batch = rep("C", 5)),
                 class = "bioLeak_unseen_batch_warning")
})

test_that("unseen = 'error' stops and unseen = 'grand_mean' is silent", {
  d <- .ba_data()
  Xnew <- data.frame(a = rnorm(5, 10), b = rnorm(5))
  expect_error(.ba_fit(d, unseen = "error")$transform(Xnew, batch = rep("C", 5)),
               class = "bioLeak_unseen_batch_error")
  gm <- .ba_fit(d, unseen = "grand_mean")
  expect_silent(out <- gm$transform(Xnew, batch = rep("C", 5)))
  expect_equal(out$a, Xnew$a)          # unseen rows pass through unadjusted
})

test_that("known batches are adjusted silently, and NA batch labels count as unseen", {
  d <- .ba_data(); gf <- .ba_fit(d)
  expect_silent(gf$transform(d$X[1:3, ], batch = d$batch[1:3]))
  expect_warning(gf$transform(d$X[1:3, ], batch = c("A", NA, "B")),
                 class = "bioLeak_unseen_batch_warning")
})

test_that("transform without batch labels warns when batch adjustment is enabled", {
  d <- .ba_data(); gf <- .ba_fit(d)
  expect_warning(gf$transform(d$X[1:3, ]), class = "bioLeak_unseen_batch_warning")
})

test_that("fit_resample gives one aggregated unseen-batch warning under batch-blocked CV", {
  set.seed(7)
  df <- data.frame(batch = rep(c("B1", "B2", "B3"), each = 16),
                   y = factor(rep(c("ctrl", "case"), 24)),
                   g1 = rnorm(48), g2 = rnorm(48))
  sp <- make_split_plan(df, outcome = "y", mode = "batch_blocked", batch = "batch", v = 3, progress = FALSE)
  w <- testthat::capture_warnings(
    fit_resample(df, outcome = "y", splits = sp, learner = "glm", metrics = "auc",
                 preprocess = list(batch_adjust = list(method = "combat_frozen", batch_col = "batch")))
  )
  hits <- grep("unseen in training", w, value = TRUE)
  expect_length(hits, 1L)
  expect_match(hits, "3 of 3 folds")
})

# ── Batch adjustment: protecting the outcome (BL-4) ───────────────────────────

.confounded <- function(seed = 3, batch_effect = 0) {
  set.seed(seed)
  y  <- rep(c(0, 1), each = 40)
  bt <- ifelse(runif(80) < ifelse(y == 1, 0.85, 0.15), "B1", "B2")
  X  <- data.frame(sig = y * 2 + rnorm(80) + batch_effect * (bt == "B1"), noise = rnorm(80))
  list(X = X, y = factor(y), batch = bt)
}
.adjust <- function(d, method = "location_scale", ...) {
  gf <- guard_fit(d$X, d$y,
                  steps = list(batch_adjust = list(method = method, batch = d$batch, ...),
                               normalize = list(method = "none")),
                  task = "binomial")
  gf$transform(d$X, batch = d$batch)
}

test_that("by default, a confounded batch warns and the outcome effect is attenuated", {
  d <- .confounded()
  expect_warning(adj <- .adjust(d), class = "bioLeak_batch_confounded_warning")
  before <- diff(tapply(d$X$sig, d$y, mean))
  expect_lt(diff(tapply(adj$sig, d$y, mean)) / before, 0.9)   # conservative, documented
})

test_that("protect = TRUE keeps the outcome effect when batch is confounded with it", {
  d <- .confounded()
  before <- diff(tapply(d$X$sig, d$y, mean))
  for (m in c("location_scale", "combat_frozen", "mean_center")) {
    w <- testthat::capture_warnings(adj <- .adjust(d, m, protect = TRUE))
    expect_true(any(grepl("bias cross-validated performance upward", w)), label = m)
    expect_gt(diff(tapply(adj$sig, d$y, mean)) / before, 0.9, label = m)
  }
})

test_that("protect = TRUE still removes a real batch effect", {
  d <- .confounded(batch_effect = 3)
  adj <- suppressWarnings(.adjust(d, protect = TRUE))
  gap <- function(x) mean(sapply(split(seq_along(x), d$y), function(i)
    abs(diff(tapply(x[i], d$batch[i], mean)))))
  expect_lt(gap(adj$sig), 0.25 * gap(d$X$sig))
  expect_gt(diff(tapply(adj$sig, d$y, mean)) / 2, 0.8)   # true effect is 2
})

test_that("no confounding warning when batch and outcome are balanced", {
  d <- .ba_data()
  expect_silent(.ba_fit(d))
  expect_silent(.ba_fit(d, protect = TRUE))
})

test_that("protect = TRUE with a fully confounded batch applies no adjustment, with a warning", {
  set.seed(4)
  y  <- factor(rep(0:1, each = 20))
  bt <- ifelse(y == "1", "B1", "B2")
  X  <- data.frame(sig = as.numeric(y) * 2 + rnorm(40))
  w <- testthat::capture_warnings(
    gf <- guard_fit(X, y, steps = list(batch_adjust = list(method = "location_scale", batch = bt, protect = TRUE),
                                       impute = list(winsor = FALSE), normalize = list(method = "none")),
                    task = "binomial"))
  expect_true(any(grepl("fully confounded", w)))
  expect_equal(gf$transform(X, batch = bt)$sig, X$sig)
})

test_that("protect = TRUE without a usable outcome warns and falls back", {
  d <- .confounded()
  expect_warning(
    guard_fit(d$X, NULL, steps = list(batch_adjust = list(method = "location_scale", batch = d$batch, protect = TRUE),
                                      normalize = list(method = "none")), task = "binomial"),
    "needs a non-survival outcome")
})

test_that("batch adjustment works without y and with a numeric outcome", {
  d <- .confounded()
  expect_silent(guard_fit(d$X, NULL,
    steps = list(batch_adjust = list(method = "location_scale", batch = d$batch),
                 normalize = list(method = "none")), task = "binomial"))
  set.seed(5)
  yn <- rnorm(80)
  bt <- ifelse(yn + rnorm(80, sd = 0.5) > 0, "B1", "B2")   # batch confounded with yn
  Xn <- data.frame(sig = 2 * yn + rnorm(80), noise = rnorm(80))
  gf <- suppressWarnings(guard_fit(Xn, yn,
    steps = list(batch_adjust = list(method = "location_scale", batch = bt, protect = TRUE),
                 normalize = list(method = "none")), task = "gaussian"))
  adj <- gf$transform(Xn, batch = bt)
  expect_gt(stats::coef(stats::lm(adj$sig ~ yn))[2] / 2, 0.9)   # true slope 2 kept
})

test_that("a single training batch is not adjusted, with a warning instead of an error", {
  set.seed(6)
  X <- data.frame(a = rnorm(20), b = rnorm(20))
  expect_warning(
    gf <- guard_fit(X, factor(rep(0:1, 10)),
                    steps = list(batch_adjust = list(method = "location_scale", batch = rep("A", 20)),
                                 impute = list(winsor = FALSE), normalize = list(method = "none")),
                    task = "binomial"),
    class = "bioLeak_batch_single_warning")
  expect_equal(gf$transform(X, batch = rep("A", 20))$a, X$a)
})

test_that("protect must be logical, and NA outcomes under protect = TRUE warn", {
  d <- .confounded()
  expect_error(.adjust(d, protect = "yes"), "TRUE or FALSE")
  y_na <- d$y; y_na[2] <- NA
  w <- testthat::capture_warnings(
    guard_fit(d$X, y_na, steps = list(batch_adjust = list(method = "location_scale", batch = d$batch, protect = TRUE),
                                      normalize = list(method = "none")), task = "binomial"))
  expect_true(any(grepl("missing values", w)))
  # the association warning must describe the fallback (unprotected), not protection
  expect_true(any(grepl("underestimated", w)))
  expect_false(any(grepl("signal is kept", w)))
})

# ── Review round (batch adjustment), cycle 2 ──────────────────────────────────

test_that("protect = TRUE handles a class confined to one batch without disabling everything", {
  set.seed(11)
  bt <- rep(c("A", "B", "C", "D"), each = 30)
  y  <- factor(c(rep(c("u", "v"), 45), rep("w", 30)))   # w only in batch D
  shift <- c(A = 0, B = 3, C = -3, D = 0)[bt]
  X  <- data.frame(V1 = rnorm(120) + shift + 2 * (y == "v"), V2 = rnorm(120))
  w <- testthat::capture_warnings(
    gf <- guard_fit(X, y, steps = list(batch_adjust = list(method = "location_scale", batch = bt, protect = TRUE),
                                       impute = list(winsor = FALSE), normalize = list(method = "none")),
                    task = "multiclass"))
  expect_true(any(grepl("cannot be separated from batch", w)))
  expect_false(any(grepl("no batch adjustment is applied", w)))
  adj <- gf$transform(X, batch = bt)
  abc <- bt %in% c("A", "B", "C")
  # within class u, the A/B/C batch means are equalised; the v effect is kept
  bm <- tapply(adj$V1[abc & y == "u"], bt[abc & y == "u"], mean)
  expect_lt(diff(range(bm)), 0.6)
  expect_gt(mean(adj$V1[abc & y == "v"]) - mean(adj$V1[abc & y == "u"]), 1.6)
})

test_that("protect = TRUE with single-class batches removes only within-class batch differences", {
  set.seed(12)
  bt <- rep(c("A", "B", "C", "D"), each = 20)
  y  <- factor(ifelse(bt %in% c("A", "B"), "0", "1"))
  X  <- data.frame(V1 = rnorm(80) + c(A = 0, B = 2, C = 0, D = 2)[bt] + 3 * (y == "1"))
  w <- testthat::capture_warnings(
    gf <- guard_fit(X, y, steps = list(batch_adjust = list(method = "location_scale", batch = bt, protect = TRUE),
                                       impute = list(winsor = FALSE), normalize = list(method = "none")),
                    task = "binomial"))
  expect_true(any(grepl("single outcome class", w)))
  adj <- gf$transform(X, batch = bt)
  m <- tapply(adj$V1, bt, mean)
  expect_lt(abs(m["A"] - m["B"]), 0.5)          # within-class batch difference removed
  expect_lt(abs(m["C"] - m["D"]), 0.5)
  expect_gt(mean(m[c("C", "D")]) - mean(m[c("A", "B")]), 2.5)   # class difference kept (true 3)
})

test_that("keep_batch_col = TRUE passes the batch column through unadjusted", {
  set.seed(13)
  X <- data.frame(g = rnorm(30) + rep(c(0, 3, -3), each = 10),
                  batch = rep(c("A", "B", "C"), each = 10), stringsAsFactors = FALSE)
  for (m in c("mean_center", "location_scale", "combat_frozen")) {
    gf <- guard_fit(X, NULL, steps = list(batch_adjust = list(method = m, batch_col = "batch", keep_batch_col = TRUE),
                                          impute = list(winsor = FALSE), normalize = list(method = "none")),
                    task = "gaussian")
    out <- gf$transform(X)
    bcols <- grep("^batch", names(out), value = TRUE)
    expect_true(length(bcols) >= 2L, label = m)
    expect_true(all(unlist(out[bcols]) %in% c(0, 1)), label = m)   # one-hot indicators intact
    expect_lt(diff(range(tapply(out$g, X$batch, mean))), 0.5, label = m)  # g still adjusted
  }
})

test_that("unseen = 'error' gives a consistent message", {
  d <- .ba_data()
  e <- tryCatch(.ba_fit(d, unseen = "error")$transform(d$X[1:2, ], batch = c("A", "Z")),
                error = function(e) conditionMessage(e))
  expect_match(e, "not seen in training \\(Z\\)")
  expect_false(grepl("left unadjusted|to stop", e))
})

test_that("missing batch labels warn once in guard_fit and once in fit_resample", {
  d <- .ba_data(); b_na <- d$batch; b_na[c(3, 25)] <- NA
  w <- testthat::capture_warnings(
    guard_fit(d$X, d$y, steps = list(batch_adjust = list(method = "location_scale", batch = b_na),
                                     normalize = list(method = "none")), task = "binomial"))
  expect_equal(sum(grepl("missing batch label", w)), 1L)

  set.seed(14)
  df <- data.frame(id = 1:60, batch = rep(c("B1", "B2", "B3"), 20), y = factor(rep(c("a", "b"), 30)),
                   g1 = rnorm(60), g2 = rnorm(60), stringsAsFactors = FALSE)
  df$batch[c(5, 17, 40)] <- NA
  sp <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "id", v = 5, progress = FALSE)
  w <- testthat::capture_warnings(
    fit_resample(df, outcome = "y", splits = sp, learner = "glm", metrics = "auc",
                 preprocess = list(batch_adjust = list(method = "location_scale", batch_col = "batch"))))
  expect_equal(sum(grepl("missing batch label", w)), 1L)
  expect_false(any(grepl("unseen in training", w)))
})

# ── Review round (batch adjustment), cycle 3 ──────────────────────────────────

test_that("protect = TRUE does not rescale the outcome when batches differ only in spread", {
  set.seed(21)
  n  <- 200
  y  <- rbinom(n, 1, 0.5)
  bt <- ifelse(runif(n) < ifelse(y == 1, 0.8, 0.2), "B1", "B2")
  sdv <- ifelse(bt == "B1", 1, 3)
  X  <- data.frame(v = 2 * y + rnorm(n, sd = sdv))           # no location batch effect
  for (m in c("location_scale", "combat_frozen")) {
    gf <- suppressWarnings(guard_fit(X, factor(y),
      steps = list(batch_adjust = list(method = m, batch = bt, protect = TRUE),
                   impute = list(winsor = FALSE), normalize = list(method = "none")),
      task = "binomial"))
    adj <- gf$transform(X, batch = bt)
    gap <- function(x, b) mean(x[bt == b & y == 1]) - mean(x[bt == b & y == 0])
    for (b in c("B1", "B2"))
      expect_equal(gap(adj$v, b), gap(X$v, b), tolerance = 1e-8, label = paste(m, b))
    ctrl <- tapply(adj$v[y == 0], bt[y == 0], mean)
    expect_lt(abs(diff(ctrl)), 0.6, label = m)                 # no new batch shift
  }
})

test_that("a numeric kept batch column does not exempt features sharing its name prefix", {
  set.seed(22)
  b <- rep(1:3, each = 10)
  X <- data.frame(g = b, g1 = rnorm(30) + c(0, 2, -2)[b], g2 = rnorm(30) + c(0, -2, 2)[b])
  gf <- guard_fit(X, NULL, steps = list(batch_adjust = list(method = "mean_center", batch_col = "g", keep_batch_col = TRUE),
                                        impute = list(winsor = FALSE), normalize = list(method = "none")),
                  task = "gaussian")
  out <- gf$transform(X)
  expect_equal(out$g, X$g)                                     # kept column untouched
  for (f in c("g1", "g2")) expect_lt(diff(range(tapply(out[[f]], b, mean))), 0.3, label = f)
})

test_that("fit_resample names folds with no adjustment and reports a single batch once", {
  set.seed(23)
  df <- data.frame(batch = rep(c("B1", "B2", "B3"), each = 20),
                   y = factor(c(rep("ctrl", 20), rep("case", 20), rep(c("ctrl", "case"), 10))),
                   g1 = rnorm(60), g2 = rnorm(60), stringsAsFactors = FALSE)
  sp <- make_split_plan(df, outcome = "y", mode = "batch_blocked", batch = "batch", v = 3, progress = FALSE)
  w <- testthat::capture_warnings(
    fit_resample(df, outcome = "y", splits = sp, learner = "glm", metrics = "auc",
                 preprocess = list(batch_adjust = list(method = "location_scale", batch_col = "batch",
                                                       protect = TRUE, unseen = "grand_mean"))))
  conf <- grep("batch_adjust: in", w, value = TRUE)
  conf <- grep("outcome|adjustment was applied", conf, value = TRUE)
  expect_length(conf, 1L)
  expect_match(conf, "no batch adjustment was applied")

  df2 <- df[df$batch != "B3", ]; df2$y <- factor(rep(c("ctrl", "case"), 20))
  sp2 <- make_split_plan(df2, outcome = "y", mode = "batch_blocked", batch = "batch", v = 2, progress = FALSE)
  w2 <- testthat::capture_warnings(
    fit_resample(df2, outcome = "y", splits = sp2, learner = "glm", metrics = "auc",
                 preprocess = list(batch_adjust = list(method = "location_scale", batch_col = "batch",
                                                       unseen = "grand_mean"))))
  expect_equal(sum(grepl("no batch adjustment was applied", w2)), 1L)
  expect_false(any(grepl("contain 1 batch level", w2)))
})

# ── Verification round: batch adjustment usability ───────────────────────────

test_that("unseen = 'error' under batch-blocked CV names the cause in the final error", {
  set.seed(31)
  df <- data.frame(batch = rep(c("B1", "B2", "B3"), each = 16), y = factor(rep(c("a", "b"), 24)),
                   g1 = rnorm(48), g2 = rnorm(48), stringsAsFactors = FALSE)
  sp <- make_split_plan(df, outcome = "y", mode = "batch_blocked", batch = "batch", v = 3, progress = FALSE)
  err <- tryCatch(suppressWarnings(fit_resample(df, outcome = "y", splits = sp, learner = "glm",
                    metrics = "auc", preprocess = list(batch_adjust = list(
                      method = "location_scale", batch_col = "batch", unseen = "error")))),
                  error = function(e) conditionMessage(e))
  expect_match(err, "No successful folds")
  expect_match(err, "not seen in training")
})

test_that("keep_batch_col on a split column warns that it has no effect", {
  set.seed(32)
  df <- data.frame(batch = rep(c("B1", "B2", "B3"), each = 16), y = factor(rep(c("a", "b"), 24)),
                   g1 = rnorm(48), g2 = rnorm(48), stringsAsFactors = FALSE)
  sp <- make_split_plan(df, outcome = "y", mode = "batch_blocked", batch = "batch", v = 3, progress = FALSE)
  w <- testthat::capture_warnings(fit_resample(df, outcome = "y", splits = sp, learner = "glm",
         metrics = "auc", preprocess = list(batch_adjust = list(method = "location_scale",
           batch_col = "batch", keep_batch_col = TRUE, unseen = "grand_mean"))))
  expect_true(any(grepl("keep_batch_col = TRUE has no effect", w)))
})

test_that("a batch vector that disagrees with batch_col warns", {
  set.seed(33)
  X <- data.frame(g = rnorm(20), bc = rep(c("A", "B"), 10), stringsAsFactors = FALSE)
  expect_warning(guard_fit(X, NULL, steps = list(batch_adjust = list(
                   method = "location_scale", batch = rep(c("x", "y"), each = 10), batch_col = "bc"),
                   normalize = list(method = "none")), task = "gaussian"),
                 "disagree", class = "bioLeak_input_warning")
})

test_that("imputation is fitted even when the training rows have no missing values", {
  set.seed(1)
  tr <- data.frame(a = rnorm(30), b = rnorm(30))
  te <- data.frame(a = c(NA, 1, 2), b = c(0.5, NA, 1))
  methods <- c("median", "knn", if (requireNamespace("randomForest", quietly = TRUE)) "missForest")
  for (m in methods) {
    g <- guard_fit(tr, steps = list(impute = list(method = m)))
    out <- if (identical(m, "missForest")) {
      # no per-column model without training NAs: documented median fallback
      expect_warning(g$transform(te), "falling back to median")
    } else {
      g$transform(te)
    }
    expect_equal(sum(is.na(out)), 0L, label = m)
  }
  # a held-out study that never measured a feature is still scored
  set.seed(2); n <- 160
  df <- data.frame(study = rep(paste0("S", 1:4), each = 40), x1 = rnorm(n), x2 = rnorm(n))
  df$y <- factor(rbinom(n, 1, plogis(1.5 * df$x1 + df$x2)))
  df$x1[df$study == "S4"] <- NA
  sp <- make_split_plan(df, outcome = "y", mode = "study_loocv", study = "study", progress = FALSE)
  fit <- suppressWarnings(fit_resample(df, outcome = "y", splits = sp, learner = "glm",
                                       custom_learners = make_custom_learners(), metrics = "auc",
                                       preprocess = list(impute = list(method = "median"))))
  expect_equal(nrow(fit@metrics), 4L)
})
