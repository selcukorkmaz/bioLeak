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
    sim_threshold = 0.999
  )

  expect_true(nrow(audit@batch_assoc) >= 1)
  expect_true(is.finite(audit@batch_assoc$pval[1]))
  expect_true(nrow(audit@duplicates) >= 1)
})
