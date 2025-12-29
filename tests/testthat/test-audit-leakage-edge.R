test_that("audit_leakage validates predictions and learner selection", {
  df <- make_class_df(12)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()
  custom$glm2 <- custom$glm

  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = c("glm", "glm2"),
                            custom_learners = custom,
                            metrics = "auc", refit = FALSE, seed = 1)

  expect_error(audit_leakage(fit, metric = "auc", B = 5),
               "Multiple learners")
  expect_error(audit_leakage(fit, metric = "auc", learner = "missing", B = 5),
               "not found")

  fit@predictions <- list()
  expect_error(audit_leakage(fit, metric = "auc", B = 5),
               "No predictions available")
})

test_that("audit_leakage warns when learner IDs are missing", {
  df <- make_class_df(12)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()
  custom$glm2 <- custom$glm

  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = c("glm", "glm2"),
                            custom_learners = custom,
                            metrics = "auc", refit = FALSE, seed = 1)
  fit@predictions <- lapply(fit@predictions, function(df_pred) {
    df_pred$learner <- NULL
    df_pred
  })

  expect_warning({
    audit <- audit_leakage(fit, metric = "auc", B = 5)
  }, "predictions lack learner IDs")

  expect_warning({
    audit <- audit_leakage(fit, metric = "auc", learner = "glm", B = 5)
  }, "ignored")
})

test_that("audit_leakage handles duplicate detection and permutation toggles", {
  df <- make_class_df(12)
  df$x1[6] <- df$x1[1]
  df$x2[6] <- df$x2[1]
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()

  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = custom,
                            metrics = "auc", refit = FALSE, seed = 1)

  audit <- audit_leakage(fit, metric = "auc", B = 5, return_perm = FALSE,
                         include_z = FALSE, X_ref = df[, c("x1", "x2")],
                         sim_method = "pearson", feature_space = "rank")
  expect_equal(length(audit@perm_values), 0)
  expect_true(is.na(audit@permutation_gap$z[1]))
  expect_true(nrow(audit@duplicates) >= 1)
})

test_that("audit_leakage warns when observed metric is non-finite", {
  indices <- list(list(train = 1:2, test = 3:5, fold = 1, repeat_id = 1))
  splits <- LeakSplits(mode = "custom", indices = indices, info = list(outcome = "outcome"))
  preds <- list(data.frame(
    id = 1:3,
    truth = factor(c(0, 0, 0), levels = c(0, 1)),
    pred = c(0.1, 0.2, 0.3),
    fold = 1,
    learner = "glm",
    stringsAsFactors = FALSE
  ))
  fit <- LeakFit(
    splits = splits,
    metrics = data.frame(fold = 1, learner = "glm", auc = NA_real_),
    metric_summary = data.frame(),
    audit = data.frame(),
    predictions = preds,
    preprocess = list(),
    learners = list(),
    outcome = "outcome",
    task = "binomial",
    feature_names = character(),
    info = list()
  )
  expect_warning({
    audit <- audit_leakage(fit, metric = "auc", B = 3)
  }, "Observed metric is NA")
  expect_true(is.na(audit@permutation_gap$metric_obs[1]))
})

test_that("audit_leakage warns on misaligned coldata", {
  df <- make_class_df(12)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = custom,
                            metrics = "auc", refit = FALSE, seed = 1)

  fit@predictions <- lapply(fit@predictions, function(df_pred) {
    df_pred$id <- paste0("id_", df_pred$id)
    df_pred
  })
  bad_cd <- df
  rownames(bad_cd) <- seq_len(nrow(bad_cd))

  expect_warning({
    audit <- audit_leakage(fit, metric = "auc", B = 5, coldata = bad_cd)
  }, "not aligned")
  expect_equal(nrow(audit@batch_assoc), 0)
})

test_that("audit_leakage_by_learner validates inputs", {
  df <- make_class_df(12)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()
  custom$glm2 <- custom$glm
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = c("glm", "glm2"), custom_learners = custom,
                            metrics = "auc", refit = FALSE, seed = 1)

  audits <- audit_leakage_by_learner(fit, metric = "auc", B = 3)
  expect_true(inherits(audits, "LeakAuditList"))
  expect_equal(length(audits), 2)

  fit@predictions <- lapply(fit@predictions, function(df_pred) {
    df_pred$learner <- NULL
    df_pred
  })
  expect_error(audit_leakage_by_learner(fit, metric = "auc", learners = c("a", "b")),
               "cannot audit per learner")
})

test_that("print.LeakAuditList returns an invisible object", {
  df <- make_class_df(10)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 2, seed = 1)
  custom <- make_custom_learners()
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = custom,
                            metrics = "auc", refit = FALSE, seed = 1)
  audits <- audit_leakage_by_learner(fit, metric = "auc", B = 2)
  output <- capture.output(ret <- print(audits))
  expect_true(inherits(ret, "LeakAuditList"))
  expect_true(any(grepl("LeakAuditList", output)))
})
