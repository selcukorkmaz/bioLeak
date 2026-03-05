test_that("audit_leakage warns when observed metric is non-finite", {
  indices <- list(list(train = 1:2, test = 3:4, fold = 1, repeat_id = 1))
  splits <- bioLeak:::LeakSplits(mode = "custom", indices = indices,
                                 info = list(outcome = "outcome"))
  preds <- list(data.frame(
    id = 1:2,
    truth = c(1, 2),
    pred = c(NA_real_, 1),
    fold = 1,
    learner = "glm",
    stringsAsFactors = FALSE
  ))
  fit <- bioLeak:::LeakFit(
    splits = splits,
    metrics = data.frame(fold = 1, learner = "glm", rmse = NA_real_),
    metric_summary = data.frame(),
    audit = data.frame(),
    predictions = preds,
    preprocess = list(),
    learners = list(),
    outcome = "outcome",
    task = "gaussian",
    feature_names = character(),
    info = list()
  )
  audit <- expect_warning_match(
    audit_leakage(fit, metric = "rmse", B = 3),
    "Observed metric is NA"
  )
  expect_true(is.na(audit@permutation_gap$metric_obs[1]))
})
