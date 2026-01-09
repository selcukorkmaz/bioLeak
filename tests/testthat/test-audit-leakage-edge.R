# test_that("audit_leakage validates predictions and learner selection", {
#   df <- make_class_df(12)
#   splits <- make_split_plan_quiet(df, outcome = "outcome",
#                               mode = "subject_grouped", group = "subject",
#                               v = 3, seed = 1)
#   custom <- make_custom_learners()
#   custom$glm2 <- custom$glm
#
#   fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
#                             learner = c("glm", "glm2"),
#                             custom_learners = custom,
#                             metrics = "auc", refit = FALSE, seed = 1)
#
#   expect_error(audit_leakage(fit, metric = "auc", B = 5),
#                "Multiple learners")
#   expect_error(audit_leakage(fit, metric = "auc", learner = "missing", B = 5),
#                "not found")
#
#   fit@predictions <- list()
#   expect_error(audit_leakage(fit, metric = "auc", B = 5),
#                "No predictions available")
# })

# test_that("audit_leakage warns when learner IDs are missing", {
#   df <- make_class_df(12)
#   splits <- make_split_plan_quiet(df, outcome = "outcome",
#                               mode = "subject_grouped", group = "subject",
#                               v = 3, seed = 1)
#   custom <- make_custom_learners()
#   custom$glm2 <- custom$glm
#
#   fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
#                             learner = c("glm", "glm2"),
#                             custom_learners = custom,
#                             metrics = "auc", refit = FALSE, seed = 1)
#   fit@predictions <- lapply(fit@predictions, function(df_pred) {
#     df_pred[df_pred$learner == "glm", , drop = FALSE]
#   })
#   fit@predictions <- Filter(function(df_pred) nrow(df_pred) > 0, fit@predictions)
#   fit@predictions <- lapply(fit@predictions, function(df_pred) {
#     df_pred$learner <- NULL
#     df_pred
#   })
#
#   audit <- expect_warning_match(
#     audit_leakage(fit, metric = "auc", B = 5),
#     "predictions lack learner IDs"
#   )
#
#   audit <- expect_warning_match(
#     audit_leakage(fit, metric = "auc", learner = "glm", B = 5),
#     "ignored"
#   )
# })

# test_that("audit_leakage handles duplicate detection and permutation toggles", {
#   df <- make_class_df(12)
#   df$x1[6] <- df$x1[1]
#   df$x2[6] <- df$x2[1]
#   splits <- make_split_plan_quiet(df, outcome = "outcome",
#                               mode = "subject_grouped", group = "subject",
#                               v = 3, seed = 1)
#   custom <- make_custom_learners()
#
#   fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
#                             learner = "glm", custom_learners = custom,
#                             metrics = "auc", refit = FALSE, seed = 1)
#
#   audit <- audit_leakage(fit, metric = "auc", B = 5, return_perm = FALSE,
#                          include_z = FALSE, X_ref = df[, c("x1", "x2")],
#                          sim_method = "pearson", feature_space = "rank",
#                          perm_stratify = FALSE)
#   expect_equal(length(audit@perm_values), 0)
#   expect_true(is.na(audit@permutation_gap$z[1]))
#   expect_true(nrow(audit@duplicates) >= 1)
# })

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

# test_that("audit_leakage warns on misaligned coldata for permutations", {
#   df <- make_class_df(12)
#   splits <- make_split_plan_quiet(df, outcome = "outcome",
#                               mode = "subject_grouped", group = "subject",
#                               v = 3, seed = 1)
#   fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
#                             learner = "glm", custom_learners = make_custom_learners(),
#                             metrics = "auc", refit = FALSE, seed = 1)
#
#   bad_cd <- df
#   rownames(bad_cd) <- paste0("id_", seq_len(nrow(bad_cd)))
#   expect_warning_match(
#     audit_leakage(fit, metric = "auc", B = 5, coldata = bad_cd,
#                   perm_stratify = TRUE),
#     "restricted permutations disabled"
#   )
# })

# test_that("audit_leakage duplicate_scope filters within-fold duplicates", {
#   df <- make_class_df(12)
#   splits <- make_split_plan_quiet(df, outcome = "outcome",
#                               mode = "subject_grouped", group = "subject",
#                               v = 3, seed = 1)
#   fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
#                             learner = "glm", custom_learners = make_custom_learners(),
#                             metrics = "auc", refit = FALSE, seed = 1)
#
#   X_ref <- df[, c("x1", "x2")]
#   for (sid in unique(df$subject)) {
#     idx <- which(df$subject == sid)
#     if (length(idx) > 1) {
#       X_ref[idx[2], ] <- X_ref[idx[1], ]
#     }
#   }
#
#   audit_tt <- audit_leakage(fit, metric = "auc", B = 3,
#                             X_ref = X_ref, sim_threshold = 0.999,
#                             perm_stratify = FALSE)
#   expect_equal(nrow(audit_tt@duplicates), 0)
#
#   audit_all <- audit_leakage(fit, metric = "auc", B = 3,
#                              X_ref = X_ref, sim_threshold = 0.999,
#                              perm_stratify = FALSE, duplicate_scope = "all")
#   expect_true(nrow(audit_all@duplicates) > 0)
#   expect_true("cross_fold" %in% names(audit_all@duplicates))
#   expect_true(all(audit_all@duplicates$cross_fold == FALSE))
# })

# test_that("audit_leakage warns on misaligned coldata", {
#   df <- make_class_df(12)
#   splits <- make_split_plan_quiet(df, outcome = "outcome",
#                               mode = "subject_grouped", group = "subject",
#                               v = 3, seed = 1)
#   custom <- make_custom_learners()
#   fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
#                             learner = "glm", custom_learners = custom,
#                             metrics = "auc", refit = FALSE, seed = 1)
#
#   fit@predictions <- lapply(fit@predictions, function(df_pred) {
#     df_pred$id <- paste0("id_", df_pred$id)
#     df_pred
#   })
#   bad_cd <- df
#   rownames(bad_cd) <- seq_len(nrow(bad_cd))
#
#   audit <- expect_warning_match(
#     audit_leakage(fit, metric = "auc", B = 5, coldata = bad_cd,
#                   perm_stratify = FALSE),
#     "not aligned"
#   )
#   expect_equal(nrow(audit@batch_assoc), 0)
# })

# test_that("audit_leakage_by_learner validates inputs", {
#   df <- make_class_df(12)
#   splits <- make_split_plan_quiet(df, outcome = "outcome",
#                               mode = "subject_grouped", group = "subject",
#                               v = 3, seed = 1)
#   custom <- make_custom_learners()
#   custom$glm2 <- custom$glm
#   fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
#                             learner = c("glm", "glm2"), custom_learners = custom,
#                             metrics = "auc", refit = FALSE, seed = 1)
#
#   audits <- audit_leakage_by_learner(fit, metric = "auc", B = 3,
#                                      perm_stratify = FALSE)
#   expect_true(inherits(audits, "LeakAuditList"))
#   expect_equal(length(audits), 2)
#
#   fit@predictions <- lapply(fit@predictions, function(df_pred) {
#     df_pred$learner <- NULL
#     df_pred
#   })
#   expect_error(audit_leakage_by_learner(fit, metric = "auc", learners = c("a", "b")),
#                "cannot audit per learner")
# })

# test_that("print.LeakAuditList returns an invisible object", {
#   df <- make_class_df(10)
#   splits <- make_split_plan_quiet(df, outcome = "outcome",
#                               mode = "subject_grouped", group = "subject",
#                               v = 2, seed = 1)
#   custom <- make_custom_learners()
#   fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
#                             learner = "glm", custom_learners = custom,
#                             metrics = "auc", refit = FALSE, seed = 1)
#   audits <- audit_leakage_by_learner(fit, metric = "auc", B = 2,
#                                      perm_stratify = FALSE)
#   output <- capture.output(ret <- print(audits))
#   expect_true(inherits(ret, "LeakAuditList"))
#   expect_true(any(grepl("LeakAuditList", output)))
# })

