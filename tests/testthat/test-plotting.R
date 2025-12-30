# test_that("plotting helpers validate ggplot2 availability", {
#   df <- make_class_df(10)
#   splits <- make_splits_quiet(df, outcome = "outcome",
#                               mode = "subject_grouped", group = "subject",
#                               v = 2, seed = 1)
#   custom <- make_custom_learners()
#   fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
#                             learner = "glm", custom_learners = custom,
#                             metrics = "accuracy", refit = FALSE, seed = 1)
#   audit <- audit_leakage(fit, metric = "auc", B = 3, perm_stratify = FALSE)
#
#   if (requireNamespace("ggplot2", quietly = TRUE)) {
#     expect_true(TRUE)
#   } else {
#     expect_error(plot_perm_distribution(audit), "ggplot2")
#     expect_error(plot_fold_balance(fit), "ggplot2")
#     expect_error(plot_overlap_checks(fit, column = "subject"), "ggplot2")
#     expect_error(plot_time_acf(fit), "ggplot2")
#   }
# })

# test_that("plotting helpers return objects when ggplot2 is installed", {
#   skip_if_not_installed("ggplot2")
#   df <- make_class_df(12)
#   splits <- make_splits_quiet(df, outcome = "outcome",
#                               mode = "subject_grouped", group = "subject",
#                               v = 3, seed = 1)
#   custom <- make_custom_learners()
#   fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
#                             learner = "glm", custom_learners = custom,
#                             metrics = "auc", refit = FALSE, seed = 1)
#   audit <- audit_leakage(fit, metric = "auc", B = 5, perm_stratify = FALSE)
#
#   with_temp_plot_device({
#     res1 <- plot_perm_distribution(audit)
#     expect_true(all(c("observed", "permuted_mean", "plot") %in% names(res1)))
#
#     res2 <- plot_fold_balance(fit)
#     expect_true(all(c("fold_summary", "plot") %in% names(res2)))
#
#     res3 <- plot_overlap_checks(fit, column = "subject")
#     expect_true(all(c("overlap_counts", "plot") %in% names(res3)))
#
#     res4 <- plot_time_acf(fit, lag.max = 3)
#     expect_true(all(c("acf", "plot") %in% names(res4)))
#   })
# })

# test_that("plot_overlap_checks and plot_time_acf validate inputs", {
#   skip_if_not_installed("ggplot2")
#   df <- make_class_df(8)
#   df <- df[, c("subject", "outcome", "x1", "x2")]
#   splits <- make_splits_quiet(df, outcome = "outcome",
#                               mode = "subject_grouped", group = "subject",
#                               v = 2, seed = 1)
#   custom <- make_custom_learners()
#   fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
#                             learner = "glm", custom_learners = custom,
#                             metrics = "accuracy", refit = FALSE, seed = 1)
#   expect_error(plot_overlap_checks(fit, column = "missing"), "Column not available")
#
#   fit@predictions <- list(data.frame(pred = 1, truth = factor(1), id = 1, fold = 1, learner = "glm"))
#   expect_error(plot_time_acf(fit, lag.max = 2), "Not enough finite predictions")
# })

# test_that("plot_perm_distribution errors without permutation values", {
#   skip_if_not_installed("ggplot2")
#   df <- make_class_df(10)
#   splits <- make_splits_quiet(df, outcome = "outcome",
#                               mode = "subject_grouped", group = "subject",
#                               v = 2, seed = 1)
#   custom <- make_custom_learners()
#   fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
#                             learner = "glm", custom_learners = custom,
#                             metrics = "auc", refit = FALSE, seed = 1)
#   audit <- audit_leakage(fit, metric = "auc", B = 3, return_perm = FALSE,
#                          perm_stratify = FALSE)
#   expect_error(plot_perm_distribution(audit), "No finite permutation values")
# })
