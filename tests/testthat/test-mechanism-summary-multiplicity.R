# The confounding_alignment rule (JSS review, L1.3) uses the batch-outcome
# test, run once per batch column (so repeats cannot inflate it) and
# Holm-adjusted across columns, and requires batches shared by train and
# test (span_frac > 0). The fold x batch test does not enter the rule.

batch_rows <- function(outcome_p, effect = 0.3, span = 0.8, fold_p = 0.5,
                       cols = paste0("b", seq_along(outcome_p)), repeats = 1L) {
  base <- data.frame(batch_col = cols, outcome_test = "chisq", outcome_pval = outcome_p,
                     outcome_pval_adj = stats::p.adjust(outcome_p, "holm"),
                     outcome_effect = effect, outcome_effect_min = 0.1,
                     span_frac = span, stringsAsFactors = FALSE)
  out <- base[rep(seq_len(nrow(base)), each = repeats), , drop = FALSE]
  out$repeat_id <- rep(seq_len(repeats), times = nrow(base))
  out$stat <- 1; out$df <- 4; out$pval <- fold_p; out$cramer_v <- 1
  out
}

confounding_row <- function(batch_df) {
  ms <- bioLeak:::.bio_mechanism_summary(
    perm_df = data.frame(), batch_df = batch_df,
    target_df = data.frame(), dup_df = data.frame()
  )
  ms[ms$mechanism_class == "confounding_alignment", ]
}

test_that("a nominal p < 0.05 for one of several batch columns is Holm-corrected", {
  p <- c(0.03, 0.5, 0.7)
  row <- confounding_row(batch_rows(p))
  expect_false(row$flagged)
  expect_equal(row$p_value, min(stats::p.adjust(p, "holm")))
})

test_that("repeats do not change the batch-outcome evidence", {
  row1 <- confounding_row(batch_rows(0.03, repeats = 1))
  row10 <- confounding_row(batch_rows(0.03, repeats = 10))
  expect_true(row1$flagged)
  expect_identical(row1$flagged, row10$flagged)
  expect_equal(row1$p_value, row10$p_value)
})

test_that("a strong association with shared batches is flagged", {
  expect_true(confounding_row(batch_rows(c(0.001, 0.5)))$flagged)
})

test_that("the effect-size floor applies", {
  expect_false(confounding_row(batch_rows(0.0001, effect = 0.05))$flagged)
})

test_that("batches never shared by train and test are not flagged", {
  expect_false(confounding_row(batch_rows(0.0001, span = 0))$flagged)
})

test_that("fold-batch alignment alone is not flagged", {
  expect_false(confounding_row(batch_rows(0.6, fold_p = 1e-10))$flagged)
})

test_that("NA p-values do not enlarge the family", {
  expect_equal(bioLeak:::.holm_adjust(c(0.01, NA, 0.04)), c(0.02, NA, 0.04))
})

test_that("audit_leakage reports Holm-adjusted fold-batch and batch-outcome p-values", {
  set.seed(5)
  df <- make_class_df(60)
  sp <- make_split_plan(df, outcome = "outcome", mode = "subject_grouped",
                        group = "subject", v = 3, repeats = 4, progress = FALSE)
  fit <- fit_resample(df, outcome = "outcome", splits = sp, learner = "glm",
                      custom_learners = make_custom_learners(), metrics = "auc",
                      refit = FALSE, seed = 1)
  aud <- suppressMessages(audit_leakage(fit, metric = "auc", B = 5,
                                        coldata = df, batch_cols = "batch",
                                        target_scan = FALSE,
                                        target_scan_multivariate = FALSE))
  ba <- aud@batch_assoc
  expect_true("pval_adj" %in% names(ba))
  expect_equal(nrow(ba), 4L)
  expect_equal(ba$pval_adj, stats::p.adjust(ba$pval, "holm"))
  expect_equal(length(unique(ba$outcome_pval)), 1L)   # once per column
  ms <- aud@info$mechanism_summary
  expect_equal(ms$p_value[ms$mechanism_class == "confounding_alignment"],
               min(ba$outcome_pval_adj[ba$span_frac > 0]))
})
