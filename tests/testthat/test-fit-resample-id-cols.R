# Regression tests: identifier columns must not silently become predictors.
# Before 0.3.8 fit_resample() dropped only the outcome and declared split
# columns, so a character sample_id was one-hot encoded into n predictors.

make_id_df <- function(n = 40, seed = 1) {
  set.seed(seed)
  df <- make_class_df(n)[, c("subject", "outcome", "x1", "x2")]
  df$sample_id <- sprintf("sample_%03d", seq_len(n))
  df
}

fit_id <- function(df, ...) {
  sp <- make_split_plan(df, outcome = "outcome", mode = "subject_grouped",
                        group = "subject", v = 4, progress = FALSE)
  fit_resample(df, outcome = "outcome", splits = sp, learner = "glm",
               custom_learners = make_custom_learners(), metrics = "auc",
               refit = FALSE, seed = 1, ...)
}

test_that("an identifier-like character column triggers a warning", {
  df <- make_id_df()
  # With one dummy per row the glm folds may fail outright; the warning must
  # be raised before fitting either way.
  w <- tryCatch(fit_id(df), bioLeak_input_warning = function(w) w)
  expect_s3_class(w, "bioLeak_input_warning")
  expect_match(conditionMessage(w), "'sample_id'.*look like identifiers")
  expect_match(conditionMessage(w), 'id_cols = c\\("sample_id"\\)')
})

test_that("id_cols excludes identifier columns from the predictors", {
  df <- make_id_df()
  expect_no_warning(fit <- fit_id(df, id_cols = "sample_id"))
  expect_false(any(grepl("sample_id", fit@feature_names)))
  expect_setequal(fit@feature_names, c("x1", "x2"))
  expect_identical(fit@info$id_cols, "sample_id")
  expect_identical(fit@info$perm_refit_spec$id_cols, "sample_id")
})

test_that("low-cardinality categorical predictors do not warn", {
  df <- make_id_df()
  df$sample_id <- NULL
  df$site <- rep(c("north", "south"), length.out = nrow(df))
  expect_no_warning(fit_id(df))
})

test_that("id_cols is validated", {
  df <- make_id_df()
  expect_error(fit_id(df, id_cols = "missing_col"), "id_cols not found")
  expect_error(fit_id(df, id_cols = 1), "character vector")
})

test_that("refit permutations reuse id_cols", {
  df <- make_id_df()
  fit <- fit_id(df, id_cols = "sample_id")
  spec <- list(x = df, outcome = "outcome", learner = "glm",
               custom_learners = make_custom_learners())
  expect_no_warning(
    aud <- audit_leakage(fit, metric = "auc", B = 3, perm_refit = TRUE,
                         perm_refit_spec = spec,
                         target_scan = FALSE, target_scan_multivariate = FALSE)
  )
  expect_true(all(is.finite(aud@perm_values)))
})
