skip_if_tune_deps <- function() {
  skip_if_not_installed("tune")
  skip_if_not_installed("dials")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("recipes")
  skip_if_not_installed("rsample")
  skip_if_not_installed("yardstick")
  skip_if_not_installed("workflows")
  skip_if_not_installed("parsnip")
}

test_that("tune_resample selects a deterministic simple model with one_std_err", {
  skip_if_tune_deps()

  df <- make_class_df(24)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                                  mode = "subject_grouped", group = "subject",
                                  v = 2, nested = TRUE, stratify = FALSE, seed = 11)

  spec <- parsnip::logistic_reg(penalty = tune::tune(), mixture = 1) |>
    parsnip::set_engine("glmnet")
  rec <- recipes::recipe(outcome ~ x1 + x2, data = df)

  local_mocked_bindings(
    tune_grid = function(...) structure(list(), class = "mock_tune_results"),
    collect_metrics = function(...) {
      data.frame(
        .metric = rep("accuracy", 4),
        .estimator = rep("binary", 4),
        .estimate = c(0.80, 0.76, 0.77, 0.77),
        .config = c("cfg_a", "cfg_a", "cfg_b", "cfg_b"),
        penalty = c(0.01, 0.01, 1.0, 1.0),
        stringsAsFactors = FALSE
      )
    },
    finalize_workflow = function(x, parameters, ...) x,
    .package = "tune"
  )
  local_mocked_bindings(
    fit_resample = function(x, outcome, splits, ...) {
      truth <- factor(rep(c(0, 1), length.out = 4), levels = c(0, 1))
      methods::new(
        "LeakFit",
        splits = splits,
        metrics = data.frame(fold = 1, learner = "mock", accuracy = 0.75),
        metric_summary = data.frame(learner = "mock", accuracy_mean = 0.75, accuracy_sd = 0),
        audit = data.frame(),
        predictions = list(data.frame(
          id = as.character(seq_len(4)),
          truth = truth,
          pred = c(0.45, 0.55, 0.40, 0.60),
          fold = 1,
          learner = "mock",
          stringsAsFactors = FALSE
        )),
        preprocess = list(),
        learners = list(),
        outcome = outcome,
        task = "binomial",
        feature_names = c("x1", "x2"),
        info = list(sample_ids = as.character(seq_len(nrow(x))))
      )
    },
    .package = "bioLeak"
  )

  tuned_best <- suppressWarnings(
    tune_resample(df, outcome = "outcome", splits = splits,
                  learner = spec, preprocess = rec, grid = 2,
                  metrics = "accuracy", selection = "best", seed = 11)
  )
  tuned_ose <- suppressWarnings(
    tune_resample(df, outcome = "outcome", splits = splits,
                  learner = spec, preprocess = rec, grid = 2,
                  metrics = "accuracy", selection = "one_std_err", seed = 11)
  )

  expect_true(nrow(tuned_best$best_params) > 0)
  expect_true(nrow(tuned_ose$best_params) > 0)
  expect_true(all(tuned_best$best_params$penalty == 0.01))
  expect_true(all(tuned_ose$best_params$penalty == 1.0))
})

test_that("tune_resample supports final refit and stores fold status", {
  skip_if_tune_deps()

  df <- make_class_df(40)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                                  mode = "subject_grouped", group = "subject",
                                  v = 2, nested = TRUE, stratify = FALSE, seed = 2)

  spec <- parsnip::logistic_reg(penalty = tune::tune(), mixture = 1) |>
    parsnip::set_engine("glmnet")
  rec <- recipes::recipe(outcome ~ x1 + x2, data = df)

  tuned <- suppressWarnings(
    tune_resample(df, outcome = "outcome", splits = splits,
                  learner = spec, preprocess = rec, grid = 2,
                  metrics = "accuracy", selection = "best",
                  refit = TRUE, seed = 2)
  )

  expect_true(isTRUE(tuned$info$refit))
  expect_true(!is.null(tuned$final_model))
  expect_s3_class(tuned$final_workflow, "workflow")
  expect_true(nrow(tuned$final_params) > 0)
  expect_true(is.data.frame(tuned$fold_status))
  expect_true(nrow(tuned$fold_status) == length(splits@indices))
  expect_true(all(tuned$fold_status$status %in% c("success", "skipped", "failed")))
})

test_that("LeakTune summary and audit handle skipped outer folds", {
  skip_if_tune_deps()

  set.seed(1)
  df <- data.frame(
    subject = paste0("s", seq_len(20)),
    outcome = factor(c(rep(0, 12), rep(1, 8)), levels = c(0, 1)),
    x1 = rnorm(20),
    x2 = rnorm(20),
    stringsAsFactors = FALSE
  )

  indices <- list(
    list(train = 1:12, test = 13:20, fold = 1L, repeat_id = 1L),
    list(train = c(1:10, 13:20), test = 11:12, fold = 2L, repeat_id = 1L)
  )
  split_info <- list(
    outcome = "outcome",
    v = 2L,
    repeats = 1L,
    seed = 1L,
    mode = "subject_grouped",
    perm_mode = "subject_grouped",
    group = "subject",
    batch = NULL,
    study = NULL,
    time = NULL,
    stratify = FALSE,
    nested = FALSE,
    horizon = 0,
    summary = data.frame(
      fold = c(1L, 2L),
      repeat_id = c(1L, 1L),
      train_n = c(12L, 18L),
      test_n = c(8L, 2L)
    ),
    hash = "manual",
    inner = NULL,
    compact = FALSE,
    fold_assignments = NULL,
    coldata = df
  )
  splits <- bioLeak:::LeakSplits(mode = "subject_grouped", indices = indices, info = split_info)

  spec <- parsnip::logistic_reg(penalty = tune::tune(), mixture = 1) |>
    parsnip::set_engine("glmnet")
  rec <- recipes::recipe(outcome ~ x1 + x2, data = df)

  tuned <- suppressWarnings(
    tune_resample(df, outcome = "outcome", splits = splits,
                  learner = spec, preprocess = rec, grid = 2,
                  metrics = "accuracy", seed = 1)
  )

  expect_true(is.null(tuned$outer_fits[[1]]))
  expect_true(!is.null(tuned$outer_fits[[2]]))
  expect_true(is.data.frame(tuned$fold_status))
  expect_equal(nrow(tuned$fold_status), 2)
  expect_true(any(tuned$fold_status$status == "skipped"))
  expect_true(any(tuned$fold_status$status == "success"))

  sum_out <- paste(capture.output(summary(tuned)), collapse = "\n")
  expect_match(sum_out, "Outer Folds: 1 successful / 2 total", fixed = TRUE)

  audit <- suppressWarnings(expect_no_error(
    audit_leakage(
      tuned,
      metric = "accuracy",
      B = 3,
      boot_B = 20,
      perm_stratify = FALSE,
      return_perm = FALSE,
      target_scan = FALSE,
      target_scan_multivariate = FALSE
    )
  ))
  expect_s4_class(audit, "LeakAudit")
})
