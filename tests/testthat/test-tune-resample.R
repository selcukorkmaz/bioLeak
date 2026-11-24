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

tune_resample_quiet <- function(...) {
  out <- NULL
  capture.output({
    out <- suppressWarnings(tune_resample(...))
  })
  out
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

  tuned_best <- tune_resample_quiet(df, outcome = "outcome", splits = splits,
                                    learner = spec, preprocess = rec, grid = 2,
                                    metrics = "accuracy", selection = "best", seed = 11)
  tuned_ose <- tune_resample_quiet(df, outcome = "outcome", splits = splits,
                                   learner = spec, preprocess = rec, grid = 2,
                                   metrics = "accuracy", selection = "one_std_err", seed = 11)

  expect_true(nrow(tuned_best$best_params) > 0)
  expect_true(nrow(tuned_ose$best_params) > 0)
  expect_true(all(tuned_best$best_params$penalty == 0.01))
  expect_true(all(tuned_ose$best_params$penalty == 1.0))
})

test_that("tune_resample supports final refit and stores fold status", {
  skip_if_tune_deps()

  df <- make_class_df(80)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                                  mode = "subject_grouped", group = "subject",
                                  v = 2, nested = TRUE, stratify = FALSE, seed = 2)

  spec <- parsnip::logistic_reg(penalty = tune::tune(), mixture = 1) |>
    parsnip::set_engine("glmnet")
  rec <- recipes::recipe(outcome ~ x1 + x2, data = df)

  tuned <- tune_resample_quiet(df, outcome = "outcome", splits = splits,
                               learner = spec, preprocess = rec, grid = 2,
                               metrics = "accuracy", selection = "best",
                               refit = TRUE, seed = 2)

  expect_true(isTRUE(tuned$info$refit))
  expect_true(!is.null(tuned$final_model))
  expect_s3_class(tuned$final_workflow, "workflow")
  expect_true(nrow(tuned$final_params) > 0)
  expect_true(is.data.frame(tuned$fold_status))
  expect_true(nrow(tuned$fold_status) == length(splits@indices))
  expect_true(all(tuned$fold_status$status %in% c("success", "skipped", "failed")))
})

test_that("final refit uses aggregated params, not single best outer fold", {
  skip_if_tune_deps()

  df <- make_class_df(80)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                                  mode = "subject_grouped", group = "subject",
                                  v = 2, nested = TRUE, stratify = FALSE, seed = 2)

  spec <- parsnip::logistic_reg(penalty = tune::tune(), mixture = 1) |>
    parsnip::set_engine("glmnet")
  rec <- recipes::recipe(outcome ~ x1 + x2, data = df)

  tuned <- tune_resample_quiet(df, outcome = "outcome", splits = splits,
                               learner = spec, preprocess = rec, grid = 3,
                               metrics = "accuracy", selection = "best",
                               refit = TRUE, seed = 2)

  # refit_method must be "aggregate", not a single-fold selection
  expect_identical(tuned$info$refit_method, "aggregate")

  # refit_fold must be NA — no single fold was selected based on outer metrics
  expect_true(is.na(tuned$info$refit_fold))

  # final_params must have exactly the hyperparameter columns, no fold column
  expect_false("fold" %in% names(tuned$final_params))

  # --- Behavioral proof that outer test metrics are not consulted ---
  # final_params must exactly equal the coordinate-wise aggregate of ALL folds'
  # best_params. This is the core independence guarantee: the old (leaky) code
  # did which.max(metrics_df[[metric]]) to pick a single fold, so final_params
  # would equal that fold's row. The new code aggregates all folds equally.
  bp <- tuned$best_params
  param_cols <- setdiff(names(bp), "fold")
  expect_true(nrow(bp) >= 2,
              info = "Need at least 2 outer folds to test aggregation")
  expect_true(length(param_cols) >= 1,
              info = "Need at least 1 hyperparameter column")

  for (col in param_cols) {
    if (is.numeric(bp[[col]])) {
      expected <- stats::median(bp[[col]], na.rm = TRUE)
      expect_equal(tuned$final_params[[col]], expected,
                   info = paste("Param", col,
                                "must equal median across ALL folds, not a single fold's value"))
    } else {
      tbl <- table(bp[[col]])
      expected <- names(tbl)[which.max(tbl)]
      expect_equal(as.character(tuned$final_params[[col]]), expected,
                   info = paste("Param", col,
                                "must equal majority vote across ALL folds"))
    }
  }

  # Verify final_params is not simply copied from the best-metric fold.
  # Identify the fold that had the best outer metric — if it differs from the
  # aggregate, the old (leaky) logic would have returned that fold's params.
  metric_col <- tuned$info$selection_metric
  if (!is.null(metric_col) && metric_col %in% names(tuned$metrics)) {
    metric_vals <- tuned$metrics[[metric_col]]
    minimize <- metric_col %in% c("rmse", "mae", "log_loss", "mn_log_loss")
    best_idx <- if (minimize) which.min(metric_vals) else which.max(metric_vals)
    best_fold <- tuned$metrics$fold[best_idx]
    best_fold_params <- bp[bp$fold == best_fold, param_cols, drop = FALSE]

    # Check if folds actually disagree on params
    folds_disagree <- FALSE
    for (col in param_cols) {
      if (length(unique(bp[[col]])) > 1L) {
        folds_disagree <- TRUE
        break
      }
    }

    if (folds_disagree && nrow(best_fold_params) == 1L) {
      # When folds disagree, the aggregate should NOT equal the best-metric
      # fold's params. This is the direct proof: the old code would have
      # returned best_fold_params, the new code returns the aggregate.
      params_match_best <- all(vapply(param_cols, function(col) {
        identical(tuned$final_params[[col]], best_fold_params[[col]])
      }, logical(1)))
      expect_false(params_match_best,
                   info = paste("final_params must not equal fold", best_fold,
                                "params (the best-metric fold) — that would be",
                                "nested-CV leakage"))
    }
  }
})

test_that("LeakTune summary and audit handle skipped outer folds", {
  skip_if_tune_deps()

  set.seed(1)
  df <- data.frame(
    subject = paste0("s", seq_len(60)),
    outcome = factor(c(rep(0, 30), rep(1, 30)), levels = c(0, 1)),
    x1 = rnorm(60),
    x2 = rnorm(60),
    stringsAsFactors = FALSE
  )

  indices <- list(
    list(train = 1:20, test = 21:60, fold = 1L, repeat_id = 1L),
    list(train = c(1:20, 31:60), test = 21:30, fold = 2L, repeat_id = 1L)
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
      train_n = c(20L, 50L),
      test_n = c(40L, 10L)
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

  tuned <- tune_resample_quiet(df, outcome = "outcome", splits = splits,
                               learner = spec, preprocess = rec, grid = 2,
                               metrics = "accuracy", seed = 1)

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

test_that("tune_resample tunes binomial thresholds from inner predictions", {
  skip_if_tune_deps()

  df <- make_class_df(24)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                                  mode = "subject_grouped", group = "subject",
                                  v = 2, nested = TRUE, stratify = FALSE, seed = 21)

  spec <- parsnip::logistic_reg(penalty = tune::tune(), mixture = 1) |>
    parsnip::set_engine("glmnet")
  rec <- recipes::recipe(outcome ~ x1 + x2, data = df)

  seen_save_pred <- logical()
  seen_thresholds <- numeric()

  local_mocked_bindings(
    tune_grid = function(..., control) {
      seen_save_pred <<- c(seen_save_pred, isTRUE(control$save_pred))
      structure(list(), class = "mock_tune_results")
    },
    collect_metrics = function(...) {
      data.frame(
        .metric = rep("accuracy", 4),
        .estimator = rep("binary", 4),
        .estimate = c(0.8, 0.8, 0.7, 0.7),
        .config = c("cfg_a", "cfg_a", "cfg_b", "cfg_b"),
        penalty = c(0.01, 0.01, 1.0, 1.0),
        stringsAsFactors = FALSE
      )
    },
    collect_predictions = function(...) {
      data.frame(
        .config = rep("cfg_a", 4),
        outcome = factor(c(0, 1, 0, 1), levels = c(0, 1)),
        .pred_0 = c(0.9, 0.6, 0.55, 0.4),
        .pred_1 = c(0.1, 0.4, 0.45, 0.6),
        stringsAsFactors = FALSE
      )
    },
    finalize_workflow = function(x, parameters, ...) x,
    .package = "tune"
  )
  local_mocked_bindings(
    fit_resample = function(x, outcome, splits, classification_threshold = 0.5, ...) {
      seen_thresholds <<- c(seen_thresholds, classification_threshold)
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
          pred = c(0.1, 0.4, 0.45, 0.6),
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

  tuned <- tune_resample_quiet(df, outcome = "outcome", splits = splits,
                               learner = spec, preprocess = rec, grid = 2,
                               metrics = "accuracy", selection = "best", seed = 21,
                               tune_threshold = TRUE, threshold_grid = c(0.2, 0.8))

  expect_true(all(seen_save_pred))
  expect_equal(length(seen_thresholds), length(splits@indices))
  expect_true(all(abs(seen_thresholds - 0.2) < 1e-12))
  expect_true(is.data.frame(tuned$thresholds))
  expect_equal(nrow(tuned$thresholds), length(splits@indices))
  expect_true(all(abs(tuned$thresholds$threshold - 0.2) < 1e-12))
  expect_true(isTRUE(tuned$info$threshold_tuned))
  expect_equal(tuned$info$threshold_metric, "accuracy")
})

test_that("tune_resample keeps the split columns out of the predictors", {
  skip_if_not_installed("tune")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("recipes")
  set.seed(3); n <- 120
  df <- data.frame(subj = paste0("P", rep(1:40, each = 3)))
  for (j in 1:5) df[[paste0("x", j)]] <- rnorm(n)
  df$y <- factor(rbinom(n, 1, plogis(df$x1)))
  sp <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subj", v = 4,
                        progress = FALSE)
  spec <- parsnip::logistic_reg(penalty = tune::tune(), mixture = 1) |>
    parsnip::set_engine("glmnet")
  tr <- suppressWarnings(tune_resample(df, outcome = "y", splits = sp, learner = spec, grid = 2,
                                       metrics = "auc", inner_v = 3, seed = 1, refit = TRUE))
  expect_equal(nrow(tr$metrics), length(sp@indices))
  mold <- workflows::extract_mold(tr$final_model)
  expect_false("subj" %in% names(mold$predictors))
  rec <- recipes::recipe(y ~ ., data = df)
  tr2 <- suppressWarnings(tune_resample(df, outcome = "y", splits = sp, learner = spec, grid = 2,
                                        preprocess = rec, metrics = "auc", inner_v = 3, seed = 1))
  expect_equal(nrow(tr2$metrics), length(sp@indices))
})

test_that("tune_resample accepts a workflow whose recipe gives a split column the id role", {
  skip_if_not_installed("tune")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("recipes")
  skip_if_not_installed("workflows")
  set.seed(5); n <- 120
  dat <- data.frame(y = factor(sample(c("a", "b"), n, TRUE)), subject = paste0("P", rep(1:30, each = 4)),
                    x1 = rnorm(n), x2 = rnorm(n))
  sp <- make_split_plan(dat, outcome = "y", mode = "subject_grouped", group = "subject", v = 3,
                        progress = FALSE)
  rec <- recipes::recipe(y ~ ., data = dat) |> recipes::update_role(subject, new_role = "id")
  wf <- workflows::workflow() |> workflows::add_recipe(rec) |>
    workflows::add_model(parsnip::logistic_reg(penalty = tune::tune(), mixture = 1) |>
                           parsnip::set_engine("glmnet"))
  tr <- suppressWarnings(tune_resample(dat, outcome = "y", splits = sp, learner = wf, grid = 2,
                                       metrics = "auc", inner_v = 3, seed = 1))
  expect_equal(nrow(tr$metrics), length(sp@indices))
})

test_that("tune_resample ranks tuning settings by a correctly oriented AUC", {
  skip_if_not_installed("tune")
  skip_if_not_installed("glmnet")
  set.seed(1); n <- 240
  df <- data.frame(id = 1:n, x1 = rnorm(n), x2 = rnorm(n))
  df$y <- factor(rbinom(n, 1, plogis(2 * df$x1)))
  sp <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "id", v = 3,
                        progress = FALSE)
  spec <- parsnip::logistic_reg(penalty = tune::tune(), mixture = 1) |>
    parsnip::set_engine("glmnet")
  for (pc in c("1", "0")) {
    tr <- suppressWarnings(tune_resample(df, outcome = "y", splits = sp, learner = spec,
                                         grid = data.frame(penalty = c(0.001, 0.3)),
                                         metrics = "auc", positive_class = pc, seed = 1))
    inner <- tune::collect_metrics(tr$inner_results[[1]])
    expect_gt(max(inner$mean[inner$.metric == "roc_auc"]), 0.7, label = pc)
    expect_true(all(tr$best_params$penalty == 0.001), label = pc)
    expect_gt(tr$metric_summary$roc_auc_mean, 0.75, label = pc)
  }
})

test_that("tune_resample accepts numeric 0/1 outcomes, as_leaksplits plans and rsample input", {
  skip_if_not_installed("tune")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("rsample")
  set.seed(2); n <- 120
  df <- data.frame(subject = rep(1:40, each = 3), x1 = rnorm(n), x2 = rnorm(n))
  df$y <- rbinom(n, 1, plogis(1.5 * df$x1))  # numeric 0/1
  spec <- parsnip::logistic_reg(penalty = tune::tune(), mixture = 1) |>
    parsnip::set_engine("glmnet")
  g <- data.frame(penalty = c(0.001, 0.05))
  sp <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject", v = 3,
                        progress = FALSE)
  tr <- suppressWarnings(tune_resample(df, "y", sp, learner = spec, grid = g, metrics = "auc", seed = 1))
  expect_equal(nrow(tr$metrics), 3L)

  df$y <- factor(df$y)
  df$sample_id <- paste0("S", seq_len(n))
  sd <- data.frame(sample_id = df$sample_id, group_id = paste0("P", df$subject), stringsAsFactors = FALSE)
  spec_sg <- structure(list(sample_data = sd, group_var = "group_id", constraint_mode = "subject",
                            ordering_required = FALSE), class = "split_spec")
  dat <- df[, c("sample_id", "y", "x1", "x2")]
  sp_al <- suppressWarnings(as_leaksplits(spec_sg, dat, "y", v = 3, progress = FALSE))
  tr2 <- suppressWarnings(tune_resample(dat, "y", sp_al, learner = spec, grid = g, metrics = "auc", seed = 1))
  expect_equal(nrow(tr2$metrics), 3L)

  df$sample_id <- NULL  # a per-row id would be a predictor below
  gv <- rsample::group_vfold_cv(df, group = subject, v = 3)  # design inferred, no attribute
  tr3 <- suppressWarnings(tune_resample(df, "y", gv, learner = spec, grid = g, metrics = "auc"))
  expect_equal(nrow(tr3$metrics), 3L)

  nc <- rsample::nested_cv(df, outside = rsample::group_vfold_cv(group = subject, v = 3),
                           inside = rsample::group_vfold_cv(group = subject, v = 3))
  sp_nc <- bioLeak:::.bio_as_leaksplits_from_rsample(nc, n = n, coldata = df)
  for (k in seq_along(sp_nc@indices)) {
    for (f in sp_nc@info$inner[[k]]) {
      expect_true(all(c(f$train, f$test) %in% sp_nc@indices[[k]]$train))
      expect_length(intersect(df$subject[f$train], df$subject[f$test]), 0L)
    }
  }
  tr4 <- suppressWarnings(tune_resample(df, "y", nc, learner = spec, grid = g, metrics = "auc"))
  expect_equal(nrow(tr4$metrics), 3L)
})

test_that("tune_resample handles fine-grained rsample resamples and design-free nested_cv", {
  skip_if_not_installed("tune")
  skip_if_not_installed("glmnet")
  skip_if_not_installed("rsample")
  set.seed(5); n <- 120
  d <- data.frame(subject = sprintf("s%02d", rep(1:30, each = 4)), x1 = rnorm(n), x2 = rnorm(n))
  d$y <- factor(rbinom(n, 1, plogis(1.5 * d$x1)))
  spec <- parsnip::logistic_reg(penalty = tune::tune(), mixture = 1) |>
    parsnip::set_engine("glmnet")
  # 25 bootstrap resamples -> 25 tiny inner folds; AUC is dropped from tuning
  # with a warning instead of failing the whole call
  tr <- suppressWarnings(tune_resample(d, "y", rsample::bootstraps(d, times = 25), learner = spec,
                                       grid = 2, seed = 1))
  expect_equal(nrow(tr$metrics), 25L)
  d0 <- d[, c("x1", "x2", "y")]
  nc <- rsample::nested_cv(d0, outside = rsample::vfold_cv(v = 3), inside = rsample::vfold_cv(v = 3))
  tr2 <- suppressWarnings(tune_resample(d0, "y", nc, learner = spec, grid = 2, metrics = "auc"))
  expect_equal(nrow(tr2$metrics), 3L)
})
