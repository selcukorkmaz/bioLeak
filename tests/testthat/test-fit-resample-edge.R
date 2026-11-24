test_that("fit_resample validates custom learner specifications", {
  df <- make_class_df(10)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 2, seed = 1)

  expect_error(fit_resample_quiet(df, outcome = "outcome", splits = splits,
                                  learner = "glm", custom_learners = "bad"),
               "custom_learners must be a named list")

  expect_error(fit_resample_quiet(df, outcome = "outcome", splits = splits,
                                  learner = "glm", custom_learners = list(list())),
               "custom_learners must be a named list")

  bad_custom <- list(glm = list(fit = function(...) NULL))
  expect_error(fit_resample_quiet(df, outcome = "outcome", splits = splits,
                                  learner = "glm", custom_learners = bad_custom),
               "fit.*predict")
})

test_that("fit_resample validates class weights and positive class", {
  df <- make_class_df(12)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()

  expect_error(fit_resample_quiet(df, outcome = "outcome", splits = splits,
                                  learner = "glm", custom_learners = custom,
                                  class_weights = c(a = 1),
                                  metrics = "auc", refit = FALSE),
               "missing levels")

  expect_error(fit_resample_quiet(df, outcome = "outcome", splits = splits,
                                  learner = "glm", custom_learners = custom,
                                  positive_class = c(0, 1),
                                  metrics = "auc", refit = FALSE),
               "single value")
})

test_that("fit_resample supports class weights with parsnip learners", {
  skip_if_not_installed("parsnip")
  df <- make_class_df(12)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  spec <- parsnip::logistic_reg() |> parsnip::set_engine("glm")
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = spec,
                            class_weights = c("0" = 1, "1" = 1),
                            metrics = "auc", refit = FALSE)
  expect_true(nrow(fit@metrics) > 0)
})

test_that("fit_resample drops invalid metrics with warnings", {
  df <- make_class_df(10)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 2, seed = 1)
  custom <- make_custom_learners()

  fit <- expect_warning_match(
    fit_resample_quiet(df, outcome = "outcome", splits = splits,
                       learner = "glm", custom_learners = custom,
                       metrics = c("auc", "rmse"), refit = FALSE),
    "Dropping metrics"
  )
  expect_true(nrow(fit@metrics) > 0)
})

test_that("fit_resample summarizes metrics with all-NA columns", {
  df <- make_class_df(10)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 2, seed = 1)
  custom <- make_custom_learners()
  na_metric <- function(y, pred) NA_real_

  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = custom,
                            metrics = list(auc = "auc", na_metric = na_metric),
                            refit = FALSE)
  expect_true(nrow(fit@metric_summary) > 0)
  expect_true("na_metric_mean" %in% colnames(fit@metric_summary))
})

test_that("fit_resample handles gaussian tasks and ignores classification options", {
  df <- make_regression_df(12)
  splits <- make_split_plan_quiet(df, outcome = "y",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()

  fit <- expect_warning_match(
    fit_resample_quiet(df, outcome = "y", splits = splits,
                       learner = "glm", custom_learners = custom,
                       class_weights = c(a = 1, b = 2),
                       positive_class = "a",
                       metrics = c("rmse", "cindex"), refit = FALSE),
    "ignored",
    all = TRUE
  )
  expect_equal(fit@task, "gaussian")
})

# test_that("fit_resample errors on unsupported outcomes", {
#   df <- make_class_df(10)
#   df$outcome <- factor(c("a", "b", "c", "a", "b", "c", "a", "b", "c", "a"))
#   splits <- make_split_plan_quiet(df, outcome = "outcome",
#                               mode = "subject_grouped", group = "subject",
#                               v = 2, seed = 1)
#   custom <- make_custom_learners()
#   expect_error(fit_resample_quiet(df, outcome = "outcome", splits = splits,
#                                   learner = "glm", custom_learners = custom),
#                "No successful folds were completed. Check learner and preprocessing settings.")
# })

test_that("fit_resample respects positive_class releveling", {
  df <- make_class_df(12)
  df$outcome <- factor(ifelse(df$outcome == 1, "yes", "no"), levels = c("no", "yes"))
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  custom <- make_custom_learners()
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = custom,
                            positive_class = "yes", metrics = "auc",
                            refit = FALSE)
  expect_equal(fit@info$positive_class, "yes")
})

test_that("fit_resample applies classification_threshold for binomial tasks", {
  df <- make_class_df(12)
  df$x1 <- seq(0.05, 0.95, length.out = nrow(df))
  custom <- list(
    fixed_prob = list(
      fit = function(x, y, task, weights, ...) list(),
      predict = function(object, newdata, task, ...) as.numeric(newdata$x1)
    )
  )
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 2)

  fit_low <- fit_resample_quiet(
    df, outcome = "outcome", splits = splits,
    learner = "fixed_prob", custom_learners = custom,
    metrics = "accuracy", classification_threshold = 0.2,
    refit = FALSE
  )
  fit_high <- fit_resample_quiet(
    df, outcome = "outcome", splits = splits,
    learner = "fixed_prob", custom_learners = custom,
    metrics = "accuracy", classification_threshold = 0.8,
    refit = FALSE
  )

  pred_low <- do.call(rbind, fit_low@predictions)
  pred_high <- do.call(rbind, fit_high@predictions)
  key_low <- order(pred_low$fold, pred_low$id)
  key_high <- order(pred_high$fold, pred_high$id)
  expect_false(identical(as.character(pred_low$pred_class[key_low]),
                         as.character(pred_high$pred_class[key_high])))
  expect_false(isTRUE(all.equal(fit_low@metrics$accuracy, fit_high@metrics$accuracy)))
  expect_equal(fit_low@info$classification_threshold, 0.2)
  expect_equal(fit_high@info$classification_threshold, 0.8)

  expect_error(
    fit_resample_quiet(df, outcome = "outcome", splits = splits,
                       learner = "fixed_prob", custom_learners = custom,
                       classification_threshold = 1.2, refit = FALSE),
    "classification_threshold must be"
  )
})

test_that("fit_resample surfaces custom learner prediction length errors", {
  df <- make_class_df(10)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 2, seed = 1)
  custom <- list(
    bad = list(
      fit = function(x, y, task, weights, ...) suppressWarnings(stats::glm(
        y ~ ., data = data.frame(y = y, x), family = stats::binomial()
      )),
      predict = function(object, newdata, task, ...) rep(0.5, nrow(newdata) + 1)
    )
  )
  expect_warning_match(
    expect_error(fit_resample_quiet(df, outcome = "outcome", splits = splits,
                                    learner = "bad", custom_learners = custom,
                                    metrics = "auc", refit = FALSE),
                 "No successful folds were completed"),
    "Custom learner 'bad' returned"
  )
})

test_that("fit_resample warns when a fold lacks both classes", {
  df <- make_class_df(8)
  df$outcome <- c(rep(0, 4), rep(1, 4))
  indices <- list(
    list(train = 1:4, test = 5:6, fold = 1, repeat_id = 1),
    list(train = c(1:3, 5:8), test = 4, fold = 2, repeat_id = 1)
  )
  splits <- bioLeak:::LeakSplits(mode = "custom", indices = indices,
                                 info = list(outcome = "outcome", coldata = df))
  custom <- make_custom_learners()
  fit <- expect_warning_match(
    fit_resample_quiet(df, outcome = "outcome", splits = splits,
                       learner = "glm", custom_learners = custom,
                       metrics = "accuracy", refit = FALSE),
    "only one class"
  )
  expect_true(nrow(fit@metrics) > 0)
  expect_true(is.data.frame(fit@info$fold_status))
  expect_equal(nrow(fit@info$fold_status), length(splits@indices))
  expect_true(any(fit@info$fold_status$status == "success"))
  expect_true(any(fit@info$fold_status$status == "skipped"))
})

test_that("fit_resample errors for compact time_series without time metadata", {
  df <- make_class_df(10)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "time_series", time = "time",
                              v = 2, seed = 1, compact = TRUE)
  splits@info$time <- NULL
  custom <- make_custom_learners()
  expect_error(fit_resample_quiet(df, outcome = "outcome", splits = splits,
                                  learner = "glm", custom_learners = custom,
                                  metrics = "auc", refit = FALSE),
               "time_series compact splits")
})

test_that("fit_resample handles compact splits with repeats", {
  df <- make_class_df(30)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "subject_grouped",
                              group = "subject",
                              v = 3, repeats = 2,
                              compact = TRUE, seed = 1)
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = make_custom_learners(),
                            metrics = "auc", refit = FALSE, seed = 1)
  pred_df <- do.call(rbind, fit@predictions)
  fold_counts <- table(pred_df$fold)
  expect_equal(length(fold_counts), length(splits@indices))
  expect_true(all(fold_counts > 0))
})

test_that("fit_resample resolves compact time_series splits with purge and embargo", {
  df <- make_class_df(12)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "time_series", time = "time",
                              v = 3, compact = TRUE, seed = 1,
                              purge = 2, embargo = 5)
  fit <- fit_resample_quiet(df, outcome = "outcome", splits = splits,
                            learner = "glm", custom_learners = make_custom_learners(),
                            metrics = "accuracy", refit = FALSE, seed = 1)

  expect_equal(fit@audit$n_train, c(3, 7))
  expect_equal(fit@splits@info$purge, 2)
  expect_equal(fit@splits@info$embargo, 5)
})

test_that("fit_resample supports parsnip learners when available", {
  skip_if_not_installed("parsnip")
  df <- make_class_df(12)
  splits <- make_split_plan_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 3, seed = 1)
  spec <- parsnip::logistic_reg() |>
    parsnip::set_engine("glm")
  fit <- expect_warning_match(
    fit_resample_quiet(df, outcome = "outcome", splits = splits,
                       learner = spec, custom_learners = make_custom_learners(),
                       learner_args = list(alpha = 1),
                       metrics = "auc", refit = FALSE, seed = 1),
    "ignored"
  )
  expect_true(nrow(fit@metrics) > 0)
  expect_true(any(nzchar(fit@metrics$learner)))
})

test_that("'No successful folds' says why when every metric is NA", {
  set.seed(1)
  df <- data.frame(id = 1:40, g = rep(1:4, each = 10), x = rnorm(40))
  df$y <- factor(ifelse(df$g <= 2, "a", "b"))  # every test fold is single-class
  sp <- make_split_plan_quiet(df, outcome = "y", mode = "subject_grouped", group = "g", v = 4)
  expect_error(suppressWarnings(fit_resample(df, outcome = "y", splits = sp, learner = "glm",
                                             metrics = "auc", refit = FALSE)),
               "Skipped 4 folds: every requested metric was NA")
})

test_that("'No successful folds' flattens multi-line errors and counts the rest", {
  set.seed(1)
  df <- data.frame(id = 1:40, x = rnorm(40), y = factor(rep(c("a", "b"), 20)))
  sp <- make_split_plan_quiet(df, outcome = "y", mode = "subject_grouped", group = "id", v = 5)
  k <- 0L
  custom <- list(bad = list(
    fit = function(x, y, task, weights, ...) { k <<- k + 1L; stop(sprintf("err %d\n  detail", k)) },
    predict = function(object, newdata, task, ...) rep(0.5, nrow(newdata))))
  e <- tryCatch(suppressWarnings(fit_resample(df, outcome = "y", splits = sp, learner = "bad",
                                              custom_learners = custom, metrics = "auc",
                                              refit = FALSE)),
                error = function(e) conditionMessage(e))
  expect_match(e, "err 1 detail | err 2 detail | err 3 detail (and 2 more).", fixed = TRUE)
})

test_that("split and stratify columns never reach the learner through a recipe", {
  skip_if_not_installed("recipes")
  set.seed(7); n <- 200
  df <- data.frame(subject = paste0("s", rep(1:50, each = 4)))
  df$outcome <- factor(rep(rbinom(50, 1, 0.5), each = 4), labels = c("ctrl", "case"))
  df$subtype <- ifelse(df$outcome == "case", sample(c("A", "B"), n, TRUE), "none")
  for (j in 1:6) df[[paste0("g", j)]] <- rnorm(n)
  sp <- make_split_plan(df, outcome = "outcome", mode = "subject_grouped", group = "subject",
                        v = 5, stratify = TRUE, stratify_by = "subtype", progress = FALSE)
  seen <- list()
  cl <- list(rec = list(
    fit = function(x, y, task, weights, ...) {
      seen[[length(seen) + 1L]] <<- colnames(x)
      stats::glm(y ~ ., data = data.frame(y = y, x), family = stats::binomial())
    },
    predict = function(object, newdata, task, ...) {
      stats::predict(object, newdata = data.frame(newdata), type = "response")
    }))
  rec <- recipes::recipe(outcome ~ ., data = df) |>
    recipes::step_dummy(recipes::all_nominal_predictors())
  fit <- suppressWarnings(fit_resample(df, outcome = "outcome", splits = sp, learner = "rec",
                                       custom_learners = cl, metrics = "auc", preprocess = rec))
  cols <- unique(unlist(seen))
  expect_setequal(cols, paste0("g", 1:6))
  expect_lt(mean(fit@metrics$auc), 0.75)
})

test_that("recipes that derive predictors from split columns by type selectors stop loudly", {
  skip_if_not_installed("recipes")
  set.seed(11); n <- 120
  df <- data.frame(subject = paste0("s", rep(1:30, each = 4)))
  df$outcome <- factor(rep(rbinom(30, 1, 0.5), each = 4), labels = c("ctrl", "case"))
  df$strat <- paste(df$outcome, rep(sample(c("x", "y"), 30, TRUE), each = 4))
  for (j in 1:4) df[[paste0("g", j)]] <- rnorm(n)
  sp <- make_split_plan(df, outcome = "outcome", mode = "subject_grouped", group = "subject",
                        v = 4, stratify = TRUE, stratify_by = "strat", progress = FALSE)
  base <- recipes::recipe(outcome ~ ., data = df)
  legacy <- base |> recipes::step_dummy(recipes::all_nominal(), -recipes::all_outcomes())
  expect_error(suppressWarnings(fit_resample(df, outcome = "outcome", splits = sp, learner = "glm",
                                             custom_learners = make_custom_learners(),
                                             preprocess = legacy, metrics = "auc")),
               "builds predictors from split columns")
  ok <- base |> recipes::step_rm(subject, strat) |>
    recipes::step_dummy(recipes::all_nominal(), -recipes::all_outcomes())
  fit <- suppressWarnings(fit_resample(df, outcome = "outcome", splits = sp, learner = "glm",
                                       custom_learners = make_custom_learners(),
                                       preprocess = ok, metrics = "auc"))
  expect_equal(nrow(fit@metrics), length(sp@indices))
})

test_that("recipes and workflows accept positive_class set to the first level", {
  skip_if_not_installed("recipes")
  skip_if_not_installed("workflows")
  set.seed(1); n <- 120
  df <- data.frame(id = 1:n, x1 = rnorm(n), x2 = rnorm(n))
  df$y <- factor(ifelse(df$x1 + rnorm(n) > 0, "b", "a"))
  sp <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "id", v = 4,
                        progress = FALSE)
  rec <- recipes::recipe(y ~ x1 + x2, data = df)
  spec <- parsnip::logistic_reg() |> parsnip::set_engine("glm")
  wf <- workflows::workflow() |> workflows::add_recipe(rec) |> workflows::add_model(spec)
  f_rec <- suppressWarnings(fit_resample(df, "y", sp, learner = spec, preprocess = rec,
                                         metrics = "auc", positive_class = "a"))
  f_wf <- suppressWarnings(fit_resample(df, "y", sp, learner = wf, metrics = "auc",
                                        positive_class = "a"))
  f_g <- suppressWarnings(fit_resample(df, "y", sp, learner = spec, metrics = "auc",
                                       positive_class = "a"))
  expect_equal(f_rec@metric_summary$auc_mean, f_g@metric_summary$auc_mean, tolerance = 1e-6)
  expect_equal(f_wf@metric_summary$auc_mean, f_g@metric_summary$auc_mean, tolerance = 1e-6)
  expect_gt(f_rec@metric_summary$auc_mean, 0.7)
})

test_that("multiclass yardstick metric sets with class probabilities work", {
  skip_if_not_installed("yardstick")
  skip_if_not_installed("glmnet")
  set.seed(2); n <- 150
  df <- data.frame(id = 1:n, x1 = rnorm(n), x2 = rnorm(n))
  df$y <- factor(c("a", "b", "c")[1 + (df$x1 > -0.4) + (df$x1 > 0.4)])
  sp <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "id", v = 3,
                        progress = FALSE)
  ms <- yardstick::metric_set(yardstick::roc_auc, yardstick::accuracy)
  f <- suppressWarnings(fit_resample(df, "y", sp,
                                     learner = parsnip::multinom_reg(penalty = 0.01) |>
                                       parsnip::set_engine("glmnet"),
                                     metrics = ms))
  expect_true(all(c("roc_auc_mean", "accuracy_mean") %in% names(f@metric_summary)))
  expect_gt(f@metric_summary$roc_auc_mean, 0.8)
})

test_that("positive_class works with non-syntactic level names in workflows", {
  skip_if_not_installed("workflows")
  set.seed(1); n <- 120
  df <- data.frame(subject = rep(1:60, each = 2), x1 = rnorm(n))
  df$y <- factor(ifelse(runif(n) < plogis(3 * df$x1), "1", "2"), levels = c("1", "2"))
  sp <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject", v = 5,
                        seed = 1, progress = FALSE)
  wf <- workflows::workflow() |> workflows::add_formula(y ~ x1) |>
    workflows::add_model(parsnip::logistic_reg())
  auc <- vapply(c("1", "2"), function(pc) {
    f <- suppressWarnings(fit_resample(df, outcome = "y", splits = sp, learner = wf,
                                       metrics = "auc", positive_class = pc, refit = FALSE))
    f@metric_summary$auc_mean
  }, numeric(1))
  expect_gt(min(auc), 0.8)            # never inverted
  expect_equal(auc[[1]], auc[[2]], tolerance = 1e-8)
})
