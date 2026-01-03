# Nested tuning ---------------------------------------------------------------

#' Leakage-aware nested tuning with tidymodels
#'
#' Runs nested cross-validation for hyperparameter tuning using leakage-aware
#' splits. Inner resamples are constructed from each outer training fold to
#' avoid information leakage during tuning. Requires tidymodels tuning
#' packages and a workflow or recipe-based preprocessing. Survival tasks are
#' not yet supported.
#'
#' @param x SummarizedExperiment or matrix/data.frame.
#' @param outcome Outcome column name (if x is SE or data.frame).
#' @param splits LeakSplits object defining the outer resamples. If the splits
#'   do not already include inner folds, they are created from each outer
#'   training fold using the same split metadata. rsample splits must already
#'   include inner folds.
#' @param split_cols Optional named list/character vector or `"auto"` (default)
#'   overriding group/batch/study/time column names when `splits` is an rsample
#'   object and its attributes are missing. `"auto"` falls back to common
#'   metadata column names (e.g., `group`, `subject`, `batch`, `study`, `time`).
#'   Supported names are `group`, `batch`, `study`, and `time`.
#' @param learner A parsnip model_spec with tunable parameters, or a workflows
#'   workflow. When a model_spec is provided, a workflow is built using
#'   `preprocess` or a formula.
#' @param preprocess Optional `recipes::recipe`. Required when you need
#'   preprocessing for tuning. Ignored when `learner` is already a workflow.
#' @param grid Tuning grid passed to `tune::tune_grid()`. Can be a data.frame or
#'   an integer size.
#' @param metrics Character vector of metric names (`auc`, `pr_auc`, `accuracy`,
#'   `macro_f1`, `log_loss`, `rmse`) or a yardstick metric set/list. Metrics are
#'   computed with yardstick; unsupported metrics are dropped with a warning.
#' @param positive_class Optional value indicating the positive class for
#'   binomial outcomes. When set, the outcome levels are reordered so the
#'   positive class is second.
#' @param selection Selection rule for tuning, either `"best"` or `"one_std_err"`.
#' @param selection_metric Metric name used for selecting hyperparameters.
#'   Defaults to the first metric in `metrics`.
#' @param inner_v Optional number of folds for inner CV when inner splits are
#'   not precomputed. Defaults to the outer `v`.
#' @param inner_repeats Optional number of repeats for inner CV when inner
#'   splits are not precomputed. Defaults to 1.
#' @param inner_seed Optional seed for inner split generation when inner splits
#'   are not precomputed. Defaults to the outer split seed.
#' @param control Optional `tune::control_grid()` settings for tuning.
#' @param parallel Logical; passed to [fit_resample()] when evaluating outer
#'   folds (single-fold, no refit).
#' @param seed Integer seed for reproducibility.
#' @return A list of class `"LeakTune"` with components:
#'   \item{metrics}{Outer-fold metrics.}
#'   \item{metric_summary}{Mean/SD metrics across outer folds.}
#'   \item{best_params}{Best hyperparameters per outer fold.}
#'   \item{inner_results}{List of inner tuning results.}
#'   \item{outer_fits}{List of outer LeakFit objects.}
#'   \item{info}{Metadata about the tuning run.}
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   subject = rep(1:10, each = 2),
#'   outcome = factor(rep(c(0, 1), each = 10)),
#'   x1 = rnorm(20),
#'   x2 = rnorm(20)
#' )
#' splits <- make_splits(df, outcome = "outcome",
#'                       mode = "subject_grouped", group = "subject",
#'                       v = 3, nested = TRUE)
#' spec <- parsnip::logistic_reg(penalty = tune::tune(), mixture = 1) |>
#'   parsnip::set_engine("glmnet")
#' rec <- recipes::recipe(outcome ~ x1 + x2, data = df)
#' tuned <- tune_resample(df, outcome = "outcome", splits = splits,
#'                        learner = spec, preprocess = rec, grid = 5)
#' tuned$metric_summary
#' }
#' @export
tune_resample <- function(x, outcome, splits,
                          learner,
                          preprocess = NULL,
                          grid = 10,
                          metrics = NULL,
                          positive_class = NULL,
                          selection = c("best", "one_std_err"),
                          selection_metric = NULL,
                          inner_v = NULL,
                          inner_repeats = 1,
                          inner_seed = NULL,
                          control = NULL,
                          parallel = FALSE,
                          seed = 1,
                          split_cols = "auto") {

  selection <- match.arg(selection)

  if (!requireNamespace("tune", quietly = TRUE)) {
    stop("Package 'tune' is required for tune_resample().", call. = FALSE)
  }
  if (!requireNamespace("dials", quietly = TRUE)) {
    stop("Package 'dials' is required for tune_resample().", call. = FALSE)
  }
  if (!requireNamespace("yardstick", quietly = TRUE)) {
    stop("Package 'yardstick' is required for tune_resample().", call. = FALSE)
  }
  if (!requireNamespace("rsample", quietly = TRUE)) {
    stop("Package 'rsample' is required for tune_resample().", call. = FALSE)
  }
  if (!requireNamespace("workflows", quietly = TRUE)) {
    stop("Package 'workflows' is required for tune_resample().", call. = FALSE)
  }

  is_parsnip_spec <- function(obj) inherits(obj, "model_spec")
  is_workflow <- function(obj) inherits(obj, "workflow")

  if (!is_workflow(learner) && !is_parsnip_spec(learner)) {
    stop("learner must be a parsnip model_spec or a workflows::workflow.", call. = FALSE)
  }

  if (is_workflow(learner) && !is.null(preprocess)) {
    warning("preprocess ignored when learner is a workflow.")
  }
  if (!is.null(preprocess) && !inherits(preprocess, "recipe") && !is_workflow(learner)) {
    stop("tune_resample requires a recipe (or a workflow) for preprocessing.", call. = FALSE)
  }
  if (!is_workflow(learner) && inherits(preprocess, "recipe") &&
      !requireNamespace("recipes", quietly = TRUE)) {
    stop("Package 'recipes' is required when preprocess is a recipe.", call. = FALSE)
  }

  if (!inherits(splits, "LeakSplits")) {
    if (.bio_is_rsample(splits)) {
      coldata <- if (.bio_is_se(x)) {
        as.data.frame(SummarizedExperiment::colData(x))
      } else if (is.data.frame(x)) {
        x
      } else if (is.matrix(x)) {
        data.frame(row_id = seq_len(nrow(x)))
      } else {
        NULL
      }
      splits <- .bio_as_leaksplits_from_rsample(splits, n = nrow(.bio_get_x(x)), coldata = coldata,
                                                split_cols = split_cols)
    } else {
      stop("splits must be a LeakSplits or rsample rset/rsplit.", call. = FALSE)
    }
  }
  if (identical(splits@mode, "rsample") && is.null(splits@info$inner)) {
    stop("rsample splits require precomputed inner folds for tune_resample().", call. = FALSE)
  }

  Xall <- .bio_get_x(x)
  yall <- .bio_get_y(x, outcome)

  if (.bio_is_survival(yall)) {
    stop("tune_resample does not yet support survival tasks.", call. = FALSE)
  }

  drop_cols <- outcome
  split_info <- splits@info
  drop_cols <- unique(c(drop_cols,
                        split_info$group,
                        split_info$batch,
                        split_info$study,
                        split_info$time))
  drop_cols <- drop_cols[!is.na(drop_cols) & nzchar(drop_cols)]
  if (length(drop_cols) && !is.null(colnames(Xall))) {
    drop_cols <- intersect(colnames(Xall), drop_cols)
  }
  if (length(drop_cols)) {
    Xall <- Xall[, setdiff(colnames(Xall), drop_cols), drop = FALSE]
  }

  task <- if (.bio_is_binomial(yall)) "binomial"
  else if (.bio_is_multiclass(yall)) "multiclass"
  else if (.bio_is_regression(yall)) "gaussian"
  else if (is.factor(yall) && nlevels(yall) == 2) "binomial"
  else if (is.factor(yall) && nlevels(yall) > 2) "multiclass"
  else stop("Unsupported outcome type for tuning.", call. = FALSE)

  class_levels <- NULL
  if (task == "binomial") {
    if (!is.factor(yall)) yall <- factor(yall)
    yall <- droplevels(yall)
    if (nlevels(yall) != 2) {
      stop("Binomial task requires exactly two outcome levels.", call. = FALSE)
    }
    if (!is.null(positive_class)) {
      pos_chr <- as.character(positive_class)
      if (length(pos_chr) != 1L) {
        stop("positive_class must be a single value.", call. = FALSE)
      }
      levels_y <- levels(yall)
      if (!pos_chr %in% levels_y) {
        stop(sprintf("positive_class '%s' not found in outcome levels: %s",
                     pos_chr, paste(levels_y, collapse = ", ")), call. = FALSE)
      }
      if (!identical(pos_chr, levels_y[2])) {
        levels_y <- c(setdiff(levels_y, pos_chr), pos_chr)
        yall <- factor(yall, levels = levels_y)
      }
    }
    class_levels <- levels(yall)
  } else if (task == "multiclass") {
    if (!is.factor(yall)) yall <- factor(yall)
    yall <- droplevels(yall)
    if (nlevels(yall) < 3) {
      stop("Multiclass task requires 3 or more outcome levels.", call. = FALSE)
    }
    if (!is.null(positive_class)) {
      warning("positive_class is ignored for multiclass tasks.")
    }
    class_levels <- levels(yall)
  } else if (!is.numeric(yall)) {
    yall <- as.numeric(yall)
    if (anyNA(yall)) stop("Gaussian task requires numeric outcome values.", call. = FALSE)
  }

  if (is.null(control)) control <- tune::control_grid()

  resolve_metrics <- function(metrics, task) {
    macro_f1 <- yardstick::metric_tweak("macro_f1", yardstick::f_meas, estimator = "macro")
    auc_metric <- yardstick::metric_tweak("auc", yardstick::roc_auc, event_level = "second")
    roc_metric <- yardstick::metric_tweak("roc_auc", yardstick::roc_auc, event_level = "second")
    pr_metric <- yardstick::metric_tweak("pr_auc", yardstick::pr_auc, event_level = "second")
    log_loss_metric <- yardstick::metric_tweak("log_loss", yardstick::mn_log_loss)
    mn_log_loss_metric <- yardstick::metric_tweak("mn_log_loss", yardstick::mn_log_loss)

    defaults <- if (task == "binomial") c("auc", "pr_auc", "accuracy")
    else if (task == "multiclass") c("accuracy", "macro_f1")
    else c("rmse")

    allowed <- if (task == "binomial") c("auc", "roc_auc", "pr_auc", "accuracy")
    else if (task == "multiclass") c("accuracy", "macro_f1", "log_loss", "mn_log_loss")
    else c("rmse")

    if (is.null(metrics)) metrics <- defaults

    metric_set <- NULL
    if (inherits(metrics, "metric_set")) {
      metric_set <- metrics
    } else if (is.list(metrics) && length(metrics) &&
               all(vapply(metrics, .bio_is_yardstick_metric, logical(1)))) {
      metric_set <- do.call(yardstick::metric_set, metrics)
    } else if (is.character(metrics)) {
      invalid <- setdiff(metrics, allowed)
      if (length(invalid)) {
        warning(sprintf("Dropping metrics not applicable to %s task: %s", task,
                        paste(invalid, collapse = ", ")))
        metrics <- setdiff(metrics, invalid)
      }
      if (!length(metrics)) metrics <- defaults
    } else {
      stop("metrics must be a character vector or yardstick metric set.", call. = FALSE)
    }

    metric_map <- list(
      auc = auc_metric,
      roc_auc = roc_metric,
      pr_auc = pr_metric,
      accuracy = yardstick::accuracy,
      macro_f1 = macro_f1,
      log_loss = log_loss_metric,
      mn_log_loss = mn_log_loss_metric,
      rmse = yardstick::rmse
    )
    if (is.null(metric_set)) {
      metric_fns <- metric_map[metrics]
      metric_set <- do.call(yardstick::metric_set, metric_fns)
    }

    metric_list <- attr(metric_set, "metrics")
    metric_names <- names(metric_list)
    invalid <- setdiff(metric_names, allowed)
    if (length(invalid)) {
      warning(sprintf("Dropping metrics not applicable to %s task: %s", task,
                      paste(invalid, collapse = ", ")))
      metric_list <- metric_list[setdiff(metric_names, invalid)]
    }
    if (task == "binomial" && length(metric_list)) {
      if ("roc_auc" %in% names(metric_list)) metric_list[["roc_auc"]] <- roc_metric
      if ("auc" %in% names(metric_list)) metric_list[["auc"]] <- auc_metric
      if ("pr_auc" %in% names(metric_list)) metric_list[["pr_auc"]] <- pr_metric
    }
    if (!length(metric_list)) {
      metric_list <- metric_map[defaults]
    }
    metric_set <- do.call(yardstick::metric_set, metric_list)
    metric_names <- names(attr(metric_set, "metrics"))
    list(set = metric_set, names = metric_names)
  }

  metrics_resolved <- resolve_metrics(metrics, task)
  tune_metrics <- metrics_resolved$set
  metric_names <- metrics_resolved$names

  if (is.null(selection_metric)) {
    selection_metric <- metric_names[[1]]
  }
  if (!selection_metric %in% metric_names) {
    stop(sprintf("selection_metric '%s' not found in metrics.", selection_metric), call. = FALSE)
  }

  make_fold_df <- function(X, y) {
    df <- as.data.frame(X, check.names = FALSE)
    df[[outcome]] <- y
    df[, c(outcome, setdiff(names(df), outcome)), drop = FALSE]
  }

  build_workflow <- function() {
    if (is_workflow(learner)) return(learner)
    if (inherits(preprocess, "recipe")) {
      return(workflows::workflow() |> workflows::add_model(learner) |>
               workflows::add_recipe(preprocess))
    }
    if (length(outcome) != 1L) {
      stop("Formula workflows require a single outcome column.", call. = FALSE)
    }
    form <- stats::as.formula(paste0("`", outcome, "` ~ ."))
    workflows::workflow() |>
      workflows::add_model(learner) |>
      workflows::add_formula(form)
  }

  base_workflow <- build_workflow()

  make_inner_rset <- function(inner_idx, data) {
    if (!length(inner_idx)) stop("Inner splits are empty.", call. = FALSE)
    ids <- vapply(seq_along(inner_idx), function(i) {
      fold <- inner_idx[[i]]
      if (!is.null(fold$fold)) {
        rep_id <- fold$repeat_id %||% 1L
        paste0("Fold", fold$fold, "_Repeat", rep_id)
      } else {
        paste0("Fold", i)
      }
    }, character(1))
    split_objs <- lapply(seq_along(inner_idx), function(i) {
      fold <- inner_idx[[i]]
      if (is.null(fold$train) || is.null(fold$test)) {
        stop("Inner split is missing train/test indices.", call. = FALSE)
      }
      .bio_make_rsplit(fold$train, fold$test, data = data, id = ids[[i]])
    })
    manual_rset <- getFromNamespace("manual_rset", "rsample")
    do.call(manual_rset, list(splits = split_objs, ids = ids))
  }

  make_outer_splits <- function(train, test) {
    info <- splits@info
    info$v <- 1L
    info$repeats <- 1L
    info$nested <- FALSE
    info$inner <- NULL
    info$compact <- FALSE
    info$fold_assignments <- NULL
    info$summary <- data.frame(
      fold = 1L,
      repeat_id = 1L,
      train_n = length(train),
      test_n = length(test),
      stringsAsFactors = FALSE
    )
    indices <- list(list(train = train, test = test, fold = 1L, repeat_id = 1L))
    info$hash <- .bio_hash_indices(indices)
    new("LeakSplits", mode = splits@mode, indices = indices, info = info)
  }

  inner_seed <- inner_seed %||% (splits@info$seed %||% seed)
  inner_v <- inner_v %||% (splits@info$v %||% 5L)
  inner_repeats <- inner_repeats %||% 1L

  df_all <- make_fold_df(Xall, yall)
  coldata <- if (.bio_is_se(x)) {
    as.data.frame(SummarizedExperiment::colData(x))
  } else if (is.data.frame(x)) {
    x
  } else if (is.matrix(x)) {
    data.frame(row_id = seq_len(nrow(Xall)))
  } else {
    NULL
  }

  folds <- splits@indices
  metrics_rows <- list()
  best_rows <- list()
  inner_results <- list()
  outer_fits <- list()

  for (i in seq_along(folds)) {
    fold <- folds[[i]]
    fold <- .bio_resolve_fold_indices(splits, fold, n = nrow(df_all), data = coldata)
    tr <- fold$train
    te <- fold$test
    if (!length(tr) || !length(te)) {
      warning(sprintf("Outer fold %d skipped: empty train/test.", i), call. = FALSE)
      next
    }

    df_train <- df_all[tr, , drop = FALSE]

    inner_idx <- NULL
    if (!is.null(splits@info$inner) && length(splits@info$inner) >= i) {
      inner_idx <- splits@info$inner[[i]]
    }
    if (is.null(inner_idx)) {
      if (identical(splits@mode, "rsample")) {
        stop("Outer splits from rsample require precomputed inner splits.", call. = FALSE)
      }
      x_inner <- if (.bio_is_se(x)) x[, tr] else x[tr, , drop = FALSE]
      inner <- make_splits(
        x_inner,
        outcome = outcome,
        mode = splits@mode,
        group = splits@info$group,
        batch = splits@info$batch,
        study = splits@info$study,
        time = splits@info$time,
        v = inner_v,
        repeats = inner_repeats,
        stratify = isTRUE(splits@info$stratify),
        nested = FALSE,
        seed = inner_seed + i,
        horizon = splits@info$horizon %||% 0,
        progress = FALSE
      )
      inner_idx <- inner@indices
    }

    inner_rset <- make_inner_rset(inner_idx, df_train)
    set.seed(seed + i)
    tune_res <- tune::tune_grid(
      base_workflow,
      resamples = inner_rset,
      grid = grid,
      metrics = tune_metrics,
      control = control
    )
    best_params <- if (selection == "best") {
      tune::select_best(tune_res, metric = selection_metric)
    } else {
      tune::select_by_one_std_err(tune_res, metric = selection_metric)
    }
    if (!nrow(best_params)) {
      warning(sprintf("Outer fold %d skipped: no tuning results.", i), call. = FALSE)
      next
    }
    final_workflow <- tune::finalize_workflow(base_workflow, best_params)
    outer_splits <- make_outer_splits(tr, te)
    outer_fit <- fit_resample(
      x,
      outcome = outcome,
      splits = outer_splits,
      preprocess = list(),
      learner = final_workflow,
      metrics = tune_metrics,
      positive_class = if (task == "binomial") class_levels[2] else NULL,
      parallel = parallel,
      refit = FALSE,
      seed = seed + i
    )

    outer_fits[[i]] <- outer_fit
    inner_results[[i]] <- tune_res

    fold_metrics <- outer_fit@metrics
    fold_metrics$fold <- i
    metrics_rows[[length(metrics_rows) + 1L]] <- fold_metrics

    best_df <- as.data.frame(best_params)
    best_df$fold <- i
    best_rows[[length(best_rows) + 1L]] <- best_df
  }

  if (!length(metrics_rows)) {
    stop("No successful outer folds were completed.", call. = FALSE)
  }

  metrics_df <- do.call(rbind, metrics_rows)
  best_params_df <- if (length(best_rows)) do.call(rbind, best_rows) else data.frame()

  metric_cols <- setdiff(names(metrics_df), "fold")
  metric_summary <- aggregate(. ~ learner, data = metrics_df[, metric_cols, drop = FALSE],
                              FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                                  sd = sd(x, na.rm = TRUE)),
                              na.action = stats::na.pass)

  structure(
    list(
      splits = splits,
      metrics = metrics_df,
      metric_summary = metric_summary,
      best_params = best_params_df,
      inner_results = inner_results,
      outer_fits = outer_fits,
      info = list(
        task = task,
        metrics_used = metric_names,
        positive_class = if (task == "binomial") class_levels[2] else NULL,
        selection = selection,
        selection_metric = selection_metric,
        grid = grid,
        seed = seed
      )
    ),
    class = "LeakTune"
  )
}
