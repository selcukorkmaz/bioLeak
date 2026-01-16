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
#'   For binomial tasks, if any inner assessment fold contains a single class,
#'   probability metrics (`auc`, `roc_auc`, `pr_auc`) are dropped for tuning with
#'   a warning.
#' @param positive_class Optional value indicating the positive class for
#'   binomial outcomes. When set, the outcome levels are reordered so the
#'   positive class is second.
#' @param selection Selection rule for tuning, either `"best"` or `"one_std_err"`.
#' @param selection_metric Metric name used for selecting hyperparameters.
#'   Defaults to the first metric in `metrics`. If the chosen metric yields
#'   no valid results, the first available metric is used with a warning.
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
#'   \item{metric_summary}{Mean/SD metrics across outer folds with columns
#'     \code{learner}, and \code{<metric>_mean} and \code{<metric>_sd} for
#'     each metric.}
#'   \item{best_params}{Best hyperparameters per outer fold.}
#'   \item{inner_results}{List of inner tuning results.}
#'   \item{outer_fits}{List of outer LeakFit objects.}
#'   \item{info}{Metadata about the tuning run.}
#' @examples
#' \donttest{
#' if (requireNamespace("tune", quietly = TRUE) &&
#'     requireNamespace("recipes", quietly = TRUE) &&
#'     requireNamespace("glmnet", quietly = TRUE) &&
#'     requireNamespace("rsample", quietly = TRUE) &&
#'     requireNamespace("workflows", quietly = TRUE) &&
#'     requireNamespace("yardstick", quietly = TRUE) &&
#'     requireNamespace("dials", quietly = TRUE)) {
#'   df <- data.frame(
#'     subject = rep(1:10, each = 2),
#'     outcome = factor(rep(c(0, 1), each = 10)),
#'     x1 = rnorm(20),
#'     x2 = rnorm(20)
#'   )
#'   splits <- make_split_plan(df, outcome = "outcome",
#'                         mode = "subject_grouped", group = "subject",
#'                         v = 3, nested = TRUE, stratify = TRUE)
#'   spec <- parsnip::logistic_reg(penalty = tune::tune(), mixture = 1) |>
#'     parsnip::set_engine("glmnet")
#'   rec <- recipes::recipe(outcome ~ x1 + x2, data = df)
#'   tuned <- tune_resample(df, outcome = "outcome", splits = splits,
#'                          learner = spec, preprocess = rec, grid = 5)
#'   tuned$metric_summary
#' }
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

  # --- Learner Name Extraction ---
  derive_learner_label <- function(obj) {
    spec <- obj
    if (inherits(obj, "workflow")) {
      if (requireNamespace("workflows", quietly = TRUE)) {
        spec <- tryCatch(workflows::extract_spec_parsnip(obj), error = function(e) NULL)
      }
    }
    if (inherits(spec, "model_spec")) {
      cls <- class(spec)
      cls <- cls[!cls %in% c("model_spec", "object")]
      model_type <- if (length(cls) > 0) cls[1] else "model"
      engine <- spec$engine
      if (!is.null(engine)) {
        return(paste0(model_type, "/", engine))
      } else {
        return(model_type)
      }
    }
    return("tuned_model")
  }
  learner_label <- derive_learner_label(learner)
  # -------------------------------

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
  y_orig <- .bio_get_y(x, outcome)
  yall <- y_orig
  y_data <- y_orig

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
    if (!is.factor(yall)) {
      yall <- factor(yall)
      y_data <- factor(y_data)
    }
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
        if (!inherits(preprocess, "recipe")) {
          y_data <- factor(y_data, levels = levels_y)
        }
      }
    }
    class_levels <- levels(yall)
  } else if (task == "multiclass") {
    if (!is.factor(yall)) {
      yall <- factor(yall)
      y_data <- factor(y_data)
    }
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
    y_data <- as.numeric(y_data)
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

    if (is.null(metrics)) metrics <- defaults

    if (is.character(metrics)) {
      metrics[metrics == "auc"] <- "roc_auc"
      valid_keys <- intersect(metrics, names(metric_map))
      invalid <- setdiff(metrics, names(metric_map))

      if (length(invalid)) {
        warning(sprintf("Dropping unsupported metrics for %s task: %s", task,
                        paste(invalid, collapse = ", ")))
      }

      if (length(valid_keys) == 0) valid_keys <- defaults[defaults %in% names(metric_map)]

      fns <- metric_map[valid_keys]
      metric_set <- do.call(yardstick::metric_set, fns)
      return(list(set = metric_set, names = valid_keys))
    }

    if (inherits(metrics, "metric_set")) {
      return(list(set = metrics, names = names(attr(metrics, "metrics"))))
    }

    if (is.list(metrics)) {
      metric_set <- do.call(yardstick::metric_set, metrics)
      return(list(set = metric_set, names = names(attr(metric_set, "metrics"))))
    }

    stop("metrics must be a character vector or yardstick metric set.", call. = FALSE)
  }

  subset_metric_set <- function(metric_set, keep) {
    metric_list <- attr(metric_set, "metrics")
    metric_list <- metric_list[names(metric_list) %in% keep]
    do.call(yardstick::metric_set, metric_list)
  }

  metrics_resolved <- resolve_metrics(metrics, task)
  tune_metrics <- metrics_resolved$set
  metric_names <- metrics_resolved$names

  if (is.null(selection_metric)) {
    if ("roc_auc" %in% metric_names) selection_metric <- "roc_auc"
    else if ("auc" %in% metric_names) selection_metric <- "auc"
    else selection_metric <- metric_names[[1]]
  }

  if (!selection_metric %in% metric_names) {
    if (selection_metric == "auc" && "roc_auc" %in% metric_names) selection_metric <- "roc_auc"
    else if (selection_metric == "roc_auc" && "auc" %in% metric_names) selection_metric <- "auc"
    else stop(sprintf("selection_metric '%s' not found in metrics.", selection_metric), call. = FALSE)
  }

  make_fold_df <- function(X, y) {
    df <- as.data.frame(X, check.names = FALSE)
    df[[outcome]] <- y
    df[, c(outcome, setdiff(names(df), outcome)), drop = FALSE]
  }

  has_two_classes <- function(y) {
    y <- y[!is.na(y)]
    length(unique(y)) >= 2L
  }

  available_metrics <- function(metrics_df) {
    if (is.null(metrics_df) || !nrow(metrics_df)) return(character())
    if (!".metric" %in% names(metrics_df)) return(character())
    value_col <- NULL
    if (".estimate" %in% names(metrics_df)) {
      value_col <- ".estimate"
    } else if ("estimate" %in% names(metrics_df)) {
      value_col <- "estimate"
    } else if ("mean" %in% names(metrics_df)) {
      value_col <- "mean"
    }
    if (is.null(value_col)) return(character())
    metric_ok <- tapply(!is.na(metrics_df[[value_col]]), metrics_df$.metric, any)
    names(metric_ok)[metric_ok]
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
    manual_rset <- utils::getFromNamespace("manual_rset", "rsample")
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

  df_all <- if (is.data.frame(x)) {
    x
  } else {
    make_fold_df(Xall, y_data)
  }

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
  skipped_single_class_outer <- 0L
  skipped_single_class_inner <- 0L
  skipped_no_results <- 0L
  metrics_used <- character()
  selection_metric_requested <- selection_metric
  selection_metric_used <- selection_metric
  selection_metric_warned <- FALSE

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
    y_train_logic <- yall[tr]

    if (task %in% c("binomial", "multiclass") && !has_two_classes(y_train_logic)) {
      warning(sprintf("Outer fold %d skipped: training data has a single outcome class.", i),
              call. = FALSE)
      skipped_single_class_outer <- skipped_single_class_outer + 1L
      next
    }

    inner_idx <- NULL
    inner_precomputed <- FALSE
    if (!is.null(splits@info$inner) && length(splits@info$inner) >= i) {
      inner_idx <- splits@info$inner[[i]]
      inner_precomputed <- TRUE
    }
    if (is.null(inner_idx)) {
      if (identical(splits@mode, "rsample")) {
        inner_mode <- "subject_grouped"
        if (!is.null(splits@info$perm_mode)) inner_mode <- splits@info$perm_mode
        grp <- splits@info$group
        if(is.null(grp) && "subject" %in% names(df_train)) grp <- "subject"
        if(is.null(grp) && "id" %in% names(df_train)) grp <- "id"
        x_inner <- if (.bio_is_se(x)) x[, tr] else x[tr, , drop = FALSE]
      } else {
        x_inner <- if (.bio_is_se(x)) x[, tr] else x[tr, , drop = FALSE]
      }

      inner <- make_split_plan(
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

    inner_global <- FALSE
    if (inner_precomputed) {
      inner_positions <- unlist(lapply(inner_idx, function(fold) {
        c(fold$train, fold$test)
      }), use.names = FALSE)
      inner_positions <- inner_positions[!is.na(inner_positions)]
      if (!length(inner_positions)) {
        stop(sprintf("Outer fold %d: inner splits are empty.", i), call. = FALSE)
      }
      if (any(inner_positions < 1L)) {
        stop(sprintf("Outer fold %d: inner split indices must be positive.", i), call. = FALSE)
      }
      if (any(inner_positions > nrow(df_all))) {
        stop(sprintf("Outer fold %d: inner split indices exceed available rows.", i),
             call. = FALSE)
      }
      inner_in_train <- all(inner_positions %in% seq_len(nrow(df_train)))
      inner_in_outer <- all(inner_positions %in% tr)
      if (!inner_in_train && inner_in_outer) {
        inner_global <- TRUE
      } else if (inner_in_outer && !identical(tr, seq_len(nrow(df_train)))) {
        inner_global <- TRUE
      } else if (!inner_in_train && !inner_in_outer) {
        stop(sprintf(
          "Outer fold %d: inner split indices do not match outer training data.",
          i
        ), call. = FALSE)
      }
    }

    # Precomputed inner indices may be global; choose the matching data source.
    y_inner <- if (inner_global) yall else y_train_logic
    if (task %in% c("binomial", "multiclass")) {
      inner_has_class <- vapply(inner_idx, function(fold) {
        idx <- fold$train
        if (is.null(idx)) return(FALSE)
        has_two_classes(y_inner[idx])
      }, logical(1))
      if (!any(inner_has_class)) {
        warning(sprintf(
          "Outer fold %d skipped: inner training folds have a single outcome class. Consider `stratify = TRUE` or reducing `inner_v`.",
          i
        ), call. = FALSE)
        skipped_single_class_inner <- skipped_single_class_inner + 1L
        next
      }
      if (task == "binomial") {
        drop_count <- sum(!inner_has_class)
        if (drop_count > 0L) {
          inner_idx <- inner_idx[inner_has_class]
          warning(sprintf(
            "Outer fold %d: dropped %d inner fold(s) with single-class training data.",
            i, drop_count
          ), call. = FALSE)
        }
      }
    }

    tune_metrics_fold <- tune_metrics
    metric_names_fold <- metric_names
    if (task == "binomial") {
      inner_assess_two <- vapply(inner_idx, function(fold) {
        idx <- fold$test
        if (is.null(idx)) return(FALSE)
        has_two_classes(y_inner[idx])
      }, logical(1))
      if (!all(inner_assess_two)) {
        drop_metrics <- intersect(metric_names_fold, c("auc", "roc_auc", "pr_auc"))
        keep_metrics <- setdiff(metric_names_fold, drop_metrics)
        if (length(drop_metrics) && length(keep_metrics)) {
          tune_metrics_fold <- subset_metric_set(tune_metrics_fold, keep_metrics)
          metric_names_fold <- keep_metrics
          warning(sprintf(
            "Outer fold %d: inner assessment folds lack both classes; dropping metrics %s for tuning.",
            i, paste(drop_metrics, collapse = ", ")
          ), call. = FALSE)
        }
      }
    }

    metrics_used <- unique(c(metrics_used, metric_names_fold))

    inner_data <- if (inner_global) df_all else df_train
    inner_rset <- make_inner_rset(inner_idx, inner_data)
    set.seed(seed + i)
    tune_res <- tune::tune_grid(
      base_workflow,
      resamples = inner_rset,
      grid = grid,
      metrics = tune_metrics_fold,
      control = control
    )
    inner_results[[i]] <- tune_res
    no_results <- FALSE
    metrics_raw <- tryCatch(
      tune::collect_metrics(tune_res, summarize = FALSE),
      error = function(e) {
        if (grepl("No results are available", conditionMessage(e), fixed = TRUE)) {
          no_results <<- TRUE
          return(NULL)
        }
        stop(e)
      }
    )
    if (isTRUE(no_results) || is.null(metrics_raw) || !nrow(metrics_raw)) {
      warning(sprintf(
        "Outer fold %d skipped: inner tuning produced no results. Check `tune::collect_notes()` for failures.",
        i
      ), call. = FALSE)
      skipped_no_results <- skipped_no_results + 1L
      next
    }
    avail_metrics <- available_metrics(metrics_raw)
    metric_used <- selection_metric_used

    if (!metric_used %in% avail_metrics) {
      if ("<prb_mtrc>" %in% avail_metrics && metric_used %in% c("roc_auc", "auc", "pr_auc")) {
        metric_used <- "<prb_mtrc>"
      } else if ("<clss_mtr>" %in% avail_metrics && metric_used %in% c("accuracy", "kap")) {
        metric_used <- "<clss_mtr>"
      } else {
        if (length(avail_metrics)) {
          metric_used <- avail_metrics[1]
          if (!selection_metric_warned) {
            warning(sprintf("Selected metric not found; falling back to '%s'.", metric_used), call. = FALSE)
            selection_metric_warned <- TRUE
          }
        } else {
          warning(sprintf("Outer fold %d skipped: no usable metrics found.", i), call. = FALSE)
          skipped_no_results <- skipped_no_results + 1L
          next
        }
      }
    }

    res_sub <- metrics_raw[metrics_raw$.metric == metric_used, ]
    res_agg <- aggregate(.estimate ~ .config, data = res_sub, FUN = mean, na.rm = TRUE)

    minimize <- metric_used %in% c("rmse", "mae", "log_loss", "mn_log_loss")
    best_idx <- if(minimize) which.min(res_agg$.estimate) else which.max(res_agg$.estimate)

    if (length(best_idx) == 0) {
      warning(sprintf("Outer fold %d skipped: could not determine best parameters.", i), call. = FALSE)
      next
    }

    best_config <- res_agg$.config[best_idx]
    param_cols <- setdiff(names(metrics_raw), c(".metric", ".estimator", ".estimate", "n", "std_err", ".config"))
    best_params_row <- metrics_raw[metrics_raw$.config == best_config, param_cols, drop = FALSE]
    best_params <- best_params_row[1, , drop = FALSE]

    final_workflow <- tune::finalize_workflow(base_workflow, best_params)
    outer_splits <- make_outer_splits(tr, te)

    # CRITICAL: Force the correct name in fit_resample by passing a named list
    learner_for_fit <- list()
    learner_for_fit[[learner_label]] <- final_workflow

    outer_fit <- fit_resample(
      x,
      outcome = outcome,
      splits = outer_splits,
      preprocess = list(),
      learner = learner_for_fit, # Passed as named list
      metrics = tune_metrics_fold,
      positive_class = if (task == "binomial") class_levels[2] else NULL,
      parallel = parallel,
      refit = FALSE,
      seed = seed + i
    )

    outer_fits[[i]] <- outer_fit
    fold_metrics <- outer_fit@metrics
    fold_metrics$fold <- i
    metrics_rows[[length(metrics_rows) + 1L]] <- fold_metrics

    best_df <- as.data.frame(best_params)
    best_df$fold <- i
    best_rows[[length(best_rows) + 1L]] <- best_df
  }

  if (!length(metrics_rows)) {
    msg <- "No successful outer folds were completed."
    reasons <- character()
    if (skipped_single_class_outer) {
      reasons <- c(reasons, "Some outer training folds had a single outcome class.")
    }
    if (skipped_single_class_inner) {
      reasons <- c(reasons, "Some inner training folds had a single outcome class.")
    }
    if (skipped_no_results) {
      reasons <- c(reasons, "Inner tuning produced no results for some folds.")
    }
    if (length(reasons)) {
      msg <- paste(msg, paste(reasons, collapse = " "),
                   "Consider `stratify = TRUE` in make_split_plan() or reducing `inner_v`.")
    }
    stop(msg, call. = FALSE)
  }

  metrics_df <- do.call(rbind, metrics_rows)
  best_params_df <- if (length(best_rows)) do.call(rbind, best_rows) else data.frame()

  metric_cols <- setdiff(names(metrics_df), "fold")
  metric_summary_raw <- aggregate(. ~ learner, data = metrics_df[, metric_cols, drop = FALSE],
                                  FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                                      sd = sd(x, na.rm = TRUE)),
                                  na.action = stats::na.pass)
  # flatten embedded matrices from aggregate() into separate columns
  metric_summary <- data.frame(learner = metric_summary_raw$learner,
                               stringsAsFactors = FALSE)
  for (col in setdiff(colnames(metric_summary_raw), "learner")) {
    mat <- metric_summary_raw[[col]]
    if (is.matrix(mat)) {
      metric_summary[[paste0(col, "_mean")]] <- mat[, "mean"]
      metric_summary[[paste0(col, "_sd")]] <- mat[, "sd"]
    } else {
      metric_summary[[col]] <- mat
    }
  }

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
        metrics_requested = metric_names,
        metrics_used = if (length(metrics_used)) metrics_used else metric_names,
        positive_class = if (task == "binomial") class_levels[2] else NULL,
        selection = selection,
        selection_metric = selection_metric_used,
        selection_metric_requested = selection_metric_requested,
        grid = grid,
        seed = seed
      )
    ),
    class = "LeakTune"
  )
}
