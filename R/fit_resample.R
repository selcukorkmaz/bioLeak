#' Fit and evaluate with leakage guards over predefined splits
#'
#' Performs cross-validated model training and evaluation using
#' leakage-protected preprocessing (.guard_fit) and user-specified learners.
#'
#' @param x SummarizedExperiment or matrix/data.frame
#' @param outcome outcome column name (if x is SE or data.frame), or a length-2
#'   character vector of time/event column names for survival outcomes.
#' @param splits LeakSplits object from make_split_plan(), or an `rsample` rset/rsplit.
#' @param split_cols Optional named list/character vector or `"auto"` (default)
#'   overriding group/batch/study/time column names when `splits` is an rsample
#'   object and its attributes are missing. `"auto"` falls back to common
#'   metadata column names (e.g., `group`, `subject`, `batch`, `study`, `time`).
#'   Supported names are `group`, `batch`, `study`, and `time`.
#' @param store_refit_data Logical; when TRUE (default), stores the original
#'   data and learner configuration inside the fit to enable refit-based
#'   permutation tests without manual `perm_refit_spec` setup.
#' @param preprocess list(impute, normalize, filter=list(...), fs) or a
#'   `recipes::recipe` object. When a recipe is supplied, the guarded preprocessing
#'   pipeline is bypassed and the recipe is prepped on training data only.
#' @param learner parsnip model_spec (or list of model_spec objects) describing
#'   the model(s) to fit, or a `workflows::workflow`. For legacy use, a character
#'   vector of learner names (e.g., "glmnet", "ranger") or custom learner IDs is
#'   still supported.
#' @param learner_args list of additional arguments passed to legacy learners
#'   (ignored when `learner` is a parsnip model_spec).
#' @param custom_learners named list of custom learner definitions used only
#'   with legacy character learners. Each entry
#'   must contain \code{fit} and \code{predict} functions. The \code{fit} function
#'   should accept \code{x}, \code{y}, \code{task}, and \code{weights}, and return
#'   a model object. The \code{predict} function should accept \code{object},
#'   \code{newdata}, and \code{task}. For binomial/regression/survival tasks it
#'   should return a numeric vector; for multiclass tasks it should return either
#'   class labels or a matrix/data.frame of class probabilities.
#' @param metrics named list of metric functions, vector of metric names, or a
#'   `yardstick::metric_set`. When a yardstick metric set (or list of yardstick
#'   metric functions) is supplied, metrics are computed using yardstick with the
#'   positive class set to the second factor level.
#' @param class_weights optional named numeric vector of weights for binomial or
#'   multiclass outcomes
#' @param positive_class optional value indicating the positive class for binomial outcomes.
#'   When set, the outcome levels are reordered so that \code{positive_class} is treated
#'   as the positive class (level 2). If NULL, the second factor level is used.
#' @param parallel logical, use future.apply for multicore execution
#' @param refit logical, if TRUE retrain final model on full data
#' @param seed integer, for reproducibility
#' @return A \code{\linkS4class{LeakFit}} S4 object containing:
#'   \describe{
#'     \item{\code{splits}}{The \code{LeakSplits} object used for resampling.}
#'     \item{\code{metrics}}{Data.frame of per-fold, per-learner performance
#'       metrics with columns \code{fold}, \code{learner}, and one column per
#'       requested metric.}
#'     \item{\code{metric_summary}}{Data.frame summarizing metrics across folds
#'       for each learner with columns \code{learner}, and \code{<metric>_mean}
#'       and \code{<metric>_sd} for each requested metric.}
#'     \item{\code{audit}}{Data.frame with per-fold audit information including
#'       \code{fold}, \code{n_train}, \code{n_test}, \code{learner}, and
#'       \code{features_final} (number of features after preprocessing).}
#'     \item{\code{predictions}}{List of data.frames containing out-of-fold
#'       predictions with columns \code{id} (sample identifier), \code{truth}
#'       (true outcome), \code{pred} (predicted value or probability), \code{fold},
#'       and \code{learner}. For classification tasks, includes \code{pred_class}.
#'       For multiclass, includes per-class probability columns.}
#'     \item{\code{preprocess}}{List of preprocessing state objects from each fold,
#'       storing imputation parameters, normalization statistics, and feature
#'       selection results.}
#'     \item{\code{learners}}{List of fitted model objects from each fold.}
#'     \item{\code{outcome}}{Character string naming the outcome variable.}
#'     \item{\code{task}}{Character string indicating the task type
#'       (\code{"binomial"}, \code{"multiclass"}, \code{"gaussian"}, or
#'       \code{"survival"}).}
#'     \item{\code{feature_names}}{Character vector of feature names after
#'       preprocessing.}
#'     \item{\code{info}}{List of additional metadata including \code{hash},
#'       \code{metrics_used}, \code{class_weights}, \code{positive_class},
#'       \code{sample_ids}, \code{fold_status}, \code{refit}, \code{final_model} (refitted model if
#'       \code{refit = TRUE}), \code{final_preprocess}, \code{learner_names},
#'       and \code{perm_refit_spec} (for permutation-based audits).}
#'   }
#'   Use \code{summary()} to print a formatted report, or access slots directly
#'   with \code{@}.
#' @details
#' Preprocessing is fit on the training fold and applied to the test fold,
#' preventing leakage from global imputation, scaling, or feature selection.
#' When a `recipes::recipe` or `workflows::workflow` is supplied, the recipe is
#' prepped on the training fold and baked on the test fold.
#' For data.frame or matrix inputs, columns used to define splits
#' (outcome, group, batch, study, time) are excluded from the predictor matrix.
#' Use \code{learner_args} to pass model-specific arguments, either as a named
#' list keyed by learner or a single list applied to all learners. For custom
#' learners, \code{learner_args[[name]]} may be a list with \code{fit} and
#' \code{predict} sublists to pass distinct arguments to each stage. For binomial
#' tasks, predictions and metrics assume the positive class is the second factor
#' level; use \code{positive_class} to control this. Parsnip learners must support
#' probability predictions for binomial metrics (AUC/PR-AUC/accuracy) and
#' multiclass log-loss when requested.
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   subject = rep(1:10, each = 2),
#'   outcome = rbinom(20, 1, 0.5),
#'   x1 = rnorm(20),
#'   x2 = rnorm(20)
#' )
#' splits <- make_split_plan(df, outcome = "outcome",
#'                       mode = "subject_grouped", group = "subject", v = 5)
#'
#' # glmnet learner (requires glmnet package)
#' fit <- fit_resample(df, outcome = "outcome", splits = splits,
#'                       learner = "glmnet", metrics = "auc")
#' summary(fit)
#'
#' # Custom learner (logistic regression) - no extra packages needed
#' custom <- list(
#'   glm = list(
#'     fit = function(x, y, task, weights, ...) {
#'       stats::glm(y ~ ., data = as.data.frame(x),
#'                  family = stats::binomial(), weights = weights)
#'     },
#'     predict = function(object, newdata, task, ...) {
#'       as.numeric(stats::predict(object, newdata = as.data.frame(newdata), type = "response"))
#'     }
#'   )
#' )
#' fit2 <- fit_resample(df, outcome = "outcome", splits = splits,
#'                      learner = "glm", custom_learners = custom,
#'                      metrics = "accuracy")
#'
#' summary(fit2)
#' @export
fit_resample <- function(x, outcome, splits,
                         preprocess = list(
                           impute = list(method = "median"),
                           normalize = list(method = "zscore"),
                           filter = list(var_thresh = 0, iqr_thresh = 0),
                           fs = list(method = "none")
                         ),
                         learner = c("glmnet", "ranger"),
                         learner_args = list(),
                         custom_learners = list(),
                         metrics = c("auc", "pr_auc", "accuracy"),
                         class_weights = NULL,
                         positive_class = NULL,
                         parallel = FALSE,
                         refit = TRUE,
                         seed = 1,
                         split_cols = "auto",
                         store_refit_data = TRUE) {

  set.seed(seed)
  learner_input <- learner
  if (is.null(custom_learners)) custom_learners <- list()
  if (!is.list(custom_learners)) {
    stop("custom_learners must be a named list of learner definitions.")
  }
  if (length(custom_learners) && is.null(names(custom_learners))) {
    stop("custom_learners must be a named list.")
  }
  if (length(custom_learners)) {
    dup_names <- intersect(names(custom_learners), c("glmnet", "ranger"))
    if (length(dup_names)) {
      stop(sprintf("custom_learners cannot override built-in learners: %s",
                   paste(dup_names, collapse = ", ")))
    }
    bad <- vapply(custom_learners, function(def) {
      !is.list(def) || !is.function(def$fit) || !is.function(def$predict)
    }, logical(1))
    if (any(bad)) {
      stop("Each custom learner must be a list with `fit` and `predict` functions.")
    }
  }

  is_parsnip_spec <- function(obj) inherits(obj, "model_spec")
  is_workflow <- function(obj) inherits(obj, "workflow")
  use_parsnip <- FALSE
  use_workflow <- FALSE
  learner_specs <- NULL
  learner_names <- NULL

  if (is_workflow(learner)) {
    use_workflow <- TRUE
    learner_specs <- list(learner)
  } else if (is.list(learner) && length(learner) &&
             all(vapply(learner, is_workflow, logical(1)))) {
    use_workflow <- TRUE
    learner_specs <- learner
  } else if (is_parsnip_spec(learner)) {
    use_parsnip <- TRUE
    learner_specs <- list(learner)
  } else if (is.list(learner) && length(learner) &&
             all(vapply(learner, is_parsnip_spec, logical(1)))) {
    use_parsnip <- TRUE
    learner_specs <- learner
  }

  if (use_workflow) {
    if (!requireNamespace("workflows", quietly = TRUE)) {
      stop("Package 'workflows' is required when passing a workflow to learner.")
    }
    if (length(custom_learners)) {
      warning("custom_learners ignored when learner is a workflow.")
    }
    if (length(learner_args)) {
      warning("learner_args ignored when learner is a workflow.")
    }
    learner_names <- names(learner_specs)
    if (is.null(learner_names)) learner_names <- rep("", length(learner_specs))
    missing_names <- !nzchar(learner_names)
    if (any(missing_names)) {
      fallback <- paste0("workflow_", seq_along(learner_specs))
      learner_names[missing_names] <- fallback[missing_names]
    }
  } else if (use_parsnip) {
    if (!requireNamespace("parsnip", quietly = TRUE)) {
      stop("Package 'parsnip' is required when passing a model_spec to learner.")
    }
    if (length(custom_learners)) {
      warning("custom_learners ignored when learner is a parsnip model_spec.")
    }
    if (length(learner_args)) {
      warning("learner_args ignored when learner is a parsnip model_spec.")
    }
    parsnip_label <- function(spec, fallback) {
      model_class <- class(spec)
      model_class <- model_class[model_class != "model_spec"]
      model_class <- model_class[1]
      label <- if (!is.null(model_class) && nzchar(model_class)) model_class else fallback
      engine <- NULL
      if (!is.null(spec$engine) && nzchar(spec$engine)) {
        engine <- spec$engine
      } else if (!is.null(spec$method) && !is.null(spec$method$engine) &&
                 nzchar(spec$method$engine)) {
        engine <- spec$method$engine
      }
      if (!is.null(engine)) label <- paste0(label, "/", engine)
      label
    }
    learner_names <- names(learner_specs)
    if (is.null(learner_names)) learner_names <- rep("", length(learner_specs))
    missing_names <- !nzchar(learner_names)
    if (any(missing_names)) {
      fallback <- paste0("spec_", seq_along(learner_specs))
      spec_labels <- vapply(seq_along(learner_specs), function(i) {
        parsnip_label(learner_specs[[i]], fallback[[i]])
      }, character(1))
      learner_names[missing_names] <- spec_labels[missing_names]
    }
  } else {
    builtin_learners <- c("glmnet", "ranger")
    all_learners <- c(builtin_learners, names(custom_learners))
    if (!is.character(learner)) {
      stop("learner must be a character vector of legacy learners, a parsnip model_spec, or a workflow.")
    }
    learner <- match.arg(learner, choices = all_learners, several.ok = TRUE)
    learner_names <- learner
  }
  use_recipe <- .bio_is_recipe(preprocess)
  if (use_recipe && !requireNamespace("recipes", quietly = TRUE)) {
    stop("Package 'recipes' is required when preprocess is a recipe.")
  }
  if (use_workflow && isTRUE(use_recipe)) {
    warning("Recipe preprocess ignored when learner is a workflow.")
    use_recipe <- FALSE
  }
  preprocess_mode <- if (use_workflow) "workflow" else if (use_recipe) "recipe" else "guard"

  Xall <- .bio_get_x(x)
  yall <- .bio_get_y(x, outcome)

  if (!inherits(splits, "LeakSplits")) {
    if (.bio_is_rsample(splits)) {
      coldata <- if (.bio_is_se(x)) {
        as.data.frame(SummarizedExperiment::colData(x))
      } else if (is.data.frame(x)) {
        x
      } else if (is.matrix(x)) {
        data.frame(row_id = seq_len(nrow(Xall)))
      } else {
        NULL
      }
      splits <- .bio_as_leaksplits_from_rsample(splits, n = nrow(Xall), coldata = coldata,
                                                split_cols = split_cols)
    } else {
      stop("splits must be a LeakSplits or rsample rset/rsplit.")
    }
  }
  drop_cols <- outcome
  if (inherits(splits, "LeakSplits")) {
    split_info <- splits@info
    drop_cols <- unique(c(drop_cols,
                          split_info$group,
                          split_info$batch,
                          split_info$study,
                          split_info$time))
  }
  drop_cols <- drop_cols[!is.na(drop_cols) & nzchar(drop_cols)]
  if (length(drop_cols) && !is.null(colnames(Xall))) {
    drop_cols <- intersect(colnames(Xall), drop_cols)
  }
  if (length(drop_cols)) {
    Xall <- Xall[, setdiff(colnames(Xall), drop_cols), drop = FALSE]
  }
  sample_ids <- NULL
  if (inherits(splits, "LeakSplits") && !is.null(splits@info$coldata)) {
    cd <- splits@info$coldata
    rn_cd <- rownames(cd)
    if (!is.null(rn_cd) && !anyNA(rn_cd) && all(nzchar(rn_cd)) && !anyDuplicated(rn_cd)) {
      sample_ids <- rn_cd
    } else if ("row_id" %in% names(cd)) {
      rid <- as.character(cd[["row_id"]])
      if (length(rid) == nrow(cd) && !anyNA(rid) && !anyDuplicated(rid) && all(nzchar(rid))) {
        sample_ids <- rid
      }
    }
  }
  if (is.null(sample_ids)) {
    rn <- rownames(Xall)
    if (!is.null(rn) && !anyNA(rn) && all(nzchar(rn)) && !anyDuplicated(rn)) {
      sample_ids <- rn
    }
  }
  if (is.null(sample_ids) || length(sample_ids) != nrow(Xall)) {
    sample_ids <- as.character(seq_len(nrow(Xall)))
  }
  ids <- sample_ids
  compact <- isTRUE(splits@info$compact)
  fold_assignments <- splits@info$fold_assignments
  split_mode <- splits@mode
  split_time <- splits@info$time
  split_horizon <- splits@info$horizon %||% 0
  split_coldata <- splits@info$coldata
  time_vec <- NULL
  if (compact && identical(split_mode, "time_series")) {
    if (is.null(split_coldata) || is.null(split_time) || !split_time %in% names(split_coldata)) {
      stop("time_series compact splits require time column in coldata.")
    }
    time_vec <- split_coldata[[split_time]]
  }
  task <- if (.bio_is_survival(yall)) "survival"
  else if (.bio_is_binomial(yall)) "binomial"
  else if (.bio_is_multiclass(yall)) "multiclass"
  else if (.bio_is_regression(yall)) "gaussian"
  else if (is.factor(yall) && nlevels(yall) == 2) "binomial"
  else if (is.factor(yall) && nlevels(yall) > 2) "multiclass"
  else stop("Unsupported outcome type: require binomial/multiclass factor, numeric regression, or survival outcome.")

  if (task == "binomial") {
    if (!is.factor(yall)) yall <- factor(yall)
    yall <- droplevels(yall)
    if (nlevels(yall) != 2) {
      stop("Binomial task requires exactly two outcome levels after preprocessing.")
    }
    if (!is.null(positive_class)) {
      pos_chr <- as.character(positive_class)
      if (length(pos_chr) != 1L) {
        stop("positive_class must be a single value.")
      }
      levels_y <- levels(yall)
      if (!pos_chr %in% levels_y) {
        stop(sprintf("positive_class '%s' not found in outcome levels: %s",
                     pos_chr, paste(levels_y, collapse = ", ")))
      }
      if (!identical(pos_chr, levels_y[2])) {
        levels_y <- c(setdiff(levels_y, pos_chr), pos_chr)
        yall <- factor(yall, levels = levels_y)
      }
    }
    class_levels <- levels(yall)
    if (!is.null(class_weights)) {
      if (!is.numeric(class_weights)) stop("class_weights must be numeric.")
      if (is.null(names(class_weights))) {
        if (length(class_weights) != length(class_levels)) {
          stop("class_weights must align with outcome levels.")
        }
        names(class_weights) <- class_levels[seq_along(class_weights)]
      }
      missing_cw <- setdiff(class_levels, names(class_weights))
      if (length(missing_cw)) {
        stop(sprintf("class_weights missing levels: %s", paste(missing_cw, collapse = ", ")))
      }
      class_weights <- class_weights[class_levels]
    }
  } else if (task == "multiclass") {
    if (!is.factor(yall)) yall <- factor(yall)
    yall <- droplevels(yall)
    if (nlevels(yall) < 3) {
      stop("Multiclass task requires 3 or more outcome levels after preprocessing.")
    }
    class_levels <- levels(yall)
    if (!is.null(class_weights)) {
      if (!is.numeric(class_weights)) stop("class_weights must be numeric.")
      if (is.null(names(class_weights))) {
        if (length(class_weights) != length(class_levels)) {
          stop("class_weights must align with outcome levels.")
        }
        names(class_weights) <- class_levels[seq_along(class_weights)]
      }
      missing_cw <- setdiff(class_levels, names(class_weights))
      if (length(missing_cw)) {
        stop(sprintf("class_weights missing levels: %s", paste(missing_cw, collapse = ", ")))
      }
      class_weights <- class_weights[class_levels]
    }
    if (!is.null(positive_class)) {
      warning("positive_class is ignored for multiclass tasks.")
    }
  } else if (task == "gaussian") {
    if (!is.numeric(yall)) {
      yall <- as.numeric(yall)
      if (anyNA(yall)) stop("Gaussian task requires numeric outcome values.")
    }
    class_levels <- NULL
    if (!is.null(class_weights)) {
      warning("class_weights is ignored for gaussian tasks.")
    }
    if (!is.null(positive_class)) {
      warning("positive_class is ignored for gaussian tasks.")
    }
  } else {
    class_levels <- NULL
    if (!inherits(yall, "Surv")) {
      stop("Survival task requires a Surv outcome.")
    }
    if (!is.null(class_weights)) {
      warning("class_weights is ignored for survival tasks.")
    }
    if (!is.null(positive_class)) {
      warning("positive_class is ignored for survival tasks.")
    }
  }

  if (use_workflow && task == "binomial" && !is.null(class_weights)) {
    warning("class_weights are ignored for workflow learners unless explicitly handled in the workflow.")
  }
  if (use_workflow && task == "multiclass" && !is.null(class_weights)) {
    warning("class_weights are ignored for workflow learners unless explicitly handled in the workflow.")
  }

  metrics_input <- metrics
  metric_mode <- "legacy"
  yardstick_set <- NULL
  yardstick_metrics <- NULL

  if (!is.null(metrics) && inherits(metrics, "metric_set")) {
    metric_mode <- "yardstick"
    yardstick_set <- metrics
  } else if (!is.null(metrics) && .bio_is_yardstick_metric(metrics)) {
    metric_mode <- "yardstick"
    yardstick_metrics <- list(metrics)
  } else if (is.list(metrics) && length(metrics) &&
             all(vapply(metrics, .bio_is_yardstick_metric, logical(1)))) {
    metric_mode <- "yardstick"
    yardstick_metrics <- metrics
  }

  if (identical(metric_mode, "yardstick")) {
    if (!requireNamespace("yardstick", quietly = TRUE)) {
      stop("Package 'yardstick' is required when metrics are a yardstick set.", call. = FALSE)
    }
    if (is.null(yardstick_set)) {
      yardstick_set <- do.call(yardstick::metric_set, yardstick_metrics)
    }
    metric_labels <- character(0)
  } else {
    if (is.null(metrics)) {
      metrics <- if (task == "binomial") c("auc", "pr_auc", "accuracy")
      else if (task == "multiclass") c("accuracy", "macro_f1")
      else if (task == "survival") c("cindex")
      else c("rmse")
    }

    if (is.character(metrics)) {
      allowed <- if (task == "binomial") c("auc", "pr_auc", "accuracy")
      else if (task == "multiclass") c("accuracy", "macro_f1", "log_loss")
      else if (task == "survival") c("cindex")
      else c("rmse", "cindex")
      invalid <- setdiff(metrics, allowed)
      if (length(invalid)) {
        warning(sprintf("Dropping metrics not applicable to %s task: %s", task,
                        paste(invalid, collapse = ", ")))
        metrics <- setdiff(metrics, invalid)
      }
      if (!length(metrics)) {
        metrics <- if (task == "binomial") c("auc", "pr_auc", "accuracy")
        else if (task == "multiclass") c("accuracy", "macro_f1")
        else if (task == "survival") c("cindex")
        else c("rmse")
      }
    }

    metrics <- if (is.list(metrics)) metrics else as.list(metrics)
    metric_labels <- vapply(seq_along(metrics), function(i) {
      nm <- names(metrics)[i]
      if (!is.null(nm) && !is.na(nm) && nzchar(nm)) return(nm)
      mi <- metrics[[i]]
      if (is.character(mi) && length(mi) == 1) return(mi)
      paste0("metric_", i)
    }, character(1))
  }

  learner_objs <- if (use_parsnip || use_workflow) learner_specs else as.list(learner_names)

  # helper: safe metric computation -------------------------------------------
  compute_metric <- function(name, y, pred) {
    if (is.function(name)) return(name(y, pred))
    yb <- NULL
    if (task == "binomial") {
      yb <- if (is.factor(y)) as.numeric(y) - 1 else as.numeric(y)
    }
    if (name == "auc" && task == "binomial") {
      if (requireNamespace("pROC", quietly = TRUE))
        return(as.numeric(pROC::auc(pROC::roc(y, pred, quiet = TRUE))))
      pos <- pred[yb == 1]
      neg <- pred[yb == 0]
      comp <- outer(pos, neg, function(a, b) (a > b) + 0.5 * (a == b))
      return(mean(comp))
    }
    if (name == "pr_auc" && task == "binomial") {
      if (requireNamespace("PRROC", quietly = TRUE)) {
        pr <- PRROC::pr.curve(scores.class0 = pred[yb == 1],
                              scores.class1 = pred[yb == 0],
                              curve = FALSE)
        return(pr$auc.integral)
      }
      return(NA_real_)
    }
    if (name == "accuracy" && task == "binomial") {
      return(mean((pred >= 0.5) == as.logical(yb)))
    }
    if (name == "rmse" && task == "gaussian")
      return(sqrt(mean((as.numeric(y) - pred)^2)))
    if (name == "cindex" && task == "gaussian")
      return(.cindex_pairwise(pred, y))
    if (name == "cindex" && task == "survival")
      return(.cindex_survival(pred, y))
    NA_real_
  }

  compute_yardstick <- function(y, pred, pred_class, prob = NULL) {
    if (task == "survival") {
      stop("Yardstick metrics are not supported for survival tasks.", call. = FALSE)
    }
    if (task == "multiclass") {
      df <- data.frame(truth = y, pred_class = pred_class, stringsAsFactors = FALSE)
      if (!is.null(prob)) {
        prob <- as.data.frame(prob, check.names = FALSE)
        prob_cols <- paste0(".pred_", make.names(class_levels))
        if (ncol(prob) == length(class_levels)) {
          names(prob) <- prob_cols
        }
        df <- cbind(df, prob)
      }
      pred_cols <- grep("^\\.pred_", names(df), value = TRUE)
      args <- list(df, truth = quote(truth), estimate = quote(pred_class))
      for (col in pred_cols) {
        args[[col]] <- as.name(col)
      }
      res <- try(do.call(yardstick_set, args), silent = TRUE)
    } else if (task == "binomial") {
      df <- data.frame(truth = y, .pred = as.numeric(pred), stringsAsFactors = FALSE)
      if (!is.null(pred_class)) df$.pred_class <- pred_class
      res <- try(yardstick_set(df, truth = truth, estimate = .pred_class,
                               .pred, event_level = "second"),
                 silent = TRUE)
    } else {
      df <- data.frame(truth = y, .pred = as.numeric(pred), stringsAsFactors = FALSE)
      res <- try(yardstick_set(df, truth = truth, estimate = .pred),
                 silent = TRUE)
    }
    if (inherits(res, "try-error")) {
      err_msg <- attr(res, "condition")$message
      stop(sprintf("Yardstick metrics failed: %s", err_msg), call. = FALSE)
    }
    stats::setNames(res$.estimate, res$.metric)
  }

  align_probabilities <- function(prob, class_levels) {
    prob_mat <- if (is.data.frame(prob)) as.matrix(prob) else as.matrix(prob)
    if (is.null(colnames(prob_mat))) {
      if (ncol(prob_mat) != length(class_levels)) {
        stop("Probability predictions do not match class levels.", call. = FALSE)
      }
      colnames(prob_mat) <- class_levels
      return(prob_mat)
    }
    exp_cols <- paste0(".pred_", make.names(class_levels))
    if (all(exp_cols %in% colnames(prob_mat))) {
      prob_mat <- prob_mat[, exp_cols, drop = FALSE]
    } else if (all(class_levels %in% colnames(prob_mat))) {
      prob_mat <- prob_mat[, class_levels, drop = FALSE]
    } else if (all(make.names(class_levels) %in% colnames(prob_mat))) {
      prob_mat <- prob_mat[, make.names(class_levels), drop = FALSE]
    } else if (ncol(prob_mat) == length(class_levels)) {
      prob_mat <- prob_mat[, seq_len(ncol(prob_mat)), drop = FALSE]
    } else {
      stop("Probability predictions do not align with class levels.", call. = FALSE)
    }
    colnames(prob_mat) <- class_levels
    prob_mat
  }

  # --- Robust design matrix builder ------------------------------------------
  make_design_matrix <- function(X, ref_cols = NULL) {
    X <- as.data.frame(X)
    X <- X[, !names(X) %in% c("y", "outcome"), drop = FALSE]

    is_num <- vapply(X, is.numeric, logical(1))
    if (all(is_num)) {
      mm <- as.matrix(X)
    } else {
      mf <- stats::model.frame(~ ., data = X, na.action = stats::na.pass)
      mm <- stats::model.matrix(~ . - 1, data = mf)
      mm <- as.matrix(mm)
    }

    if (!is.null(ref_cols)) {
      missing_cols <- setdiff(ref_cols, colnames(mm))
      if (length(missing_cols)) {
        mm <- cbind(
          mm,
          matrix(0, nrow = nrow(mm), ncol = length(missing_cols),
                 dimnames = list(NULL, missing_cols))
        )
      }
      extra_cols <- setdiff(colnames(mm), ref_cols)
      if (length(extra_cols)) {
        mm <- mm[, setdiff(colnames(mm), extra_cols), drop = FALSE]
      }
      mm <- mm[, ref_cols, drop = FALSE]
      return(list(matrix = mm, columns = ref_cols))
    }

    if (!ncol(mm)) {
      stop("All predictors have zero variance after preprocessing.")
    }

    keep <- apply(mm, 2, sd, na.rm = TRUE) > 0
    if (!any(keep)) {
      stop("All predictors have zero variance after preprocessing.")
    }
    mm <- mm[, keep, drop = FALSE]

    list(matrix = mm, columns = colnames(mm))
  }

  resolve_weights <- function(y, weights_spec) {
    if (is.null(weights_spec) || !task %in% c("binomial", "multiclass")) return(NULL)
    if (!is.numeric(weights_spec)) stop("class_weights must be numeric.")
    cw <- weights_spec
    if (is.null(names(cw))) {
      if (length(cw) != length(class_levels)) {
        stop("Provide class_weights as a named vector matching outcome levels.")
      }
      names(cw) <- class_levels[seq_along(cw)]
    }
    missing_levels <- setdiff(class_levels, names(cw))
    if (length(missing_levels)) {
      stop(sprintf("class_weights missing levels: %s", paste(missing_levels, collapse = ", ")))
    }
    cw[as.character(y)]
  }

  resolve_args <- function(name, defaults = list()) {
    if (!length(learner_args)) return(modifyList(defaults, list()))
    if (!is.null(names(learner_args)) && all(names(learner_args) %in% learner)) {
      extras <- learner_args[[name]] %||% list()
    } else {
      extras <- learner_args
    }
    modifyList(defaults, extras)
  }

  resolve_custom_args <- function(name) {
    if (!length(learner_args)) return(list(fit = list(), predict = list()))
    if (!is.null(names(learner_args)) && all(names(learner_args) %in% learner)) {
      extras <- learner_args[[name]] %||% list()
    } else {
      extras <- learner_args
    }
    if (is.list(extras) && (("fit" %in% names(extras)) || ("predict" %in% names(extras)))) {
      return(list(fit = extras$fit %||% list(),
                  predict = extras$predict %||% list()))
    }
    list(fit = extras, predict = extras)
  }

  # single learner wrapper ----------------------------------------------------
  train_one_learner <- function(learner_obj, learner_label, Xtrg, ytr, Xteg, yte, weights = NULL) {
    if (inherits(learner_obj, "model_spec")) {
      if (!requireNamespace("parsnip", quietly = TRUE)) {
        stop("Package 'parsnip' is required when learner is a model_spec.")
      }
      if (!is.null(weights) && !inherits(weights, "hardhat_case_weights")) {
        if (!requireNamespace("hardhat", quietly = TRUE)) {
          stop("Package 'hardhat' is required for case weights with parsnip learners.")
        }
        weights <- hardhat::frequency_weights(weights)
      }
      y_for_fit <- if (task %in% c("binomial", "multiclass")) {
        factor(ytr, levels = class_levels)
      } else if (task == "survival") {
        ytr
      } else {
        as.numeric(ytr)
      }
      fit <- if (is.null(weights)) {
        parsnip::fit_xy(learner_obj, x = Xtrg, y = y_for_fit)
      } else {
        parsnip::fit_xy(learner_obj, x = Xtrg, y = y_for_fit, case_weights = weights)
      }
      if (task == "binomial") {
        prob <- try(stats::predict(fit, new_data = Xteg, type = "prob"), silent = TRUE)
        if (inherits(prob, "try-error")) {
          # Extract the actual error text from parsnip/xgboost
          err_msg <- attr(prob, "condition")$message
          stop(sprintf("Parsnip learner '%s' failed to predict: %s",
                       learner_label, err_msg))
        }
        prob_df <- as.data.frame(prob)
        pos_col <- paste0(".pred_", make.names(class_levels[2]))
        if (!pos_col %in% names(prob_df)) {
          if (ncol(prob_df) >= 2L) {
            pos_col <- names(prob_df)[2]
          } else {
            stop(sprintf("Parsnip learner '%s' did not return class probabilities.",
                         learner_label))
          }
        }
        pred <- as.numeric(prob_df[[pos_col]])
        pred_class <- factor(ifelse(pred >= 0.5, class_levels[2], class_levels[1]),
                             levels = class_levels)
        return(list(pred = pred, pred_class = pred_class, fit = fit))
      }
      if (task == "multiclass") {
        prob <- try(stats::predict(fit, new_data = Xteg, type = "prob"), silent = TRUE)
        if (inherits(prob, "try-error")) {
          err_msg <- attr(prob, "condition")$message
          stop(sprintf("Parsnip learner '%s' failed to predict: %s",
                       learner_label, err_msg))
        }
        prob_mat <- align_probabilities(prob, class_levels)
        class_pred <- try(stats::predict(fit, new_data = Xteg, type = "class"), silent = TRUE)
        if (inherits(class_pred, "try-error")) {
          err_msg <- attr(class_pred, "condition")$message
          stop(sprintf("Parsnip learner '%s' failed to predict classes: %s",
                       learner_label, err_msg))
        }
        class_df <- as.data.frame(class_pred)
        pred_class <- factor(as.character(class_df[[1]]), levels = class_levels)
        return(list(pred = pred_class, pred_class = pred_class, prob = prob_mat, fit = fit))
      }
      if (task == "survival") {
        pred_df <- try(stats::predict(fit, new_data = Xteg, type = "numeric"), silent = TRUE)
        if (inherits(pred_df, "try-error")) {
          pred_df <- try(stats::predict(fit, new_data = Xteg, type = "risk"), silent = TRUE)
        }
        if (inherits(pred_df, "try-error")) {
          err_msg <- attr(pred_df, "condition")$message
          stop(sprintf("Parsnip learner '%s' failed to predict: %s",
                       learner_label, err_msg))
        }
        pred <- as.numeric(as.data.frame(pred_df)[[1]])
        return(list(pred = pred, fit = fit))
      } else {
        pred_df <- stats::predict(fit, new_data = Xteg, type = "numeric")
        pred <- as.numeric(pred_df[[1]])
        return(list(pred = pred, fit = fit))
      }
    }
    if (learner_obj %in% names(custom_learners)) {
      def <- custom_learners[[learner_obj]]
      args <- resolve_custom_args(learner_obj)
      fit_args <- c(list(x = Xtrg, y = ytr, task = task, weights = weights), args$fit)
      model <- do.call(def$fit, fit_args)
      pred_args <- c(list(object = model, newdata = Xteg, task = task), args$predict)
      pred <- do.call(def$predict, pred_args)
      if (task == "multiclass") {
        if (is.data.frame(pred) || is.matrix(pred)) {
          prob_mat <- align_probabilities(pred, class_levels)
          pred_class <- factor(class_levels[max.col(prob_mat, ties.method = "first")],
                               levels = class_levels)
          return(list(pred = pred_class, pred_class = pred_class, prob = prob_mat, fit = model))
        }
        if (is.factor(pred) || is.character(pred)) {
          pred_class <- factor(as.character(pred), levels = class_levels)
          if (length(pred_class) != nrow(Xteg)) {
            stop(sprintf("Custom learner '%s' returned %d predictions for %d rows.",
                         learner_obj, length(pred_class), nrow(Xteg)))
          }
          return(list(pred = pred_class, pred_class = pred_class, fit = model))
        }
        stop(sprintf("Custom learner '%s' must return class labels or class probabilities for multiclass tasks.",
                     learner_obj))
      }
      pred <- as.numeric(pred)
      if (length(pred) != nrow(Xteg)) {
        stop(sprintf("Custom learner '%s' returned %d predictions for %d rows.",
                     learner_obj, length(pred), nrow(Xteg)))
      }
      return(list(pred = pred, fit = model))
    }
    if (learner_obj == "glmnet") {
      if (!requireNamespace("glmnet", quietly = TRUE)) stop("Install 'glmnet'.")
      fam <- if (task == "binomial") "binomial"
      else if (task == "multiclass") "multinomial"
      else if (task == "survival") "cox"
      else "gaussian"
      la  <- resolve_args("glmnet", list(alpha = 0.9, standardize = FALSE))
      Xtr_design <- make_design_matrix(Xtrg)
      Xte_design <- make_design_matrix(Xteg, ref_cols = Xtr_design$columns)
      y_for_fit <- if (task == "binomial") {
        as.numeric(factor(ytr, levels = class_levels)) - 1
      } else if (task == "multiclass") {
        factor(ytr, levels = class_levels)
      } else if (task == "survival") {
        ytr
      } else {
        as.numeric(ytr)
      }
      cv_args <- c(list(x = Xtr_design$matrix,
                        y = y_for_fit,
                        family = fam,
                        alpha = la$alpha,
                        standardize = la$standardize %||% FALSE,
                        weights = weights),
                   la[setdiff(names(la), c("alpha", "standardize"))])
      cvfit <- do.call(glmnet::cv.glmnet, cv_args)
      if (task == "multiclass") {
        pred_arr <- predict(cvfit, Xte_design$matrix, s = "lambda.min", type = "response")
        if (length(dim(pred_arr)) == 3L) {
          prob_mat <- pred_arr[, , 1, drop = FALSE][,,1]
        } else {
          prob_mat <- pred_arr
        }
        prob_mat <- align_probabilities(prob_mat, class_levels)
        pred_class <- factor(class_levels[max.col(prob_mat, ties.method = "first")],
                             levels = class_levels)
        return(list(pred = pred_class, pred_class = pred_class, prob = prob_mat, fit = cvfit))
      }
      pred_type <- if (task == "survival") "link" else "response"
      pred  <- as.numeric(predict(cvfit, Xte_design$matrix, s = "lambda.min",
                                  type = pred_type))

      return(list(pred = pred, fit = cvfit))
    }
    if (learner_obj == "ranger") {
      if (!requireNamespace("ranger", quietly = TRUE)) stop("Install 'ranger'.")
      if (task == "survival") {
        stop("Learner 'ranger' does not support survival tasks in bioLeak; use parsnip/workflow or a custom learner.")
      }
      y_for_fit <- if (task %in% c("binomial", "multiclass")) {
        factor(ytr, levels = class_levels)
      } else {
        as.numeric(ytr)
      }
      dftr <- data.frame(y = y_for_fit, Xtrg, check.names = FALSE)
      frm  <- stats::as.formula("y ~ .")
      rng_args <- resolve_args("ranger", list())
      if (task %in% c("binomial", "multiclass")) {
        rng_args$class.weights <- rng_args$class.weights %||% class_weights
      }
      rg   <- do.call(ranger::ranger, c(list(formula = frm, data = dftr,
                                             probability = (task %in% c("binomial", "multiclass"))),
                                        rng_args))
      dfte <- data.frame(Xteg, check.names = FALSE)
      pr   <- predict(rg, dfte)
      if (task == "binomial") {
        pred <- pr$predictions[, 2]
        pred_class <- factor(ifelse(pred >= 0.5, class_levels[2], class_levels[1]),
                             levels = class_levels)
        return(list(pred = pred, pred_class = pred_class, fit = rg))
      }
      if (task == "multiclass") {
        prob_mat <- pr$predictions
        prob_mat <- align_probabilities(prob_mat, class_levels)
        pred_class <- factor(class_levels[max.col(prob_mat, ties.method = "first")],
                             levels = class_levels)
        return(list(pred = pred_class, pred_class = pred_class, prob = prob_mat, fit = rg))
      }
      pred <- pr$predictions
      return(list(pred = pred, fit = rg))
    }
    stop("Unsupported learner.")
  }

  train_one_workflow <- function(learner_obj, learner_label, dftr, dfte, weights = NULL) {
    if (!requireNamespace("workflows", quietly = TRUE)) {
      stop("Package 'workflows' is required when learner is a workflow.")
    }
    fit <- try(workflows::fit(learner_obj, data = dftr), silent = TRUE)
    if (inherits(fit, "try-error")) {
      err_msg <- attr(fit, "condition")$message
      stop(sprintf("Workflow learner '%s' failed to fit: %s", learner_label, err_msg))
    }
    if (task == "binomial") {
      prob <- try(stats::predict(fit, new_data = dfte, type = "prob"), silent = TRUE)
      if (inherits(prob, "try-error")) {
        err_msg <- attr(prob, "condition")$message
        stop(sprintf("Workflow learner '%s' failed to predict: %s", learner_label, err_msg))
      }
      prob_df <- as.data.frame(prob)
      pos_col <- paste0(".pred_", make.names(class_levels[2]))
      if (!pos_col %in% names(prob_df)) {
        if (ncol(prob_df) >= 2L) {
          pos_col <- names(prob_df)[2]
        } else {
          stop(sprintf("Workflow learner '%s' did not return class probabilities.",
                       learner_label))
        }
      }
      pred <- as.numeric(prob_df[[pos_col]])
      pred_class <- factor(ifelse(pred >= 0.5, class_levels[2], class_levels[1]),
                           levels = class_levels)
      return(list(pred = pred, pred_class = pred_class, fit = fit))
    }
    if (task == "multiclass") {
      prob <- try(stats::predict(fit, new_data = dfte, type = "prob"), silent = TRUE)
      if (inherits(prob, "try-error")) {
        err_msg <- attr(prob, "condition")$message
        stop(sprintf("Workflow learner '%s' failed to predict: %s", learner_label, err_msg))
      }
      prob_mat <- align_probabilities(prob, class_levels)
      class_pred <- try(stats::predict(fit, new_data = dfte, type = "class"), silent = TRUE)
      if (inherits(class_pred, "try-error")) {
        err_msg <- attr(class_pred, "condition")$message
        stop(sprintf("Workflow learner '%s' failed to predict classes: %s",
                     learner_label, err_msg))
      }
      class_df <- as.data.frame(class_pred)
      pred_class <- factor(as.character(class_df[[1]]), levels = class_levels)
      return(list(pred = pred_class, pred_class = pred_class, prob = prob_mat, fit = fit))
    }
    if (task == "survival") {
      pred_df <- try(stats::predict(fit, new_data = dfte, type = "numeric"), silent = TRUE)
      if (inherits(pred_df, "try-error")) {
        pred_df <- try(stats::predict(fit, new_data = dfte, type = "risk"), silent = TRUE)
      }
      if (inherits(pred_df, "try-error")) {
        err_msg <- attr(pred_df, "condition")$message
        stop(sprintf("Workflow learner '%s' failed to predict: %s", learner_label, err_msg))
      }
      pred <- as.numeric(as.data.frame(pred_df)[[1]])
      return(list(pred = pred, fit = fit))
    } else {
      pred_df <- stats::predict(fit, new_data = dfte, type = "numeric")
      pred <- as.numeric(pred_df[[1]])
      return(list(pred = pred, fit = fit))
    }
  }

  resolve_fold_indices <- function(fold) {
    if (!isTRUE(compact) || !is.null(fold$train)) return(fold)
    if (is.null(fold_assignments) || !length(fold_assignments)) {
      stop("Compact splits require fold assignments to compute indices.")
    }
    r <- fold$repeat_id
    if (is.null(r) || !is.finite(r)) r <- 1L
    assign_vec <- fold_assignments[[r]]
    if (is.null(assign_vec)) {
      stop(sprintf("Missing fold assignments for repeat %s.", r))
    }
    test <- which(assign_vec == fold$fold)
    if (identical(split_mode, "time_series")) {
      if (is.null(time_vec) || !length(time_vec)) {
        stop("time_series compact splits require time column values.")
      }
      if (!length(test)) {
        train <- integer(0)
      } else {
        tmin <- min(time_vec[test])
        if (split_horizon == 0) {
          train <- which(time_vec < tmin)
        } else {
          train <- which(time_vec <= (tmin - split_horizon))
        }
      }
    } else {
      train <- setdiff(seq_len(nrow(Xall)), test)
    }
    fold_seq <- fold$fold_seq %||% fold$fold
    list(train = train, test = test, fold = fold$fold,
         repeat_id = fold$repeat_id, fold_seq = fold_seq)
  }

  make_fold_df <- function(X, y) {
    df <- as.data.frame(X, check.names = FALSE)
    if (length(outcome) == 2L) {
      if (!inherits(y, "Surv")) {
        stop("Survival tasks require a Surv outcome when building fold data.")
      }
      y_mat <- as.matrix(y)
      if (ncol(y_mat) < 2L) {
        stop("Survival outcome must include time and event columns.")
      }
      df[[outcome[[1]]]] <- y_mat[, 1]
      df[[outcome[[2]]]] <- y_mat[, ncol(y_mat)]
      df <- df[, c(outcome, setdiff(names(df), outcome)), drop = FALSE]
      return(df)
    }
    df[[outcome]] <- y
    df[, c(outcome, setdiff(names(df), outcome)), drop = FALSE]
  }

  # fold-level function -------------------------------------------------------
  do_fold <- function(fold) {
    fold_full <- resolve_fold_indices(fold)
    fold_id <- fold_full$fold_seq %||% fold_full$fold
    set.seed(seed + fold_id)
    tr <- fold_full$train
    te <- fold_full$test

    Xtr <- Xall[tr, , drop = FALSE]
    Xte <- Xall[te, , drop = FALSE]
    ytr <- yall[tr]
    yte <- yall[te]

    if (task %in% c("binomial", "multiclass")) {
      ytr <- factor(ytr, levels = class_levels)
      yte <- factor(yte, levels = class_levels)
      if (nlevels(droplevels(ytr)) < 2) {
        warning(sprintf("Fold %s skipped: only one class in training data", fold_id))
        empty_metrics <- if (identical(metric_mode, "yardstick")) {
          numeric(0)
        } else {
          setNames(rep(NA_real_, length(metrics)), metric_labels)
        }
        skipped <- lapply(learner_names, function(ln) {
          list(
            metrics = empty_metrics,
            pred = data.frame(
              id = integer(0),
              truth = factor(character(0), levels = class_levels),
              pred = numeric(0),
              fold = integer(0),
              learner = character(0),
              stringsAsFactors = FALSE
            ),
            guard = list(state = NULL),
            learner = NULL,
            feat_names = colnames(Xtr)
          )
        })
        names(skipped) <- learner_names
        return(skipped)
      }
      fold_weights <- resolve_weights(ytr, class_weights)
    } else if (task == "survival") {
      fold_weights <- NULL
    } else {
      ytr <- as.numeric(ytr)
      yte <- as.numeric(yte)
      fold_weights <- NULL
    }

    preprocess_state <- NULL
    Xtrg <- Xtr
    Xteg <- Xte
    dftr <- NULL
    dfte <- NULL

    if (identical(preprocess_mode, "guard")) {
      guard <- .guard_fit(
        X = Xtr,
        y = ytr,
        steps = if (exists("preprocess") && is.list(preprocess)) preprocess else list(),
        task  = task
      )
      Xtrg <- guard$transform(Xtr)
      Xteg <- guard$transform(Xte)
      preprocess_state <- guard$state
      colnames(Xtrg) <- make.names(colnames(Xtrg))
      colnames(Xteg) <- make.names(colnames(Xteg))
    } else if (identical(preprocess_mode, "recipe")) {
      dftr <- make_fold_df(Xtr, ytr)
      dfte <- make_fold_df(Xte, yte)
      recipe_prep <- recipes::prep(preprocess, training = dftr, retain = TRUE)
      Xtrg <- as.data.frame(recipes::juice(recipe_prep, recipes::all_predictors()),
                            check.names = FALSE)
      Xteg <- as.data.frame(recipes::bake(recipe_prep, new_data = dfte,
                                          recipes::all_predictors()),
                            check.names = FALSE)
      if (!ncol(Xtrg)) stop("All predictors have zero variance after preprocessing.")
      preprocess_state <- list(type = "recipe", recipe = recipe_prep)
      colnames(Xtrg) <- make.names(colnames(Xtrg))
      colnames(Xteg) <- make.names(colnames(Xteg))
    } else {
      dftr <- make_fold_df(Xtr, ytr)
      dfte <- make_fold_df(Xte, yte)
      preprocess_state <- list(type = "workflow")
    }

    results <- list()
    for (i in seq_along(learner_names)) {
      ln <- learner_names[[i]]
      learner_obj <- learner_objs[[i]]
      # Train one learner; ensure consistent level matching
      if (identical(preprocess_mode, "workflow")) {
        dfte_pred <- dfte[, setdiff(names(dfte), outcome), drop = FALSE]
        model <- train_one_workflow(learner_obj, ln, dftr, dfte_pred, weights = fold_weights)
        feat_names <- setdiff(names(dftr), outcome)
      } else {
        model <- train_one_learner(learner_obj, ln, Xtrg, ytr, Xteg, yte, weights = fold_weights)
        feat_names <- colnames(Xtrg)
      }

      pred_class <- if (task == "binomial") {
        if (!is.null(model$pred_class)) {
          model$pred_class
        } else if (is.factor(model$pred) || is.character(model$pred)) {
          factor(as.character(model$pred), levels = class_levels)
        } else {
          factor(ifelse(model$pred >= 0.5, class_levels[2], class_levels[1]),
                 levels = class_levels)
        }
      } else if (task == "multiclass") {
        model$pred_class %||% {
          if (is.factor(model$pred) || is.character(model$pred)) {
            factor(as.character(model$pred), levels = class_levels)
          } else {
            NULL
          }
        }
      } else {
        NULL
      }

      prob_mat <- model$prob %||% NULL

      ms <- if (identical(metric_mode, "yardstick")) {
        compute_yardstick(yte, model$pred, pred_class, prob = prob_mat)
      } else if (task == "multiclass") {
        vals <- vapply(seq_along(metrics), function(idx) {
          mname <- metrics[[idx]]
          if (is.function(mname)) return(mname(yte, pred_class))
          if (identical(mname, "accuracy")) return(.multiclass_accuracy(yte, pred_class))
          if (identical(mname, "macro_f1")) return(.multiclass_macro_f1(yte, pred_class))
          if (identical(mname, "log_loss")) {
            if (is.null(prob_mat)) {
              stop("log_loss requires class probability predictions for multiclass tasks.")
            }
            return(.multiclass_log_loss(yte, prob_mat))
          }
          NA_real_
        }, numeric(1))
        names(vals) <- metric_labels
        vals
      } else {
        vals <- vapply(seq_along(metrics), function(idx) {
          compute_metric(metrics[[idx]], yte, model$pred)
        }, numeric(1))
        names(vals) <- metric_labels
        vals
      }
      pred_tbl <- if (task == "survival") {
        yte_mat <- as.matrix(yte)
        time_col <- if ("time" %in% colnames(yte_mat)) "time" else colnames(yte_mat)[1]
        status_col <- if ("status" %in% colnames(yte_mat)) "status" else colnames(yte_mat)[ncol(yte_mat)]
        data.frame(
          id = ids[te],
          truth_time = yte_mat[, time_col],
          truth_event = yte_mat[, status_col],
          pred = model$pred,
          fold = fold_id,
          learner = ln,
          stringsAsFactors = FALSE
        )
      } else {
        data.frame(
          id = ids[te],
          truth = yte,
          pred = model$pred,
          fold = fold_id,
          learner = ln,
          stringsAsFactors = FALSE
        )
      }
      if (!is.null(pred_class)) pred_tbl$pred_class <- pred_class
      if (!is.null(prob_mat) && task == "multiclass") {
        prob_df <- as.data.frame(prob_mat, check.names = FALSE)
        names(prob_df) <- paste0(".pred_", make.names(class_levels))
        pred_tbl <- cbind(pred_tbl, prob_df)
      }
      results[[ln]] <- list(
        metrics = ms,
        pred = pred_tbl,
        guard = preprocess_state,
        learner = model$fit,
        feat_names = feat_names
      )
    }

    results
  }


  folds <- splits@indices
  nfold <- length(folds)
  fold_errors <- rep(NA_character_, nfold)

  # progress bar --------------------------------------------------------------
  pb <- utils::txtProgressBar(min = 0, max = nfold, style = 3)
  pb_counter <- 0
  progress_wrap <- function(f) {
    fold_id <- f$fold_seq %||% f$fold
    res <- tryCatch(do_fold(f), error = function(e) {
      fold_errors[[fold_id]] <<- conditionMessage(e)
      warning(sprintf("Fold %s failed: %s", fold_id, e$message)); NULL
    })
    pb_counter <<- pb_counter + 1
    utils::setTxtProgressBar(pb, pb_counter)
    res
  }

  # parallel or sequential execution -----------------------------------------
  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    out <- future.apply::future_lapply(seq_along(folds), function(i) {
      fold <- folds[[i]]
      fold$fold_seq <- i
      progress_wrap(fold)
    }, future.seed = TRUE)
  } else {
    out <- lapply(seq_along(folds), function(i) {
      fold <- folds[[i]]
      fold$fold_seq <- i
      progress_wrap(fold)
    })
  }

  close(pb)

  # collect results -----------------------------------------------------------
  met_rows <- list()
  preds <- list()
  guards <- list()
  lears <- list()
  featn <- NULL
  audit_rows <- list()
  fold_status_rows <- vector("list", length(out))

  for (fold_idx in seq_along(out)) {
    fold_res <- out[[fold_idx]]
    learner_total <- length(learner_names)
    learner_success <- 0L
    if (is.null(fold_res)) next
    fold_info <- resolve_fold_indices(folds[[fold_idx]])
    fold_id <- fold_idx
    for (ln in names(fold_res)) {
      res <- fold_res[[ln]]
      if (is.null(res)) next
      m <- res$metrics
      if (all(is.na(m))) {
        next
      }
      learner_success <- learner_success + 1L
      metric_row <- c(list(fold = fold_id, learner = ln), as.list(m))
      met_rows[[length(met_rows) + 1]] <- as.data.frame(metric_row,
                                                        row.names = NULL,
                                                        check.names = FALSE)
      if (is.data.frame(res$pred) || is.matrix(res$pred)) {
        preds[[length(preds) + 1]] <- res$pred
      }
      if (!is.null(res$guard) && is.list(res$guard)) {
        guards[[length(guards) + 1]] <- res$guard
        if (is.null(featn) && !is.null(res$feat_names)) featn <- res$feat_names
      } else {
        guards[[length(guards) + 1]] <- NULL
      }
      lears[[length(lears) + 1]] <- res$learner
      audit_rows[[length(audit_rows) + 1]] <- data.frame(
        fold = fold_id,
        n_train = length(fold_info$train),
        n_test = length(fold_info$test),
        learner = ln,
        features_final = if (!is.null(res$guard) &&
                              is.list(res$guard) &&
                              !is.null(res$guard$filter) &&
                              !is.null(res$guard$filter$keep)) {
          sum(res$guard$filter$keep)
        } else if (!is.null(res$feat_names)) {
          length(res$feat_names)
        } else NA_integer_,
        row.names = NULL
      )
    }
    all_preds_empty <- length(fold_res) > 0 && all(vapply(fold_res, function(res) {
      is.null(res$pred) || nrow(res$pred) == 0L
    }, logical(1)))
    fold_status_rows[[fold_idx]] <- data.frame(
      fold = fold_idx,
      stage = if (learner_success > 0L) "fold_run" else "fold_precheck",
      status = if (learner_success > 0L) "success" else "skipped",
      reason = if (learner_success > 0L) NA_character_ else if (all_preds_empty) {
        "single_class_training"
      } else {
        "no_valid_metrics"
      },
      notes = sprintf("Successful learners: %d/%d", learner_success, learner_total),
      stringsAsFactors = FALSE
    )
  }

  for (fold_idx in seq_along(out)) {
    if (!is.null(fold_status_rows[[fold_idx]])) next
    err_note <- fold_errors[[fold_idx]]
    if (is.na(err_note)) err_note <- "Fold failed before producing any learner results."
    fold_status_rows[[fold_idx]] <- data.frame(
      fold = fold_idx,
      stage = "fold_run",
      status = "failed",
      reason = "fold_error",
      notes = err_note,
      stringsAsFactors = FALSE
    )
  }
  fold_status_df <- do.call(rbind, fold_status_rows)

  if (!length(met_rows)) {
    stop("No successful folds were completed. Check learner and preprocessing settings.")
  }

  met_df <- do.call(rbind, met_rows)
  audit_df <- do.call(rbind, audit_rows)
  metrics_used <- setdiff(colnames(met_df), c("fold", "learner"))

  # summarize metrics ---------------------------------------------------------
  metric_summary_raw <- aggregate(. ~ learner, data = met_df[, -1, drop = FALSE],
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

  # optional refit ------------------------------------------------------------
  final_model <- NULL
  final_guard <- NULL
  if (refit) {
    ln <- learner_names[[1]]
    learner_obj <- learner_objs[[1]]
    full_weights <- if (task %in% c("binomial", "multiclass")) resolve_weights(yall, class_weights) else NULL
    if (identical(preprocess_mode, "guard")) {
      guard_full <- .guard_fit(Xall, yall, preprocess, task)
      Xfullg <- guard_full$transform(Xall)
      colnames(Xfullg) <- make.names(colnames(Xfullg))
      final_model <- train_one_learner(learner_obj, ln, Xfullg, yall, Xfullg, yall,
                                       weights = full_weights)$fit
      final_guard <- guard_full$state
    } else if (identical(preprocess_mode, "recipe")) {
      df_full <- make_fold_df(Xall, yall)
      recipe_prep <- recipes::prep(preprocess, training = df_full, retain = TRUE)
      Xfullg <- as.data.frame(recipes::juice(recipe_prep, recipes::all_predictors()),
                              check.names = FALSE)
      colnames(Xfullg) <- make.names(colnames(Xfullg))
      final_model <- train_one_learner(learner_obj, ln, Xfullg, yall, Xfullg, yall,
                                       weights = full_weights)$fit
      final_guard <- list(type = "recipe", recipe = recipe_prep)
    } else {
      df_full <- make_fold_df(Xall, yall)
      final_model <- workflows::fit(learner_obj, data = df_full)
      final_guard <- list(type = "workflow")
    }
  }

  # assemble LeakFit object ---------------------------------------------------
  perm_refit_spec <- NULL
  if (isTRUE(store_refit_data)) {
    perm_refit_spec <- list(
      x = x,
      outcome = outcome,
      preprocess = preprocess,
      learner = learner_input,
      learner_args = learner_args,
      custom_learners = custom_learners,
      class_weights = class_weights,
      positive_class = positive_class,
      parallel = parallel
    )
  }

  new("LeakFit",
      splits = splits,
      metrics = met_df,
      metric_summary = metric_summary,
      audit = audit_df,
      predictions = preds,
      preprocess = guards,
      learners = lears,
      outcome = outcome,
      task = task,
      feature_names = featn,
      info = list(hash = .bio_hash_indices(folds),
                  metrics_requested = metrics_input,
                  metrics_used = metrics_used,
                  class_weights = class_weights,
                  positive_class = if (task == "binomial") class_levels[2] else NULL,
                  sample_ids = ids,
                  fold_seeds = setNames(seed + seq_along(folds),
                                        paste0("fold", seq_along(folds))),
                  fold_status = fold_status_df,
                  refit = refit,
                  final_model = final_model,
                  final_preprocess = final_guard,
                  learner_names = learner_names,
                  perm_refit_spec = perm_refit_spec))
}
