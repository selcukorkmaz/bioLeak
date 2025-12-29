#' Fit and evaluate with leakage guards over predefined splits
#'
#' Performs cross-validated model training and evaluation using
#' leakage-protected preprocessing (.guard_fit) and user-specified learners.
#'
#' @param x SummarizedExperiment or matrix/data.frame
#' @param outcome outcome column (if x is SE)
#' @param splits LeakSplits object from make_splits()
#' @param preprocess list(impute, normalize, filter=list(...), fs)
#' @param learner parsnip model_spec (or list of model_spec objects) describing
#'   the model(s) to fit. For legacy use, a character vector of learner names
#'   (e.g., "glmnet", "ranger") or custom learner IDs is still supported.
#' @param learner_args list of additional arguments passed to legacy learners
#'   (ignored when `learner` is a parsnip model_spec).
#' @param custom_learners named list of custom learner definitions used only
#'   with legacy character learners. Each entry
#'   must contain \code{fit} and \code{predict} functions. The \code{fit} function
#'   should accept \code{x}, \code{y}, \code{task}, and \code{weights}, and return
#'   a model object. The \code{predict} function should accept \code{object},
#'   \code{newdata}, and \code{task}, and return numeric predictions.
#' @param metrics named list of metric functions or vector of metric names
#' @param class_weights optional named numeric vector of weights for binomial outcomes
#' @param positive_class optional value indicating the positive class for binomial outcomes.
#'   When set, the outcome levels are reordered so that \code{positive_class} is treated
#'   as the positive class (level 2). If NULL, the second factor level is used.
#' @param parallel logical, use future.apply for multicore execution
#' @param refit logical, if TRUE retrain final model on full data
#' @param seed integer, for reproducibility
#' @return LeakFit object
#' @details
#' Preprocessing is fit on the training fold and applied to the test fold,
#' preventing leakage from global imputation, scaling, or feature selection.
#' For data.frame or matrix inputs, columns used to define splits
#' (outcome, group, batch, study, time) are excluded from the predictor matrix.
#' Use \code{learner_args} to pass model-specific arguments, either as a named
#' list keyed by learner or a single list applied to all learners. For custom
#' learners, \code{learner_args[[name]]} may be a list with \code{fit} and
#' \code{predict} sublists to pass distinct arguments to each stage. For binomial
#' tasks, predictions and metrics assume the positive class is the second factor
#' level; use \code{positive_class} to control this. Parsnip learners must support
#' probability predictions for binomial metrics (AUC/PR-AUC/accuracy).
#' @examples
#' \dontrun{
#' set.seed(1)
#' df <- data.frame(
#'   subject = rep(1:10, each = 2),
#'   outcome = rbinom(20, 1, 0.5),
#'   x1 = rnorm(20),
#'   x2 = rnorm(20)
#' )
#' splits <- make_splits(df, outcome = "outcome",
#'                       mode = "subject_grouped", group = "subject", v = 5)
#' fit <- fit_resample(df, outcome = "outcome", splits = splits,
#'                     learner = "glmnet", metrics = "auc")
#' summary(fit)
#'
#' # Parsnip model_spec
#' \dontrun{
#' if (requireNamespace("parsnip", quietly = TRUE)) {
#'   spec <- parsnip::boost_tree(mode = "classification", trees = 200) |>
#'     parsnip::set_engine("xgboost")
#'   fit3 <- fit_resample(df, outcome = "outcome", splits = splits,
#'                        learner = spec, metrics = "auc")
#' }
#' }
#'
#' # Custom learner (logistic regression)
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
#' }
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
                         seed = 1) {

  set.seed(seed)
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
  use_parsnip <- FALSE
  learner_specs <- NULL
  learner_names <- NULL

  if (is_parsnip_spec(learner)) {
    use_parsnip <- TRUE
    learner_specs <- list(learner)
  } else if (is.list(learner) && length(learner) &&
             all(vapply(learner, is_parsnip_spec, logical(1)))) {
    use_parsnip <- TRUE
    learner_specs <- learner
  }

  if (use_parsnip) {
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
      stop("learner must be a character vector of legacy learners or a parsnip model_spec.")
    }
    learner <- match.arg(learner, choices = all_learners, several.ok = TRUE)
    learner_names <- learner
  }
  Xall <- .bio_get_x(x)
  yall <- .bio_get_y(x, outcome)
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
  ids  <- seq_len(nrow(Xall))
  task <- if (.bio_is_binomial(yall)) "binomial"
  else if (.bio_is_regression(yall)) "gaussian"
  else if (is.factor(yall) && nlevels(yall) == 2) "binomial"
  else stop("Unsupported outcome type: require binomial (2-class) or numeric regression outcome.")

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
  } else {
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
  }

  metrics_input <- metrics
  if (is.null(metrics)) {
    metrics <- if (task == "binomial") c("auc", "pr_auc", "accuracy") else c("rmse")
  }

  if (is.character(metrics)) {
    allowed <- if (task == "binomial") c("auc", "pr_auc", "accuracy") else c("rmse", "cindex")
    invalid <- setdiff(metrics, allowed)
    if (length(invalid)) {
      warning(sprintf("Dropping metrics not applicable to %s task: %s", task,
                      paste(invalid, collapse = ", ")))
      metrics <- setdiff(metrics, invalid)
    }
    if (!length(metrics)) {
      metrics <- if (task == "binomial") c("auc", "pr_auc", "accuracy") else c("rmse")
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

  learner_objs <- if (use_parsnip) learner_specs else as.list(learner_names)

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
    NA_real_
  }

  # --- Robust design matrix builder ------------------------------------------
  make_design_matrix <- function(X, ref_cols = NULL) {
    X <- as.data.frame(X)
    X <- X[, !names(X) %in% c("y", "outcome"), drop = FALSE]

    mf <- stats::model.frame(~ ., data = X, na.action = stats::na.pass)
    mm <- stats::model.matrix(~ . - 1, data = mf)
    mm <- as.matrix(mm)

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
    if (is.null(weights_spec) || task != "binomial") return(NULL)
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
      y_for_fit <- if (task == "binomial") factor(ytr, levels = class_levels) else as.numeric(ytr)
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
      } else {
        pred_df <- stats::predict(fit, new_data = Xteg, type = "numeric")
        pred <- as.numeric(pred_df[[1]])
      }
      return(list(pred = pred, fit = fit))
    }
    if (learner_obj %in% names(custom_learners)) {
      def <- custom_learners[[learner_obj]]
      args <- resolve_custom_args(learner_obj)
      fit_args <- c(list(x = Xtrg, y = ytr, task = task, weights = weights), args$fit)
      model <- do.call(def$fit, fit_args)
      pred_args <- c(list(object = model, newdata = Xteg, task = task), args$predict)
      pred <- do.call(def$predict, pred_args)
      pred <- as.numeric(pred)
      if (length(pred) != nrow(Xteg)) {
        stop(sprintf("Custom learner '%s' returned %d predictions for %d rows.",
                     learner_obj, length(pred), nrow(Xteg)))
      }
      return(list(pred = pred, fit = model))
    }
    if (learner_obj == "glmnet") {
      if (!requireNamespace("glmnet", quietly = TRUE)) stop("Install 'glmnet'.")
      fam <- if (task == "binomial") "binomial" else "gaussian"
      la  <- resolve_args("glmnet", list(alpha = 0.9, standardize = FALSE))
      Xtr_design <- make_design_matrix(Xtrg)
      Xte_design <- make_design_matrix(Xteg, ref_cols = Xtr_design$columns)
      y_for_fit <- if (task == "binomial") {
        as.numeric(factor(ytr, levels = class_levels)) - 1
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
      pred  <- as.numeric(predict(cvfit, Xte_design$matrix, s = "lambda.min",
                                  type = "response"))

      return(list(pred = pred, fit = cvfit))
    }
    if (learner_obj == "ranger") {
      if (!requireNamespace("ranger", quietly = TRUE)) stop("Install 'ranger'.")
      y_for_fit <- if (task == "binomial") factor(ytr, levels = class_levels) else as.numeric(ytr)
      dftr <- data.frame(y = y_for_fit, Xtrg, check.names = FALSE)
      frm  <- stats::as.formula("y ~ .")
      rng_args <- resolve_args("ranger", list())
      if (task == "binomial") {
        rng_args$class.weights <- rng_args$class.weights %||% class_weights
      }
      rg   <- do.call(ranger::ranger, c(list(formula = frm, data = dftr,
                                             probability = (task == "binomial")),
                                        rng_args))
      dfte <- data.frame(Xteg, check.names = FALSE)
      pr   <- predict(rg, dfte)
      pred <- if (task == "binomial") pr$predictions[, 2] else pr$predictions
      return(list(pred = pred, fit = rg))
    }
    stop("Unsupported learner.")
  }

  # fold-level function -------------------------------------------------------
  do_fold <- function(fold) {
    set.seed(seed + fold$fold)
    tr <- fold$train
    te <- fold$test

    Xtr <- Xall[tr, , drop = FALSE]
    Xte <- Xall[te, , drop = FALSE]
    ytr <- yall[tr]
    yte <- yall[te]

    if (task == "binomial") {
      ytr <- factor(ytr, levels = class_levels)
      yte <- factor(yte, levels = class_levels)
      if (nlevels(droplevels(ytr)) < 2) {
        warning(sprintf("Fold %s skipped: only one class in training data", fold$fold))
        empty_metrics <- setNames(rep(NA_real_, length(metrics)), metric_labels)
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
    } else {
      ytr <- as.numeric(ytr)
      yte <- as.numeric(yte)
      fold_weights <- NULL
    }

    # --- Continue preprocessing and training ---
    guard <- .guard_fit(
      X = Xtr,
      y = ytr,
      steps = if (exists("preprocess") &&
                  is.list(preprocess))
        preprocess
      else
        list(),
      task  = task
    )

    Xtrg <- guard$transform(Xtr)
    Xteg <- guard$transform(Xte)
    colnames(Xtrg) <- make.names(colnames(Xtrg))
    colnames(Xteg) <- make.names(colnames(Xteg))

    results <- list()
    for (i in seq_along(learner_names)) {
      ln <- learner_names[[i]]
      learner_obj <- learner_objs[[i]]
      # Train one learner; ensure consistent level matching
      model <- train_one_learner(learner_obj, ln, Xtrg, ytr, Xteg, yte, weights = fold_weights)

      ms <- vapply(seq_along(metrics), function(idx) {
        compute_metric(metrics[[idx]], yte, model$pred)
      }, numeric(1))
      names(ms) <- metric_labels
      results[[ln]] <- list(
        metrics = ms,
        pred = data.frame(
          id = ids[te],
          truth = yte,
          pred = model$pred,
          fold = fold$fold,
          learner = ln,
          stringsAsFactors = FALSE
        ),
        guard = guard$state,
        learner = model$fit,
        feat_names = colnames(Xtrg)
      )
    }

    results
  }


  folds <- splits@indices
  nfold <- length(folds)

  # progress bar --------------------------------------------------------------
  pb <- utils::txtProgressBar(min = 0, max = nfold, style = 3)
  pb_counter <- 0
  progress_wrap <- function(f) {
    res <- tryCatch(do_fold(f), error = function(e) {
      warning(sprintf("Fold %s failed: %s", f$fold, e$message)); NULL
    })
    pb_counter <<- pb_counter + 1
    utils::setTxtProgressBar(pb, pb_counter)
    res
  }

  # parallel or sequential execution -----------------------------------------
  if (parallel && requireNamespace("future.apply", quietly = TRUE)) {
    out <- future.apply::future_lapply(seq_along(folds), function(i) {
      fold <- folds[[i]]
      fold$fold <- i  # fold numarasını ekle
      progress_wrap(fold)
    }, future.seed = TRUE)
  } else {
    out <- lapply(seq_along(folds), function(i) {
      fold <- folds[[i]]
      fold$fold <- i  # fold numarasını ekle
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

  for (fold_idx in seq_along(out)) {
    fold_res <- out[[fold_idx]]
    if (is.null(fold_res)) next
    fold_info <- folds[[fold_idx]]
    for (ln in names(fold_res)) {
      res <- fold_res[[ln]]
      if (is.null(res)) next
      m <- res$metrics
      if (all(is.na(m))) next
      metric_row <- c(list(fold = fold_info$fold, learner = ln), as.list(m))
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
        fold = fold_info$fold,
        n_train = length(fold_info$train),
        n_test = length(fold_info$test),
        learner = ln,
        features_final = if (!is.null(res$guard) &&
                              is.list(res$guard) &&
                              !is.null(res$guard$filter) &&
                              !is.null(res$guard$filter$keep)) {
          sum(res$guard$filter$keep)
        } else NA_integer_,
        row.names = NULL
      )
    }
  }

  if (!length(met_rows)) {
    stop("No successful folds were completed. Check learner and preprocessing settings.")
  }

  met_df <- do.call(rbind, met_rows)
  audit_df <- do.call(rbind, audit_rows)

  # summarize metrics ---------------------------------------------------------
  metric_summary <- aggregate(. ~ learner, data = met_df[, -1, drop = FALSE],
                              FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                                  sd = sd(x, na.rm = TRUE)))

  # optional refit ------------------------------------------------------------
  final_model <- NULL
  final_guard <- NULL
  if (refit) {
    guard_full <- .guard_fit(Xall, yall, preprocess, task)
    Xfullg <- guard_full$transform(Xall)
    ln <- learner_names[[1]]
    learner_obj <- learner_objs[[1]]
    full_weights <- if (task == "binomial") resolve_weights(yall, class_weights) else NULL
    final_model <- train_one_learner(learner_obj, ln, Xfullg, yall, Xfullg, yall,
                                     weights = full_weights)$fit
    final_guard <- guard_full$state
  }

  # assemble LeakFit object ---------------------------------------------------
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
                  metrics_used = metric_labels,
                  class_weights = class_weights,
                  positive_class = if (task == "binomial") class_levels[2] else NULL,
                  fold_seeds = setNames(seed + seq_along(folds),
                                        paste0("fold", seq_along(folds))),
                  refit = refit,
                  final_model = final_model,
                  final_preprocess = final_guard))
}
