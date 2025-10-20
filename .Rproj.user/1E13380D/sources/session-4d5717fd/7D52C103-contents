#' Fit and evaluate with leakage guards over predefined splits
#'
#' Performs cross-validated model training and evaluation using
#' leakage-protected preprocessing (.guard_fit) and user-specified learners.
#'
#' @param x SummarizedExperiment or matrix/data.frame
#' @param outcome outcome column (if x is SE)
#' @param splits LeakSplits object from make_splits()
#' @param preprocess list(impute, normalize, filter=list(...), fs)
#' @param learner character vector, e.g. "glmnet","ranger"
#' @param learner_args list of additional arguments passed to learner(s)
#' @param metrics named list of metric functions or vector of metric names
#' @param parallel logical, use future.apply for multicore execution
#' @param refit logical, if TRUE retrain final model on full data
#' @param seed integer, for reproducibility
#' @return LeakFit object
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
                         metrics = c("auc", "pr_auc", "accuracy"),
                         parallel = FALSE,
                         refit = TRUE,
                         seed = 1) {

  set.seed(seed)
  learner <- match.arg(learner, several.ok = TRUE)
  Xall <- .bio_get_x(x)
  yall <- as.numeric(unlist(x[[outcome]]))
  Xall <- Xall[, setdiff(colnames(Xall), outcome), drop = FALSE]
  ids  <- seq_len(nrow(Xall))
  task <- if (.bio_is_binomial(yall)) "binomial"
  else if (is.numeric(yall)) "gaussian"
  else "binomial"

  # helper: safe metric computation -------------------------------------------
  compute_metric <- function(name, y, pred) {
    if (is.function(name)) return(name(y, pred))
    if (name == "auc" && task == "binomial") {
      if (requireNamespace("pROC", quietly = TRUE))
        return(as.numeric(pROC::auc(pROC::roc(y, pred, quiet = TRUE))))
      return(mean(outer(pred[y==1], pred[y==0], ">")))
    }
    if (name == "pr_auc" && task == "binomial") {
      if (requireNamespace("PRROC", quietly = TRUE)) {
        yb <- if (is.factor(y)) as.numeric(y) - 1 else y
        pr <- PRROC::pr.curve(scores.class0 = pred[yb==1],
                              scores.class1 = pred[yb==0],
                              curve = FALSE)
        return(pr$auc.integral)
      }
      return(NA_real_)
    }
    if (name == "accuracy" && task == "binomial") {
      yb <- if (is.factor(y)) as.numeric(y) - 1 else y
      return(mean((pred >= 0.5) == as.logical(yb)))
    }
    if (name == "rmse" && task == "gaussian")
      return(sqrt(mean((as.numeric(y) - pred)^2)))
    if (name == "cindex" && requireNamespace("survival", quietly = TRUE))
      return(survival::concordance.index(pred, y)$c.index)
    NA_real_
  }

  # single learner wrapper ----------------------------------------------------
  train_one_learner <- function(learner_name, Xtrg, ytr, Xteg, yte) {
    if (learner_name == "glmnet") {
      if (!requireNamespace("glmnet", quietly = TRUE)) stop("Install 'glmnet'.")
      fam <- if (task == "binomial") "binomial" else "gaussian"
      la  <- modifyList(list(alpha = 0.9, standardize = FALSE), learner_args)
      cvfit <- glmnet::cv.glmnet(as.matrix(Xtrg), ytr, family = fam, alpha = la$alpha)
      pred  <- as.numeric(predict(cvfit, as.matrix(Xteg), s = "lambda.min",
                                  type = if (task == "binomial") "response" else "link"))

      return(list(pred = pred, fit = cvfit))
    }
    if (learner_name == "ranger") {
      if (!requireNamespace("ranger", quietly = TRUE)) stop("Install 'ranger'.")
      dftr <- data.frame(y = ytr, Xtrg, check.names = FALSE)
      frm  <- stats::as.formula("y ~ .")
      rg   <- do.call(ranger::ranger, c(list(formula = frm, data = dftr,
                                             probability = (task == "binomial")),
                                        learner_args))
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
    tr <- fold$train; te <- fold$test
    Xtr <- Xall[tr, , drop = FALSE]; ytr <- yall[tr]
    Xte <- Xall[te, , drop = FALSE]; yte <- yall[te]

    if (length(unique(ytr)) < 2)
      return(list(metrics = as.list(setNames(rep(NA_real_, length(metrics)), metrics)),
                  pred = data.frame(id = integer(0), truth = numeric(0), pred = numeric(0)),
                  guard = list(state = NULL),
                  learner = NULL,
                  feat_names = colnames(Xtr)))

    ytr <- as.numeric(as.factor(ytr))

    if (length(unique(ytr)) < 2 || anyNA(ytr)) {
      return(list(
        metrics = as.list(setNames(rep(NA_real_, length(metrics)), metrics)),
        pred = data.frame(id = integer(0), truth = numeric(0), pred = numeric(0)),
        guard = list(state = NULL),
        learner = NULL,
        feat_names = colnames(Xtr)
      ))
    }


    guard <- .guard_fit(
      Xtr, ytr,
      steps = if (exists("preprocess") && is.list(preprocess)) preprocess else list(),
      task  = task
    )

    Xtrg  <- guard$transform(Xtr)
    Xteg  <- guard$transform(Xte)
    colnames(Xtrg) <- make.names(colnames(Xtrg))
    colnames(Xteg) <- make.names(colnames(Xteg))

    results <- list()
    for (ln in learner) {
      model <- train_one_learner(ln, Xtrg, ytr, Xteg, yte)
      ms <- vapply(metrics, function(m) compute_metric(m, yte, model$pred), numeric(1))
      results[[ln]] <- list(metrics = ms,
                            pred = data.frame(id = ids[te], truth = yte, pred = model$pred),
                            guard = guard$state,
                            learner = model$fit,
                            feat_names = colnames(Xtrg))
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
      met_rows[[length(met_rows) + 1]] <- data.frame(
        fold = fold_info$fold,
        learner = ln,
        t(as.data.frame(m)),
        row.names = NULL,
        check.names = FALSE
      )
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
  if (refit) {
    guard_full <- .guard_fit(Xall, yall, preprocess, task)
    Xfullg <- guard_full$transform(Xall)
    ln <- learner[1]
    final_model <- train_one_learner(ln, Xfullg, yall, Xfullg, yall)$fit
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
                  metrics_requested = metrics,
                  refit = refit,
                  final_model = final_model))
}
