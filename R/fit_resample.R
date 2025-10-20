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
                         preprocess = list(impute = "median",
                                           normalize = "zscore",
                                           filter = list(var_thresh = 0, iqr_thresh = 0),
                                           fs = "none"),
                         learner = c("glmnet", "ranger"),
                         learner_args = list(),
                         metrics = c("auc", "pr_auc", "accuracy"),
                         parallel = FALSE,
                         refit = TRUE,
                         seed = 1) {

  set.seed(seed)
  learner <- match.arg(learner, several.ok = TRUE)
  Xall <- .bio_get_x(x)
  yall <- .bio_get_y(x, outcome)
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
      cvfit <- glmnet::cv.glmnet(Xtrg, ytr, family = fam, alpha = la$alpha)
      pred  <- as.numeric(predict(cvfit, Xteg, s = "lambda.min",
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
      return(list(metrics = NA, pred = NA, guard = NA, learner = NULL))

    guard <- .guard_fit(Xtr, ytr, preprocess, task)
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
    out <- future.apply::future_lapply(folds, progress_wrap, future.seed = TRUE)
  } else {
    out <- lapply(folds, progress_wrap)
  }
  close(pb)

  # collect results -----------------------------------------------------------
out_flat <- do.call(c, out)
out_flat <- Filter(Negate(is.null), out_flat)
learners_used <- rep(learner, each = length(splits@indices))

  met_df <- do.call(rbind, lapply(seq_along(out_flat), function(i) {
    m <- out_flat[[i]]$metrics
    data.frame(fold = ((i - 1) %% nfold) + 1,
               learner = learners_used[i],
               t(as.data.frame(m)), row.names = NULL, check.names = FALSE)
  }))
  preds <- lapply(out_flat, `[[`, "pred")
  guards <- lapply(out_flat, `[[`, "guard")
  lears <- lapply(out_flat, `[[`, "learner")
  featn <- out_flat[[1]]$feat_names

  # summarize metrics ---------------------------------------------------------
  metric_summary <- aggregate(. ~ learner, data = met_df[, -1, drop = FALSE],
                              FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                                                  sd = sd(x, na.rm = TRUE)))

  # leakage audit -------------------------------------------------------------
  audit_df <- data.frame(
    fold = seq_len(nfold),
    n_train = sapply(folds, function(f) length(f$train)),
    n_test  = sapply(folds, function(f) length(f$test)),
    learner = rep(learner, each = nfold),
    features_final = sapply(out_flat, function(o)
      if (!is.null(o$guard$filter$keep)) sum(o$guard$filter$keep) else NA)
  )

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
