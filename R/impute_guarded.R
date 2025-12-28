#' Leakage-safe data imputation (robust, enhanced version)
#'
#' @description
#' Performs simple and advanced missing value imputation methods while preventing
#' information leakage. All imputation parameters are learned *only* from training data
#' and then applied to the test data. Includes Median Absolute Deviation (MAD)–based
#' Winsorization for robust outlier handling, automatic factor level preservation,
#' and diagnostic summaries.
#' For kNN and missForest, imputation is fit on the training set and then
#' train-derived fill values are applied to test rows to avoid leakage.
#'
#' @param train data frame (training set)
#' @param test data frame (test set)
#' @param method one of "mean", "median", "mode", "constant", "knn", "mice", or "missForest"
#' @param constant_value numeric or character value if method = "constant"
#' @param k number of neighbors for kNN imputation (if method = "knn")
#' @param seed random seed for reproducibility
#' @param winsor logical; apply MAD-based Winsorization before imputation
#' @param winsor_thresh numeric; MAD cutoff (default = 3)
#' @param parallel logical; enable parallel imputation for heavy models (default = FALSE)
#' @param return_outliers logical; return logical matrix of Winsorized outlier flags
#' @param vars optional character vector; impute only selected variables
#' @return A LeakImpute object with imputed data, model, and diagnostic summary.
#' @seealso
#' [VIM::kNN()], [mice::mice()], [missForest::missForest()]
#' @references
#' Rousseeuw & Croux (1993), *Alternative to the Median Absolute Deviation.*
#' @examples
#' train <- data.frame(x = c(1, 2, NA, 4), y = c(NA, 1, 1, 0))
#' test <- data.frame(x = c(NA, 5), y = c(1, NA))
#' imp <- impute_guarded(train, test, method = "median", winsor = FALSE)
#' imp$train
#' imp$test
#' @export
impute_guarded <- function(train,
                           test,
                           method = c("median", "mean", "mode", "constant", "knn", "mice", "missForest"),
                           constant_value = 0,
                           k = 5,
                           seed = 123,
                           winsor = TRUE,
                           winsor_thresh = 3,
                           parallel = FALSE,
                           return_outliers = FALSE,
                           vars = NULL) {
  method <- match.arg(method)
  has_old_seed <- exists(".Random.seed", envir = .GlobalEnv)
  old_seed <- if (has_old_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (exists("old_seed", inherits = FALSE)) {
      if (!is.null(old_seed)) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (!has_old_seed && exists(".Random.seed", envir = .GlobalEnv)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }
  }, add = TRUE)
  set.seed(seed)

  if (!is.data.frame(train) || !is.data.frame(test))
    stop("train and test must be data frames.")

  # Subset variables if requested
  if (!is.null(vars)) {
    vars <- intersect(vars, names(train))
    train <- train[, vars, drop = FALSE]
    test  <- test[, vars, drop = FALSE]
  }

  # Ensure consistent types
  train <- type.convert(train, as.is = TRUE)
  test  <- type.convert(test, as.is = TRUE)
  factor_levels <- lapply(train, function(col) if (is.factor(col)) levels(col) else NULL)

  nmiss_train <- sum(is.na(train))
  nmiss_test  <- sum(is.na(test))
  message(sprintf("Guarded imputation using method = '%s' (%d missing in train, %d in test)",
                  method, nmiss_train, nmiss_test))

  # ---------- Robust Winsorization ----------
  winsorize_mad <- function(x, thresh = 3) {
    if (!is.numeric(x)) return(x)
    med <- median(x, na.rm = TRUE)
    mad_val <- mad(x, constant = 1.4826, na.rm = TRUE)
    if (is.na(med) || mad_val == 0) return(x)
    lower <- med - thresh * mad_val
    upper <- med + thresh * mad_val
    x[x < lower] <- lower
    x[x > upper] <- upper
    x
  }

  outlier_flags <- NULL
  if (winsor) {
    outlier_flags <- as.data.frame(lapply(train, function(x) {
      if (!is.numeric(x)) return(rep(FALSE, length(x)))
      med <- median(x, na.rm = TRUE)
      mad_val <- mad(x, constant = 1.4826, na.rm = TRUE)
      if (mad_val == 0 || is.na(med)) return(rep(FALSE, length(x)))
      (x < med - winsor_thresh * mad_val) | (x > med + winsor_thresh * mad_val)
    }))
    train <- as.data.frame(lapply(train, winsorize_mad, thresh = winsor_thresh))
    message(sprintf("Applied MAD-based Winsorization (±%.1f MADs) to numeric variables.", winsor_thresh))
  }

  train_impute_values <- function(train_imp) {
    vals <- lapply(names(train_imp), function(col) {
      x <- train_imp[[col]]
      if (is.numeric(x)) {
        if (all(is.na(x))) return(NA_real_)
        return(stats::median(x, na.rm = TRUE))
      }
      non_missing <- x[!is.na(x)]
      if (!length(non_missing)) return(NA)
      tab <- sort(table(non_missing), decreasing = TRUE)
      val <- names(tab)[1L]
      if (is.factor(x)) {
        return(factor(val, levels = levels(x)))
      }
      if (is.logical(x)) {
        return(as.logical(val))
      }
      if (is.character(x)) {
        return(val)
      }
      type.convert(val, as.is = TRUE)
    })
    names(vals) <- names(train_imp)
    vals
  }

  apply_impute_values <- function(df, impute_values) {
    for (col in names(impute_values)) {
      if (!col %in% names(df)) next
      val <- impute_values[[col]]
      if (length(val) == 1 && all(is.na(val))) next
      idx <- is.na(df[[col]])
      if (any(idx)) df[[col]][idx] <- val
    }
    df
  }

  # ---------- Imputation ----------
  impute_values <- list()
  train_imp <- test_imp <- model <- NULL

  # Simple deterministic methods
  if (method %in% c("mean", "median", "mode", "constant")) {
    impute_values <- lapply(names(train), function(col) {
      x <- train[[col]]
      if (!anyNA(x)) return(NA)
      if (is.numeric(x)) {
        switch(method,
               mean   = mean(x, na.rm = TRUE),
               median = median(x, na.rm = TRUE),
               mode   = {
                 non_missing <- x[!is.na(x)]
                 tab <- sort(table(non_missing), decreasing = TRUE)
                 type.convert(names(tab)[1L], as.is = TRUE)
               },
               constant = constant_value)
      } else {
        if (method == "constant") {
          if (is.factor(x)) {
            return(factor(constant_value, levels = levels(x)))
          }
          return(constant_value)
        }

        # For non-numeric columns, 'mean' and 'median' fall back to the mode
        non_missing <- x[!is.na(x)]
        if (length(non_missing) == 0) {
          if (is.factor(x)) return(factor(NA, levels = levels(x)))
          return(NA)
        }
        tab <- sort(table(non_missing), decreasing = TRUE)
        if (is.factor(x)) {
          factor(names(tab)[1L], levels = levels(x))
        } else {
          type.convert(names(tab)[1L], as.is = TRUE)
        }
      }
    })
    names(impute_values) <- names(train)

    train_imp <- train
    test_imp <- test

    for (col in names(train)) {
      val <- impute_values[[col]]
      if (length(val) == 1 && all(is.na(val))) next
      if (anyNA(train_imp[[col]])) {
        train_imp[[col]][is.na(train_imp[[col]])] <- val
      }
      if (col %in% names(test_imp) && anyNA(test_imp[[col]])) {
        test_imp[[col]][is.na(test_imp[[col]])] <- val
      }
    }
    model <- impute_values
  }

  # kNN
  else if (method == "knn") {
    if (!requireNamespace("VIM", quietly = TRUE))
      stop("Install 'VIM' for kNN imputation.")
    args <- list(data = train, k = k, imp_var = FALSE)
    train_imp <- suppressWarnings(do.call(VIM::kNN, args))
    impute_values <- train_impute_values(train_imp)
    model <- list(method = "knn", k = k, impute_values = impute_values)
    test_imp <- apply_impute_values(test, impute_values)
  }

  # MICE
  else if (method == "mice") {
    if (!requireNamespace("mice", quietly = TRUE))
      stop("Install 'mice' for MICE imputation.")
    imp_model <- mice::mice(train, m = 1, maxit = 5, printFlag = FALSE, seed = seed)
    train_imp <- mice::complete(imp_model)
    test_imp <- test
    for (col in names(test)) {
      if (anyNA(test[[col]]) && col %in% names(train_imp)) {
        mu <- mean(train_imp[[col]], na.rm = TRUE)
        test_imp[[col]][is.na(test_imp[[col]])] <- mu
      }
    }
    model <- imp_model
  }

  # missForest
  else if (method == "missForest") {
    if (!requireNamespace("missForest", quietly = TRUE))
      stop("Install 'missForest' for robust nonparametric imputation.")
    mf_args <- list(xmis = train, verbose = FALSE)
    if (parallel) mf_args$parallelize <- "forests"
    imp_model <- do.call(missForest::missForest, mf_args)
    train_imp <- imp_model$ximp
    impute_values <- train_impute_values(train_imp)
    model <- imp_model
    model$impute_values <- impute_values
    test_imp <- apply_impute_values(test, impute_values)
  }

  if (!is.null(test_imp)) {
    test_imp[is.na(test_imp)] <- NA
  }

  # ---------- Postprocessing ----------
  # Restore factor levels
  for (col in names(test_imp)) {
    if (!is.null(factor_levels[[col]]))
      test_imp[[col]] <- factor(test_imp[[col]], levels = factor_levels[[col]])
  }

  diagnostics <- data.frame(
    variable = names(train),
    missing_train = colSums(is.na(train)),
    missing_test  = colSums(is.na(test)),
    post_missing_train = colSums(is.na(train_imp)),
    post_missing_test  = colSums(is.na(test_imp)),
    stringsAsFactors = FALSE
  )

  diagnostics$missing_train <- as.integer(diagnostics$missing_train)
  diagnostics$missing_test  <- as.integer(diagnostics$missing_test)

  outlier_count <- if (!is.null(outlier_flags)) colSums(outlier_flags) else rep(0L, ncol(train))
  diagnostics$outliers_train <- outlier_count

  # ---------- Assemble result ----------
  structure(
    list(
      train = train_imp,
      test  = test_imp,
      model = model,
      method = method,
      summary = list(
        nmiss_train = nmiss_train,
        nmiss_test  = nmiss_test,
        winsorized = winsor,
        winsor_thresh = winsor_thresh,
        imputed_vars = names(which(colSums(is.na(train)) > 0)),
        diagnostics = diagnostics,
        model_info = if (method == "missForest") model$OOBerror else NULL
      ),
      outliers = if (return_outliers) outlier_flags else NULL
    ),
    class = "LeakImpute"
  )
}
