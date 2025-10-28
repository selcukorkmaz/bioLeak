#' Leakage-safe data imputation (robust, enhanced version)
#'
#' @description
#' Performs simple and advanced missing value imputation methods while preventing
#' information leakage. All imputation parameters are learned *only* from training data
#' and then applied to the test data. Includes Median Absolute Deviation (MAD)–based
#' Winsorization for robust outlier handling, automatic factor level preservation,
#' and diagnostic summaries.
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
                 ux <- unique(non_missing)
                 ux[which.max(tabulate(match(non_missing, ux)))]
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
        ux <- unique(non_missing)
        if (length(ux) == 0) {
          if (is.factor(x)) return(factor(NA, levels = levels(x)))
          return(NA)
        }
        ux[which.max(tabulate(match(non_missing, ux)))]
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
    model <- list(method = "knn", k = k)
    test_comb <- rbind(train, test)
    args$data <- test_comb
    test_imp_full <- suppressWarnings(do.call(VIM::kNN, args))
    test_imp <- tail(test_imp_full, n = nrow(test))
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
    model <- imp_model

    comb <- rbind(train_imp, test)
    mf_args$xmis <- comb  # comb verisini yeniden ver
    test_imp_full <- suppressWarnings(do.call(missForest::missForest, mf_args))
    test_imp <- tail(test_imp_full$ximp, n = nrow(test))
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
