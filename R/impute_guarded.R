#' Leakage-safe data imputation via guarded preprocessing
#'
#' @description
#' Fits imputation parameters on the training data only, then applies the same
#' guarded transformation to the test data. This function is a thin wrapper
#' around the guarded preprocessing used by \code{fit_resample()}.
#' Output is the transformed feature matrix used by the guarded pipeline
#' (categorical variables are one-hot encoded).
#'
#' @param train data frame (training set)
#' @param test data frame (test set)
#' @param method one of "median", "knn", "missForest", or "none"
#' @param constant_value unused; retained for backward compatibility
#' @param k number of neighbors for kNN imputation (if method = "knn")
#' @param seed random seed for reproducibility
#' @param winsor logical; apply MAD-based winsorization before imputation
#' @param winsor_thresh numeric; MAD cutoff (default = 3)
#' @param parallel logical; unused (kept for compatibility)
#' @param return_outliers logical; unused (outlier flags not returned)
#' @param vars optional character vector; impute only selected variables
#' @return A LeakImpute object with imputed data and guard state.
#' @seealso [fit_resample()], [predict_guard()]
#' @examples
#' train <- data.frame(x = c(1, 2, NA, 4), y = c(NA, 1, 1, 0))
#' test <- data.frame(x = c(NA, 5), y = c(1, NA))
#' imp <- impute_guarded(train, test, method = "median", winsor = FALSE)
#' imp$train
#' imp$test
#' @export
impute_guarded <- function(train,
                           test,
                           method = c("median", "knn", "missForest", "none"),
                           constant_value = 0,
                           k = 5,
                           seed = 123,
                           winsor = TRUE,
                           winsor_thresh = 3,
                           parallel = FALSE,
                           return_outliers = FALSE,
                           vars = NULL) {
  method <- match.arg(method)
  if (!is.data.frame(train) || !is.data.frame(test)) {
    stop("train and test must be data frames.", call. = FALSE)
  }

  if (!is.null(vars)) {
    vars <- intersect(vars, names(train))
    train <- train[, vars, drop = FALSE]
    test  <- test[, vars, drop = FALSE]
  }

  if (isTRUE(parallel)) {
    warning("parallel is unused in impute_guarded; guarded preprocessing runs sequentially.",
            call. = FALSE)
  }
  if (isTRUE(return_outliers)) {
    warning("return_outliers is not supported in impute_guarded; outlier flags are not returned.",
            call. = FALSE)
  }
  if (!identical(constant_value, 0) && method != "none") {
    warning("constant_value is unused in impute_guarded and will be ignored.",
            call. = FALSE)
  }

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

  steps <- list(
    impute = list(method = method, k = k, winsor = winsor, winsor_k = winsor_thresh),
    normalize = list(method = "none"),
    filter = list(var_thresh = 0, iqr_thresh = 0),
    fs = list(method = "none")
  )

  guard <- .guard_fit(train, y = NULL, steps = steps, task = "gaussian")
  train_imp <- guard$transform(train)
  test_imp <- guard$transform(test)

  diagnostics <- data.frame(
    variable = names(train),
    missing_train = colSums(is.na(train)),
    missing_test  = colSums(is.na(test)),
    stringsAsFactors = FALSE
  )
  diagnostics$missing_train <- as.integer(diagnostics$missing_train)
  diagnostics$missing_test  <- as.integer(diagnostics$missing_test)

  structure(
    list(
      train = train_imp,
      test  = test_imp,
      model = guard$state,
      method = method,
      summary = list(
        nmiss_train = sum(is.na(train)),
        nmiss_test  = sum(is.na(test)),
        winsorized = winsor,
        winsor_thresh = winsor_thresh,
        diagnostics = diagnostics,
        guard_audit = guard$state$audit
      ),
      outliers = NULL
    ),
    class = "LeakImpute"
  )
}
