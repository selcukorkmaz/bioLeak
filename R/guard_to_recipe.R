# Guard-to-recipe bridge -----------------------------------------------------

#' Convert guard preprocessing steps to a recipes recipe
#'
#' Maps bioLeak guard preprocessing steps (impute, normalize, filter, fs) to
#' their closest \pkg{recipes} equivalents. Requires the \pkg{recipes} package.
#' Steps that have no direct recipe equivalent are skipped with a warning.
#'
#' @param steps A named list of guard preprocessing steps, e.g.,
#'   \code{list(impute = list(method = "median"), normalize = list(method = "zscore"))}.
#' @param formula A model formula (e.g., \code{outcome ~ .}).
#' @param training_data A data.frame used to initialize the recipe.
#' @return A \code{recipes::recipe} object with the mapped steps added.
#' @details
#' Mapping:
#' \itemize{
#'   \item \code{impute$method = "median"}: \code{step_impute_median(all_numeric_predictors())}
#'   \item \code{impute$method = "knn"}: \code{step_impute_knn(all_predictors(), neighbors = k)}
#'   \item \code{impute$method = "missForest"} or \code{"mice"}: Warning + \code{step_impute_median()} fallback
#'   \item \code{normalize$method = "zscore"}: \code{step_normalize(all_numeric_predictors())}
#'   \item \code{normalize$method = "robust"}: Warning + \code{step_normalize()} fallback
#'   \item \code{normalize$method = "none"}: No step added
#'   \item \code{filter$var_thresh > 0}: \code{step_nzv(all_numeric_predictors())}
#'   \item \code{fs$method = "pca"}: \code{step_pca(all_numeric_predictors(), num_comp = ncomp)}
#'   \item \code{fs$method = "ttest"} or \code{"lasso"}: Warning, skipped (no recipe equivalent)
#' }
#' @export
guard_to_recipe <- function(steps, formula, training_data) {
  if (!requireNamespace("recipes", quietly = TRUE)) {
    stop("The 'recipes' package is required for guard_to_recipe(). ",
         "Install it with install.packages('recipes').", call. = FALSE)
  }

  stopifnot(is.list(steps))
  stopifnot(inherits(formula, "formula"))
  stopifnot(is.data.frame(training_data))

  rec <- recipes::recipe(formula, data = training_data)

  # Imputation
  imp <- steps$impute
  if (!is.null(imp) && !is.null(imp$method)) {
    if (identical(imp$method, "median")) {
      rec <- recipes::step_impute_median(rec, recipes::all_numeric_predictors())
    } else if (identical(imp$method, "knn")) {
      k <- imp$k %||% imp$neighbors %||% 5L
      rec <- recipes::step_impute_knn(rec, recipes::all_predictors(), neighbors = k)
    } else if (imp$method %in% c("missForest", "mice")) {
      warning(sprintf("Imputation method '%s' has no direct recipe equivalent; falling back to step_impute_median().",
                      imp$method), call. = FALSE)
      rec <- recipes::step_impute_median(rec, recipes::all_numeric_predictors())
    }
  }

  # Normalization
  norm <- steps$normalize
  if (!is.null(norm) && !is.null(norm$method)) {
    if (identical(norm$method, "zscore")) {
      rec <- recipes::step_normalize(rec, recipes::all_numeric_predictors())
    } else if (identical(norm$method, "robust")) {
      warning("Normalization method 'robust' has no direct recipe equivalent; falling back to step_normalize().",
              call. = FALSE)
      rec <- recipes::step_normalize(rec, recipes::all_numeric_predictors())
    } else if (!identical(norm$method, "none")) {
      # Unknown normalization method — skip with a note
      warning(sprintf("Normalization method '%s' not recognized; skipping.", norm$method), call. = FALSE)
    }
  }

  # Variance filter
  filt <- steps$filter
  if (!is.null(filt)) {
    var_thresh <- filt$var_thresh %||% 0
    if (var_thresh > 0) {
      rec <- recipes::step_nzv(rec, recipes::all_numeric_predictors())
    }
  }

  # Feature selection
  fs <- steps$fs
  if (!is.null(fs) && !is.null(fs$method)) {
    if (identical(fs$method, "pca")) {
      ncomp <- fs$ncomp %||% fs$num_comp %||% 5L
      rec <- recipes::step_pca(rec, recipes::all_numeric_predictors(), num_comp = ncomp)
    } else if (fs$method %in% c("ttest", "lasso")) {
      warning(sprintf("Feature selection method '%s' has no direct recipe equivalent; skipping.",
                      fs$method), call. = FALSE)
    }
  }

  rec
}
