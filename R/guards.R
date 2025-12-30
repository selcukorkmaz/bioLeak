# Utilities ---------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

.guard_is_binary <- function(y) {
  if (is.null(y)) return(FALSE)
  if (is.factor(y)) return(nlevels(y) == 2)
  ux <- unique(y[!is.na(y)])
  length(ux) == 2
}

.guard_majority <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

.guard_mad_winsorize <- function(x, k = 5, center = NULL, scale = NULL) {
  if (!is.numeric(x)) return(x)
  if (all(!is.finite(x))) return(x)
  m <- if (is.null(center)) stats::median(x, na.rm = TRUE) else center
  s <- if (is.null(scale)) stats::mad(x, na.rm = TRUE, constant = 1) else scale
  if (!is.finite(m) || !is.finite(s)) return(x)
  lo <- m - k * s; hi <- m + k * s
  x[x < lo] <- lo
  x[x > hi] <- hi
  x
}

.guard_align_cols <- function(X, ref_cols) {
  X <- as.data.frame(X, check.names = FALSE)
  missing_cols <- setdiff(ref_cols, names(X))
  if (length(missing_cols)) {
    for (nm in missing_cols) X[[nm]] <- NA_real_
  }
  extra_cols <- setdiff(names(X), ref_cols)
  if (length(extra_cols)) {
    X <- X[, setdiff(names(X), extra_cols), drop = FALSE]
  }
  X[, ref_cols, drop = FALSE]
}

.guard_fill_na <- function(X, fill) {
  X <- as.data.frame(X, check.names = FALSE)
  for (nm in names(fill)) {
    if (!nm %in% names(X)) next
    idx <- is.na(X[[nm]])
    if (any(idx)) X[[nm]][idx] <- fill[[nm]]
  }
  X
}

.guard_knn_search <- function(train, query, k) {
  train <- as.matrix(train)
  query <- as.matrix(query)
  if (!nrow(train) || !nrow(query) || k < 1L) {
    return(matrix(integer(0), nrow = nrow(query), ncol = 0))
  }
  k <- min(k, nrow(train))
  if (requireNamespace("RANN", quietly = TRUE)) {
    nn <- RANN::nn2(train, query, k = k)
    return(nn$nn.idx)
  }
  if (requireNamespace("FNN", quietly = TRUE)) {
    nn <- FNN::get.knnx(train, query, k = k)
    return(nn$nn.index)
  }
  idx <- matrix(NA_integer_, nrow(query), k)
  for (i in seq_len(nrow(query))) {
    diff <- sweep(train, 2, query[i, ], "-")
    d <- rowSums(diff * diff)
    ord <- order(d, na.last = NA)
    idx[i, ] <- ord[seq_len(k)]
  }
  idx
}

.guard_knn_impute_train <- function(X, k = 5) {
  X <- as.data.frame(X, check.names = FALSE)
  med <- vapply(X, function(col) {
    m <- stats::median(col, na.rm = TRUE)
    if (!is.finite(m)) 0 else m
  }, numeric(1))
  X_fill <- .guard_fill_na(X, med)
  scale <- vapply(X_fill, stats::sd, numeric(1))
  scale[!is.finite(scale) | scale == 0] <- 1
  X_scaled <- sweep(as.matrix(X_fill), 2, scale, "/")

  n <- nrow(X_fill)
  if (n < 2L) {
    return(list(imp = X_fill, med = med, scale = scale))
  }

  idx <- .guard_knn_search(X_scaled, X_scaled, k = min(k + 1L, n))
  X_imp <- X_fill
  for (i in seq_len(n)) {
    miss_cols <- which(is.na(X[i, ]))
    if (!length(miss_cols)) next
    neigh <- idx[i, ]
    neigh <- neigh[neigh != i]
    if (!length(neigh)) next
    if (length(neigh) > k) neigh <- neigh[seq_len(k)]
    for (j in miss_cols) {
      val <- mean(X_fill[neigh, j], na.rm = TRUE)
      if (!is.finite(val)) val <- med[j]
      X_imp[i, j] <- val
    }
  }
  list(imp = X_imp, med = med, scale = scale)
}

.guard_knn_impute_new <- function(Xnew, ref, k, med, scale) {
  Xnew <- .guard_align_cols(Xnew, names(med))
  Xnew_raw <- Xnew
  Xnew_fill <- .guard_fill_na(Xnew, med)
  if (!nrow(ref)) return(Xnew_fill)

  ref_scaled <- sweep(as.matrix(ref), 2, scale, "/")
  new_scaled <- sweep(as.matrix(Xnew_fill), 2, scale, "/")
  idx <- .guard_knn_search(ref_scaled, new_scaled, k = min(k, nrow(ref)))

  X_imp <- Xnew_fill
  for (i in seq_len(nrow(Xnew_fill))) {
    miss_cols <- which(is.na(Xnew_raw[i, ]))
    if (!length(miss_cols)) next
    neigh <- idx[i, ]
    if (!length(neigh)) next
    for (j in miss_cols) {
      val <- mean(ref[neigh, j], na.rm = TRUE)
      if (!is.finite(val)) val <- med[j]
      X_imp[i, j] <- val
    }
  }
  X_imp
}

.guard_rf_impute_train <- function(X, maxiter = 1, ntree = 100, seed = NULL) {
  X <- as.data.frame(X, check.names = FALSE)
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package 'randomForest' is required for impute$method='missForest'.", call. = FALSE)
  }
  med <- vapply(X, function(col) {
    m <- stats::median(col, na.rm = TRUE)
    if (!is.finite(m)) 0 else m
  }, numeric(1))
  X_fill <- .guard_fill_na(X, med)
  miss_cols <- names(which(vapply(X, function(col) any(is.na(col)), logical(1))))
  models <- list()

  if (!length(miss_cols)) {
    return(list(imp = X_fill, med = med, models = models))
  }

  for (iter in seq_len(maxiter)) {
    for (col in miss_cols) {
      obs <- !is.na(X[[col]])
      if (sum(obs) < 2L) next
      pred_cols <- setdiff(names(X_fill), col)
      if (!length(pred_cols)) next
      df_obs <- data.frame(y = X_fill[obs, col],
                           X_fill[obs, pred_cols, drop = FALSE],
                           check.names = FALSE)
      if (!is.null(seed)) set.seed(seed + iter)
      rf <- randomForest::randomForest(y ~ ., data = df_obs, ntree = ntree)
      if (iter == maxiter) {
        models[[col]] <- list(model = rf, predictors = pred_cols)
      }
      if (any(!obs)) {
        newdata <- X_fill[!obs, pred_cols, drop = FALSE]
        X_fill[!obs, col] <- stats::predict(rf, newdata = newdata)
      }
    }
  }

  med_out <- vapply(X_fill, function(col) {
    m <- stats::median(col, na.rm = TRUE)
    if (!is.finite(m)) 0 else m
  }, numeric(1))
  list(imp = X_fill, med = med_out, models = models)
}

.guard_rf_impute_new <- function(Xnew, models, med) {
  Xnew <- .guard_align_cols(Xnew, names(med))
  Xnew_raw <- Xnew
  Xnew_fill <- .guard_fill_na(Xnew, med)

  if (length(models)) {
    for (col in names(models)) {
      if (!col %in% names(Xnew_fill)) next
      idx <- is.na(Xnew_raw[[col]])
      if (!any(idx)) next
      pred_cols <- models[[col]]$predictors
      missing_pred <- setdiff(pred_cols, names(Xnew_fill))
      if (length(missing_pred)) {
        for (nm in missing_pred) Xnew_fill[[nm]] <- med[[nm]]
      }
      newdata <- Xnew_fill[idx, pred_cols, drop = FALSE]
      Xnew_fill[idx, col] <- stats::predict(models[[col]]$model, newdata = newdata)
    }
  }

  Xnew_fill <- .guard_fill_na(Xnew_fill, med)
  Xnew_fill
}

.guard_hash <- function(obj) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(obj))
  }
  paste0("h", sprintf("%08X", as.integer(stats::runif(1, 0, .Machine$integer.max))))
}

#' Ensure consistent categorical levels for guarded preprocessing
#'
#' Converts character/logical columns to factors and aligns factor levels with
#' a training-time \code{levels_map}. Adds a dummy level when a column has only
#' one observed level so that downstream one-hot encoding retains a column.
#'
#' @param df data.frame to normalize factor levels.
#' @param levels_map optional named list of factor levels learned from training data.
#' @param dummy_prefix prefix used when adding a dummy level to single-level factors.
#' @return List with elements \code{data} (data.frame) and \code{levels} (named list of levels).
#' @keywords internal
#' @export
.guard_ensure_levels <- function(df, levels_map = NULL, dummy_prefix = "__dummy__") {
  stopifnot(is.data.frame(df))
  if (!is.null(levels_map)) {
    missing_cols <- setdiff(names(levels_map), names(df))
    if (length(missing_cols)) {
      new_cols_df <- as.data.frame(lapply(levels_map[missing_cols], function(lvls) {
        factor(rep(NA_character_, nrow(df)), levels = lvls)
      }))
      df <- cbind(df, new_cols_df)
    }
  }

  cat_cols <- names(df)[!vapply(df, is.numeric, logical(1))]
  if (!length(cat_cols) && is.null(levels_map)) {
    return(list(data = df, levels = list()))
  }

  if (is.null(levels_map)) levels_map <- list()

  for (nm in cat_cols) {
    col <- df[[nm]]
    if (is.logical(col)) {
      col <- factor(col, levels = c(FALSE, TRUE))
    } else if (!is.factor(col)) {
      col <- factor(col)
    }

    lvls <- levels_map[[nm]]
    if (is.null(lvls)) {
      lvls <- levels(col)
      if (length(lvls) <= 1) {
        dummy <- paste0(dummy_prefix, nm)
        while (dummy %in% lvls) dummy <- paste0(dummy, "_")
        lvls <- unique(c(lvls, dummy))
      }
    }

    col <- factor(col, levels = lvls)
    df[[nm]] <- col
    levels_map[[nm]] <- lvls
  }

  list(data = df, levels = levels_map)
}

#' @title Fit leakage-safe preprocessing pipeline
#' @description Builds and fits a guarded preprocessing pipeline on training data,
#' then constructs a transformer for consistent application to new data.
#' @details
#' The pipeline applies, in order:
#' \itemize{
#'   \item Winsorization (optional) to limit outliers.
#'   \item Imputation learned on training data only.
#'   \item Normalization (z-score or robust).
#'   \item Variance/IQR filtering.
#'   \item Feature selection (optional; t-test, lasso, PCA).
#' }
#' All statistics are estimated on the training data and re-used for new data.
#' @param X matrix/data.frame of predictors (training).
#' @param y Optional outcome for supervised feature selection.
#' @param steps List of configuration options (see Details).
#' @param task "binomial" or "gaussian".
#' @return An object of class "GuardFit" with elements `transform`, `state`, `p_out`, and `steps`.
#' @seealso [predict_guard()]
#' @examples
#' x <- data.frame(a = c(1, 2, NA), b = c(3, 4, 5))
#' fit <- .guard_fit(x, y = c(1, 2, 3),
#'                   steps = list(impute = list(method = "median")),
#'                   task = "gaussian")
#' fit$transform(x)
#' @export
.guard_fit <- function(X, y = NULL, steps = list(), task = c("binomial","gaussian")) {
  task <- match.arg(task)
  # Coerce to data.frame for flexible handling
  if (is.matrix(X)) X <- as.data.frame(X, check.names = FALSE)
  stopifnot(is.data.frame(X))
  # --- Karakter tipli sütunları sayısala çevir (eğer sayısal görünümlüyse) ---
  for (nm in names(X)) {
    if (is.character(X[[nm]])) {
      nz <- X[[nm]][!is.na(X[[nm]])]
      if (length(nz) && all(grepl("^\\s*-?\\d+(\\.\\d+)?\\s*$", nz))) {
        X[[nm]] <- as.numeric(X[[nm]])
      } else {
        X[[nm]] <- as.factor(X[[nm]])
      }
    }
  }

  n0 <- nrow(X)

  # Keep original column names for diagnostics
  orig_colnames <- colnames(X)

  # 0) Mixed-type handling: one-hot encode non-numeric columns -----------------
  if (!ncol(X)) stop("Input X has no columns after preprocessing.")
  all_numeric <- all(vapply(X, is.numeric, logical(1)))
  if (all_numeric) {
    encoding_levels <- list()
    X <- as.data.frame(X, check.names = FALSE)
    p_after_encode <- ncol(X)
  } else {
    prep <- .guard_ensure_levels(X)
    X <- prep$data
    if (!ncol(X)) stop("Input X has no columns after preprocessing.")
    encoding_levels <- prep$levels
    mf <- stats::model.frame(~ . , data = as.data.frame(X), na.action = stats::na.pass)
    mm <- stats::model.matrix(~ . - 1, data = mf)  # satır sayısı korunur
    X <- as.data.frame(mm, check.names = FALSE)
    p_after_encode <- ncol(X)
  }

  # Replace non-finite with NA to unify missingness treatment
  for (j in seq_len(ncol(X))) {
    X[[j]][!is.finite(X[[j]])] <- NA
  }

  state <- list()
  state$encoding <- list(levels = encoding_levels, has_factors = !all_numeric)
  audit <- list()

  # 1) Robust Winsorization (train-only) --------------------------------------
  impute_cfg <- steps$impute %||% list()
  winsor_on <- impute_cfg$winsor %||% TRUE
  winsor_k  <- impute_cfg$winsor_k %||% 5

  if (isTRUE(winsor_on)) {
    winsor_med <- vapply(X, function(col) {
      if (is.numeric(col)) stats::median(col, na.rm = TRUE) else NA_real_
    }, numeric(1))
    winsor_mad <- vapply(X, function(col) {
      if (is.numeric(col)) stats::mad(col, na.rm = TRUE, constant = 1) else NA_real_
    }, numeric(1))
    X <- as.data.frame(Map(function(col, med, mad) {
      .guard_mad_winsorize(col, k = winsor_k, center = med, scale = mad)
    }, X, winsor_med, winsor_mad),
                       check.names = FALSE)
    state$winsor <- list(enabled = TRUE, k = winsor_k, med = winsor_med, mad = winsor_mad)
    audit <- c(audit, list(paste0("winsor: k=", winsor_k)))
  } else {
    state$winsor <- list(enabled = FALSE)
  }

  # 2) Imputation (train-only fitting) ----------------------------------------
  imp_method <- impute_cfg$method %||% "median"

  has_na <- any(vapply(X, function(col) any(is.na(col)), logical(1)))


  if (!has_na) {
    state$impute <- list(method = "none")
    state$impute$label <- "impute: none (no NA)"
  } else if (imp_method == "median") {
    meds <- vapply(X, function(col) stats::median(col, na.rm = TRUE), numeric(1))
    for (j in seq_len(ncol(X))) {
      idx <- is.na(X[[j]])
      if (any(idx)) X[[j]][idx] <- meds[j]
    }
    state$impute <- list(method = "median", med = meds)

  } else if (imp_method == "knn") {
    k <- impute_cfg$k %||% 5
    knn_fit <- .guard_knn_impute_train(X, k = k)
    X <- knn_fit$imp
    state$impute <- list(
      method = "knn",
      k = k,
      med = knn_fit$med,
      knn_ref = knn_fit$imp,
      knn_scale = knn_fit$scale
    )

  } else if (imp_method == "mice") {
    stop("impute$method='mice' is not supported for guarded prediction; use 'knn', 'missForest', or 'median'.",
         call. = FALSE)

  } else if (imp_method == "missForest") {
    rf_maxiter <- impute_cfg$maxiter %||% 10
    rf_ntree <- impute_cfg$ntree %||% 100
    rf_fit <- .guard_rf_impute_train(X, maxiter = rf_maxiter, ntree = rf_ntree)
    X <- rf_fit$imp
    state$impute <- list(
      method = "missForest",
      maxiter = rf_maxiter,
      ntree = rf_ntree,
      med = rf_fit$med,
      rf_models = rf_fit$models
    )
  } else if (imp_method == "none") {
    orig_cols <- colnames(X)
    med <- vapply(orig_cols, function(nm) {
      col <- X[[nm]]
      if (all(is.na(col))) return(0)
      stats::median(col, na.rm = TRUE)
    }, numeric(1))

    indicator_map <- list()
    fallback_used <- FALSE
    for (nm in orig_cols) {
      col <- X[[nm]]
      miss_idx <- is.na(col)
      if (!any(miss_idx)) next
      fallback_used <- TRUE
      ind_name <- paste0(nm, "__missing")
      indicator_map[[nm]] <- ind_name
      X[[ind_name]] <- as.numeric(miss_idx)
      fill_val <- med[[nm]]
      if (!is.finite(fill_val)) fill_val <- 0
      if (!any(!miss_idx)) {
        # Column is entirely NA; keep indicator but set baseline to 0
        X[[nm]] <- rep(fill_val, length(col))
      } else {
        X[[nm]][miss_idx] <- fill_val
      }
    }

    state$impute <- list(
      method = "none",
      fallback = if (fallback_used) "median" else NULL,
      med = med,
      indicators = indicator_map,
      base_cols = orig_cols
    )

    if (fallback_used) {
      warning("Missing values detected with impute$method='none'; applying median fallback with missingness indicators to ensure compatibility.",
              call. = FALSE)
      state$impute$label <- "impute: none (fallback median)"
    } else {
      state$impute$label <- "impute: none"
    }
  } else {
    stop("Unknown impute method. Use 'median', 'knn', 'missForest', or 'none'.")
  }

  state$impute$order <- colnames(X)
  if (is.null(state$impute$label))
    state$impute$label <- paste0("impute: ", state$impute$method)

  audit <- c(audit, list(state$impute$label))

  # 3) Normalization -----------------------------------------------------------
  norm_cfg <- steps$normalize %||% list()
  norm_method <- norm_cfg$method %||% "zscore"

  if (norm_method == "zscore") {
    mu <- vapply(X, mean, numeric(1))
    sdv <- vapply(X, stats::sd, numeric(1))
    sdv[sdv == 0 | !is.finite(sdv)] <- 1
    X <- sweep(X, 2, mu, "-")
    X <- sweep(X, 2, sdv, "/")
    state$normalize <- list(method = "zscore", mu = mu, sd = sdv)

  } else if (norm_method == "robust") {
    mu <- vapply(X, stats::median, numeric(1))
    sdv <- vapply(X, stats::mad, numeric(1))
    sdv[sdv < 1e-12 | !is.finite(sdv)] <- 1
    X <- sweep(X, 2, mu, "-")
    X <- sweep(X, 2, sdv, "/")
    state$normalize <- list(method = "robust", mu = mu, sd = sdv)

  } else if (norm_method == "none") {
    state$normalize <- list(method = "none")
  } else {
    stop("Unknown normalize method. Use 'zscore', 'robust', or 'none'.")
  }
  audit <- c(audit, list(paste0("normalize: ", state$normalize$method)))

  # 4) Filtering (variance / IQR) ---------------------------------------------
  filt_cfg <- steps$filter %||% list()
  var_th   <- filt_cfg$var_thresh %||% 0
  iqr_th   <- filt_cfg$iqr_thresh %||% 0
  min_keep <- filt_cfg$min_keep %||% NULL

  keep <- rep(TRUE, ncol(X))
  names(keep) <- colnames(X)
  if (var_th > 0) {
    v <- vapply(X, stats::var, numeric(1))
    keep <- keep & (v >= var_th)
  }
  if (iqr_th > 0) {
    i <- vapply(X, stats::IQR, numeric(1))
    keep <- keep & (i >= iqr_th)
  }

  # Adaptive: ensure at least min_keep variables remain
  if (!is.null(min_keep) && sum(keep) < min_keep) {
    # Rank by variance then IQR to recover features until min_keep
    v <- vapply(X, stats::var, numeric(1))
    i <- vapply(X, stats::IQR, numeric(1))
    ord <- order(v + i, decreasing = TRUE)
    keep[ord[seq_len(min_keep)]] <- TRUE
  }

  X <- X[, keep, drop = FALSE]
  state$filter <- list(var_thresh = var_th, iqr_thresh = iqr_th, keep = keep)
  audit <- c(audit, list(paste0("filter: var>=", var_th, ", iqr>=", iqr_th,
                                if (!is.null(min_keep)) paste0(", min_keep=", min_keep) else "")))

  # 5) Feature Selection -------------------------------------------------------
  fs_cfg <- steps$fs %||% list()
  fs_method <- fs_cfg$method %||% "none"
  selected <- seq_len(ncol(X))

  if (fs_method == "ttest") {
    if (!.guard_is_binary(y))
      stop("fs='ttest' requires a binary outcome 'y'.")
    yy <- as.factor(y)
    if (nlevels(yy) != 2) stop("fs='ttest' requires y with exactly 2 levels.")

    g1 <- X[which(yy == levels(yy)[2]), , drop = FALSE]
    g0 <- X[which(yy == levels(yy)[1]), , drop = FALSE]
    if (nrow(g1) < 2 || nrow(g0) < 2)
      stop("Not enough samples per class for t-test feature selection.")

    m1 <- vapply(g1, mean, numeric(1)); m0 <- vapply(g0, mean, numeric(1))
    s1 <- vapply(g1, stats::sd, numeric(1)); s0 <- vapply(g0, stats::sd, numeric(1))
    n1 <- nrow(g1); n0 <- nrow(g0)
    se <- sqrt(s1^2 / n1 + s0^2 / n0); se[se == 0] <- Inf
    tstat <- (m1 - m0) / se
    k <- max(50, ceiling(0.10 * length(tstat)))
    ord <- order(abs(tstat), decreasing = TRUE)
    selected <- ord[seq_len(min(k, length(ord)))]
    X <- X[, selected, drop = FALSE]
    state$fs <- list(method = "ttest", sel = selected)

  } else if (fs_method == "lasso") {
    if (!requireNamespace("glmnet", quietly = TRUE))
      stop("Package 'glmnet' is required for fs='lasso'.", call. = FALSE)
    fam <- if (task == "binomial") "binomial" else "gaussian"

    # guard against constant columns
    nonconst <- vapply(X, stats::sd, numeric(1)) > 1e-8
    if (!all(nonconst)) X <- X[, nonconst, drop = FALSE]

    fit0 <- glmnet::glmnet(as.matrix(X), y, family = fam, alpha = 1, standardize = FALSE)
    beta <- as.matrix(fit0$beta)
    # use lambda.1se if cv available; otherwise last column
    # (we didn't run cv.glmnet to keep speed/dep light)
    take <- ncol(beta)
    nz <- which(beta[, take] != 0)
    if (length(nz) == 0) {
      # fall back: take top |beta| coefficients
      nz <- order(rowSums(abs(beta)), decreasing = TRUE)[1:min(50, nrow(beta))]
    }
    selected <- nz
    X <- X[, selected, drop = FALSE]
    state$fs <- list(method = "lasso", sel = selected)

  } else if (fs_method == "pca") {
    p <- fs_cfg$pca_comp %||% min(50, ncol(X))
    pc <- stats::prcomp(X, center = FALSE, scale. = FALSE)
    X <- as.data.frame(pc$x[, seq_len(min(p, ncol(pc$x))), drop = FALSE])
    state$fs <- list(method = "pca",
                     p = p,
                     rotation = pc$rotation[, seq_len(min(p, ncol(pc$rotation))), drop = FALSE])

  } else if (fs_method == "none") {
    state$fs <- list(method = "none", sel = selected)
  } else {
    stop("Unknown fs method. Use 'none','ttest','lasso','pca'.")
  }
  audit <- c(audit, list(paste0("fs: ", state$fs$method)))

  # 6) Diagnostics, audit, hash -----------------------------------------------
  state$diagnostics <- list(
    n_train = n0,
    p_in    = length(orig_colnames),
    p_after_encode = p_after_encode,
    p_out   = ncol(X),
    removed_by_filter = sum(!state$filter$keep)
  )
  state$audit <- data.frame(
    step = seq_along(audit),
    action = as.character(unlist(audit)),
    stringsAsFactors = FALSE
  )
  impute_hash <- state$impute
  if (!is.null(impute_hash$knn_ref)) impute_hash$knn_ref <- NULL
  if (!is.null(impute_hash$rf_models)) impute_hash$rf_models <- NULL
  state$hash  <- .guard_hash(list(impute = impute_hash, normalize = state$normalize,
                                  filter = list(var_th = var_th, iqr_th = iqr_th),
                                  fs = state$fs))

  # 7) Build transformer (apply to new data) ----------------------------------
  transformer <- function(Xnew) {
    if (is.matrix(Xnew)) Xnew <- as.data.frame(Xnew, check.names = FALSE)
    stopifnot(is.data.frame(Xnew))

    # Handle mixed types with the same level structure learned during training
    has_factors <- isTRUE(state$encoding$has_factors) || length(state$encoding$levels) > 0
    if (has_factors) {
      prep_new <- .guard_ensure_levels(Xnew, state$encoding$levels)
      Xnew <- prep_new$data
      if (!ncol(Xnew)) stop("Input X has no columns after preprocessing.")
      mf <- stats::model.frame(~ ., data = as.data.frame(Xnew), na.action = stats::na.pass)
      Xnew <- stats::model.matrix(~ . - 1, data = mf)
      Xnew <- as.data.frame(Xnew, check.names = FALSE)
    } else {
      for (nm in names(Xnew)) {
        if (is.character(Xnew[[nm]])) {
          nz <- Xnew[[nm]][!is.na(Xnew[[nm]])]
          if (length(nz) && all(grepl("^\\s*-?\\d+(\\.\\d+)?\\s*$", nz))) {
            Xnew[[nm]] <- as.numeric(Xnew[[nm]])
          }
        }
      }
      if (!all(vapply(Xnew, is.numeric, logical(1)))) {
        stop("Non-numeric columns detected in new data; guarded fit was trained on numeric predictors only.")
      }
      if (!ncol(Xnew)) stop("Input X has no columns after preprocessing.")
      Xnew <- as.data.frame(Xnew, check.names = FALSE)
    }

    # Replace non-finite
    for (j in seq_len(ncol(Xnew))) {
      Xnew[[j]][!is.finite(Xnew[[j]])] <- NA
    }

    # Winsorization: apply train-fitted centers/scales to new data
    if (isTRUE(state$winsor$enabled)) {
      k <- state$winsor$k
      w_med <- state$winsor$med
      w_mad <- state$winsor$mad
      if (!is.null(w_med) && !is.null(w_mad)) {
        common <- intersect(names(w_med), colnames(Xnew))
        for (nm in common) {
          Xnew[[nm]] <- .guard_mad_winsorize(Xnew[[nm]], k = k,
                                             center = w_med[[nm]],
                                             scale = w_mad[[nm]])
        }
      } else {
        warning("Winsorization statistics missing; skipping winsorization to avoid leakage.",
                call. = FALSE)
      }
    }

    # Impute using TRAIN-FITTED parameters
    if (identical(state$impute$method, "none") && !is.null(state$impute$fallback)) {
      base_cols <- state$impute$base_cols %||% character(0)
      med <- state$impute$med %||% numeric(0)
      ind_map <- state$impute$indicators %||% list()

      # Ensure all training columns exist in new data
      missing_base <- setdiff(base_cols, colnames(Xnew))
      if (length(missing_base) > 0) {
        for (mc in missing_base) Xnew[[mc]] <- NA_real_
      }

      # Drop columns unseen during training imputation step
      expected_cols <- unique(c(base_cols, unlist(ind_map, use.names = FALSE)))
      extra_cols <- setdiff(colnames(Xnew), expected_cols)
      if (length(extra_cols) > 0)
        Xnew <- Xnew[, setdiff(colnames(Xnew), extra_cols), drop = FALSE]

      # Apply fallback imputation and create indicators
      for (nm in base_cols) {
        if (!nm %in% colnames(Xnew)) next
        idx <- is.na(Xnew[[nm]])
        ind_name <- ind_map[[nm]]
        if (!is.null(ind_name)) {
          Xnew[[ind_name]] <- as.numeric(idx)
        }
        fill_val <- med[[nm]]
        if (is.null(fill_val) || !is.finite(fill_val)) fill_val <- 0
        if (any(idx)) Xnew[[nm]][idx] <- fill_val
        if (!is.null(ind_name) && !ind_name %in% colnames(Xnew)) {
          Xnew[[ind_name]] <- 0
        }
      }

      # Ensure indicator columns exist even if no missing values in new data
      for (nm in base_cols) {
        ind_name <- ind_map[[nm]]
        if (!is.null(ind_name) && !ind_name %in% colnames(Xnew)) {
          Xnew[[ind_name]] <- 0
        }
      }

      # Reorder to match training order for downstream steps
      order_cols <- state$impute$order %||% colnames(Xnew)
      missing_ord <- setdiff(order_cols, colnames(Xnew))
      if (length(missing_ord) > 0) {
        for (mo in missing_ord) Xnew[[mo]] <- 0
      }
      Xnew <- Xnew[, order_cols, drop = FALSE]

    } else if (identical(state$impute$method, "knn")) {
      med <- state$impute$med
      ref <- state$impute$knn_ref
      scale <- state$impute$knn_scale
      k <- state$impute$k %||% 5
      if (is.null(ref) || is.null(scale) || is.null(med)) {
        warning("KNN imputation state missing; falling back to median imputation.",
                call. = FALSE)
        if (!is.null(med)) {
          Xnew <- .guard_align_cols(Xnew, names(med))
          Xnew <- .guard_fill_na(Xnew, med)
        }
      } else {
        Xnew <- .guard_knn_impute_new(Xnew, ref, k = k, med = med, scale = scale)
      }
    } else if (identical(state$impute$method, "missForest")) {
      med <- state$impute$med
      models <- state$impute$rf_models
      if (is.null(med)) {
        warning("missForest imputation state missing; skipping imputation.",
                call. = FALSE)
      } else if (is.null(models) || !length(models)) {
        warning("missForest models missing; falling back to median imputation.",
                call. = FALSE)
        Xnew <- .guard_align_cols(Xnew, names(med))
        Xnew <- .guard_fill_na(Xnew, med)
      } else {
        Xnew <- .guard_rf_impute_new(Xnew, models, med)
      }
    } else if (!identical(state$impute$method, "none")) {
      med <- state$impute$med
      # align columns if train/test encodings differ (missing columns -> fill with 0 then impute)
      missing_cols <- setdiff(names(med), colnames(Xnew))
      if (length(missing_cols) > 0) {
        for (mc in missing_cols) Xnew[[mc]] <- NA_real_
      }
      extra_cols <- setdiff(colnames(Xnew), names(med))
      if (length(extra_cols) > 0) {
        # drop unseen columns (keeps mapping to training space)
        Xnew <- Xnew[, setdiff(colnames(Xnew), extra_cols), drop = FALSE]
      }
      Xnew <- Xnew[, names(med), drop = FALSE]
      for (j in seq_along(med)) {
        idx <- is.na(Xnew[[j]])
        if (any(idx)) Xnew[[j]][idx] <- med[[j]]
      }
    }

    # Normalize with TRAIN-FITTED centers/scales
    if (identical(state$normalize$method, "zscore") || identical(state$normalize$method, "robust")) {
      mu <- state$normalize$mu; sdv <- state$normalize$sd
      # align columns
      missing_cols <- setdiff(names(mu), colnames(Xnew))
      if (length(missing_cols) > 0) for (mc in missing_cols) Xnew[[mc]] <- 0
      extra_cols <- setdiff(colnames(Xnew), names(mu))
      if (length(extra_cols) > 0) Xnew <- Xnew[, names(mu), drop = FALSE]
      Xnew <- sweep(Xnew, 2, mu, "-")
      Xnew <- sweep(Xnew, 2, sdv, "/")
    }

    # Filter: keep as per train
    keep <- state$filter$keep
    # keep may be indexed on pre-encoded space; rebuild logical selector by name
    keep_names <- names(keep)
    if (!is.null(keep_names)) {
      # ensure all expected columns exist and align order explicitly by name
      miss <- setdiff(keep_names, colnames(Xnew))
      if (length(miss) > 0) for (m in miss) Xnew[[m]] <- 0
      Xnew <- Xnew[, keep_names, drop = FALSE]
      keep <- keep[keep_names]
      Xnew <- Xnew[, keep, drop = FALSE]
    } else {
      # if names missing (rare), fall back to length match
      Xnew <- Xnew[, keep, drop = FALSE]
    }

    # Feature selection mapping
    if (state$fs$method %in% c("ttest", "lasso")) {
      sel <- state$fs$sel
      # sel indexes within filtered set at train-time; ensure bounds
      sel <- sel[sel >= 1 & sel <= ncol(Xnew)]
      Xnew <- Xnew[, sel, drop = FALSE]
    } else if (state$fs$method == "pca") {
      # Project using train PCA rotation
      R <- state$fs$rotation
      # align columns by names
      miss <- setdiff(rownames(R), colnames(Xnew))
      if (length(miss) > 0) for (m in miss) Xnew[[m]] <- 0
      Xnew <- Xnew[, rownames(R), drop = FALSE]
      Xnew <- as.matrix(Xnew) %*% as.matrix(R)
      Xnew <- as.data.frame(Xnew, check.names = FALSE)
    }

    Xnew
  }

  out <- structure(
    list(
      transform = transformer,
      state     = state,
      p_out     = ncol(X),
      steps     = list(
        impute = imp_method,
        normalize = norm_method,
        filter = list(var_thresh = var_th, iqr_thresh = iqr_th, min_keep = min_keep),
        fs = fs_method
      )
    ),
    class = "GuardFit"
  )
  out$hash <- state$hash
  out
}

# S3 helpers --------------------------------------------------------------

#' @export
print.GuardFit <- function(x, ...) {
  cat("Guarded preprocessing pipeline\n")
  cat(sprintf(" - Imputation: %s\n", x$steps$impute))
  cat(sprintf(" - Normalization: %s\n", x$steps$normalize))
  cat(sprintf(" - Filter: var>=%s, iqr>=%s\n",
              x$steps$filter$var_thresh, x$steps$filter$iqr_thresh))
  cat(sprintf(" - Feature selection: %s\n", x$steps$fs))
  cat(sprintf(" - Output features: %d\n", x$p_out))
  invisible(x)
}

#' @export
summary.GuardFit <- function(object, ...) {
  s <- object$state
  out <- list(
    diagnostics = s$diagnostics,
    audit = s$audit,
    hash  = s$hash
  )
  class(out) <- "summary.GuardFit"
  out
}

#' @export
print.summary.GuardFit <- function(x, ...) {
  cat("GuardFit summary\n")
  if (!is.null(x$diagnostics)) {
    d <- x$diagnostics
    cat(sprintf(" n_train=%s, p_out=%s, removed_by_filter=%s\n",
                d$n_train, d$p_out, d$removed_by_filter))
  }
  if (!is.null(x$audit)) {
    cat(" Steps:\n")
    print(x$audit, row.names = FALSE)
  }
  if (!is.null(x$hash))
    cat(sprintf(" hash: %s\n", x$hash))
  invisible(x)
}

# Convenience -------------------------------------------------------------

#' Apply a fitted GuardFit transformer to new data
#'
#' @description
#' Applies the preprocessing steps stored in a \code{GuardFit} object to new
#' data without refitting any statistics. This is designed to prevent
#' validation leakage that would occur if imputation, scaling, filtering, or
#' feature selection were recomputed on evaluation data. It enforces the
#' training schema by aligning columns and factor levels, and it errors when a
#' numeric-only training fit receives non-numeric predictors. It does not
#' detect label leakage, duplicate samples, or train/test contamination.
#'
#' @param fit A \code{GuardFit} object created by [.guard_fit()]. This required
#' argument (no default) contains the training-time preprocessing settings and
#' statistics. Changing \code{fit} (for example, a different imputation method
#' or feature selection step) changes the output columns and values.
#' @param newdata A matrix or data.frame of predictors with one row per sample.
#' This required argument (no default) is transformed using the training-time
#' parameters in \code{fit} only. Missing columns are added and filled, extra
#' columns are dropped, and factor levels are aligned to the training levels; if
#' the training fit was numeric-only, non-numeric columns in \code{newdata}
#' trigger an error.
#'
#' @return A data.frame of transformed predictors with the same number of rows
#' as \code{newdata}. Column order and content match the training pipeline and
#' may include derived features (one-hot encodings, missingness indicators, or
#' PCA components). This output is not a prediction; it is intended as input to
#' a downstream model and assumes the training-time preprocessing is valid for
#' the new data.
#'
#' @examples
#' x_train <- data.frame(a = c(1, 2, NA, 4), b = c(10, 11, 12, 13))
#' fit <- .guard_fit(
#'   x_train,
#'   y = c(0.1, 0.2, 0.3, 0.4),
#'   steps = list(impute = list(method = "median")),
#'   task = "gaussian"
#' )
#' x_new <- data.frame(a = c(NA, 5), b = c(9, 14))
#' out <- predict_guard(fit, x_new)
#' out
#' @export
predict_guard <- function(fit, newdata) {
  stopifnot(inherits(fit, "GuardFit"))
  fit$transform(newdata)
}
