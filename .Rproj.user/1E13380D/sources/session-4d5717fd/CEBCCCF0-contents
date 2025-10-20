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

.guard_mad_winsorize <- function(x, k = 5) {
  if (!is.numeric(x)) return(x)
  if (all(!is.finite(x))) return(x)
  m <- stats::median(x, na.rm = TRUE)
  s <- stats::mad(x, na.rm = TRUE, constant = 1)
  lo <- m - k * s; hi <- m + k * s
  x[x < lo] <- lo
  x[x > hi] <- hi
  x
}

.guard_hash <- function(obj) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(obj))
  }
  paste0("h", sprintf("%08X", as.integer(stats::runif(1, 0, .Machine$integer.max))))
}

#' Guarded preprocessing: fit on training, apply to new data (leakage-safe)
#'
#' @param X matrix/data.frame of predictors (training)
#' @param y optional outcome vector (used by some FS methods)
#' @param steps named list configuring stages:
#' \itemize{
#'   \item \code{impute = list(method = "median"|"knn"|"mice"|"none", k = 5,
#'         winsor = TRUE, winsor_k = 5)}
#'   \item \code{normalize = list(method = "zscore"|"robust"|"none")}
#'   \item \code{filter = list(var_thresh = 0, iqr_thresh = 0,
#'         min_keep = NULL)}  # adaptive to keep at least this many
#'   \item \code{fs = list(method = "none"|"ttest"|"lasso"|"pca",
#'         pca_comp = 50)}
#'   \item \code{parallel = FALSE}
#' }
#' @param task "binomial" or "gaussian" (affects FS)
#' @return S3 object of class "GuardFit" with fields:
#'   \code{transform}, \code{state}, \code{p_out}, \code{steps}
#' @export
.guard_fit <- function(X, y = NULL, steps = list(), task = c("binomial","gaussian")) {
  task <- match.arg(task)
  # Coerce to data.frame for flexible handling
  if (is.matrix(X)) X <- as.data.frame(X, check.names = FALSE)
  stopifnot(is.data.frame(X))
  n0 <- nrow(X)

  # Keep original column names for diagnostics
  orig_colnames <- colnames(X)

  # 0) Mixed-type handling: one-hot encode non-numeric columns -----------------
  num_idx <- vapply(X, is.numeric, logical(1))
  if (!all(num_idx)) {
    # Safely handle factors/characters via model.matrix
    mm <- stats::model.matrix(~ . - 1, data = X)
    X <- as.data.frame(mm, check.names = FALSE)
  }

  # Replace non-finite with NA to unify missingness treatment
  for (j in seq_len(ncol(X))) {
    X[[j]][!is.finite(X[[j]])] <- NA
  }

  state <- list()
  audit <- list()

  # 1) Robust Winsorization (train-only) --------------------------------------
  impute_cfg <- steps$impute %||% list()
  winsor_on <- impute_cfg$winsor %||% TRUE
  winsor_k  <- impute_cfg$winsor_k %||% 5

  if (isTRUE(winsor_on)) {
    X <- as.data.frame(lapply(X, .guard_mad_winsorize, k = winsor_k),
                       check.names = FALSE)
    state$winsor <- list(enabled = TRUE, k = winsor_k)
    audit <- c(audit, list(paste0("winsor: k=", winsor_k)))
  } else {
    state$winsor <- list(enabled = FALSE)
  }

  # 2) Imputation (train-only fitting) ----------------------------------------
  imp_method <- impute_cfg$method %||% "median"

  if (imp_method == "median") {
    meds <- vapply(X, function(col) stats::median(col, na.rm = TRUE), numeric(1))
    for (j in seq_len(ncol(X))) {
      idx <- is.na(X[[j]])
      if (any(idx)) X[[j]][idx] <- meds[j]
    }
    state$impute <- list(method = "median", med = meds)

  } else if (imp_method == "knn") {
    # Fit imputation on train, store column medians as safe application rule
    if (!requireNamespace("VIM", quietly = TRUE)) {
      stop("Package 'VIM' is required for impute$method='knn'. Install it or use 'median'.", call. = FALSE)
    }
    k <- impute_cfg$k %||% 5
    X_imp <- suppressWarnings(VIM::kNN(X, k = k, imp_var = FALSE))
    # After KNN on train, derive column medians to apply on newdata (leak-safe)
    meds <- vapply(X_imp, function(col) stats::median(col, na.rm = TRUE), numeric(1))
    X <- X_imp
    state$impute <- list(method = "knn", k = k, med = meds)

  } else if (imp_method == "mice") {
    if (!requireNamespace("mice", quietly = TRUE)) {
      stop("Package 'mice' is required for impute$method='mice'.", call. = FALSE)
    }
    maxit <- impute_cfg$maxit %||% 5
    mfit <- mice::mice(X, m = 1, maxit = maxit, printFlag = FALSE)
    X_imp <- mice::complete(mfit)
    # Store column medians from completed train to apply safely to new data
    meds <- vapply(X_imp, function(col) stats::median(col, na.rm = TRUE), numeric(1))
    X <- X_imp
    state$impute <- list(method = "mice", maxit = maxit, med = meds)

  } else if (imp_method == "none") {
    state$impute <- list(method = "none")
  } else {
    stop("Unknown impute method. Use 'median', 'knn', 'mice', or 'none'.")
  }
  audit <- c(audit, list(paste0("impute: ", state$impute$method)))

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
    sdv[sdv == 0 | !is.finite(sdv)] <- 1
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
    p_after_encode = ncol(as.data.frame(stats::model.matrix(~ . - 1, data = as.data.frame(setNames(rep(0, length(orig_colnames)), orig_colnames))))),
    p_out   = ncol(X),
    removed_by_filter = sum(!state$filter$keep)
  )
  state$audit <- data.frame(step = seq_along(audit), action = unlist(audit), stringsAsFactors = FALSE)
  state$hash  <- .guard_hash(list(impute = state$impute, normalize = state$normalize,
                                  filter = list(var_th = var_th, iqr_th = iqr_th),
                                  fs = state$fs))

  # 7) Build transformer (apply to new data) ----------------------------------
  transformer <- function(Xnew) {
    if (is.matrix(Xnew)) Xnew <- as.data.frame(Xnew, check.names = FALSE)
    stopifnot(is.data.frame(Xnew))

    # Handle mixed types: model.matrix with columns as in train encoding space
    if (!all(vapply(Xnew, is.numeric, logical(1)))) {
      Xnew <- stats::model.matrix(~ . - 1, data = Xnew)
      Xnew <- as.data.frame(Xnew, check.names = FALSE)
    }

    # Replace non-finite
    for (j in seq_len(ncol(Xnew))) {
      Xnew[[j]][!is.finite(Xnew[[j]])] <- NA
    }

    # Winsorization: apply same MAD k to new data (centers from new data)
    if (isTRUE(state$winsor$enabled)) {
      k <- state$winsor$k
      Xnew <- as.data.frame(lapply(Xnew, .guard_mad_winsorize, k = k), check.names = FALSE)
    }

    # Impute using TRAIN-FITTED parameters (medians from train/imputed-train)
    if (!identical(state$impute$method, "none")) {
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
      if (!all(keep_names %in% colnames(Xnew))) {
        # align: add missing as zeros
        miss <- setdiff(keep_names, colnames(Xnew))
        if (length(miss) > 0) for (m in miss) Xnew[[m]] <- 0
        Xnew <- Xnew[, keep_names, drop = FALSE]
      }
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
#' @param fit GuardFit object returned by .guard_fit()
#' @param newdata matrix/data.frame to transform
#' @export
predict_guard <- function(fit, newdata) {
  stopifnot(inherits(fit, "GuardFit"))
  fit$transform(newdata)
}
