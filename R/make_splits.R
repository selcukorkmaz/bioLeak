#' Create leakage-resistant splits
#'
#' @description
#' Generates leakage-safe cross-validation splits for common biomedical setups:
#' subject-grouped, batch-blocked, study leave-one-out, and time-series
#' rolling-origin. Supports repeats, optional stratification, nested inner CV,
#' and an optional prediction horizon for time series. Note that splits store
#' explicit indices, which can be memory-intensive for large \code{n} and many
#' repeats.
#'
#' @param x SummarizedExperiment or data.frame/matrix (samples x features).
#'   If SummarizedExperiment, metadata are taken from colData(x). If data.frame,
#'   metadata are taken from x (columns referenced by \code{group}, \code{batch}, \code{study}, \code{time}, \code{outcome}).
#' @param outcome character, outcome column name (used for stratification).
#' @param mode one of "subject_grouped","batch_blocked","study_loocv","time_series".
#' @param group optional subject/group id column (for subject_grouped). If NULL,
#'   each sample is treated as its own group (plain k-fold CV).
#' @param batch batch/plate/center column (for batch_blocked).
#' @param study study id column (for study_loocv).
#' @param time time column (numeric or POSIXct) for time_series.
#' @param v integer, number of folds (k) or rolling partitions.
#' @param repeats integer, number of repeats (>=1) for non-LOOCV modes.
#' @param stratify logical, keep outcome proportions similar across folds.
#'   For grouped modes, stratification is applied at the group level (by
#'   majority class per group) if \code{outcome} is provided; otherwise ignored.
#' @param nested logical, whether to attach inner CV splits (per outer fold)
#'   using the same \code{mode} on the outer training set (with \code{v} folds, 1 repeat).
#' @param seed integer seed.
#' @param horizon numeric (>=0), minimal time gap for time_series so that the
#'   training set only contains samples with time < min(test_time) when horizon = 0,
#'   and time <= min(test_time) - horizon otherwise.
#' @param progress logical, print progress for large jobs.
#' @param compact logical; store fold assignments instead of explicit train/test
#'   indices to reduce memory usage for large datasets. Not supported when
#'   \code{nested = TRUE}.
#'
#' @return LeakSplits S4 object
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   subject = rep(1:10, each = 2),
#'   outcome = rbinom(20, 1, 0.5),
#'   x1 = rnorm(20),
#'   x2 = rnorm(20)
#' )
#' splits <- make_splits(df, outcome = "outcome",
#'                       mode = "subject_grouped", group = "subject", v = 5)
#' @export
make_splits <- function(x, outcome = NULL,
                        mode = c("subject_grouped", "batch_blocked", "study_loocv", "time_series"),
                        group = NULL, batch = NULL, study = NULL, time = NULL,
                        v = 5, repeats = 1, stratify = FALSE, nested = FALSE,
                        seed = 1, horizon = 0, progress = TRUE, compact = FALSE) {

  mode <- match.arg(mode)
  stopifnot(is.numeric(v), v >= 2 || mode == "study_loocv")
  stopifnot(is.numeric(repeats), repeats >= 1)
  stopifnot(is.numeric(horizon), horizon >= 0)
  set.seed(seed)

  if (isTRUE(compact) && isTRUE(nested)) {
    warning("compact=TRUE is not supported with nested=TRUE; using compact=FALSE.", call. = FALSE)
    compact <- FALSE
  }

  # ---- Extract metadata ----
  X <- .bio_get_x(x)
  n <- nrow(X)

  if (.bio_is_se(x)) {
    cd <- as.data.frame(SummarizedExperiment::colData(x))
  } else if (is.data.frame(x)) {
    cd <- x
  } else if (is.matrix(x)) {
    cd <- data.frame(row_id = seq_len(n))
  } else {
    stop("x must be SummarizedExperiment, data.frame, or matrix.")
  }

  # Ensure unique row_id always exists
  if (!"row_id" %in% names(cd)) cd$row_id <- seq_len(nrow(cd))

  # Helper to ensure column existence
  .need_col <- function(col, label) {
    if (is.null(col)) stop(sprintf("Provide '%s' column name.", label))
    if (!col %in% names(cd)) stop(sprintf("Column '%s' not found in metadata.", col))
    if (anyNA(cd[[col]])) stop(sprintf("Missing values found in column '%s'.", col))
  }

  # Column checks
  if (mode == "subject_grouped") {
    if (is.null(group)) {
      group <- "row_id"  # auto-assign for plain CV
      message("No 'group' column provided -> using sample-wise CV (each sample = its own group).")
    }
    .need_col(group, "group")
  }
  if (mode == "batch_blocked")   .need_col(batch, "batch")
  if (mode == "study_loocv")     .need_col(study, "study")
  if (mode == "time_series")     .need_col(time, "time")

  if (isTRUE(stratify) && (is.null(outcome) || !outcome %in% names(cd)))
    warning("Stratification requested but 'outcome' not found; proceeding unstratified.")
  if (isTRUE(stratify) && !is.null(outcome) && length(outcome) == 2L) {
    warning("Stratification ignored for time/event outcomes; proceeding unstratified.")
    stratify <- FALSE
  }
  if (isTRUE(stratify) && !is.null(outcome) && outcome %in% names(cd)) {
    if (inherits(cd[[outcome]], "Surv")) {
      warning("Stratification ignored for survival outcomes; proceeding unstratified.")
      stratify <- FALSE
    }
  }

  # Utilities
  .as_factor <- function(z) as.factor(z)
  .majority_class <- function(y) {
    ux <- unique(y)
    if (length(ux) == 0) return(NA)
    ux[which.max(tabulate(match(y, ux)))]
  }
  .msg <- function(...) if (isTRUE(progress)) message(sprintf(...))

  indices <- list()
  inner_indices <- NULL
  summary_rows <- list()
  repeats_eff <- if (mode %in% c("study_loocv", "time_series")) 1 else repeats
  fold_assignments <- if (isTRUE(compact)) vector("list", repeats_eff) else NULL

  # ---- Subject-grouped ----
  if (mode == "subject_grouped") {
    gid <- .as_factor(cd[[group]])
    lev <- levels(gid)

    # Stratification
    if (isTRUE(stratify) && !is.null(outcome) && outcome %in% names(cd)) {
      y <- cd[[outcome]]
      g_lab <- tapply(seq_len(n), gid, function(ix) .majority_class(y[ix]))
      g_lab <- factor(g_lab)
      if (length(unique(g_lab)) < 2) {
        warning("Stratification ignored: only one class present at group level.")
        stratify <- FALSE
      }
    }

    for (r in seq_len(repeats)) {
      set.seed(seed + 1000 * r)
      if (isTRUE(compact)) {
        fold_assign <- rep(NA_integer_, n)
      }
      if (isTRUE(stratify) && exists("g_lab")) {
        gfold <- integer(length(lev)); names(gfold) <- lev
        for (cl in levels(g_lab)) {
          g_in <- names(g_lab)[g_lab == cl]
          gfold[g_in] <- sample(rep(seq_len(v), length.out = length(g_in)))
        }
      } else {
        gfold <- sample(rep(seq_len(v), length.out = length(lev)))
        names(gfold) <- lev
      }

      for (k in seq_len(v)) {
        test_g <- lev[gfold == k]
        test <- which(gid %in% test_g)
        train <- setdiff(seq_len(n), test)
        if (length(test) == 0L || length(train) == 0L) next
        if (isTRUE(compact)) {
          fold_assign[test] <- k
          indices[[length(indices) + 1L]] <- list(fold = k, repeat_id = r)
        } else {
          indices[[length(indices) + 1L]] <- list(train = train, test = test, fold = k, repeat_id = r)
        }
        summary_rows[[length(summary_rows) + 1L]] <- data.frame(
          fold = k,
          repeat_id = r,
          train_n = length(train),
          test_n = length(test),
          stringsAsFactors = FALSE
        )
      }
      if (isTRUE(compact)) fold_assignments[[r]] <- fold_assign
      .msg("subject_grouped: repeat %d/%d done.", r, repeats)
    }
  }

  # ---- Batch-blocked ----
  if (mode == "batch_blocked") {
    bid <- factor(cd[[batch]], levels = unique(cd[[batch]]))  # preserve order
    blevels <- levels(bid)

    if (length(blevels) < v) v <- length(blevels)

    if (isTRUE(stratify) && !is.null(outcome) && outcome %in% names(cd)) {
      y <- cd[[outcome]]
      b_lab <- tapply(seq_len(n), bid, function(ix) .majority_class(y[ix]))
      b_lab <- factor(b_lab)
      if (length(unique(b_lab)) < 2) {
        warning("Stratification ignored: only one class present at group level.")
        stratify <- FALSE
      }
    }

    for (r in seq_len(repeats)) {
      set.seed(seed + 1000 * r)
      if (isTRUE(compact)) {
        fold_assign <- rep(NA_integer_, n)
      }
      # if (isTRUE(stratify) && exists("b_lab")) {
      #   bfold <- integer(length(blevels)); names(bfold) <- blevels
      #   for (cl in levels(b_lab)) {
      #     b_in <- names(b_lab)[b_lab == cl]
      #     bfold[b_in] <- sample(rep(seq_len(v), length.out = length(b_in)))
      #   }
      # } else {
      #   if (v >= length(blevels)) {
      #     bfold <- seq_len(length(blevels))
      #     names(bfold) <- blevels
      #     bfold <- sample(bfold)
      #   } else {
      #     bfold <- sample(rep(seq_len(v), length.out = length(blevels)))
      #     names(bfold) <- blevels
      #   }
      # }

      if (v >= length(blevels)) {
        # One batch per fold: force unique assignment, ignore stratify
        bfold <- seq_len(length(blevels))
        names(bfold) <- blevels
        bfold <- sample(bfold)  # randomize order
      } else if (isTRUE(stratify) && exists("b_lab")) {
        bfold <- integer(length(blevels)); names(bfold) <- blevels
        for (cl in levels(b_lab)) {
          b_in <- names(b_lab)[b_lab == cl]
          bfold[b_in] <- sample(rep(seq_len(v), length.out = length(b_in)))
        }
      } else {
        bfold <- sample(rep(seq_len(v), length.out = length(blevels)))
        names(bfold) <- blevels
      }


      for (k in seq_len(v)) {
        test_b <- blevels[bfold == k]
        test <- which(bid %in% test_b)
        train <- setdiff(seq_len(n), test)
        if (length(train) == 0L || length(test) == 0L) next
        if (isTRUE(compact)) {
          fold_assign[test] <- k
          indices[[length(indices) + 1L]] <- list(fold = k, repeat_id = r)
        } else {
          indices[[length(indices) + 1L]] <- list(train = train, test = test, fold = k, repeat_id = r)
        }
        summary_rows[[length(summary_rows) + 1L]] <- data.frame(
          fold = k,
          repeat_id = r,
          train_n = length(train),
          test_n = length(test),
          stringsAsFactors = FALSE
        )
      }

      if (isTRUE(compact)) fold_assignments[[r]] <- fold_assign
      .msg("batch_blocked: repeat %d/%d done.", r, repeats)
    }
  }


  # ---- Study LOOCV ----
  if (mode == "study_loocv") {
    sid <- .as_factor(cd[[study]])
    if (isTRUE(compact)) {
      fold_assign <- rep(NA_integer_, n)
    }
    for (s in levels(sid)) {
      test <- which(sid == s)
      train <- setdiff(seq_len(n), test)
      if (length(train) == 0L || length(test) == 0L) next
      fold_id <- as.integer(which(levels(sid) == s))
      if (isTRUE(compact)) {
        fold_assign[test] <- fold_id
        indices[[length(indices) + 1L]] <- list(
          fold = fold_id,
          repeat_id = as.integer(1)
        )
      } else {
        indices[[length(indices) + 1L]] <- list(
          train = train,
          test = test,
          fold = fold_id,
          repeat_id = as.integer(1)
        )
      }
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        fold = fold_id,
        repeat_id = 1,
        train_n = length(train),
        test_n = length(test),
        stringsAsFactors = FALSE
      )
    }
    if (isTRUE(compact)) fold_assignments[[1]] <- fold_assign
  }

  # ---- Time-series ----
  if (mode == "time_series") {
    tt <- cd[[time]]
    if (!is.numeric(tt) && !inherits(tt, c("POSIXct", "Date")))
      stop("'time' column must be numeric, Date, or POSIXct.")
    ord <- order(tt)
    Xidx <- seq_len(n)[ord]

    fold_size <- max(1L, floor(n / v))
    if (isTRUE(compact)) {
      fold_assign <- rep(NA_integer_, n)
    }
    for (k in seq_len(v)) {
      test <- Xidx[max(1, (k - 1) * fold_size + 1):min(n, k * fold_size)]
      tmin <- min(tt[test])
      if (horizon == 0) {
        train <- Xidx[tt[Xidx] < tmin]
      } else {
        train <- Xidx[tt[Xidx] <= (tmin - horizon)]
      }
      if (length(train) < 1L) next
      if (length(test) < 3L) next
      if (isTRUE(compact)) {
        fold_assign[test] <- k
        indices[[length(indices) + 1L]] <- list(fold = as.integer(k), repeat_id = as.integer(1))
      } else {
        indices[[length(indices) + 1L]] <- list(train = train, test = test, fold = as.integer(k), repeat_id = as.integer(1))
      }
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        fold = as.integer(k),
        repeat_id = 1,
        train_n = length(train),
        test_n = length(test),
        stringsAsFactors = FALSE
      )
    }
    if (isTRUE(compact)) fold_assignments[[1]] <- fold_assign
  }

  # ---- Nested ----
  if (isTRUE(nested) && length(indices) > 0) {
    inner_indices <- vector("list", length(indices))
    for (i in seq_along(indices)) {
      tr_idx <- indices[[i]]$train
      x_inner <- if (.bio_is_se(x)) x[, tr_idx] else cd[tr_idx, , drop = FALSE]
      set.seed(seed + 1L)
      inner <- make_splits(
        x_inner, outcome = outcome, mode = mode, group = group,
        batch = batch, study = study, time = time, v = v,
        repeats = 1, stratify = stratify, nested = FALSE,
        seed = seed + 1L, horizon = horizon, progress = FALSE
      )
      inner_indices[[i]] <- inner@indices
    }
  }

  # ---- Summary ----
  if (length(indices) == 0L)
    stop("No valid folds generated. Check inputs (group/batch/study/time, v, repeats, horizon).")

  fold_repeat_keys <- vapply(indices, function(z) paste(z$fold, z$repeat_id), character(1))
  if (anyDuplicated(fold_repeat_keys)) {
    warning("Duplicate fold/repeat combinations detected; check v and repeats settings.")
  }

  split_summary <- if (isTRUE(compact)) {
    do.call(rbind, summary_rows)
  } else {
    data.frame(
      fold      = vapply(indices, `[[`, integer(1), "fold"),
      repeat_id = vapply(indices, `[[`, integer(1), "repeat_id"),
      train_n   = vapply(indices, function(z) length(z$train), integer(1)),
      test_n    = vapply(indices, function(z) length(z$test), integer(1)),
      stringsAsFactors = FALSE
    )
  }

  indices_hash <- .bio_hash_indices(indices)
  hash_val <- tryCatch({
    if (requireNamespace("digest", quietly = TRUE)) {
      digest::digest(list(
        mode = mode, n = n, v = v, repeats = repeats, stratify = stratify,
        seed = seed, horizon = horizon, group = group, batch = batch,
        study = study, time = time, indices_hash = indices_hash
      ))
    } else {
      indices_hash
    }
  }, error = function(e) NA_character_)

  info <- list(
    outcome = outcome, v = v, repeats = repeats, seed = seed, mode = mode,
    group = group, batch = batch, study = study, time = time,
    stratify = stratify, nested = nested, horizon = horizon,
    summary = split_summary, hash = hash_val, inner = inner_indices,
    compact = compact, fold_assignments = fold_assignments,
    coldata = cd
  )

  new("LeakSplits", mode = mode, indices = indices, info = info)
}


#' @title Display summary for LeakSplits objects
#' @description Prints fold counts, sizes, and hash metadata for quick inspection.
#' @param object LeakSplits object.
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   subject = rep(1:10, each = 2),
#'   outcome = rbinom(20, 1, 0.5),
#'   x1 = rnorm(20),
#'   x2 = rnorm(20)
#' )
#' splits <- make_splits(df, outcome = "outcome",
#'                       mode = "subject_grouped", group = "subject", v = 5)
#' show(splits)
#' }
#' @importMethodsFrom methods show
#' @export
setMethod("show", "LeakSplits", function(object) {
  cat(sprintf("LeakSplits object (mode = %s, v = %d, repeats = %d)\n",
              object@mode, object@info$v, object@info$repeats))
  cat(sprintf("Outcome: %s | Stratified: %s | Nested: %s\n",
              ifelse(is.null(object@info$outcome), "None", object@info$outcome),
              object@info$stratify, object@info$nested))
  cat("------------------------------------------------------\n")
  print(object@info$summary)
  cat("------------------------------------------------------\n")
  cat(sprintf("Total folds: %d | Hash: %s\n",
              nrow(object@info$summary), object@info$hash))
})
