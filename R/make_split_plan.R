# Internal helper shared by explicit and compact time-series split builders.
.bio_time_series_train_indices <- function(time_vec, test_idx, candidate_idx = seq_along(time_vec),
                                           horizon = 0, purge = 0, embargo = 0) {
  if (!length(test_idx)) return(integer(0))
  tmin <- min(time_vec[test_idx])
  train <- if (horizon == 0 && purge == 0) {
    candidate_idx[time_vec[candidate_idx] < tmin]
  } else {
    candidate_idx[time_vec[candidate_idx] <= (tmin - horizon - purge)]
  }
  if (embargo > 0 && length(train)) {
    tmax <- max(time_vec[test_idx])
    train <- train[time_vec[train] <= (tmax - embargo)]
  }
  train
}

# Internal helper for combined (N-axis) splits.
.bio_make_combined_splits <- function(cd, n, constraints,
                                      v, repeats, stratify, outcome, seed,
                                      compact) {
  p_col <- constraints[[1]]$col
  p_groups <- as.factor(cd[[p_col]])
  p_levels <- levels(p_groups)

  # Precompute constraint axis vectors (all axes beyond the primary)
  constraint_axes <- lapply(constraints[-1], function(ax) cd[[ax$col]])

  indices <- list()
  summary_rows <- list()
  fold_assignments <- if (isTRUE(compact)) vector("list", repeats) else NULL

  .majority_class <- function(y) {
    ux <- unique(y)
    if (length(ux) == 0) return(NA)
    ux[which.max(tabulate(match(y, ux)))]
  }

  # Optional stratification at primary group level
  g_lab <- NULL
  if (isTRUE(stratify) && !is.null(outcome) && outcome %in% names(cd)) {
    y <- cd[[outcome]]
    g_lab <- tapply(seq_len(n), p_groups, function(ix) .majority_class(y[ix]))
    g_lab <- factor(g_lab)
    if (length(unique(g_lab)) < 2) g_lab <- NULL
  }

  for (r in seq_len(repeats)) {
    set.seed(seed + 1000 * r)
    if (isTRUE(compact)) fold_assign <- rep(NA_integer_, n)

    # Step 1: assign primary groups to folds
    if (!is.null(g_lab)) {
      gfold <- integer(length(p_levels)); names(gfold) <- p_levels
      for (cl in levels(g_lab)) {
        g_in <- names(g_lab)[g_lab == cl]
        gfold[g_in] <- sample(rep(seq_len(v), length.out = length(g_in)))
      }
    } else {
      gfold <- sample(rep(seq_len(v), length.out = length(p_levels)))
      names(gfold) <- p_levels
    }

    for (k in seq_len(v)) {
      # Test = all samples belonging to primary groups assigned to fold k
      test_p <- p_levels[gfold == k]
      test <- which(p_groups %in% test_p)
      if (length(test) == 0L) next

      # Step 2-3: for each constraint axis, remove from remaining any samples
      # whose constraint-axis level appears in test
      remaining <- setdiff(seq_len(n), test)
      for (ax_vec in constraint_axes) {
        test_levels <- unique(ax_vec[test])
        remaining <- remaining[!ax_vec[remaining] %in% test_levels]
      }
      train <- remaining

      if (length(train) == 0L) next

      if (isTRUE(compact)) {
        fold_assign[test] <- k
        indices[[length(indices) + 1L]] <- list(fold = k, repeat_id = r)
      } else {
        indices[[length(indices) + 1L]] <- list(
          train = train, test = test, fold = k, repeat_id = r
        )
      }
      summary_rows[[length(summary_rows) + 1L]] <- data.frame(
        fold = k, repeat_id = r,
        train_n = length(train), test_n = length(test),
        stringsAsFactors = FALSE
      )
    }
    if (isTRUE(compact)) fold_assignments[[r]] <- fold_assign
  }

  list(indices = indices, summary_rows = summary_rows,
       fold_assignments = fold_assignments)
}

#' Create leakage-resistant splits
#'
#' @description
#' Generates leakage-safe cross-validation splits for common biomedical setups:
#' subject-grouped, batch-blocked, study leave-one-out, and time-series
#' rolling-origin. Supports repeats, optional stratification, nested inner CV,
#' and optional prediction horizon/purge/embargo gaps for time series. Note that splits store
#' explicit indices, which can be memory-intensive for large \code{n} and many
#' repeats.
#'
#' @param x SummarizedExperiment or data.frame/matrix (samples x features).
#'   If SummarizedExperiment, metadata are taken from colData(x). If data.frame,
#'   metadata are taken from x (columns referenced by \code{group}, \code{batch}, \code{study}, \code{time}, \code{outcome}).
#' @param outcome character, outcome column name (used for stratification).
#' @param mode one of "subject_grouped","batch_blocked","study_loocv","time_series","combined".
#' @param group subject/group id column (for subject_grouped). Required when
#'   mode is `subject_grouped`; use `group = "row_id"` to explicitly request
#'   sample-wise CV.
#' @param batch batch/plate/center column (for batch_blocked).
#' @param study study id column (for study_loocv).
#' @param time time column (numeric or POSIXct) for time_series.
#' @param primary_axis List with elements \code{type} (one of \code{"subject"},
#'   \code{"batch"}, \code{"study"}) and \code{col} (column name). Used only
#'   when \code{mode = "combined"} to define the primary grouping axis.
#'   Deprecated in favor of \code{constraints}; still supported for backward
#'   compatibility.
#' @param secondary_axis List with elements \code{type} and \code{col}. Used
#'   only when \code{mode = "combined"} to define the secondary constraint axis.
#'   Training sets exclude samples whose secondary-axis levels appear in the
#'   test set. Deprecated in favor of \code{constraints}; still supported for
#'   backward compatibility.
#' @param constraints A list of constraint specifications for \code{mode = "combined"}.
#'   Each element is a list with \code{type} (one of \code{"subject"}, \code{"batch"},
#'   \code{"study"}) and \code{col} (column name). The first element defines the
#'   primary grouping axis (fold driver); subsequent elements define exclusion
#'   constraints (training samples sharing constraint-axis levels with the test
#'   set are removed). Requires at least 2 elements. Cannot be used together with
#'   \code{primary_axis}/\code{secondary_axis}.
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
#' @param purge numeric (>=0), additional gap removed immediately before each
#'   time-series test block.
#' @param embargo numeric (>=0), additional exclusion window anchored at the end
#'   of each time-series test block. Training rows with
#'   \code{time > max(test_time) - embargo} are removed.
#' @param progress logical, print progress for large jobs.
#' @param compact logical; store fold assignments instead of explicit train/test
#'   indices to reduce memory usage for large datasets. Not supported when
#'   \code{nested = TRUE}.
#' @param strict logical; deprecated and ignored. `subject_grouped` always
#'   requires a non-NULL `group`.
#'
#' @return A \code{\linkS4class{LeakSplits}} S4 object containing:
#'   \describe{
#'     \item{\code{mode}}{Character string indicating the splitting mode
#'       (\code{"subject_grouped"}, \code{"batch_blocked"}, \code{"study_loocv"},
#'       or \code{"time_series"}).}
#'     \item{\code{indices}}{List of fold descriptors, each containing
#'       \code{train} (integer vector of training indices), \code{test}
#'       (integer vector of test indices), \code{fold} (fold number), and
#'       \code{repeat_id} (repeat identifier). When \code{compact = TRUE},
#'       indices are stored as fold assignments instead.}
#'     \item{\code{info}}{List of metadata including \code{outcome}, \code{v},
#'       \code{repeats}, \code{seed}, grouping columns (\code{group},
#'       \code{batch}, \code{study}, \code{time}), \code{stratify},
#'       \code{nested}, \code{horizon}, \code{purge}, \code{embargo},
#'       \code{summary} (data.frame of fold
#'       sizes), \code{hash} (reproducibility checksum), \code{inner}
#'       (nested inner splits if \code{nested = TRUE}), and \code{coldata}
#'       (sample metadata).}
#'   }
#'   Use the \code{show} method to print a summary, or access slots directly
#'   with \code{@}.
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   subject = rep(1:10, each = 2),
#'   outcome = rbinom(20, 1, 0.5),
#'   x1 = rnorm(20),
#'   x2 = rnorm(20)
#' )
#' splits <- make_split_plan(df, outcome = "outcome",
#'                       mode = "subject_grouped", group = "subject", v = 5)
#' @export
make_split_plan <- function(x, outcome = NULL,
                        mode = c("subject_grouped", "batch_blocked", "study_loocv", "time_series", "combined"),
                        group = NULL, batch = NULL, study = NULL, time = NULL,
                        primary_axis = NULL, secondary_axis = NULL,
                        constraints = NULL,
                        v = 5, repeats = 1, stratify = FALSE, nested = FALSE,
                        seed = 1, horizon = 0, purge = 0, embargo = 0,
                        progress = TRUE, compact = FALSE,
                        strict = TRUE) {

  mode <- match.arg(mode)
  .bio_strict_checks(context = "make_split_plan", seed = seed,
                     nested = nested, mode = mode)
  stopifnot(is.numeric(v), v >= 2 || mode == "study_loocv")
  stopifnot(is.numeric(repeats), repeats >= 1)
  stopifnot(is.numeric(horizon), horizon >= 0)
  stopifnot(is.numeric(purge), purge >= 0)
  stopifnot(is.numeric(embargo), embargo >= 0)
  set.seed(seed)
  # Seed lineage:
  #   make_split_plan:  set.seed(seed), repeats: seed+1000*r, nested: seed+1
  #   fit_resample:     set.seed(seed), per-fold: seed+fold_id
  #   audit_leakage:    set.seed(seed), per-perm: seed+b
  # These are NOT isolated streams. For practical v/repeats/B values,
  # collisions do not occur.

  if (!identical(mode, "time_series") && (purge > 0 || embargo > 0)) {
    .bio_warn("purge and embargo are only used when mode = 'time_series'; ignoring both.",
              "bioLeak_fold_warning")
  }

  if (isTRUE(compact) && isTRUE(nested)) {
    .bio_warn("compact=TRUE is not supported with nested=TRUE; using compact=FALSE.",
              "bioLeak_fold_warning")
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
    .bio_stop("x must be SummarizedExperiment, data.frame, or matrix.",
              "bioLeak_input_error")
  }

  # Ensure unique row_id always exists
  if (!"row_id" %in% names(cd)) cd$row_id <- seq_len(nrow(cd))

  # Helper to ensure column existence
  .need_col <- function(col, label) {
    if (is.null(col)) .bio_stop(sprintf("Provide '%s' column name.", label),
                                "bioLeak_column_error")
    if (!col %in% names(cd)) .bio_stop(sprintf("Column '%s' not found in metadata.", col),
                                       "bioLeak_column_error")
    if (anyNA(cd[[col]])) .bio_stop(sprintf("Missing values found in column '%s'.", col),
                                    "bioLeak_column_error")
  }

  # Column checks
  if (mode == "subject_grouped") {
    if (is.null(group)) {
      .bio_stop("subject_grouped mode requires a non-NULL 'group' column.",
                "bioLeak_input_error")
    }
    .need_col(group, "group")
  }
  if (mode == "batch_blocked")   .need_col(batch, "batch")
  if (mode == "study_loocv")     .need_col(study, "study")
  if (mode == "time_series")     .need_col(time, "time")
  if (mode == "combined") {
    valid_axis_types <- c("subject", "batch", "study")
    .validate_axis <- function(ax, label) {
      if (is.null(ax) || !is.list(ax))
        .bio_stop(sprintf("'%s' must be a list with 'type' and 'col' elements.", label),
                  "bioLeak_input_error")
      if (!all(c("type", "col") %in% names(ax)))
        .bio_stop(sprintf("'%s' must have 'type' and 'col' elements.", label),
                  "bioLeak_input_error")
      if (!ax$type %in% valid_axis_types)
        .bio_stop(sprintf("'%s$type' must be one of: %s.", label, paste(valid_axis_types, collapse = ", ")),
                  "bioLeak_input_error")
      .need_col(ax$col, paste0(label, "$col"))
    }

    # Resolve constraints vs primary_axis/secondary_axis
    has_legacy <- !is.null(primary_axis) || !is.null(secondary_axis)
    has_constraints <- !is.null(constraints)

    if (has_constraints && has_legacy) {
      .bio_stop("Cannot specify both 'constraints' and 'primary_axis'/'secondary_axis'. Use one or the other.",
                "bioLeak_input_error")
    }

    if (has_constraints) {
      # Validate constraints list
      if (!is.list(constraints) || length(constraints) < 2L)
        .bio_stop("'constraints' must be a list with at least 2 elements (primary + constraint).",
                  "bioLeak_input_error")
      for (i in seq_along(constraints)) {
        .validate_axis(constraints[[i]], sprintf("constraints[[%d]]", i))
      }
      # Derive primary_axis/secondary_axis for backward compat in info slot
      primary_axis <- constraints[[1]]
      secondary_axis <- constraints[[2]]
    } else {
      # Legacy path: validate and convert to constraints
      .validate_axis(primary_axis, "primary_axis")
      .validate_axis(secondary_axis, "secondary_axis")
      constraints <- list(primary_axis, secondary_axis)
    }
  }

  if (isTRUE(stratify) && (is.null(outcome) || !outcome %in% names(cd)))
    .bio_warn("Stratification requested but 'outcome' not found; proceeding unstratified.",
              "bioLeak_fold_warning")
  if (isTRUE(stratify) && !is.null(outcome) && length(outcome) == 2L) {
    .bio_warn("Stratification ignored for time/event outcomes; proceeding unstratified.",
              "bioLeak_fold_warning")
    stratify <- FALSE
  }
  if (isTRUE(stratify) && !is.null(outcome) && outcome %in% names(cd)) {
    if (inherits(cd[[outcome]], "Surv")) {
      .bio_warn("Stratification ignored for survival outcomes; proceeding unstratified.",
                "bioLeak_fold_warning")
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
      .bio_stop("'time' column must be numeric, Date, or POSIXct.",
                "bioLeak_input_error")
    ord <- order(tt)
    Xidx <- seq_len(n)[ord]

    block_ids <- cut(seq_along(Xidx), breaks = v, include.lowest = TRUE, labels = FALSE)
    blocks <- split(Xidx, block_ids)
    if (isTRUE(compact)) {
      fold_assign <- rep(NA_integer_, n)
    }
    for (k in seq_len(v)) {
      test <- blocks[[k]]
      if (!length(test)) next
      train <- .bio_time_series_train_indices(
        time_vec = tt,
        test_idx = test,
        candidate_idx = Xidx,
        horizon = horizon,
        purge = purge,
        embargo = embargo
      )
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

  # ---- Combined ----
  if (mode == "combined") {
    combined <- .bio_make_combined_splits(
      cd = cd, n = n,
      constraints = constraints,
      v = v, repeats = repeats, stratify = stratify, outcome = outcome,
      seed = seed, compact = compact
    )
    indices <- combined$indices
    summary_rows <- combined$summary_rows
    fold_assignments <- combined$fold_assignments
  }

  # ---- Nested ----
  if (isTRUE(nested) && length(indices) > 0) {
    inner_indices <- vector("list", length(indices))
    for (i in seq_along(indices)) {
      tr_idx <- indices[[i]]$train
      x_inner <- if (.bio_is_se(x)) x[, tr_idx] else cd[tr_idx, , drop = FALSE]
      set.seed(seed + 1L)
      inner <- make_split_plan(
        x_inner, outcome = outcome, mode = mode, group = group,
        batch = batch, study = study, time = time, v = v,
        repeats = 1, stratify = stratify, nested = FALSE,
        seed = seed + 1L, horizon = horizon, purge = purge, embargo = embargo,
        progress = FALSE
      )
      inner_indices[[i]] <- inner@indices
    }
  }

  # ---- Summary ----
  if (length(indices) == 0L)
    .bio_stop("No valid folds generated. Check inputs (group/batch/study/time, v, repeats, horizon).",
              "bioLeak_split_error")

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
        seed = seed, horizon = horizon, purge = purge, embargo = embargo,
        group = group, batch = batch,
        study = study, time = time, indices_hash = indices_hash
      ))
    } else {
      indices_hash
    }
  }, error = function(e) NA_character_)

  info <- list(
    outcome = outcome, v = v, repeats = repeats_eff, seed = seed, mode = mode,
    group = group, batch = batch, study = study, time = time,
    primary_axis = primary_axis, secondary_axis = secondary_axis,
    constraints = constraints,
    stratify = stratify, nested = nested, horizon = horizon,
    purge = purge, embargo = embargo,
    summary = split_summary, hash = hash_val, inner = inner_indices,
    compact = compact, fold_assignments = fold_assignments,
    coldata = cd
  )

  result <- new("LeakSplits", mode = mode, indices = indices, info = info)

  if (.bio_is_strict() && !isTRUE(compact)) {
    check_split_overlap(result)
  }

  result
}


#' @title Display summary for LeakSplits objects
#' @description Prints fold counts, sizes, and hash metadata for quick inspection.
#' @param object LeakSplits object.
#' @return No return value, called for side effects (prints a summary to the
#'   console showing mode, fold count, repeats, outcome, stratification status,
#'   nested status, per-fold train/test sizes, and the reproducibility hash).
#' @examples
#' df <- data.frame(
#'   subject = rep(1:10, each = 2),
#'   outcome = rbinom(20, 1, 0.5),
#'   x1 = rnorm(20),
#'   x2 = rnorm(20)
#' )
#' splits <- make_split_plan(df, outcome = "outcome",
#'                       mode = "subject_grouped", group = "subject", v = 5)
#' show(splits)
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


#' Check split overlap invariants
#'
#' Verifies that a \code{\linkS4class{LeakSplits}} object satisfies the
#' expected no-overlap constraints for one or more grouping columns. For each
#' fold, the function checks that no group-level value appearing in the test
#' set is also present in the training set.
#'
#' @param splits A \code{LeakSplits} object from \code{\link{make_split_plan}}.
#' @param coldata A data.frame of sample metadata. When \code{NULL} (default),
#'   the function uses \code{splits@@info$coldata} if available.
#' @param cols Character vector of column names to check for overlap. When
#'   \code{NULL} (default), the function infers columns from the split mode
#'   (e.g., \code{group} for \code{subject_grouped}, \code{batch} for
#'   \code{batch_blocked}, both axes for \code{combined}).
#' @return A data.frame with one row per (fold × column) combination and
#'   columns \code{fold}, \code{repeat_id}, \code{col}, \code{n_overlap}
#'   (number of overlapping group values), and \code{pass} (logical).
#'   Invisible. Raises an error if any fold fails and \code{stop_on_fail = TRUE}.
#' @param stop_on_fail Logical; if \code{TRUE} (default), raises an error when
#'   any overlap is detected.
#' @export
check_split_overlap <- function(splits, coldata = NULL, cols = NULL,
                                stop_on_fail = TRUE) {
  stopifnot(inherits(splits, "LeakSplits"))

  cd <- coldata %||% splits@info$coldata
  if (is.null(cd) || !is.data.frame(cd)) {
    .bio_stop("coldata must be provided or present in splits@info$coldata.",
              "bioLeak_input_error")
  }

  # Infer columns to check from split mode when not supplied
  if (is.null(cols)) {
    mode <- splits@mode
    cols <- character(0)
    if (mode == "subject_grouped" && !is.null(splits@info$group))
      cols <- c(cols, splits@info$group)
    if (mode == "batch_blocked" && !is.null(splits@info$batch))
      cols <- c(cols, splits@info$batch)
    if (mode == "study_loocv" && !is.null(splits@info$study))
      cols <- c(cols, splits@info$study)
    if (mode == "combined") {
      if (!is.null(splits@info$constraints)) {
        for (ax in splits@info$constraints) {
          if (!is.null(ax$col)) cols <- c(cols, ax$col)
        }
      } else {
        if (!is.null(splits@info$primary_axis$col))
          cols <- c(cols, splits@info$primary_axis$col)
        if (!is.null(splits@info$secondary_axis$col))
          cols <- c(cols, splits@info$secondary_axis$col)
      }
    }
    if (!length(cols)) {
      message("No grouping columns detected; nothing to check.")
      return(invisible(data.frame()))
    }
  }

  missing_cols <- setdiff(cols, names(cd))
  if (length(missing_cols)) {
    .bio_stop(sprintf("Column(s) not found in coldata: %s",
                      paste(missing_cols, collapse = ", ")),
              "bioLeak_column_error")
  }

  # Resolve explicit indices only (compact mode not supported)
  folds <- splits@indices
  has_explicit <- all(vapply(folds, function(f) !is.null(f$train), logical(1)))
  if (!has_explicit) {
    message("compact splits (no explicit train/test indices) cannot be checked; skipping.")
    return(invisible(data.frame()))
  }

  result_rows <- list()
  for (fi in seq_along(folds)) {
    f <- folds[[fi]]
    tr <- f$train
    te <- f$test
    for (col in cols) {
      train_vals <- unique(cd[[col]][tr])
      test_vals  <- unique(cd[[col]][te])
      overlap    <- intersect(train_vals, test_vals)
      result_rows[[length(result_rows) + 1L]] <- data.frame(
        fold       = f$fold %||% fi,
        repeat_id  = f$repeat_id %||% 1L,
        col        = col,
        n_overlap  = length(overlap),
        pass       = length(overlap) == 0L,
        stringsAsFactors = FALSE
      )
    }
  }

  result <- do.call(rbind, result_rows)
  failures <- result[!result$pass, , drop = FALSE]

  if (nrow(failures) > 0L && isTRUE(stop_on_fail)) {
    .bio_stop(sprintf(
      "Overlap detected in %d fold(s). First failure: fold=%s, col=%s, n_overlap=%d.",
      nrow(failures),
      failures$fold[[1]], failures$col[[1]], failures$n_overlap[[1]]
    ), "bioLeak_overlap_error")
  } else if (nrow(failures) > 0L) {
    .bio_warn(sprintf("Overlap detected in %d fold(s).", nrow(failures)),
              "bioLeak_overlap_error")
  }

  invisible(result)
}
