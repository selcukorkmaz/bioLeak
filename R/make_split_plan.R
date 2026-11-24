# Stratified assignment of units (groups, batches) to v folds. Units are
# shuffled within each class and dealt round-robin as one sequence, class by
# class, continuing the fold counter across classes from a random start. Each
# fold gets floor or ceiling of (units / v) units and each class is spread as
# evenly as possible. Dealing every class from fold 1 (the previous approach)
# left high-numbered folds empty when classes had fewer units than v (e.g. 6
# of 10 requested folds) and could worsen balance.
.bio_stratified_deal <- function(labels, v) {
  labels <- droplevels(as.factor(labels))
  ids <- names(labels)
  if (is.null(ids)) ids <- as.character(seq_along(labels))
  ord <- unlist(lapply(levels(labels), function(cl) {
    u <- ids[labels == cl]
    if (length(u) > 1L) u[sample.int(length(u))] else u
  }), use.names = FALSE)
  start <- sample.int(v, 1L) - 1L
  fold <- ((start + seq_along(ord) - 1L) %% v) + 1L
  stats::setNames(as.integer(fold), ord)
}

# Row-balanced stratified assignment of units (groups, batches) to v folds,
# following the idea of scikit-learn's StratifiedGroupKFold. `unit` and `y`
# are per-row. Units are shuffled, then ordered largest first, and each is
# placed in the fold that minimises the mean (over classes) SD across folds of
# the share of that class's rows; ties go to the fold with fewer rows, then at
# random. When the units left are only enough to fill the empty folds, only
# empty folds are candidates, so every fold gets a unit when units >= v.
# Mixed-outcome groups are handled through their row counts per class.
.bio_row_balanced_assign <- function(unit, y, v) {
  unit <- droplevels(as.factor(unit))
  keep <- !is.na(y)
  y <- factor(as.character(y))
  cnt <- matrix(0, nlevels(unit), nlevels(y),
                dimnames = list(levels(unit), levels(y)))
  if (any(keep)) {
    tab <- table(unit[keep], droplevels(y[keep]))
    cnt[rownames(tab), colnames(tab)] <- unclass(tab)
  }
  size <- as.numeric(table(unit))
  tot <- colSums(cnt)
  tot[tot == 0] <- 1
  n_units <- nlevels(unit)
  ord <- sample.int(n_units)
  ord <- ord[order(-size[ord])]
  share_u <- sweep(cnt, 2L, tot, "/")  # unit's contribution to each class share
  share <- matrix(0, v, ncol(cnt))      # per-fold class shares
  fold_n <- numeric(v)
  fold <- integer(n_units)
  for (i in seq_along(ord)) {
    u <- ord[i]
    cand <- seq_len(v)
    empty <- which(fold_n == 0)
    if (length(empty) && (n_units - i + 1L) <= length(empty)) cand <- empty
    # SD across folds of each class share if u joins fold k (only row k
    # changes): from column sums and sums of squares.
    d <- share_u[u, ]
    s1 <- colSums(share)
    s2 <- colSums(share^2)
    old <- share[cand, , drop = FALSE]
    s2k <- sweep(-old^2 + sweep(old, 2L, d, "+")^2, 2L, s2, "+")
    s1k <- s1 + d
    var_k <- sweep(s2k, 2L, s1k^2 / v, "-") / (v - 1)
    score <- rowMeans(sqrt(pmax(var_k, 0)))
    best <- cand[score <= min(score) + 1e-12]
    best <- best[fold_n[best] == min(fold_n[best])]
    k <- if (length(best) == 1L) best else best[sample.int(length(best), 1L)]
    fold[u] <- k
    share[k, ] <- share[k, ] + d
    fold_n[k] <- fold_n[k] + size[u]
  }
  stats::setNames(as.integer(fold), levels(unit))
}

# TRUE when per-row labels have at least two classes (row-balanced method).
.bio_has_two_classes <- function(y) length(unique(y[!is.na(y)])) >= 2L

# Internal helper shared by explicit and compact time-series split builders.
.bio_time_series_train_indices <- function(time_vec, test_idx, candidate_idx = seq_along(time_vec),
                                           horizon = 0, purge = 0, embargo = 0,
                                           group_vec = NULL) {
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
  if (!is.null(group_vec) && length(train)) {
    test_groups <- unique(group_vec[test_idx])
    train <- train[!group_vec[train] %in% test_groups]
  }
  train
}

# Training indices for a compact fold, mirroring how fit_resample() builds
# them: time_series trains on earlier rows only (minus horizon / purge /
# embargo, and excluding test groups when a group column is set); combined
# drops rows sharing a secondary-constraint level with the test rows; all
# other modes train on every non-test row.
.bio_compact_train_indices <- function(splits, test_idx, n) {
  info <- splits@info
  cd   <- info$coldata
  if (identical(splits@mode, "time_series") && is.data.frame(cd) &&
      !is.null(info$time) && info$time %in% names(cd)) {
    gv <- if (!is.null(info$group) && info$group %in% names(cd)) cd[[info$group]] else NULL
    return(.bio_time_series_train_indices(
      time_vec = cd[[info$time]], test_idx = test_idx, candidate_idx = order(cd[[info$time]]),
      horizon = info$horizon %||% 0, purge = info$purge %||% 0,
      embargo = info$embargo %||% 0, group_vec = gv
    ))
  }
  remaining <- setdiff(seq_len(n), test_idx)
  if (identical(splits@mode, "combined") && is.data.frame(cd)) {
    cst <- info$constraints
    if (length(cst) > 1L) {
      for (ax in cst[-1]) {
        v <- cd[[ax$col]]
        if (!is.null(v)) remaining <- remaining[!v[remaining] %in% unique(v[test_idx])]
      }
    }
  }
  remaining
}

# Internal helper for combined (N-axis) splits.
.bio_make_combined_splits <- function(cd, n, constraints,
                                      v, repeats, stratify, outcome, seed,
                                      compact, stratify_method = "groups") {
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
  rows_y <- NULL
  if (isTRUE(stratify) && !is.null(outcome) && outcome %in% names(cd)) {
    y <- cd[[outcome]]
    if (identical(stratify_method, "rows")) {
      if (.bio_has_two_classes(y)) rows_y <- y
    } else {
      g_lab <- tapply(seq_len(n), p_groups, function(ix) .majority_class(y[ix]))
      g_lab <- factor(g_lab)
      if (length(unique(g_lab)) < 2) g_lab <- NULL
    }
  }

  # Build constraint-aware components of primary groups.
  # Subjects sharing any constraint-axis level form a component and should
  # be assigned to the same fold to avoid excessive exclusion.
  comp_map <- setNames(seq_along(p_levels), p_levels)
  for (ax_vec in constraint_axes) {
    ax_to_p <- split(as.character(p_groups), ax_vec)
    for (members in ax_to_p) {
      members <- unique(members)
      if (length(members) > 1L) {
        target_comp <- comp_map[members[1L]]
        for (m in members[-1L]) {
          old_comp <- comp_map[m]
          if (old_comp != target_comp) comp_map[comp_map == old_comp] <- target_comp
        }
      }
    }
  }
  components <- split(p_levels, comp_map)
  n_comp <- length(components)
  use_components <- n_comp >= v
  if (!use_components && length(constraint_axes) > 0L) {
    warning(sprintf(
      paste0("Only %d constraint-linked component(s) of primary groups ",
             "but v=%d folds; falling back to random group assignment."),
      n_comp, v), call. = FALSE)
  }

  for (r in seq_len(repeats)) {
    set.seed(seed + 1000 * r)
    if (isTRUE(compact)) fold_assign <- rep(NA_integer_, n)

    # Assign primary groups to folds
    gfold <- integer(length(p_levels))
    names(gfold) <- p_levels
    if (!is.null(rows_y)) {
      # Row-balanced: components (or primary groups) are the units
      unit_of_group <- if (use_components && length(constraint_axes) > 0L) {
        comp_map[p_levels]
      } else {
        stats::setNames(seq_along(p_levels), p_levels)
      }
      ufold <- .bio_row_balanced_assign(unit_of_group[as.character(p_groups)], rows_y, v)
      gfold[] <- ufold[as.character(unit_of_group[p_levels])]
    } else if (use_components && length(constraint_axes) > 0L) {
      # Constraint-aware: assign components via greedy bin packing
      comp_order <- sample(n_comp)
      if (!is.null(g_lab)) {
        comp_labels <- vapply(components, function(comp) {
          as.character(.majority_class(g_lab[comp]))
        }, character(1))
        g_levels <- levels(factor(g_lab))
        fold_counts <- matrix(0L, nrow = v, ncol = length(g_levels),
                              dimnames = list(NULL, g_levels))
        for (ci in comp_order) {
          cl <- comp_labels[ci]
          candidates <- which(fold_counts[, cl] == min(fold_counts[, cl]))
          target_fold <- if (length(candidates) == 1L) candidates else sample(candidates, 1L)
          for (pg in components[[ci]]) gfold[pg] <- target_fold
          fold_counts[target_fold, cl] <- fold_counts[target_fold, cl] + length(components[[ci]])
        }
      } else {
        fold_counts <- rep(0L, v)
        for (ci in comp_order) {
          candidates <- which(fold_counts == min(fold_counts))
          target_fold <- if (length(candidates) == 1L) candidates else sample(candidates, 1L)
          for (pg in components[[ci]]) gfold[pg] <- target_fold
          fold_counts[target_fold] <- fold_counts[target_fold] + length(components[[ci]])
        }
      }
    } else {
      # Random assignment (original behavior / fallback for highly connected data)
      if (!is.null(g_lab)) {
        gfold[names(g_lab)] <- .bio_stratified_deal(g_lab, v)[names(g_lab)]
      } else {
        gfold <- sample(rep(seq_len(v), length.out = length(p_levels)))
        names(gfold) <- p_levels
      }
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
       fold_assignments = fold_assignments, stratified = !is.null(g_lab) || !is.null(rows_y))
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
#' @param group subject/group id column. Required when mode is
#'   `subject_grouped` (use `group = "row_id"` to explicitly request sample-wise
#'   CV). With `mode = "time_series"` it makes the folds group-disjoint: no
#'   group appears in both the training and test rows of a fold. When regular
#'   rolling-origin folds lose time blocks (typically because groups recur
#'   across the time range), a panel fallback (each fold tests its groups' rows
#'   after the earliest timepoint and trains on other groups' earlier rows) is
#'   built too, and the plan with the larger sum over folds of (test rows x
#'   training rows) is used, so a plan that tests more rows but trains on very
#'   few does not win. A warning reports the lost blocks and their reasons; when
#'   the fallback is used, it also gives the fallback's coverage and training
#'   share.
#' @param batch batch/plate/center column (for batch_blocked).
#' @param study study id column (for study_loocv).
#' @param time time column (numeric, Date or POSIXct) for time_series. Rows
#'   are split into \code{v} contiguous blocks in time order; block 1 is never
#'   tested (nothing earlier to train on). A later block whose training set
#'   comes out empty is dropped with a warning giving the reason (tied times,
#'   group separation, or horizon/purge/embargo) and the number of rows tested.
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
#'   Each element is a list with \code{type} (one of \code{"subject"},
#'   \code{"batch"}, \code{"study"}, \code{"site"}, \code{"region"},
#'   \code{"platform"}, \code{"assay"}, \code{"relatedness"},
#'   \code{"spatial"}, \code{"group"}) and \code{col} (column name). The first element defines the
#'   primary grouping axis (fold driver); subsequent elements define exclusion
#'   constraints (training samples sharing constraint-axis levels with the test
#'   set are removed). Requires at least 2 elements. Cannot be used together with
#'   \code{primary_axis}/\code{secondary_axis}.
#' @param v integer, number of folds (k) or rolling partitions.
#' @param repeats integer, number of repeats (>=1) for non-LOOCV modes.
#' @param stratify logical, keep outcome (or \code{stratify_by}) proportions
#'   similar across folds. Groups are never split; how they are assigned is
#'   set by \code{stratify_method}.
#'   For \code{mode = "time_series"}, rolling-origin folds are contiguous
#'   time blocks and are never stratified (a warning says so); only the
#'   grouped panel fallback stratifies its groups. \code{splits@info$stratify}
#'   is \code{TRUE} only when stratification was applied.
#' @param stratify_by optional character column name in metadata to stratify on
#'   instead of \code{outcome} when \code{stratify = TRUE}.
#' @param stratify_method how grouped modes stratify when
#'   \code{stratify = TRUE}.
#'   \code{"groups"} (default) labels each group with its majority class and
#'   deals the groups of each class to folds round-robin in random order:
#'   folds get almost equal numbers of groups and each class's groups are
#'   spread evenly, but rows are not counted, so with unequal group sizes or
#'   mixed-outcome groups the per-fold class proportions in rows can be
#'   poorly balanced.
#'   \code{"rows"} balances rows: groups are shuffled, ordered largest first,
#'   and each is placed in the fold that keeps every class's share of rows
#'   most even across folds (the smallest mean standard deviation over
#'   classes; ties go to the fold with fewer rows), using the group's row
#'   count in each class. This is the approach of scikit-learn's
#'   \code{StratifiedGroupKFold}; it gives closer class proportions and fold
#'   sizes in rows, while the number of groups per fold may differ more.
#'   Applies to \code{subject_grouped}, \code{batch_blocked},
#'   \code{combined} (constraint-linked components are the units) and the
#'   grouped fallback of \code{time_series}.
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
#'       \code{stratify_method} (when stratified), \code{nested}, \code{horizon}, \code{purge}, \code{embargo},
#'       \code{summary} (data.frame of fold
#'       sizes), \code{hash} (reproducibility checksum), \code{inner}
#'       (nested inner splits if \code{nested = TRUE}), and \code{coldata}
#'       (sample metadata).}
#'   }
#'   Use the \code{show} method to print a summary; downstream access
#'   to the indices and metadata is normally done through the
#'   functions that consume a \code{LeakSplits} (for example
#'   \code{\link{fit_resample}}) rather than by reading slots directly.
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
                        v = 5, repeats = 1, stratify = FALSE, stratify_by = NULL,
                        stratify_method = c("groups", "rows"),
                        nested = FALSE,
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
    cd <- .bio_coldata_df(x)
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
  if (mode == "time_series") {
    .need_col(time, "time")
    if (!is.null(group)) .need_col(group, "group")
  }
  if (mode == "combined") {
    valid_axis_types <- c("subject", "batch", "study", "site", "region",
                          "platform", "assay", "relatedness", "spatial", "group")
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

  if (is.character(stratify) && length(stratify) == 1L) {
    stratify_by <- stratify
    stratify <- TRUE
  }
  strat_col <- stratify_by %||% outcome
  stratify_method_given <- !missing(stratify_method)
  stratify_method <- match.arg(stratify_method)
  if (stratify_method_given && identical(stratify_method, "rows") && !isTRUE(stratify)) {
    .bio_warn("stratify_method = \"rows\" has no effect because stratify = FALSE.",
              "bioLeak_fold_warning")
  }

  if (isTRUE(stratify) && (is.null(strat_col) || !all(strat_col %in% names(cd)))) {
    .bio_warn(sprintf("Stratification requested but the %s column %s was not found; proceeding unstratified.",
                      if (!is.null(stratify_by)) "stratify_by" else "outcome",
                      if (is.null(strat_col)) "(none given)" else sprintf("'%s'", paste(strat_col, collapse = "', '"))),
              "bioLeak_fold_warning")
    stratify <- FALSE
  }
  if (isTRUE(stratify) && !is.null(strat_col) && length(strat_col) == 2L) {
    .bio_warn("Stratification ignored for time/event outcomes; proceeding unstratified.",
              "bioLeak_fold_warning")
    stratify <- FALSE
  }
  if (isTRUE(stratify) && !is.null(strat_col) && strat_col %in% names(cd)) {
    if (inherits(cd[[strat_col]], "Surv")) {
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
    if (isTRUE(stratify) && !is.null(strat_col) && strat_col %in% names(cd)) {
      y <- cd[[strat_col]]
      if (identical(stratify_method, "rows")) {
        if (!.bio_has_two_classes(y)) {
          warning("Stratification ignored: only one class present.")
          stratify <- FALSE
        }
      } else {
        g_lab <- tapply(seq_len(n), gid, function(ix) .majority_class(y[ix]))
        g_lab <- factor(g_lab)
        if (length(unique(g_lab)) < 2) {
          warning("Stratification ignored: only one class present at group level.")
          stratify <- FALSE
        }
      }
    }

    for (r in seq_len(repeats)) {
      set.seed(seed + 1000 * r)
      if (isTRUE(compact)) {
        fold_assign <- rep(NA_integer_, n)
      }
      if (isTRUE(stratify) && identical(stratify_method, "rows")) {
        gfold <- .bio_row_balanced_assign(gid, y, v)[lev]
      } else if (isTRUE(stratify) && exists("g_lab")) {
        gfold <- integer(length(lev)); names(gfold) <- lev
        gfold[names(g_lab)] <- .bio_stratified_deal(g_lab, v)[names(g_lab)]
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

    if (isTRUE(stratify) && !is.null(strat_col) && strat_col %in% names(cd)) {
      y <- cd[[strat_col]]
      if (identical(stratify_method, "rows")) {
        if (!.bio_has_two_classes(y)) {
          warning("Stratification ignored: only one class present.")
          stratify <- FALSE
        }
      } else {
        b_lab <- tapply(seq_len(n), bid, function(ix) .majority_class(y[ix]))
        b_lab <- factor(b_lab)
        if (length(unique(b_lab)) < 2) {
          warning("Stratification ignored: only one class present at group level.")
          stratify <- FALSE
        }
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
      } else if (isTRUE(stratify) && identical(stratify_method, "rows")) {
        bfold <- .bio_row_balanced_assign(bid, y, v)[blevels]
      } else if (isTRUE(stratify) && exists("b_lab")) {
        bfold <- integer(length(blevels)); names(bfold) <- blevels
        bfold[names(b_lab)] <- .bio_stratified_deal(b_lab, v)[names(b_lab)]
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
    ts_group_vec <- if (!is.null(group)) cd[[group]] else NULL

    block_ids <- cut(seq_along(Xidx), breaks = v, include.lowest = TRUE, labels = FALSE)
    blocks <- split(Xidx, block_ids)
    if (isTRUE(compact)) {
      fold_assign <- rep(NA_integer_, n)
    }
    ts_dropped <- character()   # reasons for later blocks that could not be used
    for (k in seq_len(v)) {
      test <- blocks[[k]]
      if (!length(test)) next
      train <- .bio_time_series_train_indices(
        time_vec = tt,
        test_idx = test,
        candidate_idx = Xidx,
        horizon = horizon,
        purge = purge,
        embargo = embargo,
        group_vec = ts_group_vec
      )
      if (length(train) < 1L) {
        # The first block has nothing earlier to train on, by design. A later
        # block losing its training data is data loss and is reported.
        if (k > 1L) {
          no_grp <- if (!is.null(ts_group_vec)) {
            .bio_time_series_train_indices(time_vec = tt, test_idx = test, candidate_idx = Xidx,
                                           horizon = horizon, purge = purge, embargo = embargo)
          } else integer(0)
          reason <- if (!any(tt[Xidx] < min(tt[test]))) {
            "tied times (block starts at the earliest time)"
          } else if (length(no_grp)) {
            "group separation"
          } else {
            "horizon/purge/embargo"
          }
          ts_dropped <- c(ts_dropped, reason)
        }
        next
      }
      if (length(test) < 3L) {
        warning(sprintf("time_series fold %d skipped: only %d test sample(s) (minimum 3).", k, length(test)),
                call. = FALSE)
        next
      }
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

    # Longitudinal panel fallback: when recurring groups span all time blocks so
    # excluding test groups from earlier windows empties `train`, partition the
    # distinct groups across `v` folds; fold k tests all of its groups' rows
    # after the earliest timepoint and trains on the other groups' rows before
    # them. Every group is tested, but training is limited to the earliest
    # timepoint(s), which is reported. (Staggered rolling-origin cutoffs train
    # on more history but leave rows, and on unbalanced panels whole groups,
    # never tested, so they are not used.)
    warn_dropped <- function() {
      tested <- length(unique(unlist(lapply(blocks[seq_len(v)][-1L], function(b) b))))
      kept_tested <- sum(vapply(summary_rows, function(r) r$test_n, numeric(1)))
      reasons <- table(ts_dropped)
      .bio_warn(sprintf(paste0(
        "time_series: %d of %d later time block(s) could not be used because no earlier ",
        "training rows remained (%s), so only %d of %d rows are tested (rolling-origin CV never ",
        "tests the first block; here %d further rows are lost).%s"),
        length(ts_dropped), v - 1L,
        paste(sprintf("%s: %d", names(reasons), as.integer(reasons)), collapse = ", "),
        as.integer(kept_tested), n, as.integer(tested - kept_tested),
        if (!is.null(ts_group_vec) && any(grepl("group", ts_dropped))) paste0(
          " With group = , this happens when the later blocks hold few groups that also make ",
          "up all earlier data; consider a smaller v or a split that does not separate groups ",
          "over time.") else if (any(grepl("tied", ts_dropped))) paste0(
          " Many rows share the earliest time; consider a smaller v.") else ""),
        "bioLeak_fold_warning")
    }

    # Grouped series whose regular folds lost any block: also build the panel
    # fallback and keep the plan with the larger sum(test_n * train_n).
    try_fallback <- !is.null(ts_group_vec) && (length(indices) == 0L || length(ts_dropped) > 0L)
    reg_state <- list(indices = indices, summary_rows = summary_rows,
                      fold_assign = if (isTRUE(compact)) fold_assign else NULL)
    fb_small <- 0L
    if (try_fallback) {
      indices <- list(); summary_rows <- list()
      if (isTRUE(compact)) fold_assign <- rep(NA_integer_, n)
    }

    panel_used <- FALSE
    ts_strat_applied <- FALSE
    if (try_fallback) {
      panel_used <- TRUE
      tmin_global <- min(tt)
      future_idx <- if (horizon == 0 && purge == 0) {
        Xidx[tt[Xidx] > tmin_global]
      } else {
        Xidx[tt[Xidx] >= (tmin_global + horizon + purge)]
      }
      if (length(future_idx) > 0L) {
        gid <- .as_factor(ts_group_vec)
        lev <- levels(gid)
        set.seed(seed + 1000L)
        if (isTRUE(stratify) && !is.null(strat_col) && strat_col %in% names(cd)) {
          y_ts <- cd[[strat_col]]
          g_lab_ts <- if (identical(stratify_method, "groups")) {
            factor(tapply(seq_len(n), gid, function(ix) .majority_class(y_ts[ix])))
          }
          if (identical(stratify_method, "rows") && .bio_has_two_classes(y_ts)) {
            ts_strat_applied <- TRUE
            gfold_ts <- .bio_row_balanced_assign(gid, y_ts, v)[lev]
          } else if (length(unique(g_lab_ts)) >= 2L) {
            ts_strat_applied <- TRUE
            gfold_ts <- integer(length(lev)); names(gfold_ts) <- lev
            gfold_ts[names(g_lab_ts)] <- .bio_stratified_deal(g_lab_ts, v)[names(g_lab_ts)]
          } else {
            gfold_ts <- sample(rep(seq_len(v), length.out = length(lev)))
            names(gfold_ts) <- lev
          }
        } else {
          gfold_ts <- sample(rep(seq_len(v), length.out = length(lev)))
          names(gfold_ts) <- lev
        }

        if (isTRUE(compact)) fold_assign <- rep(NA_integer_, n)
        for (k in seq_len(v)) {
          test_g <- lev[gfold_ts == k]
          test <- future_idx[gid[future_idx] %in% test_g]
          if (!length(test)) next
          if (length(test) < 3L) {
            # same minimum as the regular folds
            fb_small <- fb_small + 1L
            next
          }
          train <- .bio_time_series_train_indices(
            time_vec = tt,
            test_idx = test,
            candidate_idx = Xidx,
            horizon = horizon,
            purge = purge,
            embargo = embargo,
            group_vec = ts_group_vec
          )
          if (length(train) < 1L) next
          if (isTRUE(compact)) {
            fold_assign[test] <- k
            indices[[length(indices) + 1L]] <- list(fold = as.integer(k), repeat_id = as.integer(1))
          } else {
            indices[[length(indices) + 1L]] <- list(
              train = train, test = test, fold = as.integer(k), repeat_id = as.integer(1)
            )
          }
          summary_rows[[length(summary_rows) + 1L]] <- data.frame(
            fold = as.integer(k),
            repeat_id = 1,
            train_n = length(train),
            test_n = length(test),
            stringsAsFactors = FALSE
          )
        }
      }
    }

    if (panel_used) {
      # Score a plan by the rows it tests, weighted by how much training data
      # the model that tests them saw: sum(test_n * train_n). The fallback is
      # kept only when it scores strictly higher.
      plan_score <- function(rows) sum(vapply(rows, function(r) r$test_n * r$train_n, numeric(1)))
      if (!length(summary_rows) || plan_score(summary_rows) <= plan_score(reg_state$summary_rows)) {
        # the regular folds score at least as high: keep them
        panel_used <- FALSE
        indices <- reg_state$indices; summary_rows <- reg_state$summary_rows
        if (isTRUE(compact)) fold_assign <- reg_state$fold_assign
      }
    }
    if (!panel_used && length(ts_dropped) && length(indices)) warn_dropped()
    if (isTRUE(stratify) && !(panel_used && ts_strat_applied)) {
      # Rolling-origin folds are contiguous time blocks; class balance cannot
      # be imposed on them. Record the plan as unstratified.
      .bio_warn(paste0(
        "stratify = TRUE is not applied to time_series folds: rolling-origin folds are ",
        "contiguous time blocks, so their class balance follows the data's time order. ",
        "The plan is recorded as unstratified."), "bioLeak_fold_warning")
      stratify <- FALSE
    }
    if (panel_used && length(summary_rows)) {
      sr <- do.call(rbind, summary_rows)
      frac <- mean(sr$train_n / pmax(1, n - sr$test_n))
      reg_tested <- sum(vapply(reg_state$summary_rows, function(r) r$test_n, numeric(1)))
      drop_why <- if (length(ts_dropped)) {
        rs <- table(ts_dropped)
        sprintf(" (%s)", paste(sprintf("%s: %d", names(rs), as.integer(rs)), collapse = ", "))
      } else ""
      reg_status <- if (length(reg_state$indices)) {
        sprintf("lost %d of %d later time blocks%s and tested %d of %d rows",
                length(ts_dropped), v - 1L, drop_why, as.integer(reg_tested), n)
      } else if (length(ts_dropped)) {
        paste0("were all empty", drop_why)
      } else {
        "were all skipped (fewer than 3 test rows each)"
      }
      .bio_warn(sprintf(paste0(
        "time_series with group: regular rolling-origin ",
        "folds %s. Used the panel fallback instead (it scores higher on the sum of test rows x ",
        "training rows; it tests %d of %d rows): each of %d fold(s) tests its groups' rows after ",
        "the earliest timepoint and trains on the other groups' rows before them, i.e. on %.0f%% ",
        "of the available rows on average.%s This keeps groups and time separated but trains on ",
        "little history; a single later-time holdout, or a split that does not separate groups ",
        "over time, may suit the question better."),
        reg_status, as.integer(sum(sr$test_n)), n, nrow(sr), 100 * frac,
        if (fb_small > 0L) sprintf(" %d fallback fold(s) with fewer than 3 test rows were skipped.",
                                   fb_small) else ""),
        "bioLeak_fold_warning")
    }
    if (!length(indices) && !is.null(ts_group_vec)) {
      rs <- table(ts_dropped)
      ts_empty_note <- paste0(
        " time_series with group: no regular fold was usable",
        if (length(rs)) sprintf(" (lost blocks: %s)",
                                paste(sprintf("%s: %d", names(rs), as.integer(rs)), collapse = ", ")) else "",
        if (fb_small > 0L) sprintf("; panel fallback: %d fold(s) had fewer than 3 test rows", fb_small) else "",
        ". Try a smaller v.")
    }
    if (isTRUE(compact)) fold_assignments[[1]] <- fold_assign
  }

  # ---- Combined ----
  if (mode == "combined") {
    combined <- .bio_make_combined_splits(
      cd = cd, n = n,
      constraints = constraints,
      v = v, repeats = repeats, stratify = stratify, outcome = strat_col,
      seed = seed, compact = compact, stratify_method = stratify_method
    )
    stratify <- isTRUE(stratify) && isTRUE(combined$stratified)
    indices <- combined$indices
    summary_rows <- combined$summary_rows
    fold_assignments <- combined$fold_assignments
  }

  # ---- Nested ----
  if (isTRUE(nested) && length(indices) > 0) {
    inner_indices <- vector("list", length(indices))
    inner_failed <- character()
    inner_notes <- character()
    for (i in seq_along(indices)) {
      tr_idx <- indices[[i]]$train
      x_inner <- if (.bio_is_se(x)) x[, tr_idx] else cd[tr_idx, , drop = FALSE]
      set.seed(seed + 1L)
      # An outer training set can be too small or too structured for an inner
      # split (e.g. a short panel); record that per fold instead of failing the
      # whole plan, and report inner-split warnings once, labelled as inner.
      inner <- tryCatch(withCallingHandlers(
        make_split_plan(
          x_inner, outcome = outcome, mode = mode, group = group,
          batch = batch, study = study, time = time, v = v,
          constraints = constraints,
          repeats = 1, stratify = stratify, stratify_by = stratify_by,
          stratify_method = if (isTRUE(stratify)) stratify_method else "groups",
          nested = FALSE,
          seed = seed + 1L, horizon = horizon, purge = purge, embargo = embargo,
          progress = FALSE
        ),
        warning = function(w) {
          inner_notes <<- c(inner_notes, conditionMessage(w))
          invokeRestart("muffleWarning")
        }),
        error = function(e) {
          inner_failed <<- c(inner_failed, sprintf("outer fold %d: %s", i, conditionMessage(e)))
          NULL
        })
      inner_indices[i] <- list(if (is.null(inner)) NULL else inner@indices)
    }
    if (length(inner_failed)) {
      .bio_warn(sprintf(paste0(
        "nested = TRUE: no inner split could be built for %d of %d outer fold(s); their inner ",
        "splits are NULL. First reason: %s"), length(inner_failed), length(indices),
        inner_failed[1L]), "bioLeak_fold_warning")
    }
    if (length(inner_notes)) {
      .bio_warn(sprintf("nested = TRUE: inner splits raised %d warning(s); first: %s",
                        length(inner_notes), inner_notes[1L]), "bioLeak_fold_warning")
    }
  }

  # ---- Summary ----
  if (length(indices) == 0L)
    .bio_stop(paste0("No valid folds generated. Check inputs (group/batch/study/time, v, repeats, horizon).",
                     get0("ts_empty_note", inherits = FALSE) %||% ""),
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
    # stratify_by is recorded whenever given, applied or not: fit_resample() and
    # tune_resample() keep that column out of the predictors (it is often
    # derived from the outcome).
    stratify = stratify, stratify_by = stratify_by,
    stratify_method = if (isTRUE(stratify)) stratify_method else NULL,
    nested = nested, horizon = horizon,
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
              if (isTRUE(object@info$stratify) &&
                  identical(object@info$stratify_method, "rows")) "TRUE (rows)"
              else object@info$stratify,
              object@info$nested))
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
#' @return A data.frame with one row per fold-by-column combination and
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
    if (mode == "time_series" && !is.null(splits@info$group))
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
