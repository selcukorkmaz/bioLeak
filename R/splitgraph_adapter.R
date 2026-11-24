# Adapter: consume splitGraph split_spec objects as bioLeak LeakSplits.
#
# splitGraph emits a tool-agnostic `split_spec` describing sample-level
# grouping, blocking, and ordering. This file converts that artifact into a
# bioLeak `LeakSplits` by joining the spec's sample_data to the user's
# feature/outcome frame on sample_id and delegating to `make_split_plan()`.

#' Convert a splitGraph split_spec into bioLeak splits
#'
#' Consume a `split_spec` produced by \pkg{splitGraph} and build a
#' corresponding \code{LeakSplits} object via \code{\link{make_split_plan}}.
#' The spec supplies the grouping/blocking/ordering assignments; the caller
#' supplies the observation frame (features + outcome), joined on sample id.
#'
#' The mapping from \code{spec$constraint_mode} to
#' \code{make_split_plan(mode=)} is:
#' \itemize{
#'   \item \code{"subject"}   -> \code{"subject_grouped"}
#'   \item \code{"batch"}     -> \code{"batch_blocked"}
#'   \item \code{"study"}     -> \code{"study_loocv"}
#'   \item \code{"time"}      -> \code{"time_series"}
#'   \item \code{"site"}, \code{"region"}, \code{"platform"},
#'     \code{"assay"}, \code{"relatedness"}, \code{"spatial"},
#'     \code{"composite"} -> \code{"subject_grouped"} with
#'     \code{group = spec$group_var} (\code{"group_id"})
#' }
#'
#' For the grouping modes in the last row, splitGraph has already resolved
#' the dependency structure into one group id per sample (for
#' \code{"composite"}, the merged components of all relations in \code{via}),
#' so keeping each group intact within a fold honours the constraint.
#' Specs that require ordering (\code{spec$ordering_required = TRUE}) and
#' supply a valid ordering column (\code{"order_rank"} or \code{spec$time_var})
#' map to \code{"time_series"} while also passing \code{group = spec$group_var}
#' (for non-\code{"time"} modes) so both temporal ordering and group
#' disjointness are enforced. Unrecognised future modes fall back to
#' grouped folds on \code{group_var} with a warning.
#'
#' Blocking variables declared on the spec (\code{batch_group},
#' \code{study_group}, \code{site_group}, \code{region_group},
#' \code{platform_group}, \code{assay_group}), stratification (\code{stratum}
#' or \code{spec$stratum_var}), and ordering (\code{order_rank}) are forwarded
#' automatically when relevant; stratification is not requested for ordered
#' specs, whose time-series folds are contiguous time blocks. When \code{use_block_constraints = TRUE},
#' non-missing \code{spec$block_vars} on a grouped spec are enforced as
#' secondary exclusion axes via \code{make_split_plan(mode = "combined")}.
#'
#' @param spec A \code{split_spec} object from \pkg{splitGraph}.
#' @param data A data.frame (or SummarizedExperiment) containing at least one
#'   identifier column matching \code{sample_id_col} and an \code{outcome}
#'   column.
#' @param outcome Name of the outcome column in \code{data}.
#' @param sample_id_col Name of the sample-id column in \code{data}
#'   (default \code{"sample_id"}).
#' @param v Number of CV folds to request from \code{make_split_plan()}.
#' @param use_block_constraints \code{FALSE} (default), \code{TRUE}, or a
#'   character vector of block variable names. \code{spec$block_vars} lists
#'   the blocking annotations available in the spec (batch_group, study_group,
#'   ...), not axes the split already enforces. For a grouped spec,
#'   \code{TRUE} enforces all of them and a character vector enforces the
#'   named ones, as secondary exclusion axes in
#'   \code{make_split_plan(mode = "combined")}: training rows sharing a level of
#'   any enforced axis with the test rows are excluded. Enforcing many axes at
#'   once can leave no training data (an error names the axes) or leave some
#'   rows never tested (a warning reports how many). Names are checked against
#'   the spec's usable block_vars for every spec. For
#'   non-grouped specs (ordered, batch, study) the option does not apply and a
#'   warning is issued.
#' @param ... Additional arguments forwarded to \code{\link{make_split_plan}}
#'   (e.g. \code{stratify}, \code{seed}, \code{horizon}, \code{purge}).
#'
#' @return A \code{LeakSplits} object.
#'
#' @seealso \code{\link{make_split_plan}}, \code{\link{as_rsample}}
#' @export
as_leaksplits <- function(spec, data, outcome,
                          sample_id_col = "sample_id",
                          v = 5, use_block_constraints = FALSE, ...) {
  if (!inherits(spec, "split_spec")) {
    stop("`spec` must be a `split_spec` (produced by splitGraph::as_split_spec()).",
         call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }
  if (!sample_id_col %in% names(data)) {
    stop(sprintf("Column '%s' not found in `data`.", sample_id_col),
         call. = FALSE)
  }
  if (is.null(outcome) || !outcome %in% names(data)) {
    stop("`outcome` must name a column in `data`.", call. = FALSE)
  }

  sd <- spec$sample_data
  group_var <- spec$group_var %||% "group_id"
  if (!group_var %in% names(sd)) {
    stop(sprintf("split_spec sample_data has no group column '%s'.", group_var),
         call. = FALSE)
  }

  mode_map <- c(subject     = "subject_grouped",
                batch       = "batch_blocked",
                study       = "study_loocv",
                time        = "time_series",
                site        = "subject_grouped",
                region      = "subject_grouped",
                platform    = "subject_grouped",
                assay       = "subject_grouped",
                relatedness = "subject_grouped",
                spatial     = "subject_grouped",
                composite   = "subject_grouped")
  src_mode <- as.character(spec$constraint_mode %||% "composite")[1L]
  if (is.na(src_mode) || !nzchar(src_mode)) src_mode <- "composite"
  if (src_mode %in% names(mode_map)) {
    bl_mode <- unname(mode_map[src_mode])
  } else {
    warning(sprintf(paste0(
      "split_spec constraint_mode '%s' is not recognised by as_leaksplits(); ",
      "using grouped folds on '%s'."), src_mode, group_var), call. = FALSE)
    bl_mode <- "subject_grouped"
  }

  time_col <- spec$time_var
  if (is.null(time_col) || !nzchar(as.character(time_col)[1L])) {
    if ("order_rank" %in% names(sd) && !all(is.na(sd$order_rank))) {
      time_col <- "order_rank"
    } else {
      time_col <- NULL
    }
  }
  has_valid_time <- !is.null(time_col) && time_col %in% names(sd) &&
                    !all(is.na(sd[[time_col]]))

  if (isTRUE(spec$ordering_required)) {
    if (!has_valid_time) {
      stop(paste0(
        "This split_spec requires ordering (ordering_required = TRUE), ",
        "but sample_data has no valid ordering column ('order_rank' or spec$time_var). ",
        "Derive a split_spec with time ordering metadata, or call make_split_plan() directly."),
        call. = FALSE)
    }
    bl_mode <- "time_series"
  }

  n_groups <- length(unique(stats::na.omit(sd[[group_var]])))
  if ((identical(bl_mode, "subject_grouped") ||
       (identical(bl_mode, "time_series") && !identical(src_mode, "time"))) &&
      n_groups < 2L) {
    stop(sprintf(paste0(
      "split_spec (constraint_mode '%s') places all samples in %d group(s), ",
      "so no grouped folds can separate them. Relax the constraint (e.g. ",
      "fewer relations in `via`) or use a different mode."),
      src_mode, n_groups), call. = FALSE)
  }

  extra_candidates <- unique(c(
    "batch_group", "study_group", "site_group", "region_group",
    "platform_group", "assay_group", "timepoint_id", "order_rank", "stratum",
    as.character(spec$block_vars), as.character(spec$stratum_var),
    as.character(time_col)
  ))
  extra_candidates <- intersect(extra_candidates, names(sd))
  # Keep core columns plus any extra candidate column that is not all-NA
  non_na_extra <- extra_candidates[vapply(extra_candidates, function(cn) {
    cn %in% c("batch_group", "study_group", "timepoint_id", "order_rank") ||
      !all(is.na(sd[[cn]]))
  }, logical(1))]
  keep <- unique(c("sample_id", group_var, non_na_extra))
  sd <- sd[, keep, drop = FALSE]
  clash <- setdiff(intersect(names(data), names(sd)), sample_id_col)
  if (length(clash)) {
    stop(sprintf(paste0(
      "`data` already has column(s) %s, which the split_spec also supplies; ",
      "rename or drop them before calling as_leaksplits()."),
      paste(sprintf("'%s'", clash), collapse = ", ")), call. = FALSE)
  }

  data_ids <- as.character(data[[sample_id_col]])
  spec_ids <- as.character(sd$sample_id)
  if (anyNA(data_ids) || anyNA(spec_ids) ||
      anyDuplicated(data_ids) > 0L || anyDuplicated(spec_ids) > 0L ||
      length(data_ids) != length(spec_ids) || !setequal(data_ids, spec_ids)) {
    stop(sprintf(
      "Join on '%s' failed ID alignment (data has %d unique non-NA IDs across %d rows; spec has %d); check ID alignment.",
      sample_id_col, length(unique(stats::na.omit(data_ids))), nrow(data), nrow(sd)
    ), call. = FALSE)
  }
  idx_sd <- match(data_ids, spec_ids)
  sd_extra <- sd[idx_sd, setdiff(names(sd), "sample_id"), drop = FALSE]
  row.names(sd_extra) <- NULL
  merged <- cbind(data, sd_extra)

  .fill_partial_na_block <- function(df, col_nm, fallback_col = NULL) {
    vals <- as.character(df[[col_nm]])
    na_ix <- which(is.na(vals) | !nzchar(vals))
    if (length(na_ix) > 0L) {
      if (!is.null(fallback_col) && fallback_col %in% names(df)) {
        fb <- as.character(df[[fallback_col]][na_ix])
        bad_fb <- is.na(fb) | !nzchar(fb)
        if (any(bad_fb)) {
          fb[bad_fb] <- paste0("__unlinked_", col_nm, "_", na_ix[bad_fb])
        }
        vals[na_ix] <- fb
      } else {
        vals[na_ix] <- paste0("__unlinked_", col_nm, "_", na_ix)
      }
      df[[col_nm]] <- vals
    }
    df
  }

  args <- list(
    x       = merged,
    outcome = outcome,
    mode    = bl_mode,
    v       = v
  )

  if (!identical(src_mode, "time")) {
    args$group <- group_var
  }
  if ("batch_group" %in% names(merged) &&
      !all(is.na(merged$batch_group))) {
    if (identical(bl_mode, "batch_blocked") && anyNA(merged$batch_group)) {
      merged <- .fill_partial_na_block(merged, "batch_group", fallback_col = group_var)
      args$x <- merged
    }
    args$batch <- "batch_group"
  }
  if ("study_group" %in% names(merged) &&
      !all(is.na(merged$study_group))) {
    if (identical(bl_mode, "study_loocv") && anyNA(merged$study_group)) {
      merged <- .fill_partial_na_block(merged, "study_group", fallback_col = group_var)
      args$x <- merged
    }
    args$study <- "study_group"
  }
  if (has_valid_time && time_col %in% names(merged)) {
    args$time <- time_col
  }

  # Forward custom stratum_var when present and non-NA in sample_data
  stratum_col <- spec$stratum_var
  if (!is.null(stratum_col) && length(stratum_col) == 1L &&
      nzchar(stratum_col) && stratum_col %in% names(merged) &&
      !all(is.na(merged[[stratum_col]]))) {
    args$stratify_by <- stratum_col
    # Ordered specs give rolling-origin time-series folds, which are not
    # stratified; pass stratify = TRUE through `...` to request it anyway.
    if (!identical(stratum_col, outcome) && !identical(bl_mode, "time_series")) {
      args$stratify <- TRUE
    }
  }

  # Optional multi-axis block constraints via mode = "combined".
  # spec$block_vars lists the blocking annotations available; enforcing all of
  # them at once can be infeasible, so a subset can be named.
  ubc_axes <- NULL
  if (is.character(use_block_constraints)) {
    ubc_axes <- unique(use_block_constraints)
    ubc_on <- length(ubc_axes) > 0L
  } else if (is.logical(use_block_constraints) && length(use_block_constraints) == 1L &&
             !is.na(use_block_constraints)) {
    ubc_on <- use_block_constraints
  } else {
    stop("`use_block_constraints` must be TRUE, FALSE, or a character vector of block_vars.",
         call. = FALSE)
  }
  usable_bv <- intersect(as.character(spec$block_vars), names(merged))
  usable_bv <- usable_bv[vapply(usable_bv, function(cn) !all(is.na(merged[[cn]])), logical(1))]
  if (!is.null(ubc_axes)) {
    unknown <- setdiff(ubc_axes, usable_bv)
    if (length(unknown)) {
      stop(sprintf("use_block_constraints: %s %s not among the spec's usable block_vars (available: %s).",
                   paste(unknown, collapse = ", "), if (length(unknown) == 1L) "is" else "are",
                   if (length(usable_bv)) paste(usable_bv, collapse = ", ") else "none"),
           call. = FALSE)
    }
  }
  if (ubc_on && !identical(bl_mode, "subject_grouped")) {
    warning(sprintf(paste0(
      "use_block_constraints is ignored: it applies to grouped specs only, and this spec ",
      "maps to make_split_plan(mode = \"%s\")."), bl_mode), call. = FALSE)
  }
  block_cols <- character()
  if (ubc_on && identical(bl_mode, "subject_grouped")) {
    block_cols <- if (!is.null(ubc_axes)) ubc_axes else usable_bv
    if (!length(block_cols)) {
      warning("use_block_constraints = TRUE has no effect: the spec has no usable block_vars.",
              call. = FALSE)
    }
    if (length(block_cols) > 0L) {
      block_type_map <- c(
        batch_group    = "batch",
        study_group    = "study",
        site_group     = "site",
        region_group   = "region",
        platform_group = "platform",
        assay_group    = "assay"
      )
      constraints <- list(list(type = "subject", col = group_var))
      for (bc in block_cols) {
        if (anyNA(merged[[bc]])) {
          merged <- .fill_partial_na_block(merged, bc)
        }
        btype <- if (bc %in% names(block_type_map)) unname(block_type_map[bc]) else "batch"
        constraints[[length(constraints) + 1L]] <- list(type = btype, col = bc)
      }
      args$x <- merged
      args$mode <- "combined"
      args$constraints <- constraints
    }
  }

  user <- list(...)
  if (length(block_cols) && any(c("mode", "constraints") %in% names(user))) {
    warning(paste0(
      "use_block_constraints: `mode` / `constraints` passed through `...` replace the ",
      "block constraints it built; the block_vars are not enforced as such."), call. = FALSE)
  }
  args[names(user)] <- user
  ubc_applied <- length(block_cols) > 0L && identical(args$mode, "combined") &&
    !"constraints" %in% names(user)

  res <- tryCatch(
    do.call(make_split_plan, args),
    bioLeak_split_error = function(e) {
      # Only explain failures of the outer split this option built; a nested
      # inner split or user-supplied constraints can fail for other reasons.
      if (ubc_applied && !isTRUE(args$nested) &&
          grepl("No valid folds generated", conditionMessage(e), fixed = TRUE)) {
        .bio_stop(sprintf(paste0(
          "use_block_constraints: excluding training rows that share a level of %s with the ",
          "test rows left no usable folds (every fold's training set was emptied). %s"),
          .ubc_axis_text(block_cols, merged), .ubc_advice(block_cols, merged, args$v %||% v)),
          "bioLeak_split_error")
      }
      stop(e)
    }
  )
  # Partial loss: constraints can leave fewer folds than requested, so part of
  # the data is never tested. Say so instead of returning a degenerate CV.
  if (ubc_applied) {
    # folds possible without the block constraints: at most one per group
    n_req <- min(args$v %||% v, n_groups) * (res@info$repeats %||% 1L)
    n_got <- length(res@indices)
    # test rows of every fold, explicit or compact
    tested <- length(unique(unlist(lapply(seq_along(res@indices), function(i)
      .dlsi_fold_test_idx(res, i)))))
    if (tested < nrow(merged)) {   # rows never tested; fewer folds alone lose no data
      warning(sprintf(paste0(
        "use_block_constraints: only %d of %d requested folds have usable training data, and ",
        "%d of %d rows are ever tested, after excluding training rows that share a level of %s ",
        "with the test rows. %s"),
        n_got, n_req, tested, nrow(merged), .ubc_axis_text(block_cols, merged),
        .ubc_advice(block_cols, merged, args$v %||% v)), call. = FALSE)
    }
  }
  res@info$sample_id_col <- sample_id_col
  res
}

# "batch_group (6 levels), ..." counting real levels only (not NA fill-ins).
.ubc_axis_text <- function(block_cols, merged) {
  lv <- vapply(block_cols, function(cn) {
    x <- as.character(merged[[cn]])
    length(unique(x[!is.na(x) & !startsWith(x, "__unlinked_")]))
  }, integer(1))
  paste(sprintf("%s (%d levels)", block_cols, lv), collapse = ", ")
}

# Advice that never suggests the axis set that just failed.
.ubc_advice <- function(block_cols, merged, v = NULL) {
  real_lv <- vapply(block_cols, function(cn) {
    x <- as.character(merged[[cn]])
    length(unique(x[!is.na(x) & !startsWith(x, "__unlinked_")]))
  }, integer(1))
  fewer_v <- if (!is.null(v) && is.finite(v) && min(real_lv) < v && min(real_lv) >= 2L) {
    sprintf(" Fewer folds may also work: %s has only %d levels, so try v <= %d.",
            block_cols[which.min(real_lv)], min(real_lv), min(real_lv))
  } else ""
  base <- if (length(block_cols) > 1L) {
    sprintf(paste0("Enforce fewer axes, e.g. use_block_constraints = \"%s\", or derive a ",
                   "composite split_spec over the axes you need in splitGraph."),
            block_cols[which.max(real_lv)])
  } else {
    paste0("This axis cannot be held out alongside the groups at this number of folds; use a ",
           "split_spec whose grouping includes it (e.g. splitGraph mode = \"composite\") ",
           "or drop use_block_constraints.")
  }
  paste0(base, fewer_v)
}
