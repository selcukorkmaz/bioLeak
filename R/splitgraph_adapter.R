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
#' automatically when relevant. When \code{use_block_constraints = TRUE},
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
#' @param use_block_constraints Logical (default \code{FALSE}). When
#'   \code{TRUE} and the spec uses a grouped mode with non-empty
#'   \code{spec$block_vars}, converts \code{group_var} and \code{block_vars}
#'   into multi-axis \code{constraints} in \code{make_split_plan(mode = "combined")}.
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
    if (!identical(stratum_col, outcome)) {
      args$stratify <- TRUE
    }
  }

  # Optional multi-axis block constraints via mode = "combined"
  if (isTRUE(use_block_constraints) && identical(bl_mode, "subject_grouped")) {
    block_cols <- intersect(as.character(spec$block_vars), names(merged))
    block_cols <- block_cols[vapply(block_cols, function(cn) !all(is.na(merged[[cn]])), logical(1))]
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
  args[names(user)] <- user

  res <- do.call(make_split_plan, args)
  res@info$sample_id_col <- sample_id_col
  res
}
