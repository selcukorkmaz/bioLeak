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
#'   \item \code{"composite"} -> \code{"combined"}
#' }
#'
#' Blocking variables declared on the spec (\code{batch_group},
#' \code{study_group}) and ordering (\code{order_rank}) are forwarded
#' automatically when relevant.
#'
#' @param spec A \code{split_spec} object from \pkg{splitGraph}.
#' @param data A data.frame (or SummarizedExperiment) containing at least one
#'   identifier column matching \code{sample_id_col} and an \code{outcome}
#'   column.
#' @param outcome Name of the outcome column in \code{data}.
#' @param sample_id_col Name of the sample-id column in \code{data}
#'   (default \code{"sample_id"}).
#' @param v Number of CV folds to request from \code{make_split_plan()}.
#' @param ... Additional arguments forwarded to \code{\link{make_split_plan}}
#'   (e.g. \code{stratify}, \code{seed}, \code{horizon}, \code{purge}).
#'
#' @return A \code{LeakSplits} object.
#'
#' @seealso \code{\link{make_split_plan}}, \code{\link{as_rsample}}
#' @export
as_leaksplits <- function(spec, data, outcome,
                          sample_id_col = "sample_id",
                          v = 5, ...) {
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
  keep <- c("sample_id", spec$group_var,
            intersect(c("batch_group", "study_group",
                        "timepoint_id", "order_rank"), names(sd)))
  sd <- sd[, keep, drop = FALSE]

  merged <- merge(data, sd,
                  by.x = sample_id_col, by.y = "sample_id",
                  all.x = FALSE, all.y = FALSE, sort = FALSE)
  if (nrow(merged) != nrow(sd)) {
    stop(sprintf(
      "Join on '%s' produced %d rows but spec has %d samples; check ID alignment.",
      sample_id_col, nrow(merged), nrow(sd)
    ), call. = FALSE)
  }

  mode_map <- c(subject   = "subject_grouped",
                batch     = "batch_blocked",
                study     = "study_loocv",
                time      = "time_series",
                composite = "combined")
  src_mode <- spec$constraint_mode %||% "composite"
  bl_mode <- mode_map[[src_mode]]
  if (is.null(bl_mode)) bl_mode <- "subject_grouped"

  args <- list(
    x       = merged,
    outcome = outcome,
    mode    = bl_mode,
    v       = v
  )

  args$group <- spec$group_var
  if ("batch_group" %in% names(merged) &&
      !all(is.na(merged$batch_group))) {
    args$batch <- "batch_group"
  }
  if ("study_group" %in% names(merged) &&
      !all(is.na(merged$study_group))) {
    args$study <- "study_group"
  }
  if (!is.null(spec$time_var) && spec$time_var %in% names(merged)) {
    args$time <- spec$time_var
  }

  user <- list(...)
  args[names(user)] <- user

  do.call(make_split_plan, args)
}
