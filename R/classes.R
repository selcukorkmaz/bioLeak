#' S4 Classes for bioLeak Pipeline
#'
#' These classes capture splits, model fits, and audit diagnostics produced by
#' \code{make_split_plan()}, \code{fit_resample()}, and \code{audit_leakage()}.
#'
#' @slot mode Splitting mode (e.g., "grouped_cv", "batch_blocked")
#' @slot indices List of resampling descriptors (train/test indices when available)
#' @slot info Metadata associated with split or fit
#' @seealso [make_split_plan()], [fit_resample()], [audit_leakage()]
#' @rdname LeakClasses
#' @exportClass LeakSplits
setClass("LeakSplits",
         slots = c(mode = "character",
                   indices = "list",
                   info = "list"))

#' @rdname LeakClasses
#' @slot splits A [`LeakSplits`] object used for resampling
#' @slot metrics Model performance metrics per resample
#' @slot metric_summary Summary of metrics across resamples
#' @slot audit Audit information per resample
#' @slot predictions List of prediction objects
#' @slot preprocess Preprocessing steps used during fitting
#' @slot learners Learner definitions used in the pipeline
#' @slot outcome Outcome variable name
#' @slot task Modeling task name
#' @slot feature_names Feature names included in the model
#' @slot info Additional metadata about the fit
#' @seealso [fit_resample()]
#' @exportClass LeakFit
setClass("LeakFit",
         slots = c(
           splits = "LeakSplits",
           metrics = "data.frame",
           metric_summary = "data.frame",
           audit = "data.frame",
           predictions = "list",
           preprocess = "list",
           learners = "list",
           outcome = "character",
           task = "character",
           feature_names = "character",
           info = "list"
         ))

#' @rdname LeakClasses
#' @slot fit A [`LeakFit`] object used to generate the audit
#' @slot permutation_gap Data frame summarising permutation gaps
#' @slot perm_values Numeric vector of permutation-based scores
#' @slot batch_assoc Data frame of batch associations
#' @slot target_assoc Data frame of feature-wise outcome associations
#' @slot duplicates Data frame detailing duplicate records
#' @slot trail List capturing audit trail information
#' @seealso [audit_leakage()], [audit_report()]
#' @exportClass LeakAudit
setClass("LeakAudit",
         slots = c(
           fit = "LeakFit",
           permutation_gap = "data.frame",
           perm_values = "numeric",
           batch_assoc = "data.frame",
           target_assoc = "data.frame",
           duplicates = "data.frame",
           trail = "list",
           info = "list"))

setValidity("LeakSplits", function(object) {
  if (!length(object@indices)) return("indices list cannot be empty")
  TRUE
})

setValidity("LeakFit", function(object) {
  if (!inherits(object@splits, "LeakSplits")) return("splits must be a LeakSplits object")
  TRUE
})

setValidity("LeakAudit", function(object) {
  if (!inherits(object@fit, "LeakFit")) return("fit must be a LeakFit object")
  TRUE
})

LeakSplits <- function(mode = NA_character_, indices = list(), info = list()) {
  methods::new("LeakSplits", mode = mode, indices = indices, info = info)
}

LeakFit <- function(...) methods::new("LeakFit", ...)

LeakAudit <- function(...) methods::new("LeakAudit", ...)
