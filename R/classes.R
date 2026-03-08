#' S4 Classes for bioLeak Pipeline
#'
#' These classes capture splits, model fits, and audit diagnostics produced by
#' \code{make_split_plan()}, \code{fit_resample()}, and \code{audit_leakage()}.
#'
#' @slot mode Splitting mode. One of "subject_grouped", "batch_blocked",
#'   "study_loocv", "time_series", or "combined".
#' @slot indices List of resampling descriptors (train/test indices when available)
#' @slot info Metadata associated with split or fit
#' @return An S4 object of the respective class.
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

#' @rdname LeakClasses
#' @slot metric Performance metric compared between pipelines
#' @slot exchangeability Exchangeability assumption used for the sign-flip test
#' @slot tier Inference tier label based on effective number of repeats
#' @slot strict Whether strict mode was requested
#' @slot R_eff Effective number of paired repeats available for inference
#' @slot delta_lsi Huber-robust point estimate of repeat-level metric difference
#' @slot delta_lsi_ci BCa 95\% CI for delta_lsi (NA when R_eff < 10)
#' @slot delta_metric Arithmetic mean of repeat-level metric differences
#' @slot delta_metric_ci BCa 95\% CI for delta_metric (NA when R_eff < 10)
#' @slot p_value Sign-flip randomization test p-value (NA when R_eff < 5 or unpaired)
#' @slot inference_ok TRUE when tier A (R_eff >= 20, paired, finite p and CI)
#' @slot folds_naive Per-fold data frame for the naive pipeline
#' @slot folds_guarded Per-fold data frame for the guarded pipeline
#' @slot repeats_naive Per-repeat aggregate data frame for the naive pipeline
#' @slot repeats_guarded Per-repeat aggregate data frame for the guarded pipeline
#' @slot info Additional metadata including R_naive, R_guarded, paired status
#' @seealso [delta_lsi()]
#' @exportClass LeakDeltaLSI
setClass("LeakDeltaLSI",
  slots = c(
    metric          = "character",
    exchangeability = "character",
    tier            = "character",
    strict          = "logical",
    R_eff           = "integer",
    delta_lsi       = "numeric",
    delta_lsi_ci    = "numeric",
    delta_metric    = "numeric",
    delta_metric_ci = "numeric",
    p_value         = "numeric",
    inference_ok    = "logical",
    folds_naive     = "data.frame",
    folds_guarded   = "data.frame",
    repeats_naive   = "data.frame",
    repeats_guarded = "data.frame",
    info            = "list"
  ),
  prototype = list(
    metric          = character(0),
    exchangeability = character(0),
    tier            = "D_insufficient",
    strict          = FALSE,
    R_eff           = 0L,
    delta_lsi       = NA_real_,
    delta_lsi_ci    = c(NA_real_, NA_real_),
    delta_metric    = NA_real_,
    delta_metric_ci = c(NA_real_, NA_real_),
    p_value         = NA_real_,
    inference_ok    = FALSE,
    folds_naive     = data.frame(),
    folds_guarded   = data.frame(),
    repeats_naive   = data.frame(),
    repeats_guarded = data.frame(),
    info            = list()
  )
)

LeakDeltaLSI <- function(...) methods::new("LeakDeltaLSI", ...)
