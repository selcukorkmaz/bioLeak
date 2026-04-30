#' S4 Classes for bioLeak Pipeline
#'
#' These classes capture splits, model fits, and audit diagnostics produced by
#' \code{make_split_plan()}, \code{fit_resample()}, and \code{audit_leakage()}.
#'
#' @slot mode Splitting mode. One of "subject_grouped", "batch_blocked",
#'   "study_loocv", "time_series", or "combined".
#' @slot indices List of resampling descriptors (train/test indices when available)
#' @slot info (LeakSplits) Metadata associated with the split plan (mode, coldata, hash, etc.)
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
#' @slot info (LeakFit) Metadata about the model fit (sample IDs, timings, provenance, etc.)
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
#' @slot info (LeakAudit) Metadata about the audit (mechanism summary, settings, provenance, etc.)
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
#' @slot info (LeakDeltaLSI) Metadata including R_naive, R_guarded, paired status, and block details
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

# ---- LeakAudit: show() method ------------------------------------------
# Brief auto-print representation. Detailed reporting is produced by
# summary(). Added in 0.3.7 (Comment 8 response): standard R idiom is
# that every public S4 class should have a show() method so that
# auto-printing at the console is informative.

#' @title Display summary for LeakAudit objects
#' @description Prints a brief one-screen summary of a LeakAudit, including
#'   task and outcome, the permutation-gap statistic, and counts of
#'   batch-association rows, target-leakage features, and duplicate pairs.
#'   Use \code{summary()} for the full diagnostic report.
#' @param object A \code{LeakAudit} object.
#' @return No return value, called for side effects (prints a brief summary
#'   to the console). Returns \code{object} invisibly.
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   subject = rep(1:6, each = 2),
#'   outcome = rbinom(12, 1, 0.5),
#'   x1 = rnorm(12),
#'   x2 = rnorm(12)
#' )
#' splits <- make_split_plan(df, outcome = "outcome",
#'                       mode = "subject_grouped", group = "subject", v = 3,
#'                       progress = FALSE)
#' custom <- list(
#'   glm = list(
#'     fit = function(x, y, task, weights, ...) {
#'       stats::glm(y ~ ., data = as.data.frame(x),
#'                  family = stats::binomial(), weights = weights)
#'     },
#'     predict = function(object, newdata, task, ...) {
#'       as.numeric(stats::predict(object,
#'                                 newdata = as.data.frame(newdata),
#'                                 type = "response"))
#'     }
#'   )
#' )
#' fit <- fit_resample(df, outcome = "outcome", splits = splits,
#'                     learner = "glm", custom_learners = custom,
#'                     metrics = "auc", refit = FALSE, seed = 1)
#' aud <- audit_leakage(fit, metric = "auc", B = 10,
#'                      X_ref = df[, c("x1", "x2")])
#' show(aud)
#' @importMethodsFrom methods show
#' @export
setMethod("show", "LeakAudit", function(object) {
  perm_str <- if (is.data.frame(object@permutation_gap) &&
                  nrow(object@permutation_gap) > 0L &&
                  "metric_obs" %in% names(object@permutation_gap)) {
    obs <- as.numeric(object@permutation_gap$metric_obs[1L])
    gap <- if ("gap" %in% names(object@permutation_gap))
      as.numeric(object@permutation_gap$gap[1L]) else NA_real_
    pv <- if ("p_value" %in% names(object@permutation_gap))
      as.numeric(object@permutation_gap$p_value[1L]) else NA_real_
    sprintf("metric=%.3f, gap=%.3f, p=%.3f", obs, gap, pv)
  } else {
    "(not computed)"
  }
  task_str <- if (length(object@fit@task) && nzchar(object@fit@task))
    object@fit@task else "(unknown)"
  outcome_str <- if (length(object@fit@outcome) && nzchar(object@fit@outcome))
    object@fit@outcome else "(unknown)"

  cat("A LeakAudit object\n")
  cat(sprintf("  Task:               %s\n", task_str))
  cat(sprintf("  Outcome:            %s\n", outcome_str))
  cat(sprintf("  Permutation-gap:    %s\n", perm_str))
  cat(sprintf("  Batch association:  %d row(s)\n",
              if (is.data.frame(object@batch_assoc)) nrow(object@batch_assoc) else 0L))
  cat(sprintf("  Target leakage:     %d feature(s)\n",
              if (is.data.frame(object@target_assoc)) nrow(object@target_assoc) else 0L))
  cat(sprintf("  Duplicate pairs:    %d\n",
              if (is.data.frame(object@duplicates)) nrow(object@duplicates) else 0L))
  cat("Use summary(<obj>) for the full diagnostic report.\n")
  invisible(object)
})
