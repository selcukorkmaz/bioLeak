## ---------------------------------------------------------------------
## R/accessors.R
## ---------------------------------------------------------------------
## Public accessor functions for the bioLeak S4 result classes.
##
## These wrappers expose slot data as named functions so that downstream
## code (replication scripts, vignettes, end-user analyses) can read the
## components of LeakFit, LeakAudit, and LeakDeltaLSI objects without
## reaching into S4 internals via `@`. Slot definitions are unchanged;
## these accessors are purely additive.
## ---------------------------------------------------------------------


# ---- LeakFit ---------------------------------------------------------

#' Per-fold metric data frame from a LeakFit
#'
#' Returns the per-fold metric data frame stored in a `LeakFit` object.
#' Each row is one (fold, repeat, learner) combination; columns include
#' the requested metric values such as `auc` and any task-specific
#' performance scores.
#'
#' @param fit A `LeakFit` object returned by [fit_resample()].
#' @return A `data.frame` with one row per (fold, repeat, learner)
#'   combination.
#' @seealso [fit_resample()], [audit_perm_gap()]
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   subject = rep(1:6, each = 2),
#'   outcome = factor(rep(c("a","b"), 6), levels = c("a","b")),
#'   x1 = rnorm(12), x2 = rnorm(12)
#' )
#' splits <- make_split_plan(df, outcome = "outcome",
#'                           mode = "subject_grouped", group = "subject", v = 3)
#' \dontrun{
#'   fit <- fit_resample(df, outcome = "outcome", splits = splits,
#'                       learner = parsnip::logistic_reg() |>
#'                                 parsnip::set_engine("glm"),
#'                       metrics = "auc")
#'   fit_metrics(fit)
#' }
#' @export
fit_metrics <- function(fit) {
  if (!methods::is(fit, "LeakFit")) {
    stop("`fit` must be a LeakFit object (got: ",
         paste(class(fit), collapse = "/"), ").",
         call. = FALSE)
  }
  fit@metrics
}


# ---- LeakAudit ------------------------------------------------------

#' Permutation-gap test results from a LeakAudit
#'
#' Returns the permutation-gap test data frame stored in a `LeakAudit`
#' object. Columns include the observed metric, permuted-null mean and
#' SD, gap, z-score, and permutation p-value.
#'
#' @param audit A `LeakAudit` object returned by [audit_leakage()].
#' @return A `data.frame` with one row per (mechanism class, repeat),
#'   summarising the permutation-gap test.
#' @seealso [audit_leakage()], [audit_target_assoc()]
#' @export
audit_perm_gap <- function(audit) {
  if (!methods::is(audit, "LeakAudit")) {
    stop("`audit` must be a LeakAudit object (got: ",
         paste(class(audit), collapse = "/"), ").",
         call. = FALSE)
  }
  audit@permutation_gap
}

#' Batch / study association results from a LeakAudit
#'
#' Returns the batch/study chi-squared association data frame stored in
#' a `LeakAudit` object. Columns include the metadata column, repeat,
#' chi-squared statistic, degrees of freedom, p-value, and Cramer's V.
#'
#' @param audit A `LeakAudit` object returned by [audit_leakage()].
#' @return A `data.frame` (possibly empty) with one row per
#'   (metadata column, repeat).
#' @seealso [audit_leakage()], [audit_perm_gap()]
#' @export
audit_batch_assoc <- function(audit) {
  if (!methods::is(audit, "LeakAudit")) {
    stop("`audit` must be a LeakAudit object (got: ",
         paste(class(audit), collapse = "/"), ").",
         call. = FALSE)
  }
  audit@batch_assoc
}

#' Target leakage scan results from a LeakAudit
#'
#' Returns the per-feature target-association data frame stored in a
#' `LeakAudit` object. Columns include the feature name, association
#' value, score, and a logical flag indicating whether the feature
#' exceeded the configured `target_threshold`.
#'
#' @param audit A `LeakAudit` object returned by [audit_leakage()].
#' @return A `data.frame` (possibly empty) with one row per scanned
#'   feature.
#' @seealso [audit_leakage()], [audit_perm_gap()]
#' @export
audit_target_assoc <- function(audit) {
  if (!methods::is(audit, "LeakAudit")) {
    stop("`audit` must be a LeakAudit object (got: ",
         paste(class(audit), collapse = "/"), ").",
         call. = FALSE)
  }
  audit@target_assoc
}

#' Near-duplicate sample pairs from a LeakAudit
#'
#' Returns the near-duplicate detection data frame stored in a
#' `LeakAudit` object. Each row is one (i, j) sample pair whose
#' similarity exceeds the configured threshold.
#'
#' @param audit A `LeakAudit` object returned by [audit_leakage()].
#' @return A `data.frame` (possibly empty) with one row per detected
#'   near-duplicate pair.
#' @seealso [audit_leakage()]
#' @export
audit_duplicates <- function(audit) {
  if (!methods::is(audit, "LeakAudit")) {
    stop("`audit` must be a LeakAudit object (got: ",
         paste(class(audit), collapse = "/"), ").",
         call. = FALSE)
  }
  audit@duplicates
}

#' Auxiliary information list from a LeakAudit
#'
#' Returns the auxiliary `info` list stored in a `LeakAudit` object.
#' This list carries non-tabular audit components such as
#' multivariate-target-scan results, configuration flags, and
#' permutation-test diagnostics.
#'
#' @param audit A `LeakAudit` object returned by [audit_leakage()].
#' @return A named list of supplementary audit components.
#' @seealso [audit_leakage()]
#' @export
audit_info <- function(audit) {
  if (!methods::is(audit, "LeakAudit")) {
    stop("`audit` must be a LeakAudit object (got: ",
         paste(class(audit), collapse = "/"), ").",
         call. = FALSE)
  }
  audit@info
}


# ---- LeakDeltaLSI ---------------------------------------------------

#' Raw delta-metric estimate from a LeakDeltaLSI
#'
#' Returns the arithmetic mean of the per-repeat raw metric differences
#' (leaky minus guarded) stored in a `LeakDeltaLSI` object.
#'
#' @param dlsi A `LeakDeltaLSI` object returned by [delta_lsi()].
#' @return A length-one numeric scalar.
#' @seealso [delta_lsi()], [dlsi_robust()], [dlsi_ci()]
#' @export
dlsi_metric <- function(dlsi) {
  if (!methods::is(dlsi, "LeakDeltaLSI")) {
    stop("`dlsi` must be a LeakDeltaLSI object (got: ",
         paste(class(dlsi), collapse = "/"), ").",
         call. = FALSE)
  }
  dlsi@delta_metric
}

#' Huber-robust delta_lsi point estimate from a LeakDeltaLSI
#'
#' Returns the Huber-robust point estimate of the per-repeat delta
#' values stored in a `LeakDeltaLSI` object.
#'
#' @param dlsi A `LeakDeltaLSI` object returned by [delta_lsi()].
#' @return A length-one numeric scalar.
#' @seealso [delta_lsi()], [dlsi_metric()], [dlsi_ci()]
#' @export
dlsi_robust <- function(dlsi) {
  if (!methods::is(dlsi, "LeakDeltaLSI")) {
    stop("`dlsi` must be a LeakDeltaLSI object (got: ",
         paste(class(dlsi), collapse = "/"), ").",
         call. = FALSE)
  }
  dlsi@delta_lsi
}

#' BCa confidence interval from a LeakDeltaLSI
#'
#' Returns the bias-corrected and accelerated (BCa) bootstrap confidence
#' interval stored in a `LeakDeltaLSI` object. By default returns the
#' interval for the Huber-robust delta estimate; set
#' `which = "metric"` to return the interval for the raw metric
#' difference instead.
#'
#' @param dlsi A `LeakDeltaLSI` object returned by [delta_lsi()].
#' @param which Either `"robust"` (default) for the Huber estimate's
#'   confidence interval, or `"metric"` for the raw arithmetic mean's
#'   confidence interval.
#' @return A length-two numeric vector `c(lower, upper)`. Returns
#'   `c(NA_real_, NA_real_)` when the interval is not computed (e.g.,
#'   the inference tier did not include CIs).
#' @seealso [delta_lsi()], [dlsi_robust()], [dlsi_metric()]
#' @export
dlsi_ci <- function(dlsi, which = c("robust", "metric")) {
  if (!methods::is(dlsi, "LeakDeltaLSI")) {
    stop("`dlsi` must be a LeakDeltaLSI object (got: ",
         paste(class(dlsi), collapse = "/"), ").",
         call. = FALSE)
  }
  which <- match.arg(which)
  if (identical(which, "robust"))   dlsi@delta_lsi_ci  else dlsi@delta_metric_ci
}

#' Sign-flip p-value from a LeakDeltaLSI
#'
#' Returns the paired sign-flip randomization-test p-value stored in a
#' `LeakDeltaLSI` object.
#'
#' @param dlsi A `LeakDeltaLSI` object returned by [delta_lsi()].
#' @return A length-one numeric scalar in `[0, 1]`, or `NA_real_` when
#'   the inference tier did not include hypothesis testing.
#' @seealso [delta_lsi()], [dlsi_tier()]
#' @export
dlsi_p_value <- function(dlsi) {
  if (!methods::is(dlsi, "LeakDeltaLSI")) {
    stop("`dlsi` must be a LeakDeltaLSI object (got: ",
         paste(class(dlsi), collapse = "/"), ").",
         call. = FALSE)
  }
  dlsi@p_value
}

#' Inference tier from a LeakDeltaLSI
#'
#' Returns the inference tier label stored in a `LeakDeltaLSI` object.
#' Possible values are `"A_full_inference"` (`R_eff >= 20`),
#' `"B_signflip_ci"` (`R_eff >= 10`), `"C_signflip"`
#' (`R_eff >= 5`), or `"D_insufficient"` (`R_eff < 5`).
#'
#' @param dlsi A `LeakDeltaLSI` object returned by [delta_lsi()].
#' @return A length-one character string giving the tier label.
#' @seealso [delta_lsi()], [dlsi_R_eff()]
#' @export
dlsi_tier <- function(dlsi) {
  if (!methods::is(dlsi, "LeakDeltaLSI")) {
    stop("`dlsi` must be a LeakDeltaLSI object (got: ",
         paste(class(dlsi), collapse = "/"), ").",
         call. = FALSE)
  }
  dlsi@tier
}

#' Effective number of paired repeats from a LeakDeltaLSI
#'
#' Returns the effective number of paired repeats `R_eff` stored in a
#' `LeakDeltaLSI` object. This is the count of repeats that contribute
#' to the inference; it equals the smaller of the leaky and guarded
#' fits' repeat counts when the comparison is paired.
#'
#' @param dlsi A `LeakDeltaLSI` object returned by [delta_lsi()].
#' @return A length-one integer.
#' @seealso [delta_lsi()], [dlsi_tier()]
#' @export
dlsi_R_eff <- function(dlsi) {
  if (!methods::is(dlsi, "LeakDeltaLSI")) {
    stop("`dlsi` must be a LeakDeltaLSI object (got: ",
         paste(class(dlsi), collapse = "/"), ").",
         call. = FALSE)
  }
  dlsi@R_eff
}

#' Per-repeat metric data frames from a LeakDeltaLSI
#'
#' Returns the per-repeat metric data frame for one of the two
#' pipelines stored in a `LeakDeltaLSI` object. The naive (or leaky)
#' pipeline's repeats are returned by default; set
#' `which = "guarded"` to return the guarded pipeline's repeats.
#'
#' @param dlsi A `LeakDeltaLSI` object returned by [delta_lsi()].
#' @param which Either `"naive"` (default) for the naive/leaky
#'   pipeline's per-repeat data frame, or `"guarded"` for the guarded
#'   pipeline's per-repeat data frame.
#' @return A `data.frame` with one row per repeat.
#' @seealso [delta_lsi()]
#' @export
dlsi_repeats <- function(dlsi, which = c("naive", "guarded")) {
  if (!methods::is(dlsi, "LeakDeltaLSI")) {
    stop("`dlsi` must be a LeakDeltaLSI object (got: ",
         paste(class(dlsi), collapse = "/"), ").",
         call. = FALSE)
  }
  which <- match.arg(which)
  if (identical(which, "naive"))   dlsi@repeats_naive  else dlsi@repeats_guarded
}
