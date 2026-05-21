## ---------------------------------------------------------------------
## R/accessors.R
## ---------------------------------------------------------------------
## Public S4 accessor methods for the bioLeak result classes.
##
## These accessors expose slot data through proper S4 generics so that
## downstream code can read the components of LeakFit, LeakAudit, and
## LeakDeltaLSI objects without reaching into S4 internals via `@`, and
## so that the accessors are visible through methods(class = "<Class>").
## Slot definitions are unchanged; these accessors are purely additive.
##
## In bioLeak 0.3.7 the accessors were defined as plain functions; in
## 0.3.8 they have been converted to S4 generic + method pairs.
## Function-form calls (e.g., fit_metrics(fit)) continue to work
## unchanged because the function name is the generic; only the
## dispatch mechanism has changed.
## ---------------------------------------------------------------------

#' @include classes.R
NULL


# ============================================================
# LeakFit accessors
# ============================================================

#' Per-fold metric data frame from a LeakFit
#'
#' Returns the per-fold metric data frame stored in a [`LeakFit`] object.
#' Each row is one (fold, repeat, learner) combination; columns include
#' the requested metric values such as `auc` and any task-specific
#' performance scores.
#'
#' Implemented as an S4 generic with a method for [`LeakFit`]; visible
#' via `methods(class = "LeakFit")`.
#'
#' @param fit A [`LeakFit`] object returned by [fit_resample()].
#' @return A `data.frame` with one row per (fold, repeat, learner)
#'   combination.
#' @seealso [LeakClasses], [fit_resample()], [audit_perm_gap()]
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
setGeneric("fit_metrics",
           function(fit) standardGeneric("fit_metrics"))

#' @rdname fit_metrics
#' @export
setMethod("fit_metrics", "LeakFit", function(fit) fit@metrics)


# ============================================================
# LeakAudit accessors
# ============================================================

#' Permutation-gap test results from a LeakAudit
#'
#' Returns the permutation-gap test data frame stored in a [`LeakAudit`]
#' object. Columns include the observed metric, permuted-null mean and
#' SD, gap, z-score, and permutation p-value.
#'
#' Implemented as an S4 generic with a method for [`LeakAudit`]; visible
#' via `methods(class = "LeakAudit")`.
#'
#' @param audit A [`LeakAudit`] object returned by [audit_leakage()].
#' @return A `data.frame` with one row per (mechanism class, repeat),
#'   summarising the permutation-gap test.
#' @seealso [LeakClasses], [audit_leakage()], [audit_target_assoc()]
#' @export
setGeneric("audit_perm_gap",
           function(audit) standardGeneric("audit_perm_gap"))

#' @rdname audit_perm_gap
#' @export
setMethod("audit_perm_gap", "LeakAudit",
          function(audit) audit@permutation_gap)


#' Batch / study association results from a LeakAudit
#'
#' Returns the batch/study chi-squared association data frame stored in
#' a [`LeakAudit`] object. Columns include the metadata column, repeat,
#' chi-squared statistic, degrees of freedom, p-value, and Cramer's V.
#'
#' Implemented as an S4 generic with a method for [`LeakAudit`]; visible
#' via `methods(class = "LeakAudit")`.
#'
#' @param audit A [`LeakAudit`] object returned by [audit_leakage()].
#' @return A `data.frame` with one row per (metadata column, repeat).
#' @seealso [LeakClasses], [audit_leakage()]
#' @export
setGeneric("audit_batch_assoc",
           function(audit) standardGeneric("audit_batch_assoc"))

#' @rdname audit_batch_assoc
#' @export
setMethod("audit_batch_assoc", "LeakAudit",
          function(audit) audit@batch_assoc)


#' Target-leakage scan results from a LeakAudit
#'
#' Returns the target-association data frame stored in a [`LeakAudit`]
#' object. Each row is one predictor with its association score
#' (rescaled AUC; `|AUC - 0.5| * 2`), threshold-based flag, and
#' (where applicable) p-value.
#'
#' Implemented as an S4 generic with a method for [`LeakAudit`]; visible
#' via `methods(class = "LeakAudit")`.
#'
#' @param audit A [`LeakAudit`] object returned by [audit_leakage()].
#' @return A `data.frame` with one row per predictor.
#' @seealso [LeakClasses], [audit_leakage()]
#' @export
setGeneric("audit_target_assoc",
           function(audit) standardGeneric("audit_target_assoc"))

#' @rdname audit_target_assoc
#' @export
setMethod("audit_target_assoc", "LeakAudit",
          function(audit) audit@target_assoc)


#' Near-duplicate pairs from a LeakAudit
#'
#' Returns the near-duplicate sample pairs data frame stored in a
#' [`LeakAudit`] object. Each row is a unique (row_a, row_b) pair above
#' the configured cosine-similarity threshold that crossed train/test
#' partitions in at least one fold.
#'
#' Implemented as an S4 generic with a method for [`LeakAudit`]; visible
#' via `methods(class = "LeakAudit")`.
#'
#' @param audit A [`LeakAudit`] object returned by [audit_leakage()].
#' @return A `data.frame` with one row per detected near-duplicate pair.
#' @seealso [LeakClasses], [audit_leakage()]
#' @export
setGeneric("audit_duplicates",
           function(audit) standardGeneric("audit_duplicates"))

#' @rdname audit_duplicates
#' @export
setMethod("audit_duplicates", "LeakAudit",
          function(audit) audit@duplicates)


#' Auxiliary information from a LeakAudit
#'
#' Returns the auxiliary information list stored in a [`LeakAudit`]
#' object. The list typically contains multivariate-target-scan
#' results, configuration flags, permutation-test diagnostics, and
#' provenance metadata.
#'
#' Implemented as an S4 generic with a method for [`LeakAudit`]; visible
#' via `methods(class = "LeakAudit")`.
#'
#' @param audit A [`LeakAudit`] object returned by [audit_leakage()].
#' @return A named `list`.
#' @seealso [LeakClasses], [audit_leakage()]
#' @export
setGeneric("audit_info",
           function(audit) standardGeneric("audit_info"))

#' @rdname audit_info
#' @export
setMethod("audit_info", "LeakAudit",
          function(audit) audit@info)


# ============================================================
# LeakDeltaLSI accessors
# ============================================================

#' Raw delta-metric estimate from a LeakDeltaLSI
#'
#' Returns the arithmetic mean of the per-repeat raw metric differences
#' (leaky minus guarded) stored in a [`LeakDeltaLSI`] object.
#'
#' Implemented as an S4 generic with a method for [`LeakDeltaLSI`];
#' visible via `methods(class = "LeakDeltaLSI")`.
#'
#' @param dlsi A [`LeakDeltaLSI`] object returned by [delta_lsi()].
#' @return A length-one numeric scalar.
#' @seealso [LeakClasses], [delta_lsi()], [dlsi_robust()], [dlsi_ci()]
#' @export
setGeneric("dlsi_metric",
           function(dlsi) standardGeneric("dlsi_metric"))

#' @rdname dlsi_metric
#' @export
setMethod("dlsi_metric", "LeakDeltaLSI",
          function(dlsi) dlsi@delta_metric)


#' Huber-robust delta_lsi point estimate from a LeakDeltaLSI
#'
#' Returns the Huber-robust point estimate of the per-repeat delta
#' values stored in a [`LeakDeltaLSI`] object.
#'
#' Implemented as an S4 generic with a method for [`LeakDeltaLSI`];
#' visible via `methods(class = "LeakDeltaLSI")`.
#'
#' @param dlsi A [`LeakDeltaLSI`] object returned by [delta_lsi()].
#' @return A length-one numeric scalar.
#' @seealso [LeakClasses], [delta_lsi()], [dlsi_metric()], [dlsi_ci()]
#' @export
setGeneric("dlsi_robust",
           function(dlsi) standardGeneric("dlsi_robust"))

#' @rdname dlsi_robust
#' @export
setMethod("dlsi_robust", "LeakDeltaLSI",
          function(dlsi) dlsi@delta_lsi)


#' BCa confidence interval from a LeakDeltaLSI
#'
#' Returns the bias-corrected and accelerated (BCa) bootstrap
#' confidence interval stored in a [`LeakDeltaLSI`] object. By default
#' returns the interval for the Huber-robust delta estimate; set
#' `which = "metric"` to return the interval for the raw metric
#' difference instead.
#'
#' Implemented as an S4 generic with a method for [`LeakDeltaLSI`];
#' visible via `methods(class = "LeakDeltaLSI")`.
#'
#' @param dlsi A [`LeakDeltaLSI`] object returned by [delta_lsi()].
#' @param which Either `"robust"` (default) for the Huber estimate's
#'   confidence interval, or `"metric"` for the raw arithmetic mean's
#'   confidence interval.
#' @return A length-two numeric vector `c(lower, upper)`. Returns
#'   `c(NA_real_, NA_real_)` when the interval is not computed (for
#'   example, when the inference tier did not include CIs).
#' @seealso [LeakClasses], [delta_lsi()], [dlsi_robust()],
#'   [dlsi_metric()]
#' @export
setGeneric("dlsi_ci",
           function(dlsi, which = c("robust", "metric"))
             standardGeneric("dlsi_ci"))

#' @rdname dlsi_ci
#' @export
setMethod("dlsi_ci", "LeakDeltaLSI",
          function(dlsi, which = c("robust", "metric")) {
            which <- match.arg(which)
            if (identical(which, "robust")) dlsi@delta_lsi_ci
            else                            dlsi@delta_metric_ci
          })


#' Sign-flip p-value from a LeakDeltaLSI
#'
#' Returns the paired sign-flip randomization-test p-value stored in a
#' [`LeakDeltaLSI`] object.
#'
#' Implemented as an S4 generic with a method for [`LeakDeltaLSI`];
#' visible via `methods(class = "LeakDeltaLSI")`.
#'
#' @param dlsi A [`LeakDeltaLSI`] object returned by [delta_lsi()].
#' @return A length-one numeric scalar in `[0, 1]`, or `NA_real_` when
#'   the inference tier did not include hypothesis testing.
#' @seealso [LeakClasses], [delta_lsi()], [dlsi_tier()]
#' @export
setGeneric("dlsi_p_value",
           function(dlsi) standardGeneric("dlsi_p_value"))

#' @rdname dlsi_p_value
#' @export
setMethod("dlsi_p_value", "LeakDeltaLSI",
          function(dlsi) dlsi@p_value)


#' Inference tier from a LeakDeltaLSI
#'
#' Returns the inference tier label stored in a [`LeakDeltaLSI`] object.
#' Possible values are `"A_full_inference"` (`R_eff >= 20`),
#' `"B_signflip_ci"` (`R_eff >= 10`), `"C_signflip"` (`R_eff >= 5`),
#' or `"D_insufficient"` (`R_eff < 5`).
#'
#' Implemented as an S4 generic with a method for [`LeakDeltaLSI`];
#' visible via `methods(class = "LeakDeltaLSI")`.
#'
#' @param dlsi A [`LeakDeltaLSI`] object returned by [delta_lsi()].
#' @return A length-one character string giving the tier label.
#' @seealso [LeakClasses], [delta_lsi()], [dlsi_R_eff()]
#' @export
setGeneric("dlsi_tier",
           function(dlsi) standardGeneric("dlsi_tier"))

#' @rdname dlsi_tier
#' @export
setMethod("dlsi_tier", "LeakDeltaLSI",
          function(dlsi) dlsi@tier)


#' Effective number of paired repeats from a LeakDeltaLSI
#'
#' Returns the effective number of paired repeats `R_eff` stored in a
#' [`LeakDeltaLSI`] object. This is the count of repeats that
#' contribute to the inference; it equals the smaller of the leaky and
#' guarded fits' repeat counts when the comparison is paired.
#'
#' Implemented as an S4 generic with a method for [`LeakDeltaLSI`];
#' visible via `methods(class = "LeakDeltaLSI")`.
#'
#' @param dlsi A [`LeakDeltaLSI`] object returned by [delta_lsi()].
#' @return A length-one integer.
#' @seealso [LeakClasses], [delta_lsi()], [dlsi_tier()]
#' @export
setGeneric("dlsi_R_eff",
           function(dlsi) standardGeneric("dlsi_R_eff"))

#' @rdname dlsi_R_eff
#' @export
setMethod("dlsi_R_eff", "LeakDeltaLSI",
          function(dlsi) dlsi@R_eff)


#' Per-repeat metric data frames from a LeakDeltaLSI
#'
#' Returns the per-repeat metric data frame for one of the two
#' pipelines stored in a [`LeakDeltaLSI`] object. The naive (or leaky)
#' pipeline's repeats are returned by default; set
#' `which = "guarded"` to return the guarded pipeline's repeats.
#'
#' Implemented as an S4 generic with a method for [`LeakDeltaLSI`];
#' visible via `methods(class = "LeakDeltaLSI")`.
#'
#' @param dlsi A [`LeakDeltaLSI`] object returned by [delta_lsi()].
#' @param which Either `"naive"` (default) for the naive/leaky
#'   pipeline's per-repeat data frame, or `"guarded"` for the guarded
#'   pipeline's per-repeat data frame.
#' @return A `data.frame` with one row per repeat.
#' @seealso [LeakClasses], [delta_lsi()]
#' @export
setGeneric("dlsi_repeats",
           function(dlsi, which = c("naive", "guarded"))
             standardGeneric("dlsi_repeats"))

#' @rdname dlsi_repeats
#' @export
setMethod("dlsi_repeats", "LeakDeltaLSI",
          function(dlsi, which = c("naive", "guarded")) {
            which <- match.arg(which)
            if (identical(which, "naive")) dlsi@repeats_naive
            else                           dlsi@repeats_guarded
          })
