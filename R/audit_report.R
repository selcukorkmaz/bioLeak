#' Render an HTML audit report
#'
#' Creates an HTML report that summarizes a leakage audit for a resampled model.
#' The report is built from a [LeakAudit] (or created from a [LeakFit]) and
#' presents: cross-validated metric summaries, a label-permutation association
#' test of the chosen performance metric (auto-refit when refit data are
#' available; otherwise fixed predictions), batch or study association tests
#' between metadata and predictions, confounder sensitivity plots, calibration
#' checks for binomial tasks, a target leakage scan based on feature-outcome
#' similarity (with multivariate scan enabled by default for supported tasks),
#' and duplicate
#' detection across training and test folds. The
#' output is a self-contained HTML file with tables and plots for these checks
#' plus the audit parameters used.
#'
#' The report does not refit models or reprocess data unless `perm_refit` triggers
#' refitting (`TRUE` or `"auto"` with a valid `perm_refit_spec`); it otherwise
#' inspects the predictions and metadata stored in the input.
#' Results are conditional on the
#' provided splits, selected metric, and any feature matrix supplied to
#' [audit_leakage()]. The univariate target leakage scan can miss multivariate
#' proxies, interaction leakage, or features not included in `X_ref`; the
#' multivariate scan (enabled by default for supported tasks) adds a model-based
#' check but still only uses features in `X_ref`. A
#' non-significant result does not prove the absence of leakage, especially with
#' small `B` or incomplete metadata. Rendering requires the
#' `rmarkdown` package and `ggplot2` for plots.
#'
#' @param audit A [LeakAudit] object from [audit_leakage()] or a [LeakFit] object
#'   from [fit_resample()]. If a [LeakAudit] is supplied, the report uses its
#'   stored results verbatim. If a [LeakFit] is supplied, `audit_report()` first
#'   computes a new audit via [audit_leakage(...)]; the fit must contain
#'   predictions and split metadata. When multiple learners were fit, pass a
#'   `learner` argument via `...` to select a single model.
#' @param output_file Character scalar. File name for the HTML report. Defaults
#'   to `"bioLeak_audit_report.html"`. If a relative name is provided, it is
#'   created inside `output_dir`. Changing this value only changes the file name,
#'   not the audit content.
#' @param output_dir Character scalar. Directory path where the report is
#'   written. Defaults to [tempdir()]. The directory must exist or be creatable
#'   by `rmarkdown::render()`. Changing this value only changes the output
#'   location.
#' @param quiet Logical scalar passed to `rmarkdown::render()`. Defaults to
#'   `TRUE`. When `FALSE`, knitting output and warnings are printed to the
#'   console. This does not change audit results.
#' @param open Logical scalar. Defaults to `FALSE`. When `TRUE`, opens the
#'   generated report in a browser via [utils::browseURL()]. This does not change
#'   the report contents.
#' @param ... Additional named arguments forwarded to [audit_leakage()] only
#'   when `audit` is a [LeakFit]. These control how the audit is computed and
#'   therefore change the report. Typical examples include `metric` (character),
#'   `B` (integer), `perm_stratify` (logical), `batch_cols` (character vector),
#'   `X_ref` (matrix/data.frame), `sim_method` (character), and
#'   `duplicate_scope` (character). When omitted, [audit_leakage()] defaults are
#'   used. Ignored when `audit` is already a [LeakAudit].
#' @return Path to the generated HTML report.
#' @examples
#' \dontrun{
#' set.seed(1)
#' df <- data.frame(
#'   subject = rep(1:6, each = 2),
#'   outcome = factor(rep(c(0, 1), 6)),
#'   x1 = rnorm(12),
#'   x2 = rnorm(12)
#' )
#'
#' splits <- make_splits(df, outcome = "outcome",
#'                       mode = "subject_grouped", group = "subject",
#'                       v = 3, progress = FALSE)
#'
#' custom <- list(
#'   glm = list(
#'     fit = function(x, y, task, weights, ...) {
#'       stats::glm(y ~ ., data = data.frame(y = y, x),
#'                  family = stats::binomial(), weights = weights)
#'     },
#'     predict = function(object, newdata, task, ...) {
#'       as.numeric(stats::predict(object,
#'                                 newdata = as.data.frame(newdata),
#'                                 type = "response"))
#'     }
#'   )
#' )
#'
#' fit <- fit_resample(df, outcome = "outcome", splits = splits,
#'                     learner = "glm", custom_learners = custom,
#'                     metrics = "auc", refit = FALSE, seed = 1)
#'
#' audit <- audit_leakage(fit, metric = "auc", B = 5, perm_stratify = FALSE)
#'
#' if (requireNamespace("rmarkdown", quietly = TRUE) &&
#'     requireNamespace("ggplot2", quietly = TRUE)) {
#'   out_file <- audit_report(audit, output_dir = tempdir(), quiet = TRUE)
#'   out_file
#' }
#' }
#' @export
audit_report <- function(audit,
                         output_file = "bioLeak_audit_report.html",
                         output_dir = tempdir(),
                         quiet = TRUE,
                         open = FALSE,
                         ...) {
  if (inherits(audit, "LeakFit")) {
    audit <- audit_leakage(audit, ...)
  }
  if (!inherits(audit, "LeakAudit")) {
    stop("audit must be a LeakAudit or LeakFit object.")
  }
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Install 'rmarkdown' to render the HTML report.", call. = FALSE)
  }
  template <- system.file("reports", "audit_report.Rmd", package = "bioLeak")
  if (!nzchar(template)) {
    stop("Report template not found. Reinstall bioLeak.", call. = FALSE)
  }

  out <- rmarkdown::render(
    input = template,
    output_file = output_file,
    output_dir = output_dir,
    params = list(audit = audit),
    quiet = quiet,
    envir = new.env(parent = globalenv())
  )

  if (isTRUE(open)) {
    utils::browseURL(out)
  }
  out
}
