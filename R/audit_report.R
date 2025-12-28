#' Render an HTML audit report
#'
#' Generates a one-click HTML report for a [LeakAudit] object. You can also pass
#' a [LeakFit] object, in which case [audit_leakage()] is run first.
#'
#' @param audit LeakAudit or LeakFit object.
#' @param output_file Output HTML filename.
#' @param output_dir Directory to write the report.
#' @param quiet Logical, pass through to rmarkdown render.
#' @param open Logical, open the report in a browser after rendering.
#' @param ... Additional arguments passed to [audit_leakage()] when `audit` is a LeakFit.
#' @return Path to the generated HTML report.
#' @examples
#' \dontrun{
#' # fit <- fit_resample(...)
#' audit <- audit_leakage(fit, metric = "auc", B = 50)
#' audit_report(audit, output_dir = ".")
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
