# Shared helpers for verify-change probes. Source from the package root.
# Each check prints one line: PASS / FAIL / WARN, an ID, and the evidence.
# IDs: CTRL-* are positive controls (must always pass; a failure means the
# probe itself is broken). BL-* match entries in dev/OPEN-ISSUES.md.
# INV-* are general invariants with no open issue.

.probe <- new.env()
.probe$fail <- 0L
.probe$warn <- 0L

check <- function(id, ok, detail, level = c("FAIL", "WARN")) {
  level <- match.arg(level)
  tag <- if (isTRUE(ok)) "PASS" else level
  if (identical(tag, "FAIL")) .probe$fail <- .probe$fail + 1L
  if (identical(tag, "WARN")) .probe$warn <- .probe$warn + 1L
  cat(sprintf("%-4s [%s] %s\n", tag, id, detail))
  invisible(isTRUE(ok))
}

# Capture warnings and errors from an expression without stopping the probe.
capture <- function(expr) {
  warns <- character()
  val <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(list(msg = conditionMessage(e)), class = "probe_error")),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = val, warnings = warns, error = if (inherits(val, "probe_error")) val$msg else NULL)
}

# Upper limit for an empirical type-I error from B null replicates.
type1_limit <- function(alpha, B) alpha + 3 * sqrt(alpha * (1 - alpha) / B)

probe_finish <- function() {
  cat(sprintf("\n%d FAIL, %d WARN\n", .probe$fail, .probe$warn))
  quit(save = "no", status = if (.probe$fail > 0L) 1L else 0L)
}

load_dev <- function(path = ".") {
  suppressMessages(pkgload::load_all(path, quiet = TRUE, export_all = TRUE))
}
