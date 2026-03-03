# Structured condition constructors for bioLeak ---------------------------------
#
# These use base-R condition construction (no rlang dependency).
# Users can catch specific error categories via:
#   tryCatch(..., bioLeak_input_error = function(e) ...)

#' @keywords internal
.bio_stop <- function(message, class = NULL, ..., call = NULL) {
  cond <- structure(
    class = c(class, "bioLeak_error", "error", "condition"),
    list(message = message, call = call, ...)
  )
  stop(cond)
}

#' @keywords internal
.bio_warn <- function(message, class = NULL, ..., call = NULL) {
  cond <- structure(
    class = c(class, "bioLeak_warning", "warning", "condition"),
    list(message = message, call = call, ...)
  )
  warning(cond)
}

# Condition class taxonomy:
#
# Errors:
#   bioLeak_input_error      — Bad argument type/value/missing
#   bioLeak_column_error     — Missing or invalid metadata column
#   bioLeak_split_error      — Split generation failures
#   bioLeak_overlap_error    — Overlap invariant violations
#   bioLeak_fit_error        — Model fitting failures
#   bioLeak_audit_error      — Audit-specific failures
#   bioLeak_dependency_error — Missing optional package
#
# Warnings:
#   bioLeak_validation_warning — Recipe/workflow leakage guardrail
#   bioLeak_fold_warning       — Fold-level issues
#   bioLeak_fallback_warning   — Fallback behavior triggered
