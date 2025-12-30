# Platform-safe symbols for console output.
.bio_symbol <- function(name) {
  use_unicode <- getOption("bioLeak.use_unicode", NA)
  if (is.na(use_unicode)) {
    utf8_ok <- FALSE
    if (requireNamespace("cli", quietly = TRUE)) {
      utf8_ok <- isTRUE(cli::is_utf8_output())
    } else {
      loc <- tryCatch(l10n_info(), error = function(e) NULL)
      if (!is.null(loc) && isTRUE(loc[["UTF-8"]])) utf8_ok <- TRUE
    }
    if (.Platform$OS.type == "windows") utf8_ok <- FALSE
    use_unicode <- utf8_ok
  } else {
    use_unicode <- isTRUE(use_unicode)
  }

  if (!use_unicode) {
    return(switch(name,
      pm = "+/-",
      chi_sq = "Chi^2",
      ge = ">=",
      check = "OK:",
      warn = "WARNING:",
      name
    ))
  }

  switch(name,
    pm = "\u00b1",
    chi_sq = "\u03c7\u00b2",
    ge = "\u2265",
    check = "\u2713",
    warn = "\u26a0",
    name
  )
}
