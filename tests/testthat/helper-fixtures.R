make_class_df <- function(n = 20) {
  stopifnot(n >= 4)
  data.frame(
    subject = rep(seq_len(ceiling(n / 2)), each = 2)[seq_len(n)],
    batch = rep(letters[1:3], length.out = n),
    study = rep(LETTERS[1:4], length.out = n),
    time = seq_len(n),
    outcome = factor(rep(c(0, 1), length.out = n), levels = c(0, 1)),
    x1 = rnorm(n),
    x2 = rnorm(n),
    stringsAsFactors = FALSE
  )
}

make_regression_df <- function(n = 20) {
  stopifnot(n >= 4)
  data.frame(
    subject = rep(seq_len(ceiling(n / 2)), each = 2)[seq_len(n)],
    batch = rep(letters[1:3], length.out = n),
    study = rep(LETTERS[1:4], length.out = n),
    time = seq_len(n),
    y = rnorm(n),
    x1 = rnorm(n),
    x2 = rnorm(n),
    stringsAsFactors = FALSE
  )
}

make_custom_learners <- function() {
  list(
    glm = list(
      fit = function(x, y, task, weights, ...) {
        df <- data.frame(y = y, x, check.names = FALSE)
        if (identical(task, "binomial")) {
          suppressWarnings(stats::glm(y ~ ., data = df,
                                      family = stats::binomial(), weights = weights))
        } else {
          suppressWarnings(stats::lm(y ~ ., data = df, weights = weights))
        }
      },
      predict = function(object, newdata, task, ...) {
        df <- as.data.frame(newdata, check.names = FALSE)
        if (identical(task, "binomial")) {
          as.numeric(suppressWarnings(stats::predict(object,
                                                     newdata = df,
                                                     type = "response")))
        } else {
          as.numeric(suppressWarnings(stats::predict(object, newdata = df)))
        }
      }
    )
  )
}

make_splits_quiet <- function(...) {
  make_splits(..., progress = FALSE)
}

fit_resample_quiet <- function(...) {
  out <- NULL
  suppress_patterns <- c(
    "glm.fit",
    "algorithm did not converge",
    "fitted probabilities numerically 0 or 1 occurred",
    "rank-deficient",
    "non-estim"
  )
  seen_warnings <- character()
  capture.output({
    out <- withCallingHandlers(
      fit_resample(...),
      warning = function(w) {
        seen_warnings <<- c(seen_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  })
  if (length(seen_warnings)) {
    keep <- !vapply(seen_warnings, function(msg) {
      msg <- tolower(msg)
      any(vapply(suppress_patterns, function(pat) {
        grepl(pat, msg, fixed = TRUE)
      }, logical(1)))
    }, logical(1))
    if (any(keep)) {
      for (msg in seen_warnings[keep]) {
        warning(msg, call. = FALSE)
      }
    }
  }
  out
}

with_temp_plot_device <- function(expr) {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit(grDevices::dev.off(), add = TRUE)
  eval.parent(substitute(expr))
}

expect_warning_match <- function(expr, pattern, all = FALSE) {
  warnings <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  if (!length(warnings)) {
    testthat::expect_true(FALSE, info = "Expected warning, but none were emitted.")
  } else if (isTRUE(all)) {
    testthat::expect_true(all(grepl(pattern, warnings)))
  } else {
    testthat::expect_true(any(grepl(pattern, warnings)))
  }
  value
}

