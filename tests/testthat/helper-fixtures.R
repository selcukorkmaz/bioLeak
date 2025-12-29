make_class_df <- function(n = 20) {
  stopifnot(n >= 4)
  data.frame(
    subject = rep(seq_len(ceiling(n / 2)), each = 2)[seq_len(n)],
    batch = rep(letters[1:3], length.out = n),
    study = rep(LETTERS[1:4], length.out = n),
    time = seq_len(n),
    outcome = rep(c(0, 1), length.out = n),
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
          stats::glm(y ~ ., data = df, family = stats::binomial(), weights = weights)
        } else {
          stats::lm(y ~ ., data = df, weights = weights)
        }
      },
      predict = function(object, newdata, task, ...) {
        df <- as.data.frame(newdata, check.names = FALSE)
        if (identical(task, "binomial")) {
          as.numeric(stats::predict(object, newdata = df, type = "response"))
        } else {
          as.numeric(stats::predict(object, newdata = df))
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
  capture.output(out <- fit_resample(...))
  out
}

with_temp_plot_device <- function(expr) {
  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit(grDevices::dev.off(), add = TRUE)
  eval.parent(substitute(expr))
}
