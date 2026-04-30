## ---------------------------------------------------------------------
## test-predict-guardfit.R
## ---------------------------------------------------------------------
## Tests for the 0.3.7 conversion of `predict_guard()` to a proper S3
## `predict.GuardFit()` method. These tests assert that:
##   (1) `predict()` dispatches to `predict.GuardFit()` on a GuardFit
##       object (the editor's Comment 7 requirement);
##   (2) the new method registers correctly via the standard S3
##       generic so that `methods(class = "GuardFit")` includes
##       "predict";
##   (3) `predict_guard()` is preserved as a backward-compatible alias
##       that yields bit-identical output to the S3 method;
##   (4) error messages for non-GuardFit input are informative;
##   (5) the new method silently accepts `...` (matching the `predict`
##       generic's signature) so user code passing extra arguments does
##       not error.
## ---------------------------------------------------------------------

make_simple_guardfit <- function() {
  X <- data.frame(a = c(1, 2, NA, 4), b = c(10, 11, 12, 13))
  guard_fit(
    X, y = c(0.1, 0.2, 0.3, 0.4),
    steps = list(impute = list(method = "median")),
    task = "gaussian"
  )
}

test_that("predict() dispatches to predict.GuardFit on a GuardFit object", {
  fit <- make_simple_guardfit()
  X_new <- data.frame(a = c(NA, 5), b = c(9, 14))
  out <- predict(fit, X_new)

  ## We verify dispatch via structural equivalence: the output of
  ## predict() must be bit-identical to fit$transform(newdata), the
  ## raw transform that the S3 method wraps. This avoids hardcoding
  ## guard-fit-internal normalization details that may evolve.
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L)
  expect_identical(out, fit$transform(X_new))
  ## The output must NOT contain NA values (the guard fit was
  ## configured with median imputation; the NA in X_new$a should be
  ## filled before normalization).
  expect_false(any(is.na(out$a)))
})

test_that("methods(class = 'GuardFit') includes predict (Comment 7 check)", {
  cls_methods <- as.character(methods(class = "GuardFit"))
  ## Strip dispatch labels like "predict.GuardFit" to bare generics.
  ## utils::methods returns entries like "predict.GuardFit" so we
  ## check for that exact name.
  expect_true("predict.GuardFit" %in% cls_methods)
  expect_true("print.GuardFit"   %in% cls_methods)
  expect_true("summary.GuardFit" %in% cls_methods)
})

test_that("predict_guard() yields bit-identical output to predict.GuardFit()", {
  fit <- make_simple_guardfit()
  X_new <- data.frame(a = c(NA, 5), b = c(9, 14))
  expect_identical(predict_guard(fit, X_new), predict(fit, X_new))
  expect_identical(predict_guard(fit, X_new), predict.GuardFit(fit, X_new))
})

test_that("predict.GuardFit() rejects non-GuardFit input with a clear error", {
  X_new <- data.frame(a = 1:3, b = 4:6)
  expect_error(predict.GuardFit(NULL, X_new),
               "must be a GuardFit object")
  expect_error(predict.GuardFit(list(), X_new),
               "must be a GuardFit object")
  expect_error(predict.GuardFit(data.frame(x = 1:3), X_new),
               "must be a GuardFit object")
})

test_that("predict_guard() also rejects non-GuardFit input (via the S3 method)", {
  X_new <- data.frame(a = 1:3, b = 4:6)
  expect_error(predict_guard(NULL,    X_new), "must be a GuardFit object")
  expect_error(predict_guard(list(),  X_new), "must be a GuardFit object")
})

test_that("predict.GuardFit() accepts and silently drops `...` args", {
  fit <- make_simple_guardfit()
  X_new <- data.frame(a = c(NA, 5), b = c(9, 14))
  ## Extra arguments should not cause an error (S3 method signature).
  out_no_extra <- predict(fit, X_new)
  out_extra    <- predict(fit, X_new, type = "response", se.fit = FALSE)
  expect_identical(out_no_extra, out_extra)
})

test_that("predict.GuardFit() and predict_guard() are both exported", {
  exported <- getNamespaceExports("bioLeak")
  expect_true("predict_guard" %in% exported)
  ## predict.GuardFit is registered as an S3 method via S3method() and
  ## therefore appears in the bioLeak S3 method table. Check via
  ## utils::methods() filtered for our class.
  cls_methods <- as.character(methods(class = "GuardFit"))
  expect_true("predict.GuardFit" %in% cls_methods)
})
