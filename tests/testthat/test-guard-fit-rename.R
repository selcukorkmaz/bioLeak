## ---------------------------------------------------------------------
## test-guard-fit-rename.R
## ---------------------------------------------------------------------
## Tests for the 0.3.7 rename of `.guard_fit` -> `guard_fit` and
## `.guard_ensure_levels` -> `guard_ensure_levels`. These tests ensure
## the editor's Comment 6 fix stays in place: no public bioLeak function
## should have a leading-dot name.
## ---------------------------------------------------------------------

test_that("guard_fit() (no leading dot) is exported and works", {
  x <- data.frame(a = c(1, 2, NA), b = c(3, 4, 5))
  fit <- guard_fit(
    x, y = c(1, 2, 3),
    steps = list(impute = list(method = "median")),
    task = "gaussian"
  )
  expect_s3_class(fit, "GuardFit")
  expect_true(is.function(fit$transform))
})

test_that("guard_ensure_levels() (no leading dot) is exported and works", {
  df <- data.frame(
    site = c("A", "B", "B"),
    status = c("yes", "no", "yes"),
    stringsAsFactors = FALSE
  )
  out <- guard_ensure_levels(df)
  expect_named(out, c("data", "levels"))
  expect_true(is.factor(out$data$site))
  expect_true(is.factor(out$data$status))
})

test_that("Old dot-prefixed names `.guard_fit` and `.guard_ensure_levels` are no longer exported", {
  exported <- getNamespaceExports("bioLeak")

  ## The dot-prefixed names should NOT appear in the export list.
  expect_false(".guard_fit"           %in% exported)
  expect_false(".guard_ensure_levels" %in% exported)

  ## The renamed functions SHOULD appear in the export list.
  expect_true("guard_fit"           %in% exported)
  expect_true("guard_ensure_levels" %in% exported)
})

test_that("help(package = 'bioLeak') no longer indexes any user-facing '.<name>' export", {
  ## The package-help index lists each export. We assert that no
  ## user-facing entry in the export list begins with a leading dot.
  ## Auto-generated S4 mangled names (".__C__<Class>", ".__T__<gen>:<pkg>")
  ## are excluded; those are R's own internal class/method registry
  ## and never appear in `help(package = "bioLeak")` output. The check
  ## targets only names that begin with `.` followed by an alphabetic
  ## character (the form the editor's Comment 6 flagged).
  exported <- getNamespaceExports("bioLeak")
  dot_prefixed_user <- grep("^\\.[a-zA-Z]", exported, value = TRUE)
  expect_length(dot_prefixed_user, 0L)
})
