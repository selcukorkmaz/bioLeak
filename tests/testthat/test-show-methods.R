## ---------------------------------------------------------------------
## test-show-methods.R
## ---------------------------------------------------------------------
## Tests for the 0.3.7 addition of show()/print() methods to the public
## result classes. These tests assert that:
##   (1) `methods(class = "LeakFit")`   includes "show",
##   (2) `methods(class = "LeakAudit")` includes "show",
##   (3) `methods(class = "LeakTune")`  includes "print",
##       all addressing the editor's Comment 8 complaint that these
##       classes were missing print/show methods.
## We also assert that auto-print is non-empty and contains the class
## name and at least one diagnostic field.
## ---------------------------------------------------------------------

# ---- Build minimal fixtures by direct slot construction --------------

make_leakfit_fixture <- function() {
  splits <- make_split_plan_quiet(
    make_class_df(8),
    outcome = "outcome",
    mode = "subject_grouped",
    group = "subject",
    v = 2,
    seed = 1
  )
  metrics_df <- data.frame(
    learner = "stub",
    fold = 1:2,
    repeat_id = 1L,
    auc = c(0.71, 0.65),
    stringsAsFactors = FALSE
  )
  bioLeak:::LeakFit(
    splits         = splits,
    metrics        = metrics_df,
    metric_summary = data.frame(),
    audit          = data.frame(),
    predictions    = list(),
    preprocess     = list(),
    learners       = list(),
    outcome        = "outcome",
    task           = "binomial",
    feature_names  = character(),
    info           = list()
  )
}

make_leakaudit_fixture <- function() {
  fit <- make_leakfit_fixture()
  perm_gap <- data.frame(
    mechanism_class = "non_random_signal",
    repeat_id = 1L,
    metric_obs = 0.68,
    perm_mean = 0.50,
    perm_sd = 0.05,
    gap = 0.18,
    z = 3.6,
    p_value = 0.02,
    n_perm = 50L,
    stringsAsFactors = FALSE
  )
  bioLeak:::LeakAudit(
    fit             = fit,
    permutation_gap = perm_gap,
    perm_values     = c(0.48, 0.51, 0.55),
    batch_assoc     = data.frame(),
    target_assoc    = data.frame(),
    duplicates      = data.frame(),
    trail           = list(),
    info            = list()
  )
}

make_leaktune_fixture <- function() {
  ## LeakTune is an S3-classed list with $info and $outer_fits;
  ## construct a minimal one that satisfies print.LeakTune's accesses.
  fit <- make_leakfit_fixture()
  tune_obj <- list(
    outer_fits = list(fit, fit, fit),
    outer_summary = data.frame(),
    info = list(
      grid = 5L,
      refit_metric = "roc_auc",
      selection_rule = "one_std_err",
      refit = TRUE
    )
  )
  class(tune_obj) <- "LeakTune"
  tune_obj
}


# ---- methods(class = ...) registration tests -------------------------

test_that("methods(class = 'LeakFit') includes show (Comment 8 check)", {
  cls_methods <- as.character(methods(class = "LeakFit"))
  expect_true("show,LeakFit-method" %in% cls_methods ||
              "show.LeakFit"        %in% cls_methods)
  expect_true("summary.LeakFit"     %in% cls_methods)  # pre-existing
})

test_that("methods(class = 'LeakAudit') includes show (Comment 8 check)", {
  cls_methods <- as.character(methods(class = "LeakAudit"))
  expect_true("show,LeakAudit-method" %in% cls_methods ||
              "show.LeakAudit"        %in% cls_methods)
  expect_true("summary.LeakAudit"     %in% cls_methods)  # pre-existing
})

test_that("methods(class = 'LeakTune') includes print (Comment 8 check)", {
  cls_methods <- as.character(methods(class = "LeakTune"))
  expect_true("print.LeakTune"   %in% cls_methods)
  expect_true("summary.LeakTune" %in% cls_methods)       # pre-existing
})


# ---- Auto-print correctness (S4 show methods) -----------------------

test_that("show(LeakFit) auto-prints the class name and key fields", {
  fit <- make_leakfit_fixture()
  txt <- paste(capture.output(show(fit)), collapse = "\n")
  expect_match(txt, "A LeakFit object")
  expect_match(txt, "Task:\\s+binomial")
  expect_match(txt, "Outcome:\\s+outcome")
  expect_match(txt, "Folds:\\s+2\\b")
  ## Hint about summary() should be present.
  expect_match(txt, "summary\\(<obj>\\)")
})

test_that("show(LeakAudit) auto-prints the class name and key fields", {
  audit <- make_leakaudit_fixture()
  txt <- paste(capture.output(show(audit)), collapse = "\n")
  expect_match(txt, "A LeakAudit object")
  expect_match(txt, "Permutation-gap:\\s+metric=0\\.680")
  expect_match(txt, "Batch association:\\s+0 row\\(s\\)")
  expect_match(txt, "Target leakage:\\s+0 feature\\(s\\)")
  expect_match(txt, "Duplicate pairs:\\s+0")
  expect_match(txt, "summary\\(<obj>\\)")
})

test_that("print(LeakTune) auto-prints the class name and key fields", {
  tn <- make_leaktune_fixture()
  txt <- paste(capture.output(print(tn)), collapse = "\n")
  expect_match(txt, "A LeakTune object")
  expect_match(txt, "Outer folds:\\s+3 successful / 3 total")
  expect_match(txt, "Tuning grid:\\s+5 candidates")
  expect_match(txt, "Selection rule:\\s+one_std_err")
  expect_match(txt, "Refit:\\s+yes")
  expect_match(txt, "summary\\(<obj>\\)")
})


# ---- Auto-print is the default at the console (no explicit show/print) ----

test_that("Auto-printing at the console invokes show()/print()", {
  fit  <- make_leakfit_fixture()
  audit <- make_leakaudit_fixture()
  tn   <- make_leaktune_fixture()
  ## When R auto-prints an object, it calls show() (S4) or print() (S3).
  expect_match(paste(capture.output(fit),  collapse = "\n"), "A LeakFit object")
  expect_match(paste(capture.output(audit),collapse = "\n"), "A LeakAudit object")
  expect_match(paste(capture.output(tn),   collapse = "\n"), "A LeakTune object")
})


# ---- Each method returns its argument invisibly --------------------

test_that("show()/print() methods return their argument invisibly", {
  fit  <- make_leakfit_fixture()
  audit <- make_leakaudit_fixture()
  tn   <- make_leaktune_fixture()
  capture.output({
    expect_identical(show(fit), fit)
    expect_identical(show(audit), audit)
    expect_identical(print(tn), tn)
  })
})
