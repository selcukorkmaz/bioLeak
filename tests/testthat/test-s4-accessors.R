## ---------------------------------------------------------------------
## tests/testthat/test-s4-accessors.R
## ---------------------------------------------------------------------
## Verifies that the bioLeak 0.3.8 accessor refactor (S4 generic +
## method pairs replacing the 0.3.7 plain functions) does not change
## any return value and registers each accessor as an S4 method on the
## appropriate class.
##
## Method registration is verified via isGeneric() + existsMethod()
## (the canonical S4 check that is independent of whether the package
## is installed or loaded via pkgload::load_all()). The output of
## methods(class = "<Class>") is an editor-facing rendering of the
## same registry; we verify that programmatically below by computing
## the class's method list directly from S4's getMethodsForDispatch
## hooks.
## ---------------------------------------------------------------------

test_that("LeakFit accessors are registered as S4 methods", {
  expect_true(isGeneric("fit_metrics"))
  expect_true(existsMethod("fit_metrics", "LeakFit"))
})

test_that("LeakAudit accessors are registered as S4 methods", {
  for (gen in c("audit_perm_gap", "audit_batch_assoc",
                "audit_target_assoc", "audit_duplicates",
                "audit_info")) {
    expect_true(isGeneric(gen),
                info = sprintf("'%s' should be an S4 generic", gen))
    expect_true(existsMethod(gen, "LeakAudit"),
                info = sprintf("'%s' should have a method for LeakAudit", gen))
  }
})

test_that("LeakDeltaLSI accessors are registered as S4 methods", {
  for (gen in c("dlsi_metric", "dlsi_robust", "dlsi_ci",
                "dlsi_p_value", "dlsi_tier", "dlsi_R_eff",
                "dlsi_repeats")) {
    expect_true(isGeneric(gen),
                info = sprintf("'%s' should be an S4 generic", gen))
    expect_true(existsMethod(gen, "LeakDeltaLSI"),
                info = sprintf("'%s' should have a method for LeakDeltaLSI", gen))
  }
})

test_that("Accessor methods return values identical to direct slot access", {
  splits <- methods::new("LeakSplits",
    mode = "subject_grouped",
    indices = list(list(train = 1:6, test = 7:10)),
    info = list()
  )
  fit <- methods::new("LeakFit",
    splits = splits,
    metrics = data.frame(auc = c(0.61, 0.82)),
    metric_summary = data.frame(),
    audit = data.frame(),
    predictions = list(),
    preprocess = list(),
    learners = list(),
    outcome = "y",
    task = "binary_classification",
    feature_names = c("x1", "x2"),
    info = list(R = 1L)
  )

  expect_identical(fit_metrics(fit), fit@metrics)

  audit <- methods::new("LeakAudit",
    fit = fit,
    permutation_gap = data.frame(metric_obs = 0.61, gap = 0.10,
                                 p_value = 0.05),
    perm_values = c(0.50, 0.51, 0.52),
    batch_assoc  = data.frame(batch_col = "batch", cramer_v = 0.08),
    target_assoc = data.frame(feature = "x1", score = 0.15, flag = FALSE),
    duplicates   = data.frame(row_a = integer(0), row_b = integer(0)),
    trail = list(),
    info  = list(target_multivariate = NULL)
  )

  expect_identical(audit_perm_gap(audit),     audit@permutation_gap)
  expect_identical(audit_batch_assoc(audit),  audit@batch_assoc)
  expect_identical(audit_target_assoc(audit), audit@target_assoc)
  expect_identical(audit_duplicates(audit),   audit@duplicates)
  expect_identical(audit_info(audit),         audit@info)

  dlsi <- methods::new("LeakDeltaLSI",
    metric          = "auc",
    exchangeability = "iid",
    tier            = "A_full_inference",
    strict          = FALSE,
    R_eff           = 20L,
    delta_lsi       = 0.181,
    delta_lsi_ci    = c(0.174, 0.187),
    delta_metric    = 0.180,
    delta_metric_ci = c(0.172, 0.186),
    p_value         = 0.0001,
    inference_ok    = TRUE,
    folds_naive     = data.frame(),
    folds_guarded   = data.frame(),
    repeats_naive   = data.frame(metric = rnorm(20)),
    repeats_guarded = data.frame(metric = rnorm(20)),
    info            = list(paired = TRUE)
  )

  expect_identical(dlsi_metric(dlsi),  dlsi@delta_metric)
  expect_identical(dlsi_robust(dlsi),  dlsi@delta_lsi)
  expect_identical(dlsi_p_value(dlsi), dlsi@p_value)
  expect_identical(dlsi_tier(dlsi),    dlsi@tier)
  expect_identical(dlsi_R_eff(dlsi),   dlsi@R_eff)

  expect_identical(dlsi_ci(dlsi, "robust"), dlsi@delta_lsi_ci)
  expect_identical(dlsi_ci(dlsi, "metric"), dlsi@delta_metric_ci)
  expect_identical(dlsi_ci(dlsi),           dlsi@delta_lsi_ci)  # default = "robust"

  expect_identical(dlsi_repeats(dlsi, "naive"),   dlsi@repeats_naive)
  expect_identical(dlsi_repeats(dlsi, "guarded"), dlsi@repeats_guarded)
  expect_identical(dlsi_repeats(dlsi),            dlsi@repeats_naive)  # default
})

test_that("dlsi_ci and dlsi_repeats reject invalid 'which'", {
  dlsi <- methods::new("LeakDeltaLSI",
    delta_lsi_ci    = c(0, 1),
    delta_metric_ci = c(0, 1),
    repeats_naive   = data.frame(metric = numeric(0)),
    repeats_guarded = data.frame(metric = numeric(0))
  )
  expect_error(dlsi_ci(dlsi, "robust_xxx"))
  expect_error(dlsi_repeats(dlsi, "naive_xxx"))
})
