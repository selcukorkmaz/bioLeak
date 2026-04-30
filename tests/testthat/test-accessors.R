## ---------------------------------------------------------------------
## test-accessors.R
## ---------------------------------------------------------------------
## Tests for the public accessor functions added in 0.3.7. Each accessor
## should:
##   (a) return exactly the same object as the corresponding @ slot read,
##   (b) error informatively when called on the wrong class, and
##   (c) (where applicable) honor the `which` argument that selects
##       between two related slots.
## ---------------------------------------------------------------------


# ---- Build minimal fixtures by direct slot construction --------------
# These tests avoid running fit_resample / audit_leakage / delta_lsi to
# keep the suite fast and to isolate accessor behavior from upstream
# computational changes. The fixtures populate every slot the accessors
# read, and only those slots, so an accessor regression is unambiguous.

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
  batch_assoc <- data.frame(
    mechanism_class = "confounding_alignment",
    batch_col = "batch",
    repeat_id = 1L,
    stat = 4.0,
    df = 2L,
    pval = 0.135,
    cramer_v = 0.30,
    stringsAsFactors = FALSE
  )
  target_assoc <- data.frame(
    mechanism_class = "proxy_target_leakage",
    feature = c("x1", "x2"),
    type = "numeric",
    metric = "auc",
    value = c(0.55, 0.72),
    score = c(0.10, 0.45),
    p_value = c(0.40, 0.01),
    n = 8L,
    n_levels = NA_integer_,
    flag = c(FALSE, FALSE),
    p_value_adj = NA_real_,
    flag_fdr = FALSE,
    stringsAsFactors = FALSE
  )
  duplicates <- data.frame(
    mechanism_class = character(0),
    i = integer(0),
    j = integer(0),
    sim = numeric(0),
    cross_fold = logical(0),
    cos_sim = numeric(0),
    stringsAsFactors = FALSE
  )
  info <- list(
    target_multivariate = list(score = 0.42, p_value = 0.07),
    perm_method = "fixed_predictions"
  )
  bioLeak:::LeakAudit(
    fit             = fit,
    permutation_gap = perm_gap,
    perm_values     = c(0.48, 0.51, 0.55),
    batch_assoc     = batch_assoc,
    target_assoc    = target_assoc,
    duplicates      = duplicates,
    trail           = list(),
    info            = info
  )
}

make_leakdeltalsi_fixture <- function() {
  repeats_naive <- data.frame(
    repeat_id = 1:3,
    metric = c(0.91, 0.92, 0.90),
    stringsAsFactors = FALSE
  )
  repeats_guarded <- data.frame(
    repeat_id = 1:3,
    metric = c(0.62, 0.60, 0.58),
    stringsAsFactors = FALSE
  )
  bioLeak:::LeakDeltaLSI(
    metric          = "auc",
    exchangeability = "iid",
    tier            = "B_signflip_ci",
    strict          = FALSE,
    R_eff           = 3L,
    delta_lsi       = 0.310,
    delta_lsi_ci    = c(0.295, 0.325),
    delta_metric    = 0.305,
    delta_metric_ci = c(0.290, 0.320),
    p_value         = 0.0033,
    inference_ok    = FALSE,
    folds_naive     = data.frame(),
    folds_guarded   = data.frame(),
    repeats_naive   = repeats_naive,
    repeats_guarded = repeats_guarded,
    info            = list()
  )
}


# ---- LeakFit accessors ----------------------------------------------

test_that("fit_metrics() returns the @metrics slot for a LeakFit", {
  fit <- make_leakfit_fixture()
  expect_identical(fit_metrics(fit), fit@metrics)
  expect_s3_class(fit_metrics(fit), "data.frame")
  expect_named(fit_metrics(fit), c("learner", "fold", "repeat_id", "auc"))
})

test_that("fit_metrics() rejects non-LeakFit input", {
  expect_error(fit_metrics(NULL),
               "must be a LeakFit object")
  expect_error(fit_metrics(data.frame(x = 1)),
               "must be a LeakFit object")
  audit <- make_leakaudit_fixture()
  expect_error(fit_metrics(audit),
               "must be a LeakFit object")
})


# ---- LeakAudit accessors --------------------------------------------

test_that("audit_perm_gap() returns the @permutation_gap slot", {
  audit <- make_leakaudit_fixture()
  expect_identical(audit_perm_gap(audit), audit@permutation_gap)
  expect_true("metric_obs" %in% names(audit_perm_gap(audit)))
})

test_that("audit_perm_gap() rejects non-LeakAudit input", {
  expect_error(audit_perm_gap(NULL), "must be a LeakAudit object")
  fit <- make_leakfit_fixture()
  expect_error(audit_perm_gap(fit), "must be a LeakAudit object")
})

test_that("audit_batch_assoc() returns the @batch_assoc slot", {
  audit <- make_leakaudit_fixture()
  expect_identical(audit_batch_assoc(audit), audit@batch_assoc)
  expect_true("cramer_v" %in% names(audit_batch_assoc(audit)))
})

test_that("audit_batch_assoc() rejects non-LeakAudit input", {
  expect_error(audit_batch_assoc("not an audit"),
               "must be a LeakAudit object")
})

test_that("audit_target_assoc() returns the @target_assoc slot", {
  audit <- make_leakaudit_fixture()
  expect_identical(audit_target_assoc(audit), audit@target_assoc)
  expect_true("score" %in% names(audit_target_assoc(audit)))
})

test_that("audit_target_assoc() rejects non-LeakAudit input", {
  expect_error(audit_target_assoc(list(target_assoc = data.frame())),
               "must be a LeakAudit object")
})

test_that("audit_duplicates() returns the @duplicates slot", {
  audit <- make_leakaudit_fixture()
  expect_identical(audit_duplicates(audit), audit@duplicates)
  expect_s3_class(audit_duplicates(audit), "data.frame")
})

test_that("audit_duplicates() rejects non-LeakAudit input", {
  expect_error(audit_duplicates(1L),
               "must be a LeakAudit object")
})

test_that("audit_info() returns the @info slot", {
  audit <- make_leakaudit_fixture()
  expect_identical(audit_info(audit), audit@info)
  expect_named(audit_info(audit),
               c("target_multivariate", "perm_method"))
  expect_equal(audit_info(audit)$target_multivariate$p_value, 0.07)
})

test_that("audit_info() rejects non-LeakAudit input", {
  expect_error(audit_info(list()), "must be a LeakAudit object")
})


# ---- LeakDeltaLSI accessors -----------------------------------------

test_that("dlsi_metric() returns the @delta_metric slot", {
  dlsi <- make_leakdeltalsi_fixture()
  expect_identical(dlsi_metric(dlsi), dlsi@delta_metric)
  expect_equal(dlsi_metric(dlsi), 0.305)
})

test_that("dlsi_robust() returns the @delta_lsi slot", {
  dlsi <- make_leakdeltalsi_fixture()
  expect_identical(dlsi_robust(dlsi), dlsi@delta_lsi)
  expect_equal(dlsi_robust(dlsi), 0.310)
})

test_that("dlsi_ci(which='robust') returns the @delta_lsi_ci slot (default)", {
  dlsi <- make_leakdeltalsi_fixture()
  expect_identical(dlsi_ci(dlsi), dlsi@delta_lsi_ci)
  expect_identical(dlsi_ci(dlsi, which = "robust"), dlsi@delta_lsi_ci)
  expect_equal(dlsi_ci(dlsi), c(0.295, 0.325))
})

test_that("dlsi_ci(which='metric') returns the @delta_metric_ci slot", {
  dlsi <- make_leakdeltalsi_fixture()
  expect_identical(dlsi_ci(dlsi, which = "metric"), dlsi@delta_metric_ci)
  expect_equal(dlsi_ci(dlsi, which = "metric"), c(0.290, 0.320))
})

test_that("dlsi_ci() rejects unknown `which` values", {
  dlsi <- make_leakdeltalsi_fixture()
  expect_error(dlsi_ci(dlsi, which = "raw"))      # match.arg() error
  expect_error(dlsi_ci(dlsi, which = "huber"))
})

test_that("dlsi_p_value() returns the @p_value slot", {
  dlsi <- make_leakdeltalsi_fixture()
  expect_identical(dlsi_p_value(dlsi), dlsi@p_value)
  expect_equal(dlsi_p_value(dlsi), 0.0033)
})

test_that("dlsi_tier() returns the @tier slot", {
  dlsi <- make_leakdeltalsi_fixture()
  expect_identical(dlsi_tier(dlsi), dlsi@tier)
  expect_equal(dlsi_tier(dlsi), "B_signflip_ci")
})

test_that("dlsi_R_eff() returns the @R_eff slot", {
  dlsi <- make_leakdeltalsi_fixture()
  expect_identical(dlsi_R_eff(dlsi), dlsi@R_eff)
  expect_equal(dlsi_R_eff(dlsi), 3L)
})

test_that("dlsi_repeats(which='naive') returns the @repeats_naive slot (default)", {
  dlsi <- make_leakdeltalsi_fixture()
  expect_identical(dlsi_repeats(dlsi), dlsi@repeats_naive)
  expect_identical(dlsi_repeats(dlsi, which = "naive"), dlsi@repeats_naive)
  expect_equal(dlsi_repeats(dlsi)$metric, c(0.91, 0.92, 0.90))
})

test_that("dlsi_repeats(which='guarded') returns the @repeats_guarded slot", {
  dlsi <- make_leakdeltalsi_fixture()
  expect_identical(dlsi_repeats(dlsi, which = "guarded"), dlsi@repeats_guarded)
  expect_equal(dlsi_repeats(dlsi, which = "guarded")$metric,
               c(0.62, 0.60, 0.58))
})

test_that("dlsi_repeats() rejects unknown `which` values", {
  dlsi <- make_leakdeltalsi_fixture()
  expect_error(dlsi_repeats(dlsi, which = "leaky"))   # match.arg() error
})

test_that("LeakDeltaLSI accessors all reject non-LeakDeltaLSI input", {
  not_dlsi <- list(delta_lsi = 1, delta_metric = 1)   # mimic-y but wrong class
  expect_error(dlsi_metric(not_dlsi),   "must be a LeakDeltaLSI object")
  expect_error(dlsi_robust(not_dlsi),   "must be a LeakDeltaLSI object")
  expect_error(dlsi_ci(not_dlsi),       "must be a LeakDeltaLSI object")
  expect_error(dlsi_p_value(not_dlsi),  "must be a LeakDeltaLSI object")
  expect_error(dlsi_tier(not_dlsi),     "must be a LeakDeltaLSI object")
  expect_error(dlsi_R_eff(not_dlsi),    "must be a LeakDeltaLSI object")
  expect_error(dlsi_repeats(not_dlsi),  "must be a LeakDeltaLSI object")
})


# ---- Round-trip equivalence to direct @ access ----------------------
# These tests assert that, for every accessor, calling it produces the
# bit-identical R object as direct slot access. This guards against
# any future helper that silently copies, reformats, or transforms the
# underlying data.

test_that("Every accessor is bit-identical to its corresponding @ slot read", {
  fit <- make_leakfit_fixture()
  audit <- make_leakaudit_fixture()
  dlsi <- make_leakdeltalsi_fixture()

  expect_identical(fit_metrics(fit),               fit@metrics)
  expect_identical(audit_perm_gap(audit),          audit@permutation_gap)
  expect_identical(audit_batch_assoc(audit),       audit@batch_assoc)
  expect_identical(audit_target_assoc(audit),      audit@target_assoc)
  expect_identical(audit_duplicates(audit),        audit@duplicates)
  expect_identical(audit_info(audit),              audit@info)
  expect_identical(dlsi_metric(dlsi),              dlsi@delta_metric)
  expect_identical(dlsi_robust(dlsi),              dlsi@delta_lsi)
  expect_identical(dlsi_ci(dlsi, "robust"),        dlsi@delta_lsi_ci)
  expect_identical(dlsi_ci(dlsi, "metric"),        dlsi@delta_metric_ci)
  expect_identical(dlsi_p_value(dlsi),             dlsi@p_value)
  expect_identical(dlsi_tier(dlsi),                dlsi@tier)
  expect_identical(dlsi_R_eff(dlsi),               dlsi@R_eff)
  expect_identical(dlsi_repeats(dlsi, "naive"),    dlsi@repeats_naive)
  expect_identical(dlsi_repeats(dlsi, "guarded"),  dlsi@repeats_guarded)
})
