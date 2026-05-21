## ---------------------------------------------------------------------
## tests/testthat/test-plot-methods.R
## ---------------------------------------------------------------------
## Verifies the S4 plot() methods added in bioLeak 0.3.8:
##   - plot(<LeakAudit>)       -> plot_perm_distribution()
##   - plot(<LeakFit>)         -> plot_fold_balance() by default
##   - plot(<LeakDeltaLSI>)    -> plot_dlsi_repeats()
##
## We verify (i) method registration, (ii) that calling plot() dispatches
## without error on representative objects, and (iii) that the return
## value is the list shape produced by the underlying helper.
## ---------------------------------------------------------------------

skip_if_no_ggplot <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    skip("ggplot2 not available")
  }
}

test_that("plot methods are registered on the S4 classes", {
  expect_true(existsMethod("plot",
                           signature(x = "LeakAudit",     y = "missing")))
  expect_true(existsMethod("plot",
                           signature(x = "LeakFit",       y = "missing")))
  expect_true(existsMethod("plot",
                           signature(x = "LeakDeltaLSI",  y = "missing")))
})

test_that("plot(<LeakDeltaLSI>) dispatches and returns the expected shape", {
  skip_if_no_ggplot()
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
    repeats_naive   = data.frame(metric = rnorm(20, mean = 0.78, sd = 0.02)),
    repeats_guarded = data.frame(metric = rnorm(20, mean = 0.60, sd = 0.02)),
    info            = list(paired = TRUE)
  )

  out <- plot(dlsi)

  expect_type(out, "list")
  expect_true(all(c("deltas", "delta_robust", "delta_mean",
                    "bca_ci", "R_eff", "plot") %in% names(out)))
  expect_length(out$deltas, 20L)
  expect_equal(out$R_eff, 20L)
  expect_true(inherits(out$plot, "ggplot"))
})

test_that("plot_dlsi_repeats() rejects non-LeakDeltaLSI input", {
  expect_error(plot_dlsi_repeats(list()),
               regexp = "LeakDeltaLSI", ignore.case = TRUE)
})

test_that("plot(<LeakDeltaLSI>) is the same as plot_dlsi_repeats()", {
  skip_if_no_ggplot()
  dlsi <- methods::new("LeakDeltaLSI",
    metric          = "auc",
    R_eff           = 20L,
    delta_lsi       = 0.181,
    delta_lsi_ci    = c(0.174, 0.187),
    delta_metric    = 0.180,
    delta_metric_ci = c(0.172, 0.186),
    p_value         = 0.0001,
    repeats_naive   = data.frame(metric = rep(0.78, 20)),
    repeats_guarded = data.frame(metric = rep(0.60, 20))
  )

  out_method <- plot(dlsi)
  out_func   <- plot_dlsi_repeats(dlsi)
  expect_equal(out_method$deltas, out_func$deltas)
  expect_equal(out_method$delta_robust, out_func$delta_robust)
  expect_equal(out_method$bca_ci, out_func$bca_ci)
})
