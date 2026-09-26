# Regression tests: simulate_leakage_suite() must honour `prevalence`.
# Before 0.3.8 the generator used pnorm(linpred - qnorm(prevalence)), which
# inverted it (requested 0.2 gave ~0.72; 0.8 gave ~0.27).

sim_prev <- function(prevalence, signal_strength = 1, rho = 0, n = 5000, seed = 11) {
  set.seed(seed)
  sim <- bioLeak:::.simulate_dataset(
    n = n, p = 10, prevalence = prevalence, mode = "subject_grouped",
    leakage = "none", rho = rho, signal_strength = signal_strength
  )
  mean(sim$data$y == "1")
}

test_that("simulated prevalence matches the requested prevalence", {
  for (prev in c(0.1, 0.2, 0.5, 0.8, 0.9)) {
    expect_equal(sim_prev(prev), prev, tolerance = 0.03,
                 info = sprintf("requested prevalence %.1f", prev))
  }
})

test_that("prevalence is calibrated across signal strengths", {
  for (s in c(0, 0.5, 2, 4)) {
    expect_equal(sim_prev(0.3, signal_strength = s), 0.3, tolerance = 0.03,
                 info = sprintf("signal_strength = %s", s))
  }
})

test_that("positive class probability still increases with the linear predictor", {
  set.seed(12)
  sim <- bioLeak:::.simulate_dataset(
    n = 3000, p = 5, prevalence = 0.2, mode = "subject_grouped",
    leakage = "none", rho = 0, signal_strength = 2
  )
  lp <- rowSums(as.matrix(sim$data[, sprintf("x%02d", 1:5)]))
  expect_gt(bioLeak:::.auc_binary(sim$data$y, lp), 0.75)
})

test_that("benchmark modality profiles produce their declared prevalence", {
  # Profiles used by benchmark_leakage_suite() (imaging_tabular / ehr_tabular).
  expect_equal(sim_prev(0.4, rho = 0.15, n = 5000), 0.4, tolerance = 0.03)
  expect_equal(sim_prev(0.3, rho = 0.5, n = 5000), 0.3, tolerance = 0.03)
})
