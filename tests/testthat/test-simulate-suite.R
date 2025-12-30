test_that("simulate_dataset adds leakage-specific columns", {
  sim_none <- bioLeak:::.simulate_dataset(
    n = 20, p = 5, prevalence = 0.5,
    mode = "subject_grouped", leakage = "none",
    rho = 0, signal_strength = 1
  )
  expect_false("leak_subj" %in% names(sim_none$data))

  sim_subj <- bioLeak:::.simulate_dataset(
    n = 20, p = 5, prevalence = 0.5,
    mode = "subject_grouped", leakage = "subject_overlap",
    rho = 0, signal_strength = 1
  )
  expect_true("leak_subj" %in% names(sim_subj$data))

  sim_batch <- bioLeak:::.simulate_dataset(
    n = 20, p = 5, prevalence = 0.5,
    mode = "batch_blocked", leakage = "batch_confounded",
    rho = 0, signal_strength = 1
  )
  expect_true("leak_batch" %in% names(sim_batch$data))

  sim_peek <- bioLeak:::.simulate_dataset(
    n = 20, p = 5, prevalence = 0.5,
    mode = "subject_grouped", leakage = "peek_norm",
    rho = 0, signal_strength = 1
  )
  expect_true("leak_global" %in% names(sim_peek$data))

  sim_future <- bioLeak:::.simulate_dataset(
    n = 20, p = 5, prevalence = 0.5,
    mode = "time_series", leakage = "lookahead",
    rho = 0, signal_strength = 1
  )
  expect_true("leak_future" %in% names(sim_future$data))
})

test_that("simulate_leakage_suite runs when required learners are available", {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    expect_error(simulate_leakage_suite(
      n = 120, p = 5, prevalence = 0.5,
      mode = "batch_blocked", learner = "glmnet",
      leakage = "none", rho = 0, K = 3, repeats = 1,
      B = 3, seeds = 1:2, parallel = FALSE, signal_strength = 1
    ), "glmnet")
  } else {
    res <- simulate_leakage_suite(
      n = 120, p = 5, prevalence = 0.5,
      mode = "batch_blocked", learner = "glmnet",
      leakage = "none", rho = 0, K = 3, repeats = 1,
      B = 3, seeds = 1:2, parallel = FALSE, signal_strength = 1
    )
    expect_true(inherits(res, "LeakSimResults"))
    expect_true(all(c("seed", "metric_obs", "gap", "p_value", "leakage", "mode") %in% names(res)))
  }
})
