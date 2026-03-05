test_that("signal_strength increases outcome separation", {
  set.seed(123)
  sim0 <- bioLeak:::.simulate_dataset(
    n = 500, p = 10, prevalence = 0.5,
    mode = "subject_grouped", leakage = "none",
    rho = 0, signal_strength = 0
  )
  set.seed(123)
  sim2 <- bioLeak:::.simulate_dataset(
    n = 500, p = 10, prevalence = 0.5,
    mode = "subject_grouped", leakage = "none",
    rho = 0, signal_strength = 2
  )

  feats <- paste0("x", sprintf("%02d", 1:5))
  score <- as.numeric(scale(rowSums(sim0$data[, feats, drop = FALSE])))
  y0 <- as.numeric(sim0$data$y) - 1
  y2 <- as.numeric(sim2$data$y) - 1

  cor0 <- cor(score, y0)
  cor2 <- cor(score, y2)

  expect_gt(cor2, cor0 + 0.2)
})
