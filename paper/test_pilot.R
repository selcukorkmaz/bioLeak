## Pilot test: verify simulate_leakage_suite works and measure timing
library(bioLeak)
cat("bioLeak loaded\n")

# Small pilot: 5 seeds, B=50, one config
t0 <- proc.time()
res <- simulate_leakage_suite(
  n        = 200,
  p        = 10,
  mode     = "subject_grouped",
  learner  = "glmnet",
  leakage  = "none",
  B        = 50,
  seeds    = 1:5,
  signal_strength = 1,
  verbose  = TRUE
)
t1 <- proc.time()
elapsed <- (t1 - t0)[3]
cat("Time for 5 seeds, B=50:", elapsed, "seconds\n")
cat("Per seed:", elapsed / 5, "seconds\n")
print(res)
cat("\nColumn names:", paste(names(res), collapse = ", "), "\n")
cat("Class:", class(res), "\n")
