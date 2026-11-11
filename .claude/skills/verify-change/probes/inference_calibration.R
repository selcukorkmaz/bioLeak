# Calibration of delta_lsi() inference with REAL model fits, end to end.
# Run from the bioLeak package root:
#   Rscript .claude/skills/verify-change/probes/inference_calibration.R
# Env: PROBE_B (datasets per design, default 150), PROBE_CORES (default 4).
#
# Why real fits: fold metrics from repeated CV are correlated through shared
# data (overlapping training sets, re-partitioned repeats). Simulations that
# draw fold metrics independently cannot show that, and hid a 70% type-I error.
#
# Null: pure-noise data; pipeline A = glm on features 1-10, pipeline B = glm on
# features 11-20. Both have true AUC 0.5, but their predictions differ.
# Alternative (power): feature 1 carries signal, so A is truly better.

source(".claude/skills/verify-change/probes/_probe_helpers.R")
load_dev(".")

B     <- as.integer(Sys.getenv("PROBE_B", "150"))
cores <- as.integer(Sys.getenv("PROBE_CORES", "4"))
alpha <- 0.05
lim   <- type1_limit(alpha, B)

auc_rank <- function(y, s) {
  r <- rank(s); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# One replicate: simulate data, split, fit both pipelines per fold, run delta_lsi.
one_rep <- function(b, design, signal = 0) {
  set.seed(b)
  n <- design$G * design$per
  X <- matrix(rnorm(n * 20), n)
  y <- rbinom(n, 1, stats::plogis(signal * X[, 1]))
  df <- data.frame(g = rep(sprintf("G%02d", seq_len(design$G)), each = design$per), y = y)
  sp <- switch(design$mode,
    batch   = make_split_plan(df, outcome = "y", mode = "batch_blocked", batch = "g",
                              v = design$v, repeats = design$R, seed = b, progress = FALSE),
    study   = make_split_plan(df, outcome = "y", mode = "study_loocv", study = "g", progress = FALSE),
    subject = make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "g",
                              v = design$v, repeats = design$R, seed = b, progress = FALSE))
  J <- length(sp@indices); aA <- aB <- numeric(J)
  for (j in seq_len(J)) {
    te <- sp@indices[[j]]$test; tr <- sp@indices[[j]]$train
    if (length(unique(y[te])) < 2L || length(unique(y[tr])) < 2L) { aA[j] <- aB[j] <- NA; next }
    fA <- suppressWarnings(stats::glm(y[tr] ~ X[tr, 1:10], family = stats::binomial))
    fB <- suppressWarnings(stats::glm(y[tr] ~ X[tr, 11:20], family = stats::binomial))
    aA[j] <- auc_rank(y[te], cbind(1, X[te, 1:10]) %*% stats::coef(fA))
    aB[j] <- auc_rank(y[te], cbind(1, X[te, 11:20]) %*% stats::coef(fB))
  }
  mk <- function(a) methods::new("LeakFit", splits = sp,
    metrics = data.frame(fold = seq_len(J), learner = "glm", auc = a),
    outcome = "y", task = "binomial")
  r <- suppressWarnings(delta_lsi(mk(aA), mk(aB), metric = "auc"))
  c(p = r@p_value, lo = r@delta_metric_ci[1], hi = r@delta_metric_ci[2], R_eff = r@R_eff)
}

run_design <- function(design, signal = 0) {
  out <- parallel::mclapply(seq_len(B), one_rep, design = design, signal = signal,
                            mc.cores = cores)
  m <- do.call(rbind, out)
  list(rate = mean(m[, "p"] < alpha, na.rm = TRUE), n_na = sum(is.na(m[, "p"])),
       ci_excl0 = mean(m[, "lo"] > 0 | m[, "hi"] < 0, na.rm = TRUE),
       R_eff = stats::median(m[, "R_eff"]))
}

designs <- list(
  list(id = "INV-1", label = "120 samples, 5-fold x 10 repeats",          mode = "subject", G = 120, per = 1,  v = 5, R = 10),
  list(id = "INV-1", label = "20 subjects x 6, 5-fold x 10 repeats",      mode = "subject", G = 20,  per = 6,  v = 5, R = 10),
  list(id = "BL-9",  label = "6 batches x 20, 3-fold x 20 repeats",       mode = "batch",   G = 6,   per = 20, v = 3, R = 20),
  list(id = "BL-8",  label = "6 batches x 20, leave-one-batch-out x 10",  mode = "batch",   G = 6,   per = 20, v = 6, R = 10),
  list(id = "INV-1", label = "study_loocv, 8 studies x 20",               mode = "study",   G = 8,   per = 20, v = 0, R = 1)
)

for (d in designs) {
  r <- run_design(d)
  check(d$id, r$rate <= lim && r$ci_excl0 <= lim,
        sprintf("H0, %s: type-I %.3f, CI excludes 0 %.3f (limit %.3f), R_eff %g, NA p %d/%d",
                d$label, r$rate, r$ci_excl0, lim, r$R_eff, r$n_na, B))
}

# Positive control: a clear real advantage must be detected (power).
pw <- run_design(designs[[1]], signal = 2)
check("CTRL-1", pw$rate >= 0.8,
      sprintf("power, %s, strong signal (logit slope 2): rejection %.3f (want >= 0.80)", designs[[1]]$label, pw$rate))
# Information only: power for a weak advantage, the price of calibrated df.
pw1 <- run_design(designs[[1]], signal = 1)
cat(sprintf("INFO [power] %s, weak signal (logit slope 1): rejection %.3f\n", designs[[1]]$label, pw1$rate))

probe_finish()
