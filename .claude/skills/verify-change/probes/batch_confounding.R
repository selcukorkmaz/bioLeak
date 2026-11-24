# Calibration and power of the confounding_alignment rule (JSS review L1.3):
# flag only when batch is associated with the outcome AND batches are shared
# by training and test rows. Run from the bioLeak package root:
#   Rscript .claude/skills/verify-change/probes/batch_confounding.R
# Env: PROBE_B (datasets per design, default 200), PROBE_CORES (default 4).

source(".claude/skills/verify-change/probes/_probe_helpers.R")
load_dev(".")

B     <- as.integer(Sys.getenv("PROBE_B", "200"))
cores <- as.integer(Sys.getenv("PROBE_CORES", "4"))
alpha <- 0.05
lim   <- type1_limit(alpha, B)

glm_learner <- list(glm = list(
  fit = function(x, y, task, weights, ...) {
    stats::glm(y ~ ., data = data.frame(y = y, x), family = stats::binomial())
  },
  predict = function(object, newdata, task, ...) {
    as.numeric(stats::predict(object, newdata = data.frame(newdata), type = "response"))
  }
))
lm_learner <- list(glm = list(
  fit = function(x, y, task, weights, ...) stats::lm(y ~ ., data = data.frame(y = y, x)),
  predict = function(object, newdata, task, ...) as.numeric(stats::predict(object, newdata = data.frame(newdata)))
))

sim <- function(b, p_case, gaussian = FALSE, n_batch = 6, per_batch = 25, two_cols = FALSE) {
  set.seed(b)
  batch <- rep(paste0("b", seq_len(n_batch)), each = per_batch)
  pb <- rep(p_case, length.out = n_batch)[as.integer(factor(batch))]
  df <- data.frame(subject = paste0("s", seq_along(batch)), batch = batch,
                   x1 = stats::rnorm(length(batch)), x2 = stats::rnorm(length(batch)))
  if (two_cols) df$site <- sample(paste0("t", 1:4), nrow(df), TRUE)
  if (gaussian) df$y <- stats::rnorm(nrow(df), mean = 2 * (pb - 0.5))
  else df$y <- factor(ifelse(stats::runif(nrow(df)) < pb, "case", "ctrl"), levels = c("ctrl", "case"))
  df
}

flag_of <- function(df, mode = "subject_grouped", repeats = 1) {
  gaussian <- is.numeric(df$y)
  sp <- if (mode == "batch_blocked") {
    make_split_plan(df, outcome = "y", mode = mode, batch = "batch", v = 3, seed = 1, progress = FALSE)
  } else if (mode == "study_loocv") {
    make_split_plan(df, outcome = "y", mode = mode, study = "study", group = "subject",
                    seed = 1, progress = FALSE)
  } else if (mode == "combined") {
    suppressWarnings(make_split_plan(df, outcome = "y", mode = "combined",
      constraints = list(list(type = "batch", col = "proc"), list(type = "subject", col = "subject")),
      v = 3, seed = 1, progress = FALSE))
  } else if (mode == "time_series") {
    make_split_plan(df, outcome = "y", mode = mode, time = "t", v = 5, seed = 1, progress = FALSE)
  } else {
    make_split_plan(df, outcome = "y", mode = mode, group = "subject", v = 5,
                    repeats = repeats, seed = 1, progress = FALSE)
  }
  keep <- c(switch(mode, batch_blocked = "batch", study_loocv = "study",
                   time_series = "t", combined = c("proc", "subject"), "subject"),
            "y", "x1", "x2")
  fit <- NULL
  utils::capture.output(fit <- suppressWarnings(suppressMessages(fit_resample(
    df[, keep], outcome = "y", splits = sp, learner = "glm",
    custom_learners = if (gaussian) lm_learner else glm_learner,
    metrics = if (gaussian) "rmse" else "auc", refit = FALSE, seed = 1))))
  bc <- intersect(c("batch", "site"), names(df))
  aud <- NULL
  utils::capture.output(aud <- suppressWarnings(suppressMessages(audit_leakage(
    fit, metric = if (gaussian) "rmse" else "auc", B = 5, coldata = df, batch_cols = bc,
    target_scan = FALSE, target_scan_multivariate = FALSE))))
  ms <- aud@info$mechanism_summary
  out <- isTRUE(ms$flagged[ms$mechanism_class == "confounding_alignment"])
  attr(out, "within") <- aud@batch_assoc$within_unit_effect[1]
  out
}

# Repeated measures: subject-level batch and outcome, `per` rows per subject.
# The test must count subjects, not rows (rows gave a 68% false-alarm rate).
sim_clustered <- function(b, p_case, n_subj = 40, per = 5, n_batch = 4) {
  set.seed(b)
  sb <- sample(paste0("b", seq_len(n_batch)), n_subj, TRUE)
  pb <- rep(p_case, length.out = n_batch)[as.integer(factor(sb, levels = paste0("b", seq_len(n_batch))))]
  sy <- ifelse(stats::runif(n_subj) < pb, "case", "ctrl")
  s <- rep(seq_len(n_subj), each = per)
  data.frame(subject = paste0("s", s), batch = sb[s],
             y = factor(sy[s], levels = c("ctrl", "case")),
             x1 = stats::rnorm(length(s)), x2 = stats::rnorm(length(s)))
}

# Batch varying within subject: each subject's rows are spread over batches.
# within = "random": batch random per row, outcome subject-level (null).
# within = "paired": per subject 1 normal row in batch N and 2 tumour rows in
# batch T (complete within-subject confounding; outside the L1.2 contract:
# reported as within_unit_effect, not flagged).
sim_within <- function(b, within = c("random", "paired"), n_subj = 30) {
  within <- match.arg(within)
  set.seed(b)
  if (within == "paired") {
    s <- rep(seq_len(n_subj), each = 3)
    y <- rep(c("ctrl", "case", "case"), n_subj)
    batch <- ifelse(y == "ctrl", "N", "T")
  } else {
    sizes <- sample(2:6, n_subj, TRUE)
    s <- rep(seq_len(n_subj), sizes)
    y <- rep(sample(c("ctrl", "case"), n_subj, TRUE), sizes)
    batch <- sample(paste0("b", 1:3), length(s), TRUE)
  }
  data.frame(subject = paste0("s", s), batch = batch,
             y = factor(y, levels = c("ctrl", "case")),
             x1 = stats::rnorm(length(s)), x2 = stats::rnorm(length(s)))
}

# rsample grouped folds (mode "rsample"): the unit must come from perm_mode.
flag_rsample <- function(df) {
  rs <- rsample::group_vfold_cv(df, group = "subject", v = 5)
  fit <- NULL
  utils::capture.output(fit <- suppressWarnings(suppressMessages(fit_resample(
    df[, c("subject", "batch", "y", "x1", "x2")], outcome = "y", splits = rs, learner = "glm",
    custom_learners = glm_learner, metrics = "auc",
    split_cols = list(group = "subject", batch = "batch")))))
  aud <- NULL
  utils::capture.output(aud <- suppressWarnings(suppressMessages(audit_leakage(
    fit, metric = "auc", B = 5, perm_refit = FALSE, target_scan = FALSE,
    target_scan_multivariate = FALSE, duplicate_null_reps = 0))))
  ms <- aud@info$mechanism_summary
  isTRUE(ms$flagged[ms$mechanism_class == "confounding_alignment"])
}

# Unequal subject sizes (2-40 rows), batch varying within subject.
sim_unequal <- function(b, link = 0.5, n_subj = 24) {
  set.seed(b)
  sizes <- sample(2:40, n_subj)
  s <- rep(seq_len(n_subj), sizes)
  y <- ifelse(stats::runif(length(s)) < 0.6, "case", "ctrl")
  batch <- ifelse(stats::runif(length(s)) < ifelse(y == "case", link, 1 - link), "A", "B")
  data.frame(subject = paste0("s", s), batch = batch,
             y = factor(y, levels = c("ctrl", "case")),
             x1 = stats::rnorm(length(s)), x2 = stats::rnorm(length(s)))
}
# Leave-study-out with a recorded subject column; processing batch crosses studies.
sim_loocv <- function(b) {
  df <- sim_clustered(b, 0.5)
  df$study <- paste0("st", (as.integer(sub("s", "", df$subject)) - 1L) %/% 10L + 1L)
  df
}
# Time series without a group: autocorrelated outcome, runs follow time order.
sim_ts <- function(b, ar = 0.9, n = 300) {
  set.seed(b)
  data.frame(t = seq_len(n), batch = rep(paste0("r", 1:12), each = n / 12),
             y = as.numeric(stats::arima.sim(list(ar = ar), n)),
             x1 = stats::rnorm(n), x2 = stats::rnorm(n))
}

# Visits processed in runs by date: 40 subjects x 12 visits, runs A/B/C by
# visit, stationary AR(0.9) outcome per subject (binary: path > 0), no batch
# effect.
sim_longitudinal <- function(b, n_subj = 40, visits = 12) {
  set.seed(b)
  s <- rep(seq_len(n_subj), each = visits)
  path <- unlist(lapply(seq_len(n_subj), function(i) as.numeric(stats::arima.sim(list(ar = 0.9), visits))))
  data.frame(subject = paste0("s", s), batch = rep(rep(c("A", "B", "C"), each = visits / 3), n_subj),
             y = factor(ifelse(path > 0, "case", "ctrl"), levels = c("ctrl", "case")),
             x1 = stats::rnorm(length(s)), x2 = stats::rnorm(length(s)))
}
# Combined plan: leave-processing-batch-out, subjects kept apart as a
# constraint; 60 subjects x 4 replicates; site per subject, null.
sim_combined <- function(b, n_subj = 60) {
  set.seed(b)
  s <- rep(seq_len(n_subj), each = 4)
  data.frame(subject = paste0("s", s), proc = sample(paste0("p", 1:6), length(s), TRUE),
             batch = sample(paste0("site", 1:3), n_subj, TRUE)[s],
             y = factor(sample(c("ctrl", "case"), n_subj, TRUE)[s], levels = c("ctrl", "case")),
             x1 = stats::rnorm(length(s)), x2 = stats::rnorm(length(s)))
}

# 1-2 samples per subject, batch and class drawn per sample (ties in the
# per-subject summary are common).
sim_ties <- function(b, n_subj = 300, classes = c("ctrl", "case")) {
  set.seed(b)
  s <- rep(seq_len(n_subj), sample(1:2, n_subj, TRUE))
  data.frame(subject = paste0("s", s), batch = sample(paste0("b", 1:4), length(s), TRUE),
             y = factor(sample(classes, length(s), TRUE), levels = classes),
             x1 = stats::rnorm(length(s)), x2 = stats::rnorm(length(s)))
}
# 40 tumour/normal pairs (batch = tissue) plus 60 single samples with batch
# and class independent: no association between subjects.
sim_paired_mix <- function(b) {
  set.seed(b)
  pr <- data.frame(subject = rep(paste0("p", 1:40), each = 2),
                   y = rep(c("ctrl", "case"), 40))
  pr$batch <- ifelse(pr$y == "ctrl", "N", "T")
  sg <- data.frame(subject = paste0("q", 1:60), y = sample(c("ctrl", "case"), 60, TRUE),
                   batch = sample(c("N", "T"), 60, TRUE))
  df <- rbind(pr, sg)
  df$y <- factor(df$y, levels = c("ctrl", "case"))
  df$x1 <- stats::rnorm(nrow(df)); df$x2 <- stats::rnorm(nrow(df))
  df
}

rate <- function(f) mean(unlist(parallel::mclapply(seq_len(B), f, mc.cores = cores)))

cat(sprintf("B = %d datasets per design; type-I limit %.3f\n\n", B, lim))

r <- rate(function(b) flag_of(sim(b, 0.5)))
check("L1.3-NULL-BIN", r <= lim, sprintf("binary, batch independent of outcome: flag rate %.3f", r))
r <- rate(function(b) flag_of(sim(b, 0.5, gaussian = TRUE)))
check("L1.3-NULL-GAUSS", r <= lim, sprintf("regression, batch independent of outcome: flag rate %.3f", r))
r <- rate(function(b) flag_of(sim(b, 0.5, two_cols = TRUE), repeats = 3))
check("L1.3-NULL-MULTI", r <= lim, sprintf("two batch columns, 3 repeats, null: flag rate %.3f", r))
r <- rate(function(b) flag_of(sim_clustered(b, 0.5)))
check("L1.3-NULL-CLUSTER", r <= lim, sprintf("40 subjects x 5 rows, subject-level batch and outcome, null: flag rate %.3f", r))
r <- rate(function(b) flag_of(sim_within(b, "random")))
check("L1.3-NULL-WITHIN", r <= lim, sprintf("batch varies within subject (2-6 rows), null: flag rate %.3f", r))
if (requireNamespace("rsample", quietly = TRUE)) {
  r <- rate(function(b) flag_rsample(sim_clustered(b, 0.5)))
  check("L1.3-NULL-RSAMPLE", r <= lim, sprintf("rsample group_vfold_cv, 40 subjects x 5 rows, null: flag rate %.3f", r))
}
r <- rate(function(b) flag_of(sim_unequal(b, 0.5)))
check("L1.3-NULL-UNEQUAL", r <= lim, sprintf("batch varies within subject, sizes 2-40, null: flag rate %.3f", r))
r <- rate(function(b) flag_of(sim_longitudinal(b)))
check("L1.3-NULL-LONGITUDINAL", r <= lim, sprintf("visits in runs by date, AR outcome, null: flag rate %.3f", r))
r <- rate(function(b) flag_of(sim_combined(b), mode = "combined"))
check("L1.3-NULL-COMBINED-SUBJECT", r <= lim, sprintf("combined plan, subject as constraint, site null: flag rate %.3f", r))
r <- rate(function(b) flag_of(sim_ties(b)))
check("L1.3-NULL-TIES", r <= lim, sprintf("1-2 samples per subject, per-sample batch and class, null: flag rate %.3f", r))
r <- rate(function(b) flag_of(sim_paired_mix(b)))
check("L1.3-NULL-PAIRED-MIX", r <= lim, sprintf("tumour/normal pairs + independent singletons: flag rate %.3f", r))
# Multiclass: the test itself (the flag rule after it is shared with binary).
r <- rate(function(b) {
  df <- sim_ties(b, classes = c("a", "b", "c"))
  res <- .batch_outcome_assoc(df$batch, df$y, "multiclass", unit = df$subject, seed = b)
  isTRUE(res$outcome_pval <= alpha && res$outcome_effect >= res$outcome_effect_min)
})
check("L1.3-NULL-MULTICLASS", r <= lim, sprintf("3 classes, 1-2 samples per subject, test-level null rejection %.3f", r))
r <- rate(function(b) flag_of(sim_loocv(b), mode = "study_loocv"))
check("L1.3-NULL-LOOCV-GROUP", r <= lim, sprintf("study_loocv with group, batch crossing studies, null: flag rate %.3f", r))
r <- rate(function(b) flag_of(sim_ts(b), mode = "time_series"))
check("L1.3-NULL-TS", r <= lim, sprintf("time series without group, AR 0.9 outcome, null: flag rate %.3f", r))
r <- rate(function(b) flag_of(sim(b, c(0.2, 0.8)), mode = "batch_blocked"))
check("L1.3-BLOCKED", r == 0, sprintf("batch_blocked, batch strongly linked to outcome: flag rate %.3f (must be 0)", r))
r <- rate(function(b) flag_of(sim(b, c(0.3, 0.7))))
check("L1.3-POWER-BIN", r >= 0.8, sprintf("binary, case rate 0.3 vs 0.7 by batch, shared batches: flag rate %.3f", r))
r <- rate(function(b) flag_of(sim_clustered(b, c(0.1, 0.9), n_subj = 60)))
check("L1.3-POWER-CLUSTER", r >= 0.8, sprintf("60 subjects x 5 rows, case rate 0.1 vs 0.9 by batch: flag rate %.3f", r))
# Outside the contract: within-subject confounding is reported, not flagged.
pw <- parallel::mclapply(seq_len(min(B, 50L)), function(b) flag_of(sim_within(b, "paired")), mc.cores = cores)
n_rep <- sum(vapply(pw, function(z) isTRUE(attr(z, "within") > 0.9), logical(1)))
check("L1.3-WITHIN-REPORTED", !any(unlist(pw)) && n_rep == length(pw),
      sprintf("paired tumour/normal: flagged in %d, row-level V > 0.9 reported in %d of %d datasets",
              sum(unlist(pw)), n_rep, length(pw)))
r <- rate(function(b) flag_of(sim(b, c(0.3, 0.7), gaussian = TRUE)))
check("L1.3-POWER-GAUSS", r >= 0.8, sprintf("regression, mean shift 0.8 SD by batch, shared batches: flag rate %.3f", r))

probe_finish()
