#!/usr/bin/env Rscript
## =================================================================
## Generate manuscript figures from stored RDS artifacts
##
## Produces:
##   paper/fig_sim_main.pdf       (Figure 2: main simulation)
##   paper/fig_target_scan.pdf    (Figure 4: target-scan supplementary)
##   paper/fig_casestudy.pdf      (Figure 5: case study)
##
## All data is read from stored .rds files — no models are re-fitted.
## =================================================================

cat("=== Generating Manuscript Figures ===\n\n")

base_dir <- "paper"

## ---------------------------------------------------------------
## Figure 2: Main simulation summaries (4-panel)
## ---------------------------------------------------------------
sim_file <- file.path(base_dir, "sim_results_all.rds")
if (!file.exists(sim_file)) stop("Missing: ", sim_file)

sim_all <- readRDS(sim_file)
sim <- sim_all[!is.na(sim_all$metric_obs), ]

## Friendly labels
lk_labels <- c(
  none = "Clean",
  subject_overlap = "Subject overlap",
  batch_confounded = "Batch confounding",
  peek_norm = "Peek normalization",
  lookahead = "Look-ahead"
)
lk_colors <- c(
  none = "gray50",
  subject_overlap = "#E69F00",
  batch_confounded = "#56B4E9",
  peek_norm = "#CC79A7",
  lookahead = "#009E73"
)

sim$leakage_label <- lk_labels[sim$leakage]
sim$leakage_label <- factor(sim$leakage_label,
                            levels = lk_labels)

## Helper: aggregate mean and 95% t-CI for a continuous metric
agg_ci <- function(df, group_vars, value_col) {
  keys <- do.call(paste, c(df[, group_vars, drop = FALSE], sep = "||"))
  ordered_keys <- unique(keys)
  out <- do.call(rbind, lapply(ordered_keys, function(k) {
    sub <- df[keys == k, , drop = FALSE]
    x <- sub[[value_col]]
    nseed <- length(x)
    mu <- mean(x)
    se <- sd(x) / sqrt(nseed)
    tval <- qt(0.975, nseed - 1)
    data.frame(sub[1, group_vars, drop = FALSE],
               mean = mu, lo = mu - tval * se, hi = mu + tval * se,
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

## Helper: Wilson 95% CI for a binomial proportion (k successes out of n)
wilson_ci <- function(k, n) {
  z <- 1.959964
  phat <- k / n
  denom <- 1 + z^2 / n
  center <- (phat + z^2 / (2 * n)) / denom
  half <- z * sqrt(phat * (1 - phat) / n + z^2 / (4 * n^2)) / denom
  list(rate = phat, lo = pmax(0, center - half), hi = pmin(1, center + half))
}

## Helper: draw jittered error bar + point for a single mechanism at x position
draw_err <- function(x, lo, hi, col, pch, jitter = 0) {
  xj <- x + jitter
  segments(xj, lo, xj, hi, col = col, lwd = 1.2)
  cap <- 0.04
  segments(xj - cap, lo, xj + cap, lo, col = col, lwd = 1.2)
  segments(xj - cap, hi, xj + cap, hi, col = col, lwd = 1.2)
}

pdf(file.path(base_dir, "fig_sim_main.pdf"), width = 10, height = 9)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1), oma = c(0, 0, 1, 0),
    cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)

lk_types <- names(lk_labels)
K <- length(lk_types)
jitter_step <- 0.08
jitters <- (seq_len(K) - (K + 1) / 2) * jitter_step

## Panel (a): Rejection rate at s = 0 by mechanism and feature dimension
sim_s0 <- sim[sim$s == 0, ]
rej_stats <- do.call(rbind, lapply(split(sim_s0, list(sim_s0$leakage, sim_s0$p),
                                         drop = TRUE), function(sub) {
  x <- sub$p_value
  x <- x[!is.na(x)]
  k <- sum(x < 0.05); n <- length(x)
  wci <- wilson_ci(k, n)
  data.frame(leakage = sub$leakage[1], p = sub$p[1],
             rate = wci$rate * 100, lo = wci$lo * 100, hi = wci$hi * 100,
             stringsAsFactors = FALSE)
}))
rownames(rej_stats) <- NULL

p_vals <- sort(unique(rej_stats$p))
plot(NA, xlim = c(0.5, length(p_vals) + 0.5), ylim = c(0, 105),
     xlab = "Feature dimension (p)", ylab = "Rejection rate (%)",
     main = "(a) Rejection rate at s = 0",
     xaxt = "n", las = 1)
axis(1, at = seq_along(p_vals), labels = p_vals)
abline(h = 5, lty = 2, col = "red", lwd = 1.5)
text(0.55, 1, "5% nominal", col = "red", cex = 0.7, adj = 0)

for (k in seq_along(lk_types)) {
  lk <- lk_types[k]
  sub <- rej_stats[rej_stats$leakage == lk, ]
  idx <- match(p_vals, sub$p)
  y_mean <- sub$rate[idx]; y_lo <- sub$lo[idx]; y_hi <- sub$hi[idx]
  x_pos <- seq_along(p_vals)
  draw_err(x_pos, y_lo, y_hi, col = lk_colors[lk], pch = 15 + k - 1,
           jitter = jitters[k])
  lines(x_pos + jitters[k], y_mean, col = lk_colors[lk], lwd = 1.5)
  points(x_pos + jitters[k], y_mean, pch = 15 + k - 1,
         col = lk_colors[lk], cex = 1.4)
}
legend("right", legend = lk_labels, col = lk_colors,
       pch = 15:19, lwd = 1.5, cex = 0.8, bg = "white")

## Panel (b): AUC by mechanism and signal strength (averaged over n and p)
auc_s_stats <- agg_ci(sim, c("leakage", "s"), "metric_obs")
s_vals <- sort(unique(auc_s_stats$s))

plot(NA, xlim = c(0.5, length(s_vals) + 0.5), ylim = c(0.45, 1.02),
     xlab = "Signal strength (s)", ylab = "Mean observed AUC",
     main = "(b) AUC by signal strength",
     xaxt = "n", las = 1)
axis(1, at = seq_along(s_vals), labels = s_vals)
abline(h = 0.5, lty = 2, col = "gray70")

for (k in seq_along(lk_types)) {
  lk <- lk_types[k]
  sub <- auc_s_stats[auc_s_stats$leakage == lk, ]
  idx <- match(s_vals, sub$s)
  y_mean <- sub$mean[idx]; y_lo <- sub$lo[idx]; y_hi <- sub$hi[idx]
  x_pos <- seq_along(s_vals)
  draw_err(x_pos, y_lo, y_hi, col = lk_colors[lk], pch = 15 + k - 1,
           jitter = jitters[k])
  lines(x_pos + jitters[k], y_mean, col = lk_colors[lk], lwd = 1.5)
  points(x_pos + jitters[k], y_mean, pch = 15 + k - 1,
         col = lk_colors[lk], cex = 1.4)
}
legend("bottomright", legend = lk_labels, col = lk_colors,
       pch = 15:19, lwd = 1.5, cex = 0.8, bg = "white")

## Panel (c): Observed AUC by leakage mechanism and sample size (s > 0)
sim_sig <- sim[sim$s > 0, ]
auc_stats <- agg_ci(sim_sig, c("leakage", "n"), "metric_obs")
n_vals <- sort(unique(auc_stats$n))

plot(NA, xlim = c(0.5, length(n_vals) + 0.5), ylim = c(0.6, 1.02),
     xlab = "Sample size (n)", ylab = "Mean observed AUC",
     main = "(c) Observed AUC by mechanism",
     xaxt = "n", las = 1)
axis(1, at = seq_along(n_vals), labels = n_vals)

for (k in seq_along(lk_types)) {
  lk <- lk_types[k]
  sub <- auc_stats[auc_stats$leakage == lk, ]
  idx <- match(n_vals, sub$n)
  y_mean <- sub$mean[idx]; y_lo <- sub$lo[idx]; y_hi <- sub$hi[idx]
  x_pos <- seq_along(n_vals)
  draw_err(x_pos, y_lo, y_hi, col = lk_colors[lk], pch = 15 + k - 1,
           jitter = jitters[k])
  lines(x_pos + jitters[k], y_mean, col = lk_colors[lk], lwd = 1.5)
  points(x_pos + jitters[k], y_mean, pch = 15 + k - 1,
         col = lk_colors[lk], cex = 1.4)
}
legend("bottomright", legend = lk_labels, col = lk_colors,
       pch = 15:19, lwd = 1.5, cex = 0.8, bg = "white")

## Panel (d): AUC inflation relative to clean baseline (s > 0)
## Per-seed paired differences: delta = leaky_AUC - clean_AUC at matched (n,p,s,seed)
clean_df <- subset(sim_sig, leakage == "none")[, c("seed", "n", "p", "s", "metric_obs")]
names(clean_df)[5] <- "metric_clean"
leaky_df <- subset(sim_sig, leakage != "none")
pair_df <- merge(leaky_df, clean_df, by = c("seed", "n", "p", "s"))
pair_df$delta <- pair_df$metric_obs - pair_df$metric_clean
infl_stats <- agg_ci(pair_df, c("leakage", "n"), "delta")

lk_leak <- setdiff(names(lk_labels), "none")
lk_leak_idx <- seq_along(lk_leak)
jitters_d <- (lk_leak_idx - (length(lk_leak) + 1) / 2) * jitter_step

plot(NA, xlim = c(0.5, length(n_vals) + 0.5), ylim = c(-0.02, 0.35),
     xlab = "Sample size (n)",
     ylab = expression("Mean AUC inflation  (" * Delta * "AUC)"),
     main = "(d) AUC inflation vs clean baseline",
     xaxt = "n", las = 1)
axis(1, at = seq_along(n_vals), labels = n_vals)
abline(h = 0, lty = 2, col = "gray70")

for (k in seq_along(lk_leak)) {
  lk <- lk_leak[k]
  sub <- infl_stats[infl_stats$leakage == lk, ]
  idx <- match(n_vals, sub$n)
  y_mean <- sub$mean[idx]; y_lo <- sub$lo[idx]; y_hi <- sub$hi[idx]
  x_pos <- seq_along(n_vals)
  draw_err(x_pos, y_lo, y_hi, col = lk_colors[lk], pch = 15 + k,
           jitter = jitters_d[k])
  lines(x_pos + jitters_d[k], y_mean, col = lk_colors[lk], lwd = 1.5)
  points(x_pos + jitters_d[k], y_mean, pch = 15 + k,
         col = lk_colors[lk], cex = 1.4)
}
legend("topright", legend = lk_labels[lk_leak], col = lk_colors[lk_leak],
       pch = 16:19, lwd = 1.5, cex = 0.8, bg = "white")

dev.off()
cat("Figure 2 saved: paper/fig_sim_main.pdf\n")

## ---------------------------------------------------------------
## Figure 4: Target-scan supplementary (2-panel)
## ---------------------------------------------------------------
scan_file <- file.path(base_dir, "sim_results", "supplementary_target_scan.rds")
if (!file.exists(scan_file)) {
  cat("WARNING: Missing target scan results, skipping Figure 4\n")
} else {
  scan <- readRDS(scan_file)

  ## Filter to s = 0 only — at s > 0 all mechanisms (including clean) reject
  ## because real signal creates legitimate target association
  scan_s0 <- scan[scan$s == 0, ]

  pdf(file.path(base_dir, "fig_target_scan.pdf"), width = 6, height = 5)
  par(mar = c(6.5, 4.5, 0.5, 3.2), cex.lab = 1.1)

  ## Order levels to match the rest of the paper
  scan_s0$leakage <- factor(scan_s0$leakage,
                            levels = names(lk_labels))

  boxplot(multivariate_p ~ leakage, data = scan_s0,
          xaxt = "n", xlab = "",
          ylab = "Multivariate scan p-value",
          main = "",
          col = lk_colors[levels(scan_s0$leakage)],
          las = 1, outline = FALSE, border = "gray30")
  ## Jittered stripchart overlay: keeps linear scale's honest view of the
  ## uniform null while making individual points visible for the collapsed
  ## leakage boxes (including near-threshold misses). Points share the
  ## mechanism color with a dark outline for visibility against the box fill.
  set.seed(42)
  for (i in seq_along(levels(scan_s0$leakage))) {
    lk <- levels(scan_s0$leakage)[i]
    y_i <- scan_s0$multivariate_p[scan_s0$leakage == lk]
    x_i <- jitter(rep(i, length(y_i)), amount = 0.15)
    points(x_i, y_i, pch = 21,
           bg = adjustcolor(lk_colors[lk], 0.7),
           col = "gray15", cex = 0.7)
  }
  ## Rotated x-axis labels
  axis(1, at = 1:5, labels = FALSE)
  text(1:5, par("usr")[3] - 0.03, labels = lk_labels[levels(scan_s0$leakage)],
       srt = 45, adj = 1, xpd = TRUE, cex = 0.9)
  abline(h = 0.05, lty = 2, col = "red", lwd = 1.5)
  ## Place label in the right margin adjacent to the red line
  text(par("usr")[2] + 0.12, 0.05, "p = 0.05", col = "red", cex = 0.8,
       adj = c(0, 0.5), xpd = NA)

  dev.off()
  cat("Figure 4 saved: paper/fig_target_scan.pdf\n")
}

## ---------------------------------------------------------------
## Figure 5: Case study (2-panel)
## ---------------------------------------------------------------
cs_file <- file.path(base_dir, "casestudy_results.rds")
if (!file.exists(cs_file)) {
  cat("WARNING: Missing case study results, skipping Figure 5\n")
} else {
  cs <- readRDS(cs_file)

  pdf(file.path(base_dir, "fig_casestudy.pdf"), width = 8, height = 4)
  par(mfrow = c(1, 2), mar = c(5, 5.5, 3, 1))

  ## Panel (a): Target-association scores (top features from leaky pipeline)
  ta <- cs$leaky$target_assoc
  if (!is.null(ta) && nrow(ta) > 0) {
    ta <- ta[order(-ta$score), ]
    top_n <- min(15, nrow(ta))
    ta_top <- ta[seq_len(top_n), ]

    ## Mark synthetic injected features with an asterisk
    synthetic <- grepl("^leak_", ta_top$feature)
    labels <- ifelse(synthetic,
                     paste0(ta_top$feature, " *"),
                     as.character(ta_top$feature))

    cols <- ifelse(ta_top$flag, "#CC79A7", "gray60")
    bp2 <- barplot(rev(ta_top$score),
                   horiz = TRUE, las = 1,
                   names.arg = rev(labels),
                   xlab = expression("Target-association score  " *
                                       group("|", "AUC" - 0.5, "|") %*% 2),
                   main = "",
                   col = rev(cols), border = NA,
                   xlim = c(0, 1.15), cex.names = 0.7)
    title(main = "(a)", adj = 0, line = 1, cex.main = 1.1, font.main = 2)
    abline(v = 0.9, lty = 2, col = "red", lwd = 1.5)
    mtext("threshold = 0.9", side = 3, at = 0.9, col = "red",
          cex = 0.65, line = 0.1, adj = 0.5)
    ## Numeric score annotations at bar ends
    text(x = rev(ta_top$score), y = bp2,
         labels = sprintf("%.3f", rev(ta_top$score)),
         pos = 4, cex = 0.65, col = "gray30", offset = 0.25)
    legend(x = 0.58, y = 2.3,
           legend = c("Flagged", "Not flagged", "* synthetic"),
           fill = c("#CC79A7", "gray60", NA),
           border = c(NA, NA, NA),
           cex = 0.65, bg = "white", box.col = "gray80")
  } else {
    plot.new()
    text(0.5, 0.5, "No target scan data available")
  }

  ## Panel (b): Delta LSI repeat-level dot plot
  dlsi_obj <- cs$delta_lsi$delta_lsi_obj
  if (!is.null(dlsi_obj)) {
    rn <- dlsi_obj@repeats_naive
    rg <- dlsi_obj@repeats_guarded
    deltas <- rn$metric - rg$metric
    R_eff <- length(deltas)

    delta_robust <- dlsi_obj@delta_lsi
    delta_mean   <- dlsi_obj@delta_metric
    ci <- dlsi_obj@delta_lsi_ci

    par(mar = c(5, 5.5, 3, 1))
    ylim <- c(-0.02, max(c(deltas, ci[2])) * 1.12)
    plot(seq_len(R_eff), deltas, pch = 19, cex = 1.0, col = "gray40",
         xlab = "Repeat", ylab = expression(Delta[italic(r)] ~ "(leaky " * minus * " guarded AUC)"),
         xlim = c(0.5, R_eff + 0.5), ylim = ylim,
         las = 1, xaxt = "n", main = "")
    axis(1, at = seq(5, R_eff, by = 5))
    title(main = "(b)", adj = 0, line = 1, cex.main = 1.1, font.main = 2)

    ## BCa CI band
    rect(0, ci[1], R_eff + 1, ci[2],
         col = adjustcolor("#0072B2", alpha.f = 0.12), border = NA)

    ## Robust estimate (solid) and mean (dashed)
    abline(h = delta_robust, col = "#0072B2", lwd = 2)
    abline(h = delta_mean, col = "#0072B2", lwd = 1.5, lty = 2)

    ## Zero reference
    abline(h = 0, col = "gray70", lty = 3)

    ## Re-draw points on top
    points(seq_len(R_eff), deltas, pch = 19, cex = 1.0, col = "gray40")

    legend("bottomright",
           legend = c(
             bquote(hat(Delta)[LSI] == .(sprintf("%.3f", delta_robust))),
             bquote(hat(Delta)[metric] == .(sprintf("%.3f", delta_mean))),
             "95% BCa CI"
           ),
           lty = c(1, 2, NA), lwd = c(2, 1.5, NA), col = c("#0072B2", "#0072B2", NA),
           fill = c(NA, NA, adjustcolor("#0072B2", alpha.f = 0.12)),
           border = c(NA, NA, "gray80"),
           cex = 0.65, bg = "white", box.col = "gray80")
  } else {
    par(mar = c(5, 5.5, 3, 1))
    plot.new()
    text(0.5, 0.5, "No delta LSI data available")
  }

  dev.off()
  cat("Figure 5 saved: paper/fig_casestudy.pdf\n")
}

cat("\n=== Figure generation complete ===\n")
