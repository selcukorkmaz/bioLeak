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
  batch_confounded = "Batch confounded",
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

pdf(file.path(base_dir, "fig_sim_main.pdf"), width = 10, height = 9)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1), oma = c(0, 0, 1, 0))

## Panel (a): Observed AUC by leakage and signal level
sim_sig <- sim[sim$s > 0, ]
auc_agg <- aggregate(metric_obs ~ leakage + s, data = sim_sig, FUN = mean)
auc_agg$leakage_label <- lk_labels[auc_agg$leakage]

s_vals <- sort(unique(auc_agg$s))
lk_types <- names(lk_labels)
x_pos <- seq_along(s_vals)

plot(NA, xlim = c(0.5, length(s_vals) + 0.5), ylim = c(0.45, 1.0),
     xlab = "Signal strength (s)", ylab = "Mean observed AUC",
     main = "(a) Observed AUC by leakage and signal",
     xaxt = "n", las = 1)
axis(1, at = x_pos, labels = s_vals)
abline(h = 0.5, lty = 2, col = "gray70")

bar_w <- 0.15
for (k in seq_along(lk_types)) {
  lk <- lk_types[k]
  sub <- auc_agg[auc_agg$leakage == lk, ]
  y_vals <- sub$metric_obs[match(s_vals, sub$s)]
  x_offset <- (k - (length(lk_types) + 1) / 2) * bar_w
  points(x_pos + x_offset, y_vals, pch = 15 + k - 1,
         col = lk_colors[lk], cex = 1.4)
  lines(x_pos + x_offset, y_vals, col = lk_colors[lk], lwd = 1.5)
}
legend("bottomright", legend = lk_labels, col = lk_colors,
       pch = 15:19, lwd = 1.5, cex = 0.7, bg = "white")

## Panel (b): Permutation gap by mechanism (s > 0)
gap_agg <- aggregate(gap ~ leakage + n, data = sim_sig, FUN = mean)
gap_agg$leakage_label <- lk_labels[gap_agg$leakage]

n_vals <- sort(unique(gap_agg$n))
plot(NA, xlim = c(0.5, length(n_vals) + 0.5), ylim = c(-0.02, 0.35),
     xlab = "Sample size (n)", ylab = "Mean permutation gap",
     main = "(b) Permutation gap by mechanism",
     xaxt = "n", las = 1)
axis(1, at = seq_along(n_vals), labels = n_vals)
abline(h = 0, lty = 2, col = "gray70")

for (k in seq_along(lk_types)) {
  lk <- lk_types[k]
  sub <- gap_agg[gap_agg$leakage == lk, ]
  y_vals <- sub$gap[match(n_vals, sub$n)]
  x_offset <- (k - (length(lk_types) + 1) / 2) * bar_w
  points(seq_along(n_vals) + x_offset, y_vals, pch = 15 + k - 1,
         col = lk_colors[lk], cex = 1.4)
  lines(seq_along(n_vals) + x_offset, y_vals, col = lk_colors[lk], lwd = 1.5)
}
legend("topright", legend = lk_labels, col = lk_colors,
       pch = 15:19, lwd = 1.5, cex = 0.7, bg = "white")

## Panel (c): Type I error at s = 0 (clean baseline)
clean_null <- sim[sim$leakage == "none" & sim$s == 0, ]
type1_by_n <- aggregate(p_value ~ n, data = clean_null,
                        FUN = function(x) mean(x < 0.05) * 100)
names(type1_by_n)[2] <- "type1"

barplot(type1_by_n$type1, names.arg = type1_by_n$n,
        xlab = "Sample size (n)", ylab = "Type I error (%)",
        main = "(c) Null-run Type I error (s = 0)",
        col = "gray70", ylim = c(0, max(15, max(type1_by_n$type1) * 1.3)),
        las = 1, border = NA)
abline(h = 5, lty = 2, col = "red", lwd = 1.5)
text(0.5, 5.8, "5% nominal", col = "red", cex = 0.8, adj = 0)

## Panel (d): AUC inflation relative to clean baseline (s > 0)
auc_by_cond <- aggregate(metric_obs ~ leakage + n + s, data = sim_sig, FUN = mean)
clean_auc <- auc_by_cond[auc_by_cond$leakage == "none", ]

lk_leak <- setdiff(names(lk_labels), "none")
inflation_by_n <- data.frame()
for (lk in lk_leak) {
  leaky_sub <- auc_by_cond[auc_by_cond$leakage == lk, ]
  for (nn in n_vals) {
    matched <- merge(
      leaky_sub[leaky_sub$n == nn, c("s", "metric_obs")],
      clean_auc[clean_auc$n == nn, c("s", "metric_obs")],
      by = "s", suffixes = c("_leak", "_clean")
    )
    if (nrow(matched) > 0) {
      delta <- mean(matched$metric_obs_leak - matched$metric_obs_clean)
      inflation_by_n <- rbind(inflation_by_n,
                              data.frame(leakage = lk, n = nn, delta = delta))
    }
  }
}

plot(NA, xlim = c(0.5, length(n_vals) + 0.5), ylim = c(-0.02, 0.35),
     xlab = "Sample size (n)", ylab = "Mean AUC inflation",
     main = "(d) AUC inflation vs clean baseline",
     xaxt = "n", las = 1)
axis(1, at = seq_along(n_vals), labels = n_vals)
abline(h = 0, lty = 2, col = "gray70")

for (k in seq_along(lk_leak)) {
  lk <- lk_leak[k]
  sub <- inflation_by_n[inflation_by_n$leakage == lk, ]
  y_vals <- sub$delta[match(n_vals, sub$n)]
  x_offset <- (k - (length(lk_leak) + 1) / 2) * 0.18
  points(seq_along(n_vals) + x_offset, y_vals, pch = 15 + k,
         col = lk_colors[lk], cex = 1.4)
  lines(seq_along(n_vals) + x_offset, y_vals, col = lk_colors[lk], lwd = 1.5)
}
legend("topright", legend = lk_labels[lk_leak], col = lk_colors[lk_leak],
       pch = 16:19, lwd = 1.5, cex = 0.7, bg = "white")

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

  pdf(file.path(base_dir, "fig_target_scan.pdf"), width = 9, height = 4.5)
  par(mfrow = c(1, 2), mar = c(5, 4.5, 2.5, 1))

  ## Panel (a): Univariate flag rate by mechanism
  flag_rate <- aggregate(leak_feature_flagged ~ leakage, data = scan,
                         FUN = function(x) mean(x, na.rm = TRUE) * 100)
  ## Include "none" for false positive context
  flag_rate$label <- lk_labels[flag_rate$leakage]
  flag_rate <- flag_rate[order(match(flag_rate$leakage, names(lk_labels))), ]

  bp <- barplot(flag_rate$leak_feature_flagged,
                names.arg = flag_rate$label,
                ylab = "Leak feature flag rate (%)",
                main = "(a) Univariate target-scan detection",
                col = lk_colors[flag_rate$leakage],
                ylim = c(0, 110), las = 2, border = NA)
  ## Add percentage text on bars
  text(bp, flag_rate$leak_feature_flagged + 3,
       labels = sprintf("%.0f%%", flag_rate$leak_feature_flagged),
       cex = 0.8)

  ## Panel (b): Multivariate p-values
  scan_leaky <- scan[scan$leakage != "none", ]
  scan_clean <- scan[scan$leakage == "none", ]

  boxplot(multivariate_p ~ leakage, data = scan,
          names = lk_labels[levels(factor(scan$leakage))],
          ylab = "Multivariate scan p-value",
          main = "(b) Multivariate target-scan p-values",
          col = lk_colors[levels(factor(scan$leakage))],
          las = 2, outline = TRUE, border = "gray30")
  abline(h = 0.05, lty = 2, col = "red", lwd = 1.5)
  text(0.6, 0.08, "p = 0.05", col = "red", cex = 0.75, adj = 0)

  dev.off()
  cat("Figure 4 saved: paper/fig_target_scan.pdf\n")
}

## ---------------------------------------------------------------
## Figure 5: Case study (3-panel)
## ---------------------------------------------------------------
cs_file <- file.path(base_dir, "casestudy_results.rds")
if (!file.exists(cs_file)) {
  cat("WARNING: Missing case study results, skipping Figure 5\n")
} else {
  cs <- readRDS(cs_file)

  pdf(file.path(base_dir, "fig_casestudy.pdf"), width = 11, height = 4)
  par(mfrow = c(1, 3), mar = c(5, 4.5, 2.5, 1))

  ## Panel (a): AUC comparison
  auc_vals <- c(cs$guarded$auc, cs$leaky$auc)
  gap_vals <- c(cs$guarded$gap, cs$leaky$gap)
  sd_vals  <- c(cs$guarded$auc_sd, cs$leaky$auc_sd)

  bp <- barplot(auc_vals, names.arg = c("Guarded\n(Study LOOCV)", "Naive\n(Row-wise CV)"),
                ylab = "AUC", main = "(a) Observed AUC",
                col = c("#4DAF4A", "#E41A1C"), ylim = c(0, 1),
                las = 1, border = NA)
  ## Error bars (1 SD)
  arrows(bp, auc_vals - sd_vals, bp, auc_vals + sd_vals,
         angle = 90, code = 3, length = 0.08, lwd = 1.5)
  ## Annotate gap
  text(bp, auc_vals + sd_vals + 0.04,
       labels = sprintf("gap=%.3f\np=%.4f", gap_vals,
                        c(cs$guarded$p_value, cs$leaky$p_value)),
       cex = 0.7)

  ## Panel (b): Target-association scores (top features from leaky pipeline)
  ta <- cs$leaky$target_assoc
  if (!is.null(ta) && nrow(ta) > 0) {
    ta <- ta[order(-ta$score), ]
    top_n <- min(15, nrow(ta))
    ta_top <- ta[seq_len(top_n), ]

    cols <- ifelse(ta_top$flag, "#CC79A7", "gray60")
    bp2 <- barplot(rev(ta_top$score),
                   horiz = TRUE, las = 1,
                   names.arg = rev(ta_top$feature),
                   xlab = "Target-association score",
                   main = "(b) Target scan (naive pipeline)",
                   col = rev(cols), border = NA,
                   xlim = c(0, 1.1), cex.names = 0.6)
    abline(v = 0.9, lty = 2, col = "red", lwd = 1.5)
    text(0.92, 0.5, "threshold", col = "red", cex = 0.7, srt = 90)
    legend("bottomright",
           legend = c("Flagged", "Not flagged"),
           fill = c("#CC79A7", "gray60"),
           cex = 0.7, border = NA, bg = "white")
  } else {
    plot.new()
    text(0.5, 0.5, "No target scan data available")
  }

  ## Panel (c): Duplicate counts
  n_dup_g <- if (!is.null(cs$guarded$duplicates)) nrow(cs$guarded$duplicates) else 0
  n_dup_l <- if (!is.null(cs$leaky$duplicates)) nrow(cs$leaky$duplicates) else 0

  bp3 <- barplot(c(n_dup_g, n_dup_l),
                 names.arg = c("Guarded", "Naive"),
                 ylab = "Near-duplicate pairs",
                 main = "(c) Cross-fold duplicate pairs",
                 col = c("#4DAF4A", "#E41A1C"),
                 las = 1, border = NA)
  text(bp3, c(n_dup_g, n_dup_l),
       labels = format(c(n_dup_g, n_dup_l), big.mark = ","),
       pos = 3, cex = 0.8)

  dev.off()
  cat("Figure 5 saved: paper/fig_casestudy.pdf\n")
}

cat("\n=== Figure generation complete ===\n")
