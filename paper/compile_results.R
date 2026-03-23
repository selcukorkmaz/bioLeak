#!/usr/bin/env Rscript
## =================================================================
## Compile simulation and case study results for manuscript
## Computes all summary statistics and outputs placeholder mappings
## =================================================================

cat("=== Compiling Results ===\n\n")

base_dir <- "paper"

## ---------------------------------------------------------------
## 1. Load simulation results
## ---------------------------------------------------------------
sim_all <- readRDS(file.path(base_dir, "sim_results_all.rds"))
cat(sprintf("Simulation: %d rows, %d NA\n", nrow(sim_all), sum(is.na(sim_all$metric_obs))))

## Remove NAs for analysis
sim <- sim_all[!is.na(sim_all$metric_obs), ]
cat(sprintf("Valid rows: %d\n", nrow(sim)))

## ---------------------------------------------------------------
## 2. Type I error (clean baseline, s=0 only)
## ---------------------------------------------------------------
## Type I error must be computed at s=0 (true null: no signal, no leakage).
## At s>0, "none" correctly detects real signal — that is NOT a false positive.
clean <- sim[sim$leakage == "none" & sim$s == 0, ]
type1_overall <- mean(clean$p_value < 0.05) * 100

type1_by_n <- aggregate(p_value ~ n, data = clean,
                        FUN = function(x) mean(x < 0.05) * 100)
names(type1_by_n)[2] <- "type1_rate"

cat(sprintf("\nType I error overall (s=0, none): %.1f%%\n", type1_overall))
cat("Type I error by n:\n")
print(type1_by_n)

TYPE_I_RATE <- sprintf("%.1f", type1_overall)
if (nrow(type1_by_n) > 0) {
  TYPE_I_LO <- sprintf("%.1f", min(type1_by_n$type1_rate))
  TYPE_I_HI <- sprintf("%.1f", max(type1_by_n$type1_rate))
} else {
  TYPE_I_LO <- TYPE_I_HI <- "NA"
}

## Per-n Type I error rates
get_type1_n <- function(n_val) {
  x <- type1_by_n$type1_rate[type1_by_n$n == n_val]
  if (length(x) == 0) return("NA")
  sprintf("%.1f", x)
}
TYPE_I_N100  <- get_type1_n(100)
TYPE_I_N250  <- get_type1_n(250)
TYPE_I_N500  <- get_type1_n(500)
TYPE_I_N1000 <- get_type1_n(1000)

## ---------------------------------------------------------------
## 3. Permutation-gap rejection rates (Table 2) and baselines
## ---------------------------------------------------------------
## NOTE: With perm_refit=FALSE, the permutation-gap test detects
## prediction-label association — not leakage specifically. At s>0,
## even the clean "none" condition rejects at high rates because the
## model exploits genuine signal. We report:
##   (a) s>0 rejection rates for all leakage types (Table 2)
##   (b) "none" baseline at s>0 (Table 2, for comparison)
##   (c) s=0 rejection rates (leakage-specific, since signal is absent)

## 3a. s>0 rejection rates by leakage and n
leaky <- sim[sim$leakage != "none" & sim$s > 0, ]
power_table <- aggregate(
  p_value ~ leakage + n, data = leaky,
  FUN = function(x) mean(x < 0.05)
)
names(power_table)[3] <- "power"

power_wide <- reshape(power_table, idvar = "leakage", timevar = "n",
                      direction = "wide")
cat("\nPermutation-gap rejection rates (s > 0):\n")
print(power_wide)

get_power <- function(lk, n_val) {
  x <- power_table[power_table$leakage == lk & power_table$n == n_val, "power"]
  if (length(x) == 0) return("NA")
  sprintf("%.2f", x)
}

SO_100  <- get_power("subject_overlap", 100)
SO_250  <- get_power("subject_overlap", 250)
SO_500  <- get_power("subject_overlap", 500)
SO_1000 <- get_power("subject_overlap", 1000)

BC_100  <- get_power("batch_confounded", 100)
BC_250  <- get_power("batch_confounded", 250)
BC_500  <- get_power("batch_confounded", 500)
BC_1000 <- get_power("batch_confounded", 1000)

PP_100  <- get_power("peek_norm", 100)
PP_250  <- get_power("peek_norm", 250)
PP_500  <- get_power("peek_norm", 500)
PP_1000 <- get_power("peek_norm", 1000)

TL_100  <- get_power("lookahead", 100)
TL_250  <- get_power("lookahead", 250)
TL_500  <- get_power("lookahead", 500)
TL_1000 <- get_power("lookahead", 1000)

## 3b. "none" baseline at s>0 (clean pipeline, real signal present)
none_sig <- sim[sim$leakage == "none" & sim$s > 0, ]
none_baseline <- aggregate(p_value ~ n, data = none_sig,
                           FUN = function(x) mean(x < 0.05))
names(none_baseline)[2] <- "reject"
NONE_OVERALL <- sprintf("%.1f", mean(none_sig$p_value < 0.05) * 100)
NONE_100  <- sprintf("%.2f", none_baseline$reject[none_baseline$n == 100])
NONE_250  <- sprintf("%.2f", none_baseline$reject[none_baseline$n == 250])
NONE_500  <- sprintf("%.2f", none_baseline$reject[none_baseline$n == 500])
NONE_1000 <- sprintf("%.2f", none_baseline$reject[none_baseline$n == 1000])

cat(sprintf("\nNone baseline at s>0: overall=%.1f%%, by n: %s %s %s %s\n",
            mean(none_sig$p_value < 0.05) * 100,
            NONE_100, NONE_250, NONE_500, NONE_1000))

## 3c. s=0 leakage-specific detection (no real signal present)
s0_leaky <- sim[sim$s == 0 & sim$leakage != "none", ]
s0_power <- aggregate(p_value ~ leakage + n, data = s0_leaky,
                      FUN = function(x) mean(x < 0.05))
names(s0_power)[3] <- "power"
cat("\nLeakage-specific detection at s=0:\n")
print(reshape(s0_power, idvar = "leakage", timevar = "n", direction = "wide"))

get_s0 <- function(lk, n_val) {
  x <- s0_power[s0_power$leakage == lk & s0_power$n == n_val, "power"]
  if (length(x) == 0) return("NA")
  sprintf("%.2f", x)
}

S0_SO_100  <- get_s0("subject_overlap", 100)
S0_SO_250  <- get_s0("subject_overlap", 250)
S0_SO_500  <- get_s0("subject_overlap", 500)
S0_SO_1000 <- get_s0("subject_overlap", 1000)

S0_BC_100  <- get_s0("batch_confounded", 100)
S0_BC_250  <- get_s0("batch_confounded", 250)
S0_BC_500  <- get_s0("batch_confounded", 500)
S0_BC_1000 <- get_s0("batch_confounded", 1000)

S0_PP_100  <- get_s0("peek_norm", 100)
S0_PP_250  <- get_s0("peek_norm", 250)
S0_PP_500  <- get_s0("peek_norm", 500)
S0_PP_1000 <- get_s0("peek_norm", 1000)

S0_TL_100  <- get_s0("lookahead", 100)
S0_TL_250  <- get_s0("lookahead", 250)
S0_TL_500  <- get_s0("lookahead", 500)
S0_TL_1000 <- get_s0("lookahead", 1000)

## Minimum power at s=2.0, n>=500
high_s <- sim[sim$leakage != "none" & sim$s == 2.0 & sim$n >= 500, ]
min_power_s2 <- aggregate(p_value ~ leakage, data = high_s,
                          FUN = function(x) mean(x < 0.05))
MIN_POWER_S2 <- if (nrow(min_power_s2) > 0) sprintf("%.2f", min(min_power_s2$p_value)) else "NA"

## ---------------------------------------------------------------
## 4. AUC inflation (s>0 only — at s=0 AUC ≈ 0.5 for all conditions)
## ---------------------------------------------------------------
sim_sig <- sim[sim$s > 0, ]
auc_by_cond <- aggregate(metric_obs ~ leakage + n + s, data = sim_sig, FUN = mean)

clean_auc <- auc_by_cond[auc_by_cond$leakage == "none", ]
CLEAN_AUC <- sprintf("%.2f", mean(clean_auc$metric_obs))

## Delta AUC: mean(leaky) - mean(clean) at matched (n, s)
delta_fn <- function(lk) {
  leaky_auc <- auc_by_cond[auc_by_cond$leakage == lk, ]
  deltas <- c()
  for (i in seq_len(nrow(leaky_auc))) {
    matched_clean <- clean_auc[clean_auc$n == leaky_auc$n[i] &
                               clean_auc$s == leaky_auc$s[i], "metric_obs"]
    if (length(matched_clean) > 0) {
      deltas <- c(deltas, leaky_auc$metric_obs[i] - matched_clean)
    }
  }
  mean(deltas, na.rm = TRUE)
}

DELTA_SO <- sprintf("%.3f", delta_fn("subject_overlap"))
DELTA_BC <- sprintf("%.3f", delta_fn("batch_confounded"))
DELTA_PP <- sprintf("%.3f", delta_fn("peek_norm"))
DELTA_TL <- sprintf("%.3f", delta_fn("lookahead"))

cat(sprintf("\nAUC inflation: SO=%s, BC=%s, PP=%s, TL=%s\n",
            DELTA_SO, DELTA_BC, DELTA_PP, DELTA_TL))

## ---------------------------------------------------------------
## 5. Supplementary: splitting modes (Table 3)
## ---------------------------------------------------------------
mode_results <- NULL
mode_file <- file.path(base_dir, "sim_results", "supplementary_modes.rds")
if (file.exists(mode_file)) {
  mode_results <- readRDS(mode_file)
  mode_results <- mode_results[!is.na(mode_results$metric_obs), ]

  ## Use s=0 rows if available (leakage-specific); fall back to all rows
  if ("s" %in% names(mode_results)) {
    mode_s0 <- mode_results[mode_results$s == 0, ]
    cat(sprintf("\nSplitting modes: %d rows total, %d at s=0\n",
                nrow(mode_results), nrow(mode_s0)))
  } else {
    mode_s0 <- mode_results
    cat("\nSplitting modes: legacy format (no s column), using all rows\n")
  }

  ## Exclude "none" — it's the clean baseline, not a leakage type
  mode_leaky <- mode_s0[mode_s0$leakage != "none", ]

  mode_power <- aggregate(p_value ~ mode + leakage, data = mode_leaky,
                          FUN = function(x) mean(x < 0.05))
  names(mode_power)[3] <- "power"
  cat("\nSplitting mode detection (s=0, leakage-specific):\n")
  print(reshape(mode_power, idvar = "mode", timevar = "leakage", direction = "wide"))

  ## Also check "none" baseline at s=0
  mode_none <- mode_s0[mode_s0$leakage == "none", ]
  if (nrow(mode_none) > 0) {
    none_by_mode <- aggregate(p_value ~ mode, data = mode_none,
                              FUN = function(x) mean(x < 0.05) * 100)
    names(none_by_mode)[2] <- "reject_pct"
    cat("\nNone baseline by mode (s=0):\n")
    print(none_by_mode)
  }

  get_mode_power <- function(md, lk) {
    x <- mode_power[mode_power$mode == md & mode_power$leakage == lk, "power"]
    if (length(x) == 0) return("NA")
    sprintf("%.2f", x)
  }

  SG_SO <- get_mode_power("subject_grouped", "subject_overlap")
  SG_BC <- get_mode_power("subject_grouped", "batch_confounded")
  SG_PP <- get_mode_power("subject_grouped", "peek_norm")
  SG_TL <- get_mode_power("subject_grouped", "lookahead")

  BB_SO <- get_mode_power("batch_blocked", "subject_overlap")
  BB_BC <- get_mode_power("batch_blocked", "batch_confounded")
  BB_PP <- get_mode_power("batch_blocked", "peek_norm")
  BB_TL <- get_mode_power("batch_blocked", "lookahead")

  SL_SO <- get_mode_power("study_loocv", "subject_overlap")
  SL_BC <- get_mode_power("study_loocv", "batch_confounded")
  SL_PP <- get_mode_power("study_loocv", "peek_norm")
  SL_TL <- get_mode_power("study_loocv", "lookahead")

  TS_SO <- get_mode_power("time_series", "subject_overlap")
  TS_BC <- get_mode_power("time_series", "batch_confounded")
  TS_PP <- get_mode_power("time_series", "peek_norm")
  TS_TL <- get_mode_power("time_series", "lookahead")
} else {
  cat("\nWARNING: Supplementary mode results not found. Using placeholders.\n")
  SG_SO <- SG_BC <- SG_PP <- SG_TL <- "NA"
  BB_SO <- BB_BC <- BB_PP <- BB_TL <- "NA"
  SL_SO <- SL_BC <- SL_PP <- SL_TL <- "NA"
  TS_SO <- TS_BC <- TS_PP <- TS_TL <- "NA"
}

## ---------------------------------------------------------------
## 6. Target leakage scan
## ---------------------------------------------------------------
## The target scan now runs at both s=0 and s=1.0.
## Use s=0 for leakage-detection metrics: at s=0 the real features are
## pure noise, so the ONLY source of target association is the leak feature.
## This lets us properly evaluate what the scan can and cannot catch.
## At s=1.0, the multivariate scan always rejects (real signal dominates),
## so those results are not leakage-specific.
scan_file <- file.path(base_dir, "sim_results", "supplementary_target_scan.rds")
if (file.exists(scan_file)) {
  scan_results <- readRDS(scan_file)

  ## Select s=0 rows if the s column exists (backward compatible)
  if ("s" %in% names(scan_results)) {
    scan_s0 <- scan_results[scan_results$s == 0, ]
    cat(sprintf("\nTarget scan: %d rows total, %d at s=0\n",
                nrow(scan_results), nrow(scan_s0)))
  } else {
    scan_s0 <- scan_results
    cat("\nTarget scan: legacy format (no s column), using all rows\n")
  }

  ## Univariate detection rate (leak feature flagged) across all leakage types
  scan_leaky <- scan_s0[scan_s0$leakage != "none", ]
  UNI_DETECT <- sprintf("%.0f", mean(scan_leaky$leak_feature_flagged, na.rm = TRUE) * 100)

  ## By type
  uni_by_type <- aggregate(leak_feature_flagged ~ leakage, data = scan_leaky,
                           FUN = function(x) mean(x, na.rm = TRUE) * 100)
  get_uni <- function(lk) {
    x <- uni_by_type$leak_feature_flagged[uni_by_type$leakage == lk]
    if (length(x) == 0) return("NA")
    sprintf("%.0f", x)
  }
  UNI_SO <- get_uni("subject_overlap")
  UNI_PP <- get_uni("peek_norm")

  ## Multivariate detection rate (leaky conditions at s=0)
  multi_valid <- scan_leaky$multivariate_p[!is.na(scan_leaky$multivariate_p)]
  if (length(multi_valid) > 0) {
    multi_detect <- mean(multi_valid < 0.05) * 100
    MULTI_DETECT <- sprintf("%.0f", multi_detect)
  } else {
    MULTI_DETECT <- "NA"
  }

  ## Multivariate false-positive rate ("none" at s=0 — honest null baseline)
  scan_none <- scan_s0[scan_s0$leakage == "none", ]
  multi_none_p <- scan_none$multivariate_p[!is.na(scan_none$multivariate_p)]
  if (length(multi_none_p) > 0) {
    MULTI_NONE <- sprintf("%.0f", mean(multi_none_p < 0.05) * 100)
  } else {
    MULTI_NONE <- "NA"
  }

  cat(sprintf("\nTarget scan (s=0): univariate=%s%%, multivariate=%s%% (none baseline=%s%%)\n",
              UNI_DETECT, MULTI_DETECT, MULTI_NONE))
  cat(sprintf("  SO=%s%%, PP=%s%%\n", UNI_SO, UNI_PP))
} else {
  cat("\nWARNING: Target scan results not found. Using placeholders.\n")
  UNI_DETECT <- MULTI_DETECT <- UNI_SO <- UNI_PP <- MULTI_NONE <- "NA"
}

## ---------------------------------------------------------------
## 7. Case study
## ---------------------------------------------------------------
cs_file <- file.path(base_dir, "casestudy_results.rds")
if (file.exists(cs_file)) {
  cs <- readRDS(cs_file)

  CS_N <- as.character(cs$N)
  CS_S <- as.character(cs$S)
  CS_G <- as.character(cs$G)

  GUARDED_AUC <- sprintf("%.3f", cs$guarded$auc)
  G_SD        <- sprintf("%.3f", cs$guarded$auc_sd)
  LEAKY_AUC   <- sprintf("%.3f", cs$leaky$auc)
  L_SD        <- sprintf("%.3f", cs$leaky$auc_sd)
  DELTA_AUC   <- sprintf("%.3f", cs$leaky$auc - cs$guarded$auc)

  G_GAP  <- sprintf("%.3f", cs$guarded$gap)
  G_PVAL <- sprintf("%.4f", cs$guarded$p_value)
  L_GAP  <- sprintf("%.3f", cs$leaky$gap)
  L_PVAL <- sprintf("%.4f", cs$leaky$p_value)

  ## Batch association
  if (!is.null(cs$guarded$batch_assoc) && nrow(cs$guarded$batch_assoc) > 0) {
    G_CRAMER  <- sprintf("%.3f", cs$guarded$batch_assoc$cramer_v[1])
    G_BATCH_P <- sprintf("%.4f", cs$guarded$batch_assoc$pval[1])
  } else {
    G_CRAMER <- "NA"; G_BATCH_P <- "NA"
  }
  if (!is.null(cs$leaky$batch_assoc) && nrow(cs$leaky$batch_assoc) > 0) {
    L_CRAMER  <- sprintf("%.3f", cs$leaky$batch_assoc$cramer_v[1])
    L_BATCH_P <- sprintf("%.4f", cs$leaky$batch_assoc$pval[1])
  } else {
    L_CRAMER <- "NA"; L_BATCH_P <- "NA"
  }

  ## Target scores for leaky features
  ta <- cs$leaky$target_assoc
  flagged_leaky <- ta[grepl("^leak", ta$feature) & ta$flag == TRUE, ]
  if (nrow(flagged_leaky) >= 3) {
    ## Sort by feature name for consistency
    flagged_leaky <- flagged_leaky[order(flagged_leaky$feature), ]
    SCORE_1 <- sprintf("%.3f", flagged_leaky$score[1])
    SCORE_2 <- sprintf("%.3f", flagged_leaky$score[2])
    SCORE_3 <- sprintf("%.3f", flagged_leaky$score[3])
  } else if (nrow(flagged_leaky) > 0) {
    scores <- sprintf("%.3f", flagged_leaky$score)
    SCORE_1 <- scores[1]
    SCORE_2 <- if (length(scores) >= 2) scores[2] else "< 0.9"
    SCORE_3 <- if (length(scores) >= 3) scores[3] else "< 0.9"
  } else {
    SCORE_1 <- SCORE_2 <- SCORE_3 <- "< 0.9"
  }

  ## Multivariate p-values
  MULTI_P_LEAKY   <- if (!is.null(cs$leaky$multi_p))
    sprintf("%.4f", cs$leaky$multi_p) else "NA"
  MULTI_P_GUARDED <- if (!is.null(cs$guarded$multi_p))
    sprintf("%.4f", cs$guarded$multi_p) else "NA"

  ## Duplicates
  n_dup_leaky <- if (!is.null(cs$leaky$duplicates)) nrow(cs$leaky$duplicates) else 0L
  n_dup_guarded <- if (!is.null(cs$guarded$duplicates)) nrow(cs$guarded$duplicates) else 0L
  N_DUPS <- as.character(n_dup_leaky)
  # Describe duplicates
  if (n_dup_leaky > 0) {
    DUP_DESC <- "technical replicates or near-identical expression profiles from the same patient across studies"
  } else if (n_dup_guarded > 0) {
    N_DUPS <- as.character(n_dup_guarded)
    DUP_DESC <- "near-identical expression profiles spanning study boundaries"
  } else {
    N_DUPS <- "0"
    DUP_DESC <- "no near-duplicate pairs were detected above the threshold"
  }

  cat(sprintf("\nCase study: N=%s, S=%s, G=%s\n", CS_N, CS_S, CS_G))
  cat(sprintf("Guarded: AUC=%s, gap=%s, p=%s\n", GUARDED_AUC, G_GAP, G_PVAL))
  cat(sprintf("Leaky:   AUC=%s, gap=%s, p=%s\n", LEAKY_AUC, L_GAP, L_PVAL))
} else {
  cat("\nWARNING: Case study results not found. Using placeholders.\n")
  CS_N <- CS_S <- CS_G <- "NA"
  GUARDED_AUC <- G_SD <- LEAKY_AUC <- L_SD <- DELTA_AUC <- "NA"
  G_GAP <- G_PVAL <- L_GAP <- L_PVAL <- "NA"
  G_CRAMER <- G_BATCH_P <- L_CRAMER <- L_BATCH_P <- "NA"
  SCORE_1 <- SCORE_2 <- SCORE_3 <- "NA"
  MULTI_P_LEAKY <- MULTI_P_GUARDED <- "NA"
  N_DUPS <- "0"
  DUP_DESC <- "no near-duplicate pairs were detected"
}

## ---------------------------------------------------------------
## 8. Delta LSI (case study + simulation)
## ---------------------------------------------------------------
## 8a. Case study delta_lsi (appended to casestudy_results.rds by run_delta_lsi.R)
if (file.exists(cs_file)) {
  cs <- readRDS(cs_file)
  if (!is.null(cs$delta_lsi)) {
    dl <- cs$delta_lsi
    DLSI_METRIC  <- sprintf("%.3f", dl$delta_metric)
    DLSI_ROBUST  <- sprintf("%.3f", dl$delta_lsi_robust)
    DLSI_PVAL    <- if (!is.null(dl$p_value) && is.finite(dl$p_value)) {
      if (dl$p_value < 0.001) "< 0.001" else sprintf("= %.4f", dl$p_value)
    } else "NA"
    DLSI_TIER    <- if (!is.null(dl$tier)) as.character(dl$tier) else "NA"
    DLSI_REFF    <- if (!is.null(dl$n_repeats)) as.character(dl$n_repeats) else "NA"
    DLSI_CI_LO   <- if (!is.null(dl$ci_lsi) && length(dl$ci_lsi) == 2 && all(is.finite(dl$ci_lsi)))
      sprintf("%.3f", dl$ci_lsi[1]) else "NA"
    DLSI_CI_HI   <- if (!is.null(dl$ci_lsi) && length(dl$ci_lsi) == 2 && all(is.finite(dl$ci_lsi)))
      sprintf("%.3f", dl$ci_lsi[2]) else "NA"

    cat(sprintf("\nDelta LSI (case study): metric=%.3f, robust=%.3f, p=%s, tier=%s, R_eff=%s\n",
                dl$delta_metric, dl$delta_lsi_robust, DLSI_PVAL, DLSI_TIER, DLSI_REFF))
  } else {
    cat("\nWARNING: delta_lsi not found in casestudy_results.rds. Run run_delta_lsi.R first.\n")
    DLSI_METRIC <- DLSI_ROBUST <- DLSI_PVAL <- DLSI_TIER <- DLSI_REFF <- "NA"
    DLSI_CI_LO <- DLSI_CI_HI <- "NA"
  }
} else {
  DLSI_METRIC <- DLSI_ROBUST <- DLSI_PVAL <- DLSI_TIER <- DLSI_REFF <- "NA"
  DLSI_CI_LO <- DLSI_CI_HI <- "NA"
}

## 8b. Delta LSI simulation
dlsi_sim_file <- file.path(base_dir, "sim_results", "delta_lsi_sim.rds")
if (file.exists(dlsi_sim_file)) {
  dlsi_sim <- readRDS(dlsi_sim_file)

  ## Null condition: type I error
  null_sub <- dlsi_sim[dlsi_sim$leakage == "none" & !is.na(dlsi_sim$p_value), ]
  if (nrow(null_sub) > 0) {
    DLSI_SIM_NULL_REJECT <- sprintf("%.1f", mean(null_sub$p_value < 0.05) * 100)
    DLSI_SIM_NULL_MEAN   <- sprintf("%.4f", mean(null_sub$delta_metric, na.rm = TRUE))
  } else {
    DLSI_SIM_NULL_REJECT <- DLSI_SIM_NULL_MEAN <- "NA"
  }

  ## Alternative condition: power
  alt_sub <- dlsi_sim[dlsi_sim$leakage == "peek_norm" & !is.na(dlsi_sim$p_value), ]
  if (nrow(alt_sub) > 0) {
    DLSI_SIM_ALT_REJECT <- sprintf("%.1f", mean(alt_sub$p_value < 0.05) * 100)
    DLSI_SIM_ALT_MEAN   <- sprintf("%.4f", mean(alt_sub$delta_metric, na.rm = TRUE))
  } else {
    DLSI_SIM_ALT_REJECT <- DLSI_SIM_ALT_MEAN <- "NA"
  }

  cat(sprintf("\nDelta LSI sim: null reject=%.1f%%, alt reject=%.1f%%\n",
              as.numeric(DLSI_SIM_NULL_REJECT), as.numeric(DLSI_SIM_ALT_REJECT)))
} else {
  cat("\nWARNING: delta_lsi_sim.rds not found. Run run_delta_lsi.R first.\n")
  DLSI_SIM_NULL_REJECT <- DLSI_SIM_NULL_MEAN <- "NA"
  DLSI_SIM_ALT_REJECT <- DLSI_SIM_ALT_MEAN <- "NA"
}

## ---------------------------------------------------------------
## 9. Output placeholder mapping
## ---------------------------------------------------------------
cat("\n\n======================================================\n")
cat("PLACEHOLDER -> VALUE MAPPING\n")
cat("======================================================\n\n")

placeholders <- list(
  ## Section 3
  "[PLATFORM]" = sprintf("macOS (%d parallel workers)", max(1L, min(4L, parallel::detectCores(logical = FALSE) - 1L))),

  ## Section 4.1.1
  "[TYPE_I_RATE]" = TYPE_I_RATE,
  "[TYPE_I_LO]" = TYPE_I_LO,
  "[TYPE_I_HI]" = TYPE_I_HI,
  "[TYPE_I_N100]" = TYPE_I_N100,
  "[TYPE_I_N250]" = TYPE_I_N250,
  "[TYPE_I_N500]" = TYPE_I_N500,
  "[TYPE_I_N1000]" = TYPE_I_N1000,

  ## Section 4.1.2 - Table 2 (s>0 rejection rates)
  "[NONE_OVERALL]" = NONE_OVERALL,
  "[NONE_100]" = NONE_100,
  "[NONE_250]" = NONE_250,
  "[NONE_500]" = NONE_500,
  "[NONE_1000]" = NONE_1000,
  "[SO_100]" = SO_100,
  "[SO_250]" = SO_250,
  "[SO_500]" = SO_500,
  "[SO_1000]" = SO_1000,
  "[BC_100]" = BC_100,
  "[BC_250]" = BC_250,
  "[BC_500]" = BC_500,
  "[BC_1000]" = BC_1000,
  "[PP_100]" = PP_100,
  "[PP_250]" = PP_250,
  "[PP_500]" = PP_500,
  "[PP_1000]" = PP_1000,
  "[TL_100]" = TL_100,
  "[TL_250]" = TL_250,
  "[TL_500]" = TL_500,
  "[TL_1000]" = TL_1000,
  "[MIN_POWER_S2]" = MIN_POWER_S2,

  ## Section 4.1.2 - s=0 leakage-specific detection
  "[S0_SO_100]" = S0_SO_100,
  "[S0_SO_250]" = S0_SO_250,
  "[S0_SO_500]" = S0_SO_500,
  "[S0_SO_1000]" = S0_SO_1000,
  "[S0_BC_100]" = S0_BC_100,
  "[S0_BC_250]" = S0_BC_250,
  "[S0_BC_500]" = S0_BC_500,
  "[S0_BC_1000]" = S0_BC_1000,
  "[S0_PP_100]" = S0_PP_100,
  "[S0_PP_250]" = S0_PP_250,
  "[S0_PP_500]" = S0_PP_500,
  "[S0_PP_1000]" = S0_PP_1000,
  "[S0_TL_100]" = S0_TL_100,
  "[S0_TL_250]" = S0_TL_250,
  "[S0_TL_500]" = S0_TL_500,
  "[S0_TL_1000]" = S0_TL_1000,

  ## Section 4.1.3
  "[CLEAN_AUC]" = CLEAN_AUC,
  "[DELTA_SO]" = DELTA_SO,
  "[DELTA_PP]" = DELTA_PP,
  "[DELTA_BC]" = DELTA_BC,
  "[DELTA_TL]" = DELTA_TL,

  ## Section 4.1.5 - Table 3
  "[SG_SO]" = SG_SO,
  "[SG_BC]" = SG_BC,
  "[SG_PP]" = SG_PP,
  "[SG_TL]" = SG_TL,
  "[BB_SO]" = BB_SO,
  "[BB_BC]" = BB_BC,
  "[BB_PP]" = BB_PP,
  "[BB_TL]" = BB_TL,
  "[SL_SO]" = SL_SO,
  "[SL_BC]" = SL_BC,
  "[SL_PP]" = SL_PP,
  "[SL_TL]" = SL_TL,
  "[TS_SO]" = TS_SO,
  "[TS_BC]" = TS_BC,
  "[TS_PP]" = TS_PP,
  "[TS_TL]" = TS_TL,

  ## Section 4.1.6
  "[UNI_DETECT]" = UNI_DETECT,
  "[UNI_SO]" = UNI_SO,
  "[UNI_PP]" = UNI_PP,
  "[MULTI_DETECT]" = MULTI_DETECT,
  "[MULTI_NONE]" = MULTI_NONE,

  ## Section 4.2
  "[N]" = CS_N,
  "[S]" = CS_S,
  "[G]" = CS_G,
  "[GUARDED_AUC]" = GUARDED_AUC,
  "[G_SD]" = G_SD,
  "[LEAKY_AUC]" = LEAKY_AUC,
  "[L_SD]" = L_SD,
  "[DELTA_AUC]" = DELTA_AUC,
  "[G_GAP]" = G_GAP,
  "[G_PVAL]" = G_PVAL,
  "[L_GAP]" = L_GAP,
  "[L_PVAL]" = L_PVAL,
  "[G_CRAMER]" = G_CRAMER,
  "[G_BATCH_P]" = G_BATCH_P,
  "[L_CRAMER]" = L_CRAMER,
  "[L_BATCH_P]" = L_BATCH_P,
  "[SCORE_1]" = SCORE_1,
  "[SCORE_2]" = SCORE_2,
  "[SCORE_3]" = SCORE_3,
  "[MULTI_P_LEAKY]" = MULTI_P_LEAKY,
  "[MULTI_P_GUARDED]" = MULTI_P_GUARDED,
  "[N_DUPS]" = N_DUPS,
  "[DUP_DESC -- e.g., technical replicates or near-identical\\nexpression profiles from the same patient across studies]" = DUP_DESC,

  ## Section 4.2 - Delta LSI (case study)
  "[DLSI_METRIC]" = DLSI_METRIC,
  "[DLSI_ROBUST]" = DLSI_ROBUST,
  "[DLSI_PVAL]" = DLSI_PVAL,
  "[DLSI_TIER]" = DLSI_TIER,
  "[DLSI_REFF]" = DLSI_REFF,
  "[DLSI_CI_LO]" = DLSI_CI_LO,
  "[DLSI_CI_HI]" = DLSI_CI_HI,

  ## Appendix - Delta LSI simulation
  "[DLSI_SIM_NULL_REJECT]" = DLSI_SIM_NULL_REJECT,
  "[DLSI_SIM_NULL_MEAN]" = DLSI_SIM_NULL_MEAN,
  "[DLSI_SIM_ALT_REJECT]" = DLSI_SIM_ALT_REJECT,
  "[DLSI_SIM_ALT_MEAN]" = DLSI_SIM_ALT_MEAN
)

for (ph in names(placeholders)) {
  cat(sprintf("  %s -> %s\n", ph, placeholders[[ph]]))
}

## Save for manuscript update script
saveRDS(placeholders, file.path(base_dir, "placeholder_values.rds"))
cat("\nPlaceholder values saved to paper/placeholder_values.rds\n")
cat("=== Done ===\n")
