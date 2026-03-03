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

## ---------------------------------------------------------------
## 3. Detection power (Table 2)
## ---------------------------------------------------------------
## Use s>0 data for detection power — at s=0, some leakage types
## (e.g., lookahead) are correctly undetectable, which would dilute power.
leaky <- sim[sim$leakage != "none" & sim$s > 0, ]

## Power pooled across p and s>0, by leakage and n
power_table <- aggregate(
  p_value ~ leakage + n, data = leaky,
  FUN = function(x) mean(x < 0.05)
)
names(power_table)[3] <- "power"

## Reshape to wide format
power_wide <- reshape(power_table, idvar = "leakage", timevar = "n",
                      direction = "wide")
cat("\nDetection power (proportion p < 0.05):\n")
print(power_wide)

## Extract individual cells
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

  ## Exclude "none" — it's the clean baseline, not a leakage type
  mode_leaky <- mode_results[mode_results$leakage != "none", ]

  mode_power <- aggregate(p_value ~ mode + leakage, data = mode_leaky,
                          FUN = function(x) mean(x < 0.05))
  names(mode_power)[3] <- "power"
  cat("\nSplitting mode power:\n")
  print(reshape(mode_power, idvar = "mode", timevar = "leakage", direction = "wide"))

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
scan_file <- file.path(base_dir, "sim_results", "supplementary_target_scan.rds")
if (file.exists(scan_file)) {
  scan_results <- readRDS(scan_file)

  ## Univariate detection rate (leak feature flagged) across all leakage types
  scan_leaky <- scan_results[scan_results$leakage != "none", ]
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

  ## Multivariate detection rate
  multi_valid <- scan_leaky$multivariate_p[!is.na(scan_leaky$multivariate_p)]
  if (length(multi_valid) > 0) {
    multi_detect <- mean(multi_valid < 0.05) * 100
    MULTI_DETECT <- sprintf("%.0f", multi_detect)
  } else {
    MULTI_DETECT <- "NA"
  }

  cat(sprintf("\nTarget scan: univariate=%s%%, multivariate=%s%%\n", UNI_DETECT, MULTI_DETECT))
  cat(sprintf("  SO=%s%%, PP=%s%%\n", UNI_SO, UNI_PP))
} else {
  cat("\nWARNING: Target scan results not found. Using placeholders.\n")
  UNI_DETECT <- MULTI_DETECT <- UNI_SO <- UNI_PP <- "NA"
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
## 8. Output placeholder mapping
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

  ## Section 4.1.2 - Table 2
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
  "[DUP_DESC -- e.g., technical replicates or near-identical\\nexpression profiles from the same patient across studies]" = DUP_DESC
)

for (ph in names(placeholders)) {
  cat(sprintf("  %s -> %s\n", ph, placeholders[[ph]]))
}

## Save for manuscript update script
saveRDS(placeholders, file.path(base_dir, "placeholder_values.rds"))
cat("\nPlaceholder values saved to paper/placeholder_values.rds\n")
cat("=== Done ===\n")
