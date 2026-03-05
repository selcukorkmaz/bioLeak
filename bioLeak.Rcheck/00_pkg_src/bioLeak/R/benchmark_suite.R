#' Simulation benchmark matrix for leakage diagnostics
#'
#' Runs a reproducible grid of simulation scenarios across modalities, leakage
#' mechanisms, and split modes using [simulate_leakage_suite()]. This function
#' is designed as a benchmarking harness to quantify detection rates and
#' performance inflation under controlled settings.
#'
#' @param modalities Character vector selecting predefined modality profiles.
#'   Supported values: `"omics"`, `"imaging_tabular"`, `"ehr_tabular"`.
#' @param leakages Character vector of leakage mechanisms passed to
#'   [simulate_leakage_suite()].
#' @param modes Character vector of split modes passed to [simulate_leakage_suite()].
#' @param learner Character scalar. `"glmnet"` (default) or `"ranger"`.
#' @param seeds Integer vector of Monte Carlo seeds.
#' @param B Integer scalar. Number of permutations per scenario.
#' @param alpha Numeric scalar in (0, 1). Detection threshold applied to
#'   permutation p-values.
#' @param parallel Logical scalar. If TRUE, evaluates scenarios in parallel
#'   when `future.apply` is available.
#' @return A data.frame with one row per simulation seed/scenario and columns:
#'   `modality`, `leakage`, `mode`, `seed`, observed metric, gap, p-value, and
#'   a logical `detected` flag. A scenario-level summary is attached as
#'   `attr(x, "summary")`.
#' @export
benchmark_leakage_suite <- function(modalities = c("omics", "imaging_tabular", "ehr_tabular"),
                                    leakages = c("none", "subject_overlap", "batch_confounded", "peek_norm", "lookahead"),
                                    modes = c("subject_grouped", "batch_blocked", "time_series"),
                                    learner = c("glmnet", "ranger"),
                                    seeds = 1:5,
                                    B = 200,
                                    alpha = 0.05,
                                    parallel = FALSE) {
  learner <- match.arg(learner)
  modalities <- unique(as.character(modalities))
  leakages <- unique(as.character(leakages))
  modes <- unique(as.character(modes))

  if (!length(modalities)) stop("modalities must be non-empty.", call. = FALSE)
  if (!length(leakages)) stop("leakages must be non-empty.", call. = FALSE)
  if (!length(modes)) stop("modes must be non-empty.", call. = FALSE)
  if (!is.numeric(alpha) || length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single numeric value in (0,1).", call. = FALSE)
  }

  modality_profiles <- list(
    omics = list(n = 220L, p = 800L, prevalence = 0.5, rho = 0.35),
    imaging_tabular = list(n = 500L, p = 120L, prevalence = 0.4, rho = 0.15),
    ehr_tabular = list(n = 1400L, p = 60L, prevalence = 0.3, rho = 0.5)
  )
  bad_mod <- setdiff(modalities, names(modality_profiles))
  if (length(bad_mod)) {
    stop(sprintf("Unsupported modalities: %s", paste(bad_mod, collapse = ", ")), call. = FALSE)
  }

  scenario_grid <- expand.grid(
    modality = modalities,
    leakage = leakages,
    mode = modes,
    stringsAsFactors = FALSE
  )

  eval_one <- function(i) {
    row <- scenario_grid[i, , drop = FALSE]
    prof <- modality_profiles[[row$modality]]
    horizon_use <- if (identical(row$mode, "time_series")) 1 else 0

    res <- suppressWarnings(simulate_leakage_suite(
      n = prof$n,
      p = prof$p,
      prevalence = prof$prevalence,
      mode = row$mode,
      learner = learner,
      leakage = row$leakage,
      rho = prof$rho,
      K = 5,
      repeats = 1,
      horizon = horizon_use,
      B = B,
      seeds = seeds,
      parallel = FALSE
    ))
    if (!nrow(res)) return(NULL)

    res$modality <- row$modality
    res$detected <- is.finite(res$p_value) & res$p_value <= alpha & is.finite(res$gap) & res$gap > 0
    res
  }

  rows <- if (isTRUE(parallel) && requireNamespace("future.apply", quietly = TRUE)) {
    future.apply::future_lapply(seq_len(nrow(scenario_grid)), eval_one, future.seed = TRUE)
  } else {
    lapply(seq_len(nrow(scenario_grid)), eval_one)
  }
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) {
    out <- data.frame()
    attr(out, "summary") <- data.frame()
    return(out)
  }

  out <- do.call(rbind, rows)
  out <- out[, c("modality", "leakage", "mode", "seed", "metric_obs", "gap", "p_value", "detected"), drop = FALSE]
  rownames(out) <- NULL

  out$detected_num <- as.numeric(out$detected)
  summary_df <- stats::aggregate(
    cbind(metric_obs, gap, detected_num) ~ modality + leakage + mode,
    data = out,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  out$detected_num <- NULL
  names(summary_df)[names(summary_df) == "detected_num"] <- "detection_rate"
  attr(out, "summary") <- summary_df
  out
}
