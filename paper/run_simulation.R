#!/usr/bin/env Rscript
## =================================================================
## Main simulation study for bioLeak manuscript (Section 3.1)
##
## NOTE: This script uses a CUSTOM data-generation pipeline that is
## intentionally different from the package-level simulate_leakage_suite().
## Do NOT replace this script with simulate_leakage_suite() calls —
## the two differ in the following ways:
##
##   peek_norm leakage:
##     paper  -> as.numeric(y) + rnorm(n, 0, 0.3)  [noisy continuous y]
##     package -> (y - mean(y)) / sd(y)             [z-scored binary y]
##
##   lookahead leakage:
##     paper  -> shifted continuous biomarker (linpred + noise)
##     package -> shifted binary outcome c(y[-1], y[n])
##
##   signal generation:
##     paper  -> AR(1) noise added to linpred (arima.sim, ar=0.9, sd=0.3)
##     package -> AR correlation on predictors via rho parameter (default rho=0)
##
##   audit_leakage() settings:
##     paper  -> perm_refit=FALSE, perm_stratify=TRUE, target_scan=FALSE
##     package -> perm_refit="auto", perm_stratify=FALSE (default), target_scan=TRUE
##
##   peek_norm preprocessing:
##     paper  -> zscore normalization applied uniformly across all conditions
##     package -> normalize="none" for peek_norm to preserve leakage signal
##
## Custom loop using perm_refit=FALSE (fixed-prediction permutations)
## Parameters:
##   B = 200, Replicates = 50, perm_refit = FALSE
##
## Full factorial: 5 leakage x 4 n x 3 p x 3 s = 180 configs
## Total seeds: 180 x 50 = 9,000
##
## Key design choices:
##   - s=0 included as a true null (y ~ Bernoulli(0.5), no signal)
##   - Batch is independent of outcome except in batch_confounded
##   - AR noise has fixed SD (not signal-dependent)
##   - Lookahead leakage uses shifted biomarker (linpred + noise), not raw linpred
##
## Memory-efficient strategy:
##   - Workers capped at min(4, physical_cores - 1) to avoid RAM exhaustion
##   - Tasks sorted smallest-first so large configs don't all run together
##   - gc() after every task to free intermediate objects immediately
##   - Results assembled from checkpoints at the end (not held in RAM)
##   - Chunks sized 1 task each for fine-grained load balancing
## =================================================================

cat("=== bioLeak Main Simulation Study ===\n")
cat("Start time:", format(Sys.time()), "\n\n")

library(bioLeak)
library(future)
library(future.apply)
library(progressr)

## ── Worker count: conservative to avoid memory pressure ──────────
## Each worker holds a ~n×p matrix, glmnet models, and B=200 permutation
## results. Cap at 4 workers max so the OS + R main session keep headroom.
physical_cores <- parallel::detectCores(logical = FALSE)
n_workers <- max(1L, min(4L, physical_cores - 1L))
## Limit per-worker RAM via future's global-size guard (500 MB each)
plan(multisession, workers = n_workers)
cat(sprintf("Parallel backend: %d workers (of %d physical cores)\n\n",
            n_workers, physical_cores))

## Experimental design
leakage_types <- c("none", "subject_overlap", "batch_confounded",
                   "peek_norm", "lookahead")
sample_sizes  <- c(100, 250, 500, 1000)
dimensions    <- c(10, 50, 100)
signals       <- c(0, 0.5, 2.0)
B             <- 200
n_seeds       <- 50

out_dir      <- "paper/sim_results"
ckpt_dir     <- file.path(out_dir, "checkpoints")
if (!dir.exists(out_dir))   dir.create(out_dir,   recursive = TRUE)
if (!dir.exists(ckpt_dir))  dir.create(ckpt_dir,  recursive = TRUE)

## Build full task list: every (config, seed) pair
grid <- expand.grid(
  leakage = leakage_types,
  n = sample_sizes,
  p = dimensions,
  s = signals,
  stringsAsFactors = FALSE
)
tasks <- data.frame(
  config_idx = rep(seq_len(nrow(grid)), each = n_seeds),
  seed       = rep(seq_len(n_seeds), times = nrow(grid))
)

## ── Sort tasks smallest-first by memory footprint (n * p) ────────
## This prevents multiple large configs (n=1000, p=100) from running
## simultaneously on different workers and exhausting RAM.
tasks$cost <- grid$n[tasks$config_idx] * grid$p[tasks$config_idx]
tasks <- tasks[order(tasks$cost), ]
tasks$cost <- NULL
rownames(tasks) <- NULL

n_tasks <- nrow(tasks)
cat("Configs:", nrow(grid), "| Seeds per config:", n_seeds,
    "| Total tasks:", n_tasks, "\n")

## ── Identify already-completed tasks via checkpoints ─────────────
## Content-based key: robust to grid reordering and parameter changes
tasks$ckpt_key <- sprintf("%s_n%d_p%d_s%.1f_seed%d",
                          grid$leakage[tasks$config_idx],
                          grid$n[tasks$config_idx],
                          grid$p[tasks$config_idx],
                          grid$s[tasks$config_idx],
                          tasks$seed)
completed <- vapply(tasks$ckpt_key, function(k) {
  file.exists(file.path(ckpt_dir, sprintf("task_%s.rds", k)))
}, logical(1), USE.NAMES = FALSE)
n_done <- sum(completed)
remaining_idx <- which(!completed)
n_remaining <- length(remaining_idx)

cat(sprintf("Already completed: %d | Remaining: %d\n", n_done, n_remaining))

if (n_remaining == 0L) {
  cat("All tasks already checkpointed. Assembling results...\n")
} else {
  ## ── Chunk remaining tasks into batches ─────────────────────────────────
  ## Keep total chunk count modest (≤ ~100) to avoid exhausting macOS file
  ## descriptors — each future uses inter-process connections that count
  ## against the OS open-file limit (default 256 on macOS).
  tasks_per_chunk <- max(5L, ceiling(n_remaining / 100L))
  chunk_ids <- split(remaining_idx,
                     ceiling(seq_along(remaining_idx) / tasks_per_chunk))
  cat(sprintf("Chunks: %d of ~%d task(s) each\n\n", length(chunk_ids),
              tasks_per_chunk))

  ## Worker function: process a batch of tasks
  ## Key difference: gc() after every task to prevent memory buildup
  process_chunk <- function(task_indices, tasks_df, grid_df, B_val,
                            ckpt_dir) {
    library(bioLeak)
    suppressMessages(library(glmnet))

    results <- vector("list", length(task_indices))
    for (k in seq_along(task_indices)) {
      ti      <- task_indices[k]
      cfg_idx <- tasks_df$config_idx[ti]
      seed    <- tasks_df$seed[ti]
      cfg     <- grid_df[cfg_idx, ]
      ckpt_key <- tasks_df$ckpt_key[ti]

      ## Resume from checkpoint if completed by a parallel run
      ckpt_file <- file.path(ckpt_dir, sprintf("task_%s.rds", ckpt_key))
      if (file.exists(ckpt_file)) {
        results[[k]] <- readRDS(ckpt_file)
        next
      }

      res <- tryCatch({
        set.seed(seed)
        n <- cfg$n; p_dim <- cfg$p; s <- cfg$s; leakage_type <- cfg$leakage

        X <- matrix(rnorm(n * p_dim), n, p_dim)
        colnames(X) <- sprintf("x%02d", seq_len(p_dim))
        subject <- sample(seq_len(max(5L, n %/% 5L)), n, replace = TRUE)

        ## Generate outcome: s=0 → pure noise, s>0 → signal + AR noise
        if (s > 0) {
          linpred <- rowSums(X[, seq_len(min(5L, p_dim)), drop = FALSE])
          linpred <- as.numeric(scale(linpred)) * s
          ar_noise <- as.numeric(arima.sim(model = list(ar = 0.9), n = n,
                                           sd = 0.3))
          linpred <- linpred + ar_noise
        } else {
          linpred <- rep(0, n)
        }

        y_prob  <- pnorm(linpred)
        y       <- rbinom(n, 1L, y_prob)

        ## Batch: independent of outcome by default
        batch <- sample(c("a","b","c"), n, replace = TRUE)
        study    <- sample(seq_len(max(3L, n %/% 80L)), n, replace = TRUE)
        time_var <- seq_len(n)

        ## Inject leakage features
        if (leakage_type == "subject_overlap") {
          X <- cbind(X, leak_subj = ave(y, subject, FUN = mean))
        } else if (leakage_type == "batch_confounded") {
          ## Make batch outcome-dependent, then leak the group mean
          batch <- ifelse(y == 1L,
                          sample(c("a","b","c"), n, replace = TRUE,
                                 prob = c(0.6, 0.2, 0.2)),
                          sample(c("a","b","c"), n, replace = TRUE,
                                 prob = c(0.15, 0.5, 0.35)))
          X <- cbind(X, leak_batch = ave(y, batch, FUN = mean))
        } else if (leakage_type == "peek_norm") {
          X <- cbind(X, leak_global = as.numeric(y) + rnorm(n, 0, 0.3))
        } else if (leakage_type == "lookahead") {
          ## Continuous biomarker (noisy proxy of latent process);
          ## lookahead = next time point's measurement (future information).
          ## NOTE: At s=0, linpred=0 so biomarker is pure i.i.d. noise
          ## (N(0, 0.5)). The shifted version carries no outcome information,
          ## making this mechanism undetectable — the s=0 entries in Table 2
          ## are type I error, not detection rates.
          biomarker <- linpred + rnorm(n, 0, 0.5)
          X <- cbind(X, leak_future = c(biomarker[-1], biomarker[n]))
        }

        df <- data.frame(X, y = factor(y))
        df$subject <- subject; df$batch <- batch
        df$study <- study; df$time <- time_var

        splits <- make_split_plan(
          df, outcome = "y", mode = "subject_grouped",
          group = "subject", batch = "batch", study = "study", time = "time",
          v = 5, stratify = TRUE, seed = seed,
          progress = FALSE
        )

        preprocess <- list(
          impute = list(method = "median"), normalize = list(method = "zscore"),
          filter = list(var_thresh = 0, iqr_thresh = 0),
          fs = list(method = "none")
        )

        fit <- fit_resample(
          df, outcome = "y", splits = splits,
          preprocess = preprocess,
          learner = "glmnet", metrics = "auc", seed = seed
        )

        aud <- audit_leakage(
          fit, metric = "auc", B = B_val,
          perm_refit = FALSE, perm_stratify = TRUE,
          seed = seed, target_scan = FALSE, return_perm = FALSE
        )

        data.frame(
          seed = seed, n = n, p = p_dim, s = s, leakage = leakage_type,
          metric_obs = aud@permutation_gap$metric_obs,
          perm_mean  = aud@permutation_gap$perm_mean,
          gap        = aud@permutation_gap$gap,
          p_value    = aud@permutation_gap$p_value,
          stringsAsFactors = FALSE
        )
      }, error = function(e) {
        data.frame(
          seed = seed, n = cfg$n, p = cfg$p, s = cfg$s, leakage = cfg$leakage,
          metric_obs = NA, perm_mean = NA, gap = NA, p_value = NA,
          stringsAsFactors = FALSE
        )
      })
      results[[k]] <- res
      saveRDS(res, ckpt_file)

      ## Free all large intermediates and return memory to the OS
      rm(res); gc(verbose = FALSE, full = TRUE)
    }
    do.call(rbind, results)
  }

  ## ── Run all chunks in parallel with progress bar ──────────────────
  t0 <- proc.time()

  handlers("txtprogressbar")

  cat(sprintf("Running %d tasks across %d chunks on %d workers...\n\n",
              n_remaining, length(chunk_ids), n_workers))

  suppressWarnings(with_progress({
    p <- progressor(steps = n_remaining)

    ## Wrap process_chunk to signal progress after each chunk completes
    future_lapply(
      chunk_ids,
      function(indices) {
        out <- process_chunk(indices, tasks, grid, B, ckpt_dir)
        p(amount = length(indices))
        out
      },
      future.seed = TRUE,
      future.chunk.size = NULL  # each element of chunk_ids is already a chunk
    )
  }))

  total_elapsed <- (proc.time() - t0)[3]
  cat(sprintf("\nProcessing complete. Elapsed: %.1f minutes\n",
              total_elapsed / 60))
}

## ── Assemble results from checkpoints (memory-efficient) ─────────
## Read one checkpoint at a time instead of holding all chunk_results in RAM
cat("Assembling results from checkpoints...\n")
result_list <- vector("list", n_tasks)
for (i in seq_len(n_tasks)) {
  ckpt_key <- tasks$ckpt_key[i]
  ckpt_file <- file.path(ckpt_dir, sprintf("task_%s.rds", ckpt_key))
  if (file.exists(ckpt_file)) {
    result_list[[i]] <- tryCatch(readRDS(ckpt_file), error = function(e) {
      warning(sprintf("Skipping corrupted checkpoint: %s (%s)", ckpt_file, e$message))
      NULL
    })
  }
}
result_list <- Filter(Negate(is.null), result_list)
sim_all <- do.call(rbind, result_list)
rm(result_list); gc(verbose = FALSE)

## ── Save results ──────────────────────────────────────────────────
saveRDS(sim_all, "paper/sim_results_all.rds")

## Also save per-config files for incremental access
for (i in seq_len(nrow(grid))) {
  cfg <- grid[i, ]
  sub <- sim_all[sim_all$leakage == cfg$leakage &
                 sim_all$n == cfg$n &
                 sim_all$p == cfg$p &
                 sim_all$s == cfg$s, ]
  fname <- sprintf("sim_%s_n%d_p%d_s%.1f.rds", cfg$leakage, cfg$n, cfg$p, cfg$s)
  saveRDS(sub, file.path(out_dir, fname))
}

cat("Total rows:", nrow(sim_all), "\n")
cat("NA rows:", sum(is.na(sim_all$metric_obs)), "\n")

## Quick sanity check
cat("\n=== Sanity Check ===\n")

## Table 1: True null (s=0) — "none" should be ~0.05, leakage types should be high
cat("\n--- s = 0 (no signal): type I error / leakage-only detection ---\n")
sim_null <- sim_all[sim_all$s == 0, ]
agg0 <- aggregate(p_value ~ leakage + n, data = sim_null,
                  FUN = function(x) mean(x < 0.05, na.rm = TRUE))
names(agg0)[3] <- "power"
print(reshape(agg0, idvar = "leakage", timevar = "n", direction = "wide"))

## Table 2: With signal (s>0) — all conditions should have high power
cat("\n--- s > 0 (real signal present): detection power ---\n")
sim_sig <- sim_all[sim_all$s > 0, ]
agg1 <- aggregate(p_value ~ leakage + n, data = sim_sig,
                  FUN = function(x) mean(x < 0.05, na.rm = TRUE))
names(agg1)[3] <- "power"
print(reshape(agg1, idvar = "leakage", timevar = "n", direction = "wide"))

cat("\nResults saved to paper/sim_results_all.rds\n")

## Clean up parallel workers
plan(sequential)
cat("Done.\n")
