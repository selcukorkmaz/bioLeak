# Restricted permutation engines ------------------------------------------------

#.majority_level helper
.majority_level <- function(vals) {
  vals <- vals[!is.na(vals)]
  if (!length(vals)) return(NA_character_)
  tab <- table(vals)
  names(tab)[which.max(tab)]
}

.permute_subject_grouped <- function(y, subj, strata = NULL) {
  subj <- factor(subj)
  blocks <- split(seq_along(y), subj)
  values <- lapply(blocks, function(ix) y[ix])
  ord <- seq_along(blocks)
  if (!is.null(strata)) {
    block_strata <- vapply(ord, function(i) .majority_level(strata[blocks[[i]]]), character(1))
    split_idx <- split(ord, block_strata)
    values_perm <- vector("list", length(blocks))
    for (nm in names(split_idx)) {
      idx <- split_idx[[nm]]
      perm_vals <- values[sample(idx)]
      values_perm[idx] <- perm_vals
    }
  } else {
    values_perm <- values[sample(ord)]
  }
  out <- y
  for (i in seq_along(blocks)) {
    out[blocks[[i]]] <- values_perm[[i]]
  }
  out
}

.permute_within_group <- function(y, group) {
  group <- factor(group)
  out <- y
  for (lvl in levels(group)) {
    ix <- which(group == lvl)
    if (length(ix) > 1L) {
      out[ix] <- sample(out[ix])
    }
  }
  out
}

.permute_within_batch <- function(y, batch, strata = NULL) {
  .permute_within_group(y, batch)
}

.permute_within_study <- function(y, study, strata = NULL) {
  .permute_within_group(y, study)
}

# Factory returning closure that produces permuted training labels per fold
.permute_labels_factory <- function(cd, outcome, mode, folds, perm_stratify,
                                    time_block, block_len, seed,
                                    group_col = NULL, batch_col = NULL,
                                    study_col = NULL) {
  if (is.null(cd) || !outcome %in% names(cd)) {
    stop("Metadata with outcome column required for restricted permutations.")
  }
  y_all <- cd[[outcome]]
  MIN_SAMPLES_FOR_REGRESSION_STRATIFICATION <- 20L
  strata_vec <- NULL
  should_stratify <- FALSE
  if (isTRUE(perm_stratify)) {
    should_stratify <- TRUE
  } else if (identical(perm_stratify, "auto")) {
    should_stratify <- is.factor(y_all) ||
      (is.numeric(y_all) && length(y_all) >= MIN_SAMPLES_FOR_REGRESSION_STRATIFICATION)
  }
  if (should_stratify) {
    if (is.factor(y_all)) {
      strata_vec <- y_all
    } else if (is.numeric(y_all)) {
      # for regression, bins by quantiles to maintain structure
      br <- stats::quantile(y_all, probs = seq(0, 1, length.out = 5), na.rm = TRUE)
      br <- unique(br)
      strata_vec <- cut(y_all, breaks = br, include.lowest = TRUE, labels = FALSE)
    }
  }
  if (identical(mode, "time_series")) {
    if (!exists(".stationary_bootstrap", mode = "function")) {
      stop("Missing .stationary_bootstrap() implementation.")
    }
    if (!exists(".circular_block_permute", mode = "function")) {
      stop("Missing .circular_block_permute() implementation.")
    }
  }
  set.seed(seed)
  function(b) {
    set.seed(seed + b)
    res <- vector("list", length(folds))
    for (i in seq_along(folds)) {
      te_idx <- folds[[i]]$test
      permuted <- switch(mode,
        subject_grouped = {
          subj_col <- if (!is.null(group_col) && group_col %in% names(cd)) cd[[group_col]] else NULL
          if (is.null(subj_col) && "group" %in% names(cd)) subj_col <- cd[["group"]]
          if (is.null(subj_col)) subj_col <- seq_along(y_all)
          subj <- subj_col[te_idx]
          strata <- if (!is.null(strata_vec)) strata_vec[te_idx] else NULL
          .permute_subject_grouped(y_all[te_idx], subj, strata)
        },
        batch_blocked = {
          batch_vals <- NULL
          if (!is.null(batch_col) && batch_col %in% names(cd)) batch_vals <- cd[[batch_col]]
          if (is.null(batch_vals) && "batch" %in% names(cd)) batch_vals <- cd[["batch"]]
          if (is.null(batch_vals)) stop("Batch column not found for batch_blocked mode.")
          .permute_within_batch(y_all[te_idx], batch_vals[te_idx])
        },
        study_loocv = {
          study_vals <- NULL
          if (!is.null(study_col) && study_col %in% names(cd)) study_vals <- cd[[study_col]]
          if (is.null(study_vals) && "study" %in% names(cd)) study_vals <- cd[["study"]]
          if (is.null(study_vals)) stop("Study column not found for study_loocv mode.")
          .permute_within_study(y_all[te_idx], study_vals[te_idx])
        },
        time_series = {
          idx_order <- order(te_idx)
          te_idx_sorted <- te_idx[idx_order]
          L <- block_len
          if (is.null(L) || !is.finite(L) || L <= 0) {
            L <- max(5L, floor(length(te_idx_sorted) * 0.1))
          }
          perm_idx <- if (identical(time_block, "stationary")) {
            .stationary_bootstrap(te_idx_sorted, mean_block = L)
          } else {
            .circular_block_permute(te_idx_sorted, block_len = L)
          }
          y_all[perm_idx]
        },
        {
          sample(y_all[te_idx])
        }
      )
      res[[i]] <- permuted
    }
    res
  }
}
