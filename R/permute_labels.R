# Safe sample that avoids the R sample(x) gotcha when length(x) == 1 -----------
.safe_sample <- function(x) {
  if (length(x) <= 1L) return(x)
  x[sample.int(length(x))]
}

# Restricted permutation engines ------------------------------------------------

#' Quantile break cache for permutation stratification
#'
#' Internal environment used to cache quantile breakpoints for numeric
#' outcomes during restricted permutation testing. This avoids recomputing
#' quantiles across repeated calls in \code{audit_leakage()}.
#'
#' @format An environment used to cache quantile breakpoints.
#' @return An environment (internal data object, not a function).
#' @keywords internal
#' @docType data
# Cached quantile breaks for numeric stratification ------------------------
.quantile_break_cache <- new.env(parent = emptyenv())

.get_cached_quantile_breaks <- function(vals, probs) {
  vals_clean <- vals[!is.na(vals)]
  if (requireNamespace("digest", quietly = TRUE)) {
    key <- digest::digest(list(vals = vals_clean, probs = probs))
    if (exists(key, envir = .quantile_break_cache, inherits = FALSE)) {
      return(get(key, envir = .quantile_break_cache, inherits = FALSE))
    }
    breaks <- stats::quantile(vals_clean, probs = probs, na.rm = TRUE)
    breaks <- unique(breaks)
    assign(key, breaks, envir = .quantile_break_cache)
    return(breaks)
  } else {
    # Fallback: compute directly without caching
    return(unique(stats::quantile(vals_clean, probs = probs, na.rm = TRUE)))
  }
}

#.majority_level helper
.majority_level <- function(vals) {
  vals <- vals[!is.na(vals)]
  if (!length(vals)) return(NA_character_)
  tab <- table(vals)
  names(tab)[which.max(tab)]
}

# Donor index for a label permutation that keeps each unit's labels together.
# Returns d such that key[d] (or y[d]) is the permuted labelling.
# Units (e.g. subjects) of the same size swap their whole label vectors, and
# every unit's labels are then shuffled within the unit (labels within a unit
# are exchangeable under H0; without this, paired designs such as one tumour
# and one normal sample per patient could not move any label). Exchanging only
# equal-size units keeps the row-level label multiset (class balance) exact:
# swapping labels between units of different sizes changes it whenever unit
# size is related to the outcome (type-I 1.0 for accuracy when cases had five
# samples and controls one), and shuffling rows across units made the null
# too narrow for clustered outcomes. With `strata`, units are exchanged only
# within the stratum of their majority level. Without `groups` every row is
# its own unit (a plain shuffle).
.design_perm_donor <- function(key, groups = NULL, strata = NULL) {
  n <- length(key)
  if (n <= 1L) return(seq_len(n))
  if (is.null(groups)) return(sample.int(n))
  blocks <- split(seq_len(n), factor(groups))
  donor <- seq_len(n)
  bstrata <- if (!is.null(strata)) {
    vapply(blocks, function(ix) .majority_level(strata[ix]), character(1))
  } else {
    rep("all", length(blocks))
  }
  sizes <- lengths(blocks)
  for (st in unique(bstrata)) {
    in_st <- which(bstrata == st)
    for (sz in unique(sizes[in_st])) {
      el <- in_st[sizes[in_st] == sz]
      src <- if (length(el) > 1L) el[sample.int(length(el))] else el
      for (k in seq_along(el)) {
        from <- blocks[[src[k]]]
        donor[blocks[[el[k]]]] <- from[sample.int(length(from))]
      }
    }
  }
  donor
}

# Share of rows whose unit has at least one same-size partner (within its
# stratum), i.e. rows whose labels the design permutation can move between
# units. When it is small the null is close to the identity and conservative.
.design_perm_coverage <- function(groups, strata = NULL) {
  if (is.null(groups) || !length(groups)) return(1)
  blocks <- split(seq_along(groups), factor(groups))
  bstrata <- if (!is.null(strata)) {
    vapply(blocks, function(ix) .majority_level(strata[ix]), character(1))
  } else {
    rep("all", length(blocks))
  }
  cls <- paste(bstrata, lengths(blocks))
  movable <- cls %in% cls[duplicated(cls)]
  sum(lengths(blocks)[movable]) / length(groups)
}

.permute_subject_grouped <- function(y, subj, strata = NULL) {
  y[.design_perm_donor(y, subj, strata)]
}

.permute_within_group <- function(y, group, strata = NULL) {
  group <- factor(group)
  out <- y
  if (!is.null(strata)) {
    # Permute within group x strata cells
    cells <- interaction(group, strata, drop = TRUE)
    for (cell in levels(cells)) {
      ix <- which(cells == cell)
      if (length(ix) > 1L) {
        out[ix] <- sample(out[ix])
      }
    }
  } else {
    for (lvl in levels(group)) {
      ix <- which(group == lvl)
      if (length(ix) > 1L) {
        out[ix] <- sample(out[ix])
      }
    }
  }
  out
}

.permute_within_batch <- function(y, batch, strata = NULL) {
  .permute_within_group(y, batch, strata = strata)
}

.permute_within_study <- function(y, study, strata = NULL) {
  .permute_within_group(y, study, strata = strata)
}

#' Restricted permutation label factory
#'
#' Builds a closure that generates permuted outcome vectors per fold while
#' respecting grouping/batch/study/time constraints used in
#' \code{audit_leakage()}. Numeric outcomes can be stratified by quantiles to
#' preserve outcome structure under permutation.
#'
#' @param cd data.frame of sample metadata.
#' @param outcome outcome column name.
#' @param mode resampling mode (subject_grouped, batch_blocked, study_loocv, time_series).
#' @param folds list of fold descriptors from \code{LeakSplits}. When compact
#'   splits are used, fold assignments are read from the
#'   \code{fold_assignments} attribute.
#' @param perm_stratify logical or "auto"; if TRUE, permute within strata.
#' @param time_block time-series block permutation method.
#' @param block_len block length for time-series permutations.
#' @param seed integer seed.
#' @param group_col,batch_col,study_col optional metadata columns.
#' @param time_col optional metadata column name for time-series ordering.
#' @param perm_refit logical; if TRUE model is retrained on permuted labels
#'   (block permutation preserves subject structure); if FALSE predictions are
#'   fixed and simple label shuffle is used for \code{subject_grouped} mode.
#' @param verbose logical; print progress messages.
#' @return A function that returns a list of permuted outcome vectors, one per fold.
#' @keywords internal
# Factory returning closure that produces permuted training labels per fold
.permute_labels_factory <- function(cd, outcome, mode, folds, perm_stratify,
                                    time_block, block_len, seed,
                                    group_col = NULL, batch_col = NULL,
                                    study_col = NULL, time_col = NULL,
                                    perm_refit = TRUE,
                                    verbose = FALSE) {
  if (is.null(cd) || !outcome %in% names(cd)) {
    stop("Metadata with outcome column required for restricted permutations.")
  }
  y_all <- cd[[outcome]]
  if (all(is.na(y_all))) {
    stop("Outcome column contains only NA values.")
  }

  fold_assignments <- attr(folds, "fold_assignments")
  resolve_test_idx <- function(fold) {
    if (!is.null(fold$test)) return(fold$test)
    if (is.null(fold_assignments) || !length(fold_assignments)) {
      stop("Fold assignments required to resolve test indices for compact splits.")
    }
    r <- fold$repeat_id
    if (is.null(r) || !is.finite(r)) r <- 1L
    assign_vec <- fold_assignments[[r]]
    if (is.null(assign_vec)) {
      stop(sprintf("Missing fold assignments for repeat %s.", r))
    }
    which(assign_vec == fold$fold)
  }

  MIN_SAMPLES_FOR_REGRESSION_STRATIFICATION <- 20L
  strata_vec <- NULL
  should_stratify <- FALSE
  if (isTRUE(perm_stratify)) {
    should_stratify <- TRUE
    if (is.numeric(y_all) &&
        length(stats::na.omit(y_all)) < MIN_SAMPLES_FOR_REGRESSION_STRATIFICATION) {
      warning("perm_stratify = TRUE requires at least 20 non-missing numeric outcomes; proceeding without stratification.")
      should_stratify <- FALSE
    }
  } else if (identical(perm_stratify, "auto")) {
    should_stratify <- is.factor(y_all) ||
      (is.numeric(y_all) &&
       length(stats::na.omit(y_all)) >= MIN_SAMPLES_FOR_REGRESSION_STRATIFICATION)
  }
  if (should_stratify) {
    if (is.factor(y_all)) {
      strata_vec <- y_all
      if (isTRUE(verbose)) {
        message("[permute_labels] Stratifying by factor outcome levels: ",
                paste(levels(strata_vec), collapse = ", "))
      }
    } else if (is.numeric(y_all)) {
      # for regression, bins by quantiles to maintain structure
      br <- .get_cached_quantile_breaks(y_all, probs = seq(0, 1, length.out = 5))
      strata_vec <- cut(y_all, breaks = br, include.lowest = TRUE, labels = FALSE)
      if (isTRUE(verbose)) {
        message("[permute_labels] Stratifying numeric outcome into ",
                length(unique(stats::na.omit(strata_vec))), " bins.")
      }
    }
  } else if (isTRUE(verbose)) {
    message("[permute_labels] Stratification disabled for outcome '", outcome, "'.")
  }
  # Refit nulls for batch/study designs permute labels within each batch or
  # study; stratifying that by the outcome itself permutes within
  # batch x outcome cells, which leaves every label in place (p = 1).
  if (!is.null(strata_vec) && isTRUE(perm_refit) &&
      mode %in% c("batch_blocked", "study_loocv")) {
    .bio_warn(sprintf(paste0(
      "perm_stratify is ignored for %s refit permutations: labels are permuted within ",
      "each %s, and stratifying that by the outcome would leave every label in place."),
      mode, if (identical(mode, "batch_blocked")) "batch" else "study"),
      "bioLeak_input_warning")
    strata_vec <- NULL
  }
  if (identical(mode, "time_series")) {
    time_col_use <- time_col
    if (is.null(time_col_use) && "time" %in% names(cd)) {
      time_col_use <- "time"
    }
    if (is.null(time_col_use) || !time_col_use %in% names(cd)) {
      stop("time_series permutations require a time column in metadata.", call. = FALSE)
    }
    time_vec <- cd[[time_col_use]]
    if (!is.numeric(time_vec) && !inherits(time_vec, c("POSIXct", "Date"))) {
      stop("time_series time column must be numeric, Date, or POSIXct.", call. = FALSE)
    }
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
    if (isTRUE(verbose)) {
      message("[permute_labels] Generating permuted labels for replicate ", b, ".")
    }
    res <- vector("list", length(folds))
    for (i in seq_along(folds)) {
      if (isTRUE(verbose)) {
        message(sprintf("[permute_labels] Permuting fold %d/%d using mode '%s'.",
                        i, length(folds), mode))
      }
      te_idx <- resolve_test_idx(folds[[i]])
      if (isTRUE(verbose) && !is.null(strata_vec)) {
        strata_fold <- stats::na.omit(strata_vec[te_idx])
        if (length(unique(strata_fold)) < 2L) {
          message("[permute_labels] Warning: Fewer than two non-NA strata present in this test fold; stratification has limited effect.")
        }
      }
      permuted <- switch(mode,
        subject_grouped = {
          # For perm_refit=FALSE predictions are fixed; simple shuffle gives the
          # correct null and avoids degenerate/inflated-Type-I-error behaviour
          # from block permutation constrained by subject size and outcome strata.
          if (!isTRUE(perm_refit)) {
            .safe_sample(y_all[te_idx])
          } else {
            subj_col <- if (!is.null(group_col) && group_col %in% names(cd)) cd[[group_col]] else NULL
            if (is.null(subj_col) && "group" %in% names(cd)) subj_col <- cd[["group"]]
            if (is.null(subj_col)) subj_col <- seq_along(y_all)
            subj <- subj_col[te_idx]
            strata <- if (!is.null(strata_vec)) strata_vec[te_idx] else NULL
            if (isTRUE(verbose) && !is.null(strata)) {
              message("[permute_labels] Subject-grouped strata used: ",
                      paste(sort(unique(stats::na.omit(strata))), collapse = ", "))
            }
            .permute_subject_grouped(y_all[te_idx], subj, strata)
          }
        },
        batch_blocked = {
          # For perm_refit=FALSE predictions are fixed; simple shuffle gives the
          # correct null.  Restricted within-batch permutation with stratification
          # collapses to the identity when each test fold is a single batch
          # (because shuffling within batch × outcome cells changes nothing).
          if (!isTRUE(perm_refit)) {
            .safe_sample(y_all[te_idx])
          } else {
            batch_vals <- NULL
            if (!is.null(batch_col) && batch_col %in% names(cd)) batch_vals <- cd[[batch_col]]
            if (is.null(batch_vals) && "batch" %in% names(cd)) batch_vals <- cd[["batch"]]
            if (is.null(batch_vals)) stop("Batch column not found for batch_blocked mode.")
            if (isTRUE(verbose)) {
              ktab <- table(batch_vals[te_idx])
              if (any(ktab == 1L)) {
                message("[permute_labels] Note: ", sum(ktab == 1L),
                        " batch level(s) in this fold have only one sample; permutation within those is identity.")
              }
            }
            .permute_within_batch(y_all[te_idx], batch_vals[te_idx],
                                  strata = if (!is.null(strata_vec)) strata_vec[te_idx] else NULL)
          }
        },
        study_loocv = {
          # For perm_refit=FALSE predictions are fixed; simple shuffle gives the
          # correct null.  Restricted within-study permutation with stratification
          # collapses to the identity when each test fold is a single study
          # (because shuffling within study × outcome cells changes nothing).
          if (!isTRUE(perm_refit)) {
            .safe_sample(y_all[te_idx])
          } else {
            study_vals <- NULL
            if (!is.null(study_col) && study_col %in% names(cd)) study_vals <- cd[[study_col]]
            if (is.null(study_vals) && "study" %in% names(cd)) study_vals <- cd[["study"]]
            if (is.null(study_vals)) stop("Study column not found for study_loocv mode.")
            if (isTRUE(verbose)) {
              ktab <- table(study_vals[te_idx])
              if (any(ktab == 1L)) {
                message("[permute_labels] Note: ", sum(ktab == 1L),
                        " study level(s) in this fold have only one sample; permutation within those is identity.")
              }
            }
            .permute_within_study(y_all[te_idx], study_vals[te_idx],
                                  strata = if (!is.null(strata_vec)) strata_vec[te_idx] else NULL)
          }
        },
        time_series = {
          time_vals <- time_vec[te_idx]
          idx_order <- order(time_vals, te_idx, na.last = TRUE)
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
          stopifnot(length(perm_idx) == length(te_idx_sorted))
          if (any(!perm_idx %in% te_idx_sorted)) {
            stop(".stationary_bootstrap/.circular_block_permute must return a permutation of the provided indices.")
          }
          perm_time <- y_all[perm_idx]
          pos <- match(te_idx, te_idx_sorted)
          if (anyNA(pos)) {
            stop("Failed to align permuted time-series labels to fold order.", call. = FALSE)
          }
          perm_time[pos]
        },
        {
          .safe_sample(y_all[te_idx])
        }
      )
      res[[i]] <- permuted
    }
    res
  }
}
