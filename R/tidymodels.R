# Tidymodels interoperability helpers ---------------------------------------

.bio_is_rset <- function(x) inherits(x, "rset")
.bio_is_rsplit <- function(x) inherits(x, "rsplit")
.bio_is_rsample <- function(x) .bio_is_rset(x) || .bio_is_rsplit(x)
.bio_is_recipe <- function(x) inherits(x, "recipe")
.bio_is_workflow <- function(x) inherits(x, "workflow")
.bio_is_yardstick_metric <- function(x) inherits(x, "metric") || inherits(x, "metric_function")

.bio_perm_mode <- function(splits) {
  if (!inherits(splits, "LeakSplits")) return(NA_character_)
  mode <- splits@mode %||% NA_character_
  valid_modes <- c("subject_grouped", "batch_blocked", "study_loocv", "time_series")
  if (!identical(mode, "rsample")) {
    return(mode)
  }
  perm_mode <- splits@info$perm_mode %||% NULL
  if (is.null(perm_mode)) return(mode)
  perm_mode <- as.character(perm_mode)
  perm_mode <- perm_mode[!is.na(perm_mode) & nzchar(perm_mode)]
  if (!length(perm_mode)) return(mode)
  perm_mode <- perm_mode[[1]]
  if (perm_mode %in% valid_modes) perm_mode else mode
}

.bio_rsample_col_name <- function(val) {
  if (is.null(val) || is.logical(val)) return(NULL)
  if (is.character(val)) {
    val <- val[!is.na(val) & nzchar(val)]
    return(if (length(val)) val else NULL)
  }
  if (inherits(val, "formula")) {
    vars <- all.vars(val)
    return(if (length(vars)) vars else NULL)
  }
  if (is.name(val)) return(as.character(val))
  if (is.call(val)) {
    vars <- all.vars(val)
    return(if (length(vars)) vars else NULL)
  }
  NULL
}

.bio_rsample_split_cols <- function(splits) {
  grab <- function(nm) .bio_rsample_col_name(attr(splits, nm, exact = TRUE))
  list(
    group = grab("group") %||% grab("groups"),
    batch = grab("batch"),
    study = grab("study"),
    time = grab("time") %||% grab("index") %||% grab("date")
  )
}

.bio_rsplit_indices <- function(split, n = NULL) {
  if (!inherits(split, "rsplit")) {
    stop("Expected an rsample rsplit object.", call. = FALSE)
  }
  in_id <- split$in_id %||% split$analysis
  out_id <- split$out_id %||% split$assessment
  if (is.null(in_id)) {
    stop("rsplit object does not expose analysis indices.", call. = FALSE)
  }
  if (is.null(out_id)) {
    if (is.null(n)) {
      stop("rsplit object does not expose assessment indices; provide n.", call. = FALSE)
    }
    out_id <- setdiff(seq_len(n), in_id)
  }
  if (!is.null(n)) {
    if (any(in_id < 1L | in_id > n)) {
      stop("rsplit analysis indices exceed data rows.", call. = FALSE)
    }
    if (any(out_id < 1L | out_id > n)) {
      stop("rsplit assessment indices exceed data rows.", call. = FALSE)
    }
  }
  list(train = as.integer(in_id), test = as.integer(out_id))
}

.bio_as_leaksplits_from_rsample <- function(splits, n = NULL, coldata = NULL, split_cols = NULL) {
  if (!requireNamespace("rsample", quietly = TRUE)) {
    stop("Package 'rsample' is required for rsample split objects.", call. = FALSE)
  }

  if (.bio_is_rsplit(splits)) {
    split_list <- list(splits)
    split_ids <- "Fold1"
    split_ids2 <- NULL
  } else if (.bio_is_rset(splits)) {
    split_list <- splits$splits
    split_ids <- if ("id" %in% names(splits)) as.character(splits$id) else NULL
    split_ids2 <- if ("id2" %in% names(splits)) as.character(splits$id2) else NULL
  } else {
    stop("splits must be a LeakSplits or rsample rset/rsplit.", call. = FALSE)
  }

  if (!length(split_list)) {
    stop("rsample object contains no splits.", call. = FALSE)
  }

  if (is.null(split_ids)) {
    split_ids <- paste0("Fold", seq_along(split_list))
  }
  fold_map <- as.integer(factor(split_ids, levels = unique(split_ids)))

  if (is.null(split_ids2)) {
    repeat_map <- rep(1L, length(split_list))
  } else {
    repeat_map <- as.integer(factor(split_ids2, levels = unique(split_ids2)))
  }

  rsample_cols_attr <- .bio_rsample_split_cols(splits)
  rsample_cols <- rsample_cols_attr
  split_cols_norm <- NULL
  auto_cols <- FALSE
  if (!is.null(split_cols)) {
    if (is.character(split_cols) && length(split_cols) == 1L && identical(split_cols, "auto")) {
      auto_cols <- TRUE
      split_cols <- NULL
    } else if (is.character(split_cols)) {
      if (is.null(names(split_cols)) || !length(names(split_cols))) {
        stop("split_cols must be a named character vector or list.", call. = FALSE)
      }
      split_cols <- as.list(split_cols)
    } else if (!is.list(split_cols)) {
      stop("split_cols must be a named character vector or list.", call. = FALSE)
    }
    if (!is.null(split_cols)) {
      split_cols_norm <- list(
        group = .bio_rsample_col_name(split_cols$group %||% split_cols$groups),
        batch = .bio_rsample_col_name(split_cols$batch),
        study = .bio_rsample_col_name(split_cols$study),
        time = .bio_rsample_col_name(split_cols$time)
      )
      for (nm in names(split_cols_norm)) {
        if (!is.null(split_cols_norm[[nm]])) {
          rsample_cols[[nm]] <- split_cols_norm[[nm]]
        }
      }
    }
  }
  if (isTRUE(auto_cols) && !is.null(coldata)) {
    cd_names <- names(coldata)
    if (is.null(rsample_cols$group)) {
      if ("group" %in% cd_names) {
        rsample_cols$group <- "group"
      } else if ("subject" %in% cd_names) {
        rsample_cols$group <- "subject"
      }
    }
    if (is.null(rsample_cols$batch)) {
      if ("batch" %in% cd_names) rsample_cols$batch <- "batch"
    }
    if (is.null(rsample_cols$study)) {
      if ("study" %in% cd_names) rsample_cols$study <- "study"
    }
    if (is.null(rsample_cols$time)) {
      if ("time" %in% cd_names) {
        rsample_cols$time <- "time"
      } else if ("date" %in% cd_names) {
        rsample_cols$time <- "date"
      } else if ("index" %in% cd_names) {
        rsample_cols$time <- "index"
      }
    }
  }

  valid_modes <- c("subject_grouped", "batch_blocked", "study_loocv", "time_series")
  mode_attr <- attr(splits, "bioLeak_mode", exact = TRUE)
  if (!is.null(mode_attr)) {
    mode_attr <- as.character(mode_attr)
    mode_attr <- mode_attr[!is.na(mode_attr) & nzchar(mode_attr)]
    mode_attr <- if (length(mode_attr)) mode_attr[[1]] else NULL
  }
  perm_mode_attr <- attr(splits, "bioLeak_perm_mode", exact = TRUE)
  if (!is.null(perm_mode_attr)) {
    perm_mode_attr <- as.character(perm_mode_attr)
    perm_mode_attr <- perm_mode_attr[!is.na(perm_mode_attr) & nzchar(perm_mode_attr)]
    perm_mode_attr <- if (length(perm_mode_attr)) perm_mode_attr[[1]] else NULL
  }
  mode_attr_valid <- !is.null(mode_attr) && mode_attr %in% valid_modes
  perm_mode_attr_valid <- !is.null(perm_mode_attr) && perm_mode_attr %in% valid_modes
  if (isTRUE(mode_attr_valid) && isTRUE(perm_mode_attr_valid) &&
      !identical(mode_attr, perm_mode_attr)) {
    stop("bioLeak_mode and bioLeak_perm_mode conflict; supply only one or ensure they match.",
         call. = FALSE)
  }
  split_mode <- "rsample"
  perm_mode <- NULL
  if (isTRUE(mode_attr_valid)) {
    split_mode <- mode_attr
    perm_mode <- mode_attr
  } else if (isTRUE(perm_mode_attr_valid)) {
    perm_mode <- perm_mode_attr
  } else {
    time_ok <- function(col, cd) {
      if (is.null(cd)) return(FALSE)
      if (!is.null(col)) {
        col <- col[!is.na(col) & nzchar(col)]
        col <- if (length(col)) col[[1]] else NULL
      }
      if (is.null(col)) {
        if (!"time" %in% names(cd)) return(FALSE)
        col <- "time"
      }
      if (!col %in% names(cd)) return(FALSE)
      vals <- cd[[col]]
      is.numeric(vals) || inherits(vals, c("POSIXct", "Date"))
    }
    time_col <- rsample_cols$time
    time_series_classes <- c("rolling_origin", "sliding_index", "sliding_period", "sliding_window")
    if (inherits(splits, time_series_classes) && time_ok(time_col, coldata)) {
      perm_mode <- "time_series"
    } else if (inherits(splits, "group_rset")) {
      perm_mode <- "subject_grouped"
    } else {
      explicit_cols <- rsample_cols_attr
      if (!is.null(split_cols_norm)) {
        for (nm in names(split_cols_norm)) {
          explicit_cols[[nm]] <- split_cols_norm[[nm]]
        }
      }
      has_explicit <- any(vapply(explicit_cols, function(x) {
        !is.null(x) && length(x)
      }, logical(1)))
      perm_cols <- if (isTRUE(has_explicit)) explicit_cols else rsample_cols
      candidates <- character(0)
      if (!is.null(perm_cols$time) && time_ok(perm_cols$time, coldata)) {
        candidates <- c(candidates, "time_series")
      }
      if (!is.null(perm_cols$study)) {
        candidates <- c(candidates, "study_loocv")
      }
      if (!is.null(perm_cols$batch)) {
        candidates <- c(candidates, "batch_blocked")
      }
      if (!is.null(perm_cols$group)) {
        candidates <- c(candidates, "subject_grouped")
      }
      if (length(candidates) == 1L) perm_mode <- candidates[[1]]
    }
  }
  if (is.null(perm_mode)) {
    stop(paste0("rsample splits require an explicit perm_mode (set attr(splits, ",
                "'bioLeak_perm_mode') or 'bioLeak_mode') or a single inferable ",
                "group/batch/study/time column via split_cols."),
         call. = FALSE)
  }

  indices <- vector("list", length(split_list))
  summary_rows <- vector("list", length(split_list))
  for (i in seq_along(split_list)) {
    idx <- .bio_rsplit_indices(split_list[[i]], n = n)
    indices[[i]] <- list(
      train = idx$train,
      test = idx$test,
      fold = fold_map[[i]],
      repeat_id = repeat_map[[i]]
    )
    summary_rows[[i]] <- data.frame(
      fold = fold_map[[i]],
      repeat_id = repeat_map[[i]],
      train_n = length(idx$train),
      test_n = length(idx$test),
      stringsAsFactors = FALSE
    )
  }

  split_summary <- do.call(rbind, summary_rows)
  info <- list(
    outcome = NULL,
    v = length(unique(fold_map)),
    repeats = length(unique(repeat_map)),
    seed = NA_integer_,
    mode = split_mode,
    perm_mode = perm_mode,
    group = rsample_cols$group,
    batch = rsample_cols$batch,
    study = rsample_cols$study,
    time = rsample_cols$time,
    stratify = NA,
    nested = FALSE,
    horizon = 0,
    summary = split_summary,
    hash = .bio_hash_indices(indices),
    inner = NULL,
    compact = FALSE,
    fold_assignments = NULL,
    coldata = coldata,
    source = "rsample",
    id = split_ids,
    id2 = split_ids2,
    split_cols_override = split_cols_norm
  )

  methods::new("LeakSplits", mode = split_mode, indices = indices, info = info)
}

.bio_resolve_fold_indices <- function(splits, fold, n, data = NULL) {
  if (!inherits(splits, "LeakSplits")) {
    stop("splits must be a LeakSplits object.", call. = FALSE)
  }
  if (!isTRUE(splits@info$compact) || !is.null(fold$train)) return(fold)

  fold_assignments <- splits@info$fold_assignments
  if (is.null(fold_assignments) || !length(fold_assignments)) {
    stop("Compact splits require fold assignments to compute indices.", call. = FALSE)
  }
  r <- fold$repeat_id %||% 1L
  assign_vec <- fold_assignments[[r]]
  if (is.null(assign_vec)) {
    stop(sprintf("Missing fold assignments for repeat %s.", r), call. = FALSE)
  }
  test <- which(assign_vec == fold$fold)
  if (identical(splits@mode, "time_series")) {
    time_col <- splits@info$time
    time_vec <- NULL
    if (!is.null(time_col) && !is.null(data) && time_col %in% names(data)) {
      time_vec <- data[[time_col]]
    } else if (!is.null(time_col) && !is.null(splits@info$coldata) &&
               time_col %in% names(splits@info$coldata)) {
      time_vec <- splits@info$coldata[[time_col]]
    }
    if (is.null(time_vec)) {
      stop("time_series compact splits require time column in data.", call. = FALSE)
    }
    split_horizon <- splits@info$horizon %||% 0
    if (!length(test)) {
      train <- integer(0)
    } else {
      tmin <- min(time_vec[test])
      if (split_horizon == 0) {
        train <- which(time_vec < tmin)
      } else {
        train <- which(time_vec <= (tmin - split_horizon))
      }
    }
  } else {
    train <- setdiff(seq_len(n), test)
  }

  list(train = train, test = test, fold = fold$fold, repeat_id = fold$repeat_id)
}

.bio_make_rsplit <- function(train, test, data, id = NULL) {
  rsample_make_splits <- utils::getFromNamespace("make_splits", "rsample")
  args <- list(list(analysis = train, assessment = test), data = data)
  if (!is.null(id)) args$id <- id
  out <- try(do.call(rsample_make_splits, args), silent = TRUE)
  if (inherits(out, "try-error") && !is.null(id)) {
    args$id <- NULL
    out <- do.call(rsample_make_splits, args)
  }
  out
}

#' Convert LeakSplits to an rsample resample set
#'
#' @param x LeakSplits object created by [make_split_plan()].
#' @param data Optional data.frame used to populate rsample splits. When NULL,
#'   the stored `coldata` from `x` is used (if available).
#' @param ... Additional arguments passed to methods (unused).
#' @return An rsample `rset` object.
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   subject = rep(1:10, each = 2),
#'   outcome = rbinom(20, 1, 0.5),
#'   x1 = rnorm(20),
#'   x2 = rnorm(20)
#' )
#' splits <- make_split_plan(df, outcome = "outcome",
#'                       mode = "subject_grouped", group = "subject", v = 5)
#' rset <- as_rsample(splits, data = df)
#' }
#' @export
as_rsample <- function(x, data = NULL, ...) UseMethod("as_rsample")

#' @export
as_rsample.LeakSplits <- function(x, data = NULL, ...) {
  if (!requireNamespace("rsample", quietly = TRUE)) {
    stop("Package 'rsample' is required to export rsample splits.", call. = FALSE)
  }
  if (!inherits(x, "LeakSplits")) {
    stop("x must be a LeakSplits object.", call. = FALSE)
  }
  if (is.null(data)) data <- x@info$coldata
  if (is.null(data)) {
    stop("Provide 'data' to export rsample splits.", call. = FALSE)
  }
  data <- as.data.frame(data, check.names = FALSE)
  n <- nrow(data)
  folds <- x@indices
  if (!length(folds)) stop("LeakSplits object contains no folds.", call. = FALSE)

  resolved <- lapply(folds, function(fold) {
    .bio_resolve_fold_indices(x, fold, n = n, data = data)
  })

  ids <- vapply(resolved, function(fold) {
    paste0("Fold", fold$fold)
  }, character(1))
  repeats <- vapply(resolved, function(fold) fold$repeat_id %||% 1L, integer(1))
  ids2 <- if (length(unique(repeats)) > 1L) {
    paste0("Repeat", repeats)
  } else {
    NULL
  }

  split_objs <- lapply(seq_along(resolved), function(i) {
    fold <- resolved[[i]]
    .bio_make_rsplit(fold$train, fold$test, data = data, id = ids[[i]])
  })

  manual_rset <- utils::getFromNamespace("manual_rset", "rsample")
  args <- list(splits = split_objs, ids = ids)
  if (!is.null(ids2) && "ids2" %in% names(formals(manual_rset))) {
    args$ids2 <- ids2
  }
  rset <- do.call(manual_rset, args)
  info <- x@info
  if (!is.null(info$group)) attr(rset, "group") <- info$group
  if (!is.null(info$batch)) attr(rset, "batch") <- info$batch
  if (!is.null(info$study)) attr(rset, "study") <- info$study
  if (!is.null(info$time)) attr(rset, "time") <- info$time
  if (!is.null(info$mode)) attr(rset, "bioLeak_mode") <- info$mode
  rset
}
