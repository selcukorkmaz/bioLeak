# Tidymodels interoperability helpers ---------------------------------------

.bio_is_rset <- function(x) inherits(x, "rset")
.bio_is_rsplit <- function(x) inherits(x, "rsplit")
.bio_is_rsample <- function(x) .bio_is_rset(x) || .bio_is_rsplit(x)
.bio_is_recipe <- function(x) inherits(x, "recipe")
.bio_is_workflow <- function(x) inherits(x, "workflow")
.bio_is_yardstick_metric <- function(x) inherits(x, "metric") || inherits(x, "metric_function")

.bio_rsample_col_name <- function(val) {
  if (is.null(val) || is.logical(val)) return(NULL)
  if (is.character(val)) {
    val <- val[!is.na(val) & nzchar(val)]
    return(if (length(val)) val else NULL)
  }
  if (inherits(val, "quosure")) {
    if (!requireNamespace("rlang", quietly = TRUE)) return(NULL)
    val <- rlang::quo_get_expr(val)
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

  rsample_cols <- .bio_rsample_split_cols(splits)
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
    mode = "rsample",
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

  methods::new("LeakSplits", mode = "rsample", indices = indices, info = info)
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
  make_splits <- getFromNamespace("make_splits", "rsample")
  args <- list(list(analysis = train, assessment = test))
  make_formals <- names(formals(make_splits))
  if ("data" %in% make_formals) args$data <- data
  if (!is.null(id) && "id" %in% make_formals) args$id <- id
  do.call(make_splits, args)
}

#' Convert LeakSplits to an rsample resample set
#'
#' @param x LeakSplits object created by [make_splits()].
#' @param data Optional data.frame used to populate rsample splits. When NULL,
#'   the stored `coldata` from `x` is used (if available).
#' @return An rsample `rset` object.
#' @examples
#' \dontrun{
#' df <- data.frame(
#'   subject = rep(1:10, each = 2),
#'   outcome = rbinom(20, 1, 0.5),
#'   x1 = rnorm(20),
#'   x2 = rnorm(20)
#' )
#' splits <- make_splits(df, outcome = "outcome",
#'                       mode = "subject_grouped", group = "subject", v = 5)
#' rset <- as_rsample(splits, data = df)
#' }
#' @export
as_rsample <- function(x, ...) UseMethod("as_rsample")

#' @export
as_rsample.LeakSplits <- function(x, data = NULL) {
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

  manual_rset <- getFromNamespace("manual_rset", "rsample")
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
  rset
}
