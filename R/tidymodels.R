# Tidymodels interoperability helpers ---------------------------------------

.bio_is_rset <- function(x) inherits(x, "rset")
.bio_is_rsplit <- function(x) inherits(x, "rsplit")
.bio_is_rsample <- function(x) .bio_is_rset(x) || .bio_is_rsplit(x)
.bio_is_recipe <- function(x) inherits(x, "recipe")
.bio_is_workflow <- function(x) inherits(x, "workflow")
.bio_is_yardstick_metric <- function(x) inherits(x, "metric") || inherits(x, "metric_function")

.bio_validation_mode <- function(mode = NULL) {
  if (.bio_is_strict()) return("error")
  mode <- mode %||% getOption("bioLeak.validation_mode", "warn")
  if (is.null(mode) || length(mode) != 1L || !is.character(mode)) {
    mode <- "warn"
  }
  mode <- tolower(mode)
  if (!mode %in% c("warn", "error", "off")) {
    warning("Unknown bioLeak validation mode; using 'warn'.", call. = FALSE)
    mode <- "warn"
  }
  mode
}

.bio_signal_policy <- function(message, mode = "warn") {
  mode <- .bio_validation_mode(mode)
  if (identical(mode, "off")) return(invisible(FALSE))
  if (identical(mode, "error")) {
    stop(message, call. = FALSE)
  }
  warning(message, call. = FALSE)
  invisible(TRUE)
}

.bio_recipe_step_label <- function(step) {
  cls <- class(step)
  cls <- cls[!cls %in% c("step", "list")]
  if (!length(cls)) return("step")
  cls[[1]]
}

# Give split/metadata columns that a recipe would otherwise treat as predictors
# (e.g. via `y ~ .`) a non-predictor role before prep(), so neither they nor
# anything derived from them (dummies, PCA scores) reaches the learner. The
# columns stay available to steps that name them. Trained recipes are left as
# they are (they already trigger a leakage warning).
.bio_recipe_exclude_cols <- function(recipe_obj, cols) {
  if (!.bio_is_recipe(recipe_obj) || !length(cols) || .bio_recipe_is_trained(recipe_obj)) {
    return(recipe_obj)
  }
  vi <- recipe_obj$var_info
  if (is.null(vi) || !all(c("variable", "role") %in% names(vi))) return(recipe_obj)
  hit <- unique(vi$variable[vi$variable %in% cols & vi$role %in% "predictor"])
  if (!length(hit)) return(recipe_obj)
  do.call(recipes::update_role,
          c(list(recipe_obj), lapply(hit, as.name), list(new_role = "bioLeak_split")))
}

# Direct check that a prepped recipe's predictors do not depend on the split
# columns: bake the training rows as they are and with each split column
# shifted by one row; predictors that change were derived from a split column.
# This catches steps that select by type (e.g. step_dummy(all_nominal(),
# -all_outcomes())), which ignore the non-predictor role set above.
.bio_recipe_meta_derived <- function(rec_prep, train_df, cols) {
  cols <- intersect(cols, names(train_df))
  n <- nrow(train_df)
  if (!length(cols) || n < 2L) return(character())
  bake_pred <- function(d) tryCatch(
    as.data.frame(recipes::bake(rec_prep, new_data = d, recipes::all_predictors())),
    error = function(e) NULL)
  base <- bake_pred(train_df)
  if (is.null(base)) return(character())
  pert <- train_df
  rot <- c(2:n, 1L)
  for (cc in cols) pert[[cc]] <- pert[[cc]][rot]
  alt <- bake_pred(pert)
  if (is.null(alt)) return(character())
  common <- intersect(names(base), names(alt))
  common[!vapply(common, function(nm) isTRUE(all.equal(base[[nm]], alt[[nm]],
                                                         check.attributes = FALSE)),
                 logical(1))]
}

.bio_stop_meta_derived <- function(derived, cols, what = "The recipe") {
  stop(sprintf(paste0(
    "%s builds predictors from split columns (%s): %s. These would carry ",
    "group, batch, study, time or stratum information into the model. Select ",
    "predictors with all_predictors() / all_nominal_predictors() / ",
    "all_numeric_predictors(), or remove the columns with step_rm()."),
    what, paste(cols, collapse = ", "),
    paste(utils::head(derived, 5L), collapse = ", "),
    if (length(derived) > 5L) sprintf(" (and %d more)", length(derived) - 5L) else ""),
    call. = FALSE)
}

# Probability column of one class in parsnip/workflow predictions. parsnip
# names it ".pred_<level>" with the level as is (".pred_1", ".pred_case 1");
# make.names() forms are accepted as well. Returns NULL when absent: callers
# must not guess by position (with the user's level order kept for recipes
# and workflows, the second column can be the negative class).
.bio_prob_col <- function(cols, level) {
  level <- as.character(level)
  cand <- unique(c(paste0(".pred_", level), paste0(".pred_", make.names(level)),
                   level, make.names(level)))
  hit <- cand[cand %in% cols]
  if (length(hit)) hit[[1]] else NULL
}

.bio_recipe_is_trained <- function(recipe_obj) {
  if (!.bio_is_recipe(recipe_obj)) return(FALSE)
  has_tr_info <- !is.null(recipe_obj$tr_info)
  has_orig_levels <- !is.null(recipe_obj$orig_lvls)
  has_fit_times <- !is.null(recipe_obj$fit_times)
  isTRUE(has_tr_info || has_orig_levels || has_fit_times)
}

.bio_validate_recipe_graph <- function(recipe_obj, context = "fit_resample", mode = NULL) {
  if (!.bio_is_recipe(recipe_obj)) return(invisible(TRUE))
  mode <- .bio_validation_mode(mode)
  if (identical(mode, "off")) return(invisible(TRUE))

  if (.bio_recipe_is_trained(recipe_obj)) {
    .bio_signal_policy(
      sprintf(
        "%s: received a trained recipe. Pass an untrained recipe to avoid fold-global leakage.",
        context
      ),
      mode
    )
  }

  steps <- recipe_obj$steps %||% list()
  if (length(steps)) {
    for (i in seq_along(steps)) {
      st <- steps[[i]]
      if (isTRUE(st$trained)) {
        .bio_signal_policy(
          sprintf(
            "%s: recipe step %d (%s) is already trained; pass an untrained recipe.",
            context, i, .bio_recipe_step_label(st)
          ),
          mode
        )
      }
      if (inherits(st, "step_lag")) {
        lag_vals <- suppressWarnings(as.numeric(st$lag %||% numeric(0)))
        if (length(lag_vals) && any(is.finite(lag_vals) & lag_vals < 0)) {
          .bio_signal_policy(
            sprintf(
              "%s: recipe step %d (%s) uses negative lag values; this can introduce lookahead leakage.",
              context, i, .bio_recipe_step_label(st)
            ),
            mode
          )
        }
      }
      if (inherits(st, c("step_mutate", "step_filter"))) {
        exprs <- st$inputs %||% list()
        expr_txt <- paste(vapply(exprs, function(e) paste(deparse(e), collapse = " "), character(1)), collapse = " ")
        if (grepl("\\b(mean|sd|var|median|mad|quantile|scale|min|max)\\s*\\(", expr_txt)) {
          .bio_signal_policy(
            sprintf(
              "%s: recipe step %d (%s) calls a global summary function (%s) during bake(), which leaks test-set distribution statistics.",
              context, i, .bio_recipe_step_label(st), expr_txt
            ),
            mode
          )
        }
      }
    }
  }

  invisible(TRUE)
}

.bio_workflow_is_trained <- function(workflow_obj) {
  if (!.bio_is_workflow(workflow_obj) || !requireNamespace("workflows", quietly = TRUE)) {
    return(FALSE)
  }

  is_trained_fn <- get0("is_trained_workflow", envir = asNamespace("workflows"), inherits = FALSE)
  if (is.function(is_trained_fn)) {
    out <- tryCatch(is_trained_fn(workflow_obj), error = function(e) FALSE)
    return(isTRUE(out))
  }

  has_fit <- tryCatch({
    workflows::extract_fit_parsnip(workflow_obj)
    TRUE
  }, error = function(e) FALSE)
  isTRUE(has_fit)
}

.bio_validate_workflow_graph <- function(workflow_obj, context = "fit_resample", mode = NULL) {
  if (!.bio_is_workflow(workflow_obj)) return(invisible(TRUE))
  mode <- .bio_validation_mode(mode)
  if (identical(mode, "off")) return(invisible(TRUE))

  if (.bio_workflow_is_trained(workflow_obj)) {
    .bio_signal_policy(
      sprintf(
        "%s: received a trained workflow. Pass an unfitted workflow to avoid fold-global leakage.",
        context
      ),
      mode
    )
  }

  if (requireNamespace("workflows", quietly = TRUE)) {
    pre <- tryCatch(workflows::extract_preprocessor(workflow_obj), error = function(e) NULL)
    if (.bio_is_recipe(pre)) {
      .bio_validate_recipe_graph(
        pre,
        context = paste0(context, " (workflow preprocessor)"),
        mode = mode
      )
    }
  }

  invisible(TRUE)
}

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
  # Treat NULL, zero-length, or all-NA out_id as "assessment not provided"
  if (is.null(out_id) || length(out_id) == 0L || all(is.na(out_id))) {
    if (is.null(n)) {
      stop("rsplit object does not expose assessment indices; provide n.", call. = FALSE)
    }
    out_id <- setdiff(seq_len(n), in_id[!is.na(in_id)])
  }
  # Drop any NA values from indices
  in_id  <- in_id[!is.na(in_id)]
  out_id <- out_id[!is.na(out_id)]
  if (!is.null(n)) {
    if (length(in_id) && any(in_id < 1L | in_id > n)) {
      stop("rsplit analysis indices exceed data rows.", call. = FALSE)
    }
    if (length(out_id) && any(out_id < 1L | out_id > n)) {
      stop("rsplit assessment indices exceed data rows.", call. = FALSE)
    }
  }
  list(train = as.integer(in_id), test = as.integer(out_id))
}

# Which rsample id column holds the fold and which the repeat. rsample's
# repeated designs (vfold_cv(repeats = ), group_vfold_cv(repeats = )) put the
# repeat in `id` ("Repeat1") and the fold in `id2` ("Fold1"). Older bioLeak
# exports used the reverse; recognise it by the value prefixes.
.bio_rsample_fold_repeat_ids <- function(id, id2) {
  if (is.null(id2)) return(list(fold = id, repeat_id = NULL))
  is_rep <- function(x) all(grepl("^Repeat", x))
  is_fold <- function(x) all(grepl("^Fold", x))
  if (is_fold(id) && is_rep(id2)) return(list(fold = id, repeat_id = id2))
  list(fold = id2, repeat_id = id)
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
  ids_fr <- .bio_rsample_fold_repeat_ids(split_ids, split_ids2)
  fold_map <- as.integer(factor(ids_fr$fold, levels = unique(ids_fr$fold)))
  repeat_map <- if (is.null(ids_fr$repeat_id)) {
    # no repeat column: an id seen again starts a new repeat (older as_rsample()
    # exports of repeated plans wrote only "FoldK" ids)
    as.integer(stats::ave(seq_along(split_ids), split_ids, FUN = seq_along))
  } else {
    as.integer(factor(ids_fr$repeat_id, levels = unique(ids_fr$repeat_id)))
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
  # nested_cv() supplies both outer and inner folds, so fitting and tuning need
  # no inferred design; audits without one fall back to their documented
  # unrestricted null with a warning.
  if (is.null(perm_mode) && !(.bio_is_rset(splits) && "inner_resamples" %in% names(splits))) {
    stop(paste0("rsample splits require an explicit perm_mode (set attr(splits, ",
                "'bioLeak_perm_mode') or 'bioLeak_mode') or a single inferable ",
                "group/batch/study/time column via split_cols."),
         call. = FALSE)
  }

  # rsample::nested_cv() carries inner resamples built on each outer analysis
  # set; their positions index that analysis set, so map them to data rows.
  inner_list <- NULL
  if (.bio_is_rset(splits) && "inner_resamples" %in% names(splits)) {
    inner_list <- lapply(seq_along(split_list), function(i) {
      outer_tr <- .bio_rsplit_indices(split_list[[i]], n = n)$train
      inner_rs <- splits$inner_resamples[[i]]
      n_in <- length(outer_tr)
      ids <- .bio_rsample_fold_repeat_ids(
        as.character(inner_rs$id),
        if ("id2" %in% names(inner_rs)) as.character(inner_rs$id2) else NULL)
      f_map <- as.integer(factor(ids$fold, levels = unique(ids$fold)))
      r_map <- if (is.null(ids$repeat_id)) rep(1L, nrow(inner_rs)) else
        as.integer(factor(ids$repeat_id, levels = unique(ids$repeat_id)))
      lapply(seq_len(nrow(inner_rs)), function(k) {
        ix <- .bio_rsplit_indices(inner_rs$splits[[k]], n = n_in)
        list(train = outer_tr[ix$train], test = outer_tr[ix$test],
             fold = f_map[[k]], repeat_id = r_map[[k]])
      })
    })
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
    purge = 0,
    embargo = 0,
    summary = split_summary,
    hash = .bio_hash_indices(indices),
    inner = inner_list,
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
    split_purge <- splits@info$purge %||% 0
    split_embargo <- splits@info$embargo %||% 0
    group_col <- splits@info$group
    group_vec <- NULL
    if (!is.null(group_col) && !is.null(data) && group_col %in% names(data)) {
      group_vec <- data[[group_col]]
    } else if (!is.null(group_col) && !is.null(splits@info$coldata) &&
               group_col %in% names(splits@info$coldata)) {
      group_vec <- splits@info$coldata[[group_col]]
    }
    train <- .bio_time_series_train_indices(
      time_vec = time_vec,
      test_idx = test,
      candidate_idx = order(time_vec),   # time order, as the explicit builder stores train
      horizon = split_horizon,
      purge = split_purge,
      embargo = split_embargo,
      group_vec = group_vec
    )
  } else if (identical(splits@mode, "combined")) {
    # Combined mode: exclude samples sharing constraint-axis levels with test
    remaining <- setdiff(seq_len(n), test)
    cst <- splits@info$constraints
    cd <- data %||% splits@info$coldata
    if (!is.null(cst) && length(cst) > 1L && !is.null(cd)) {
      for (ax in cst[-1]) {
        ax_vec <- cd[[ax$col]]
        test_levels <- unique(ax_vec[test])
        remaining <- remaining[!ax_vec[remaining] %in% test_levels]
      }
    }
    train <- remaining
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
#' @return An rsample \code{rset} object compatible with tidymodels workflows.
#'   The returned object is a tibble with class \code{rset} containing:
#'   \describe{
#'     \item{\code{splits}}{List-column of \code{rsplit} objects, each with
#'       \code{analysis} (training indices) and \code{assessment} (test indices).}
#'     \item{\code{id}}{Character column with fold identifiers (e.g., "Fold1");
#'       with several repeats, repeat identifiers (e.g., "Repeat1"), following
#'       rsample's convention for repeated designs.}
#'     \item{\code{id2}}{With several repeats, character column with fold
#'       identifiers (e.g., "Fold1"); otherwise absent.}
#'   }
#'   The object also carries attributes for \code{group}, \code{batch},
#'   \code{study}, \code{time} (when available from the original \code{LeakSplits}),
#'   and \code{bioLeak_mode} indicating the original splitting mode. This allows
#'   the splits to be used with \code{tune::tune_grid()}, \code{rsample::fit_resamples()},
#'   and other tidymodels functions.
#' @examples
#' if (requireNamespace("rsample", quietly = TRUE)) {
#'   df <- data.frame(
#'     subject = rep(1:10, each = 2),
#'     outcome = rbinom(20, 1, 0.5),
#'     x1 = rnorm(20),
#'     x2 = rnorm(20)
#'   )
#'   splits <- make_split_plan(df, outcome = "outcome",
#'                         mode = "subject_grouped", group = "subject", v = 5)
#'   rset <- as_rsample(splits, data = df)
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

  fold_ids <- vapply(resolved, function(fold) paste0("Fold", fold$fold), character(1))
  repeats <- vapply(resolved, function(fold) as.integer(fold$repeat_id %||% 1L), integer(1))
  multi_rep <- length(unique(repeats)) > 1L

  split_objs <- lapply(seq_along(resolved), function(i) {
    fold <- resolved[[i]]
    .bio_make_rsplit(fold$train, fold$test, data = data, id = fold_ids[[i]])
  })

  if (multi_rep) {
    # rsample's convention for repeated designs: id = "RepeatR", id2 = "FoldK"
    if (!requireNamespace("tibble", quietly = TRUE)) {
      stop("Package 'tibble' is required to export repeated splits to rsample.", call. = FALSE)
    }
    id_tbl <- tibble::tibble(id = paste0("Repeat", repeats), id2 = fold_ids)
    rset <- rsample::new_rset(split_objs, ids = id_tbl, subclass = c("manual_rset", "rset"))
  } else {
    manual_rset <- utils::getFromNamespace("manual_rset", "rsample")
    rset <- manual_rset(split_objs, ids = fold_ids)
  }
  info <- x@info
  if (!is.null(info$group)) attr(rset, "group") <- info$group
  if (!is.null(info$batch)) attr(rset, "batch") <- info$batch
  if (!is.null(info$study)) attr(rset, "study") <- info$study
  if (!is.null(info$time)) attr(rset, "time") <- info$time
  if (!is.null(info$mode)) attr(rset, "bioLeak_mode") <- info$mode
  rset
}
