# Split-constraint derivation for dependency_graph objects.

.depgraph_constraint_mode_map <- c(
  subject = "Subject",
  batch = "Batch",
  study = "Study",
  time = "Timepoint"
)

.depgraph_constraint_edge_map <- c(
  subject = "sample_belongs_to_subject",
  batch = "sample_processed_in_batch",
  study = "sample_from_study",
  time = "sample_collected_at_timepoint"
)

.depgraph_normalize_constraint_mode <- function(mode) {
  mode <- tolower(as.character(mode)[1L])
  .depgraph_assert(mode %in% c("subject", "batch", "study", "time", "composite"), paste0("Unsupported constraint mode: ", mode))
  mode
}

.depgraph_normalize_constraint_modes <- function(modes) {
  modes <- unique(tolower(as.character(modes)))
  .depgraph_assert(length(modes) > 0L, "At least one constraint mode is required.")
  .depgraph_assert(all(modes %in% names(.depgraph_constraint_mode_map)), paste0(
    "Unsupported constraint mode(s): ",
    paste(setdiff(modes, names(.depgraph_constraint_mode_map)), collapse = ", ")
  ))
  modes
}

.depgraph_constraint_samples <- function(graph, samples = NULL) {
  sample_ids <- .depgraph_resolve_sample_node_ids(graph, samples)
  sample_nodes <- graph$nodes$data[graph$nodes$data$node_type == "Sample", , drop = FALSE]
  sample_nodes[match(sample_ids, sample_nodes$node_id), , drop = FALSE]
}

.depgraph_direct_assignment <- function(graph, mode, samples = NULL) {
  mode <- .depgraph_normalize_constraint_mode(mode)
  .depgraph_assert(mode != "composite", "Composite assignments must be derived separately.")

  sample_nodes <- .depgraph_constraint_samples(graph, samples)
  node_data <- graph$nodes$data
  edge_data <- graph$edges$data
  edge_type <- .depgraph_constraint_edge_map[[mode]]
  target_type <- .depgraph_constraint_mode_map[[mode]]

  edges <- edge_data[
    edge_data$edge_type == edge_type &
      edge_data$from %in% sample_nodes$node_id,
    c("from", "to", "edge_type"),
    drop = FALSE
  ]
  edges <- edges[!duplicated(edges$from), , drop = FALSE]

  target_rows <- node_data[node_data$node_id %in% edges$to, c("node_id", "node_key", "label", "attrs"), drop = FALSE]
  target_map <- stats::setNames(split(target_rows, target_rows$node_id), names(split(target_rows, target_rows$node_id)))

  rows <- lapply(seq_len(nrow(sample_nodes)), function(i) {
    sample_row <- sample_nodes[i, , drop = FALSE]
    edge_row <- edges[edges$from == sample_row$node_id, , drop = FALSE]

    if (nrow(edge_row) == 0L) {
      return(data.frame(
        sample_id = sample_row$node_key,
        sample_node_id = sample_row$node_id,
        linked_node_id = NA_character_,
        linked_node_type = target_type,
        linked_key = NA_character_,
        linked_label = NA_character_,
        edge_type = edge_type,
        stringsAsFactors = FALSE
      ))
    }

    target_row <- target_map[[edge_row$to[[1L]]]]
    data.frame(
      sample_id = sample_row$node_key,
      sample_node_id = sample_row$node_id,
      linked_node_id = target_row$node_id[[1L]],
      linked_node_type = target_type,
      linked_key = target_row$node_key[[1L]],
      linked_label = target_row$label[[1L]],
      edge_type = edge_row$edge_type[[1L]],
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

.depgraph_warning_if_missing <- function(assignments, mode) {
  missing_samples <- assignments$sample_id[is.na(assignments$linked_node_id)]
  if (length(missing_samples) == 0L) {
    return(character())
  }

  paste0(
    "Samples missing ", mode, " assignments were retained as singleton groups: ",
    paste(missing_samples, collapse = ", ")
  )
}

.depgraph_build_sample_map <- function(assignments, mode, explanation_prefix = NULL) {
  mode <- .depgraph_normalize_constraint_mode(mode)

  rows <- lapply(seq_len(nrow(assignments)), function(i) {
    row <- assignments[i, , drop = FALSE]
    linked <- !is.na(row$linked_node_id[[1L]]) && nzchar(row$linked_node_id[[1L]])
    group_id <- if (linked) {
      paste0(mode, ":", row$linked_key[[1L]])
    } else {
      paste0(mode, ":unlinked:", row$sample_id[[1L]])
    }
    group_label <- if (linked) row$linked_key[[1L]] else paste0("unlinked_", row$sample_id[[1L]])
    explanation <- if (linked) {
      paste0(
        explanation_prefix %||% paste0("Grouped by ", mode),
        " through ", row$edge_type[[1L]],
        " -> ", row$linked_key[[1L]], "."
      )
    } else {
      paste0(
        "No ", mode, " assignment was available; sample retained as an unlinked singleton group."
      )
    }

    data.frame(
      sample_id = row$sample_id,
      sample_node_id = row$sample_node_id,
      group_id = group_id,
      constraint_type = mode,
      group_label = group_label,
      explanation = explanation,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

.depgraph_compute_time_order <- function(graph, time_node_ids) {
  time_nodes <- graph$nodes$data[match(time_node_ids, graph$nodes$data$node_id), , drop = FALSE]
  raw_time_index <- vapply(
    time_nodes$attrs,
    function(x) {
      value <- .depgraph_extract_attr(x, "time_index")
      if (length(value) == 0L || is.na(value)) {
        return(NA_real_)
      }
      suppressWarnings(as.numeric(as.character(value)))
    },
    FUN.VALUE = numeric(1)
  )

  warnings <- character()
  order_rank <- rep(NA_integer_, length(time_node_ids))
  names(order_rank) <- time_node_ids
  names(raw_time_index) <- time_node_ids

  if (length(time_node_ids) == 0L) {
    return(list(
      source = "none",
      time_index = raw_time_index,
      order_rank = order_rank,
      warnings = warnings
    ))
  }

  if (all(!is.na(raw_time_index))) {
    unique_order <- sort(unique(raw_time_index))
    rank_map <- stats::setNames(seq_along(unique_order), as.character(unique_order))
    order_rank <- as.integer(rank_map[as.character(raw_time_index)])
    names(order_rank) <- time_node_ids
    return(list(
      source = "time_index",
      time_index = raw_time_index,
      order_rank = order_rank,
      warnings = warnings
    ))
  }

  precedence_edges <- graph$edges$data[
    graph$edges$data$edge_type == "timepoint_precedes" &
      graph$edges$data$from %in% time_node_ids &
      graph$edges$data$to %in% time_node_ids,
    c("from", "to"),
    drop = FALSE
  ]

  if (nrow(precedence_edges) > 0L) {
    precedence_graph <- igraph::graph_from_data_frame(
      d = precedence_edges,
      vertices = data.frame(name = time_node_ids, stringsAsFactors = FALSE),
      directed = TRUE
    )
    if (igraph::is_dag(precedence_graph)) {
      topo <- igraph::as_ids(igraph::topo_sort(precedence_graph, mode = "out"))
      order_rank[topo] <- seq_along(topo)
      return(list(
        source = "timepoint_precedes",
        time_index = raw_time_index,
        order_rank = order_rank,
        warnings = warnings
      ))
    }

    warnings <- c(warnings, "Timepoint precedence edges are cyclic; ordering metadata could not be derived from graph structure.")
  }

  if (any(!is.na(raw_time_index))) {
    available <- sort(unique(raw_time_index[!is.na(raw_time_index)]))
    rank_map <- stats::setNames(seq_along(available), as.character(available))
    order_rank[!is.na(raw_time_index)] <- as.integer(rank_map[as.character(raw_time_index[!is.na(raw_time_index)])])
    warnings <- c(warnings, "Time ordering is only partially available because some timepoints are missing `time_index` values.")
  } else {
    warnings <- c(warnings, "Time ordering is unavailable because neither `time_index` metadata nor usable `timepoint_precedes` edges were found.")
  }

  list(
    source = "partial",
    time_index = raw_time_index,
    order_rank = order_rank,
    warnings = warnings
  )
}

.derive_subject_constraints <- function(graph, samples = NULL) {
  assignments <- .depgraph_direct_assignment(graph, "subject", samples = samples)
  warnings <- .depgraph_warning_if_missing(assignments, "subject")
  sample_map <- .depgraph_build_sample_map(assignments, "subject")

  split_constraint(
    strategy = "subject",
    sample_map = sample_map,
    recommended_bioleak_args = list(
      group_var = "group_id",
      block_var = NULL,
      time_var = NULL,
      ordering_required = FALSE
    ),
    metadata = list(
      mode = "subject",
      strategy = "subject",
      relations_used = "sample_belongs_to_subject",
      n_groups = length(unique(sample_map$group_id)),
      n_samples = nrow(sample_map),
      warnings = warnings
    )
  )
}

.derive_batch_constraints <- function(graph, samples = NULL) {
  assignments <- .depgraph_direct_assignment(graph, "batch", samples = samples)
  warnings <- .depgraph_warning_if_missing(assignments, "batch")
  sample_map <- .depgraph_build_sample_map(assignments, "batch")

  split_constraint(
    strategy = "batch",
    sample_map = sample_map,
    recommended_bioleak_args = list(
      group_var = "group_id",
      block_var = "group_id",
      time_var = NULL,
      ordering_required = FALSE
    ),
    metadata = list(
      mode = "batch",
      strategy = "batch",
      relations_used = "sample_processed_in_batch",
      n_groups = length(unique(sample_map$group_id)),
      n_samples = nrow(sample_map),
      warnings = warnings
    )
  )
}

.derive_study_constraints <- function(graph, samples = NULL) {
  assignments <- .depgraph_direct_assignment(graph, "study", samples = samples)
  warnings <- .depgraph_warning_if_missing(assignments, "study")
  sample_map <- .depgraph_build_sample_map(assignments, "study")

  split_constraint(
    strategy = "study",
    sample_map = sample_map,
    recommended_bioleak_args = list(
      group_var = "group_id",
      block_var = "group_id",
      time_var = NULL,
      ordering_required = FALSE
    ),
    metadata = list(
      mode = "study",
      strategy = "study",
      relations_used = "sample_from_study",
      n_groups = length(unique(sample_map$group_id)),
      n_samples = nrow(sample_map),
      warnings = warnings
    )
  )
}

.derive_time_constraints <- function(graph, samples = NULL) {
  assignments <- .depgraph_direct_assignment(graph, "time", samples = samples)
  warnings <- .depgraph_warning_if_missing(assignments, "timepoint")
  valid_time_nodes <- stats::na.omit(assignments$linked_node_id)
  order_info <- .depgraph_compute_time_order(graph, unique(valid_time_nodes))
  warnings <- unique(c(warnings, order_info$warnings))

  order_rank <- order_info$order_rank[assignments$linked_node_id]
  order_rank[is.na(assignments$linked_node_id)] <- NA_integer_
  time_index <- order_info$time_index[assignments$linked_node_id]
  time_index[is.na(assignments$linked_node_id)] <- NA_real_
  timepoint_id <- assignments$linked_key

  sample_map <- .depgraph_build_sample_map(assignments, "time")
  sample_map$time_index <- unname(time_index)
  sample_map$timepoint_id <- timepoint_id
  sample_map$order_rank <- unname(as.integer(order_rank))

  split_constraint(
    strategy = "time",
    sample_map = sample_map,
    recommended_bioleak_args = list(
      group_var = "group_id",
      block_var = NULL,
      time_var = if (any(!is.na(sample_map$order_rank))) "order_rank" else NULL,
      ordering_required = any(!is.na(sample_map$order_rank))
    ),
    metadata = list(
      mode = "time",
      strategy = "time",
      relations_used = c("sample_collected_at_timepoint", if (identical(order_info$source, "timepoint_precedes")) "timepoint_precedes" else NULL),
      n_groups = length(unique(sample_map$group_id)),
      n_samples = nrow(sample_map),
      time_order_source = order_info$source,
      warnings = warnings
    )
  )
}

.depgraph_normalize_via_modes <- function(via) {
  if (is.null(via)) {
    return(c("subject", "batch", "study", "time"))
  }

  via <- as.character(via)
  via_modes <- vapply(via, function(x) {
    x_lower <- tolower(x)
    if (x_lower %in% names(.depgraph_constraint_mode_map)) {
      return(x_lower)
    }

    match_idx <- match(x_lower, tolower(.depgraph_constraint_mode_map))
    .depgraph_assert(!is.na(match_idx), paste0("Unsupported composite dependency source: ", x))
    names(.depgraph_constraint_mode_map)[[match_idx]]
  }, character(1), USE.NAMES = FALSE)

  unique(via_modes)
}

.derive_composite_strict_constraints <- function(graph, samples = NULL, via = NULL) {
  via_modes <- .depgraph_normalize_via_modes(via)
  via_types <- unname(.depgraph_constraint_mode_map[via_modes])
  components <- detect_dependency_components(graph, via = via_types, min_size = 1)
  table <- components$table

  if (!is.null(samples)) {
    sample_nodes <- .depgraph_constraint_samples(graph, samples)
    table <- table[table$sample_node_id %in% sample_nodes$node_id, , drop = FALSE]
  }

  sample_map <- data.frame(
    sample_id = table$sample_id,
    sample_node_id = table$sample_node_id,
    group_id = table$component_id,
    constraint_type = "composite_strict",
    group_label = table$component_id,
    explanation = paste0(
      "Strict composite grouping via transitive closure over: ",
      paste(via_types, collapse = ", "),
      "."
    ),
    stringsAsFactors = FALSE
  )

  warnings <- character()
  if (nrow(sample_map) > 0L) {
    singleton_fraction <- mean(table$component_size == 1L)
    if (singleton_fraction > 0.5) {
      warnings <- c(warnings, "Most strict composite groups are singletons; dependency coverage may be sparse.")
    }
  }

  split_constraint(
    strategy = "strict",
    sample_map = sample_map,
    recommended_bioleak_args = list(
      group_var = "group_id",
      block_var = NULL,
      time_var = NULL,
      ordering_required = FALSE
    ),
    metadata = list(
      mode = "composite",
      strategy = "strict",
      relations_used = unname(.depgraph_constraint_edge_map[via_modes]),
      via = via_types,
      n_groups = length(unique(sample_map$group_id)),
      n_samples = nrow(sample_map),
      warnings = warnings,
      projection_edges = components$metadata$projection_edges
    )
  )
}

.derive_rule_based_constraints <- function(graph, samples = NULL, priority = NULL, via = NULL) {
  via_modes <- .depgraph_normalize_via_modes(via)
  priority <- if (is.null(priority)) via_modes else .depgraph_normalize_constraint_modes(priority)
  .depgraph_assert(all(priority %in% via_modes), "`priority` must be a subset of the selected composite dependency sources.")

  sample_nodes <- .depgraph_constraint_samples(graph, samples)
  assignments <- lapply(priority, function(mode) .depgraph_direct_assignment(graph, mode, samples = sample_nodes$node_id))
  names(assignments) <- priority
  time_info <- if ("time" %in% priority) .derive_time_constraints(graph, samples = sample_nodes$node_id)$sample_map else NULL

  rows <- lapply(seq_len(nrow(sample_nodes)), function(i) {
    sample_row <- sample_nodes[i, , drop = FALSE]
    sample_id <- sample_row$node_key[[1L]]
    sample_node_id <- sample_row$node_id[[1L]]

    available <- lapply(assignments, function(x) x[x$sample_node_id == sample_node_id, , drop = FALSE])
    available_modes <- names(available)[vapply(available, function(x) nrow(x) == 1L && !is.na(x$linked_node_id[[1L]]), logical(1))]

    if (length(available_modes) == 0L) {
      return(data.frame(
        sample_id = sample_id,
        sample_node_id = sample_node_id,
        group_id = paste0("composite:unlinked:", sample_id),
        constraint_type = "unlinked",
        group_label = paste0("unlinked_", sample_id),
        explanation = "No prioritized dependency assignment was available; sample retained as a singleton group.",
        time_index = NA_real_,
        timepoint_id = NA_character_,
        order_rank = NA_integer_,
        stringsAsFactors = FALSE
      ))
    }

    chosen_mode <- available_modes[[1L]]
    chosen_row <- available[[chosen_mode]]
    additional <- setdiff(available_modes, chosen_mode)
    explanation <- paste0(
      "Composite rule-based grouping selected ", chosen_mode,
      " based on priority order ", paste(priority, collapse = " > "),
      " -> ", chosen_row$linked_key[[1L]], "."
    )
    if (length(additional) > 0L) {
      detail <- vapply(additional, function(mode) {
        paste0(mode, "=", available[[mode]]$linked_key[[1L]])
      }, character(1), USE.NAMES = FALSE)
      explanation <- paste0(explanation, " Additional available dependencies: ", paste(detail, collapse = ", "), ".")
    }

    time_row <- if (!is.null(time_info)) time_info[time_info$sample_node_id == sample_node_id, , drop = FALSE] else NULL

    data.frame(
      sample_id = sample_id,
      sample_node_id = sample_node_id,
      group_id = paste0("composite_", chosen_mode, ":", chosen_row$linked_key[[1L]]),
      constraint_type = chosen_mode,
      group_label = chosen_row$linked_key[[1L]],
      explanation = explanation,
      time_index = if (!is.null(time_row) && nrow(time_row) == 1L) time_row$time_index[[1L]] else NA_real_,
      timepoint_id = if (!is.null(time_row) && nrow(time_row) == 1L) time_row$timepoint_id[[1L]] else NA_character_,
      order_rank = if (!is.null(time_row) && nrow(time_row) == 1L) time_row$order_rank[[1L]] else NA_integer_,
      stringsAsFactors = FALSE
    )
  })

  sample_map <- do.call(rbind, rows)
  row.names(sample_map) <- NULL
  warnings <- character()
  if (any(sample_map$constraint_type == "unlinked")) {
    warnings <- c(warnings, paste0(
      "Some samples did not match any prioritized dependency and were retained as singleton groups: ",
      paste(sample_map$sample_id[sample_map$constraint_type == "unlinked"], collapse = ", ")
    ))
  }

  split_constraint(
    strategy = "rule_based",
    sample_map = sample_map,
    recommended_bioleak_args = list(
      group_var = "group_id",
      block_var = NULL,
      time_var = if (any(!is.na(sample_map$order_rank))) "order_rank" else NULL,
      ordering_required = any(!is.na(sample_map$order_rank))
    ),
    metadata = list(
      mode = "composite",
      strategy = "rule_based",
      relations_used = unname(.depgraph_constraint_edge_map[priority]),
      priority = priority,
      n_groups = length(unique(sample_map$group_id)),
      n_samples = nrow(sample_map),
      warnings = warnings
    )
  )
}

derive_split_constraints <- function(graph, mode = c("subject", "batch", "study", "time", "composite"), samples = NULL, strategy = c("strict", "rule_based"), via = NULL, priority = NULL, include_warnings = TRUE) {
  .depgraph_assert(inherits(graph, "dependency_graph"), "`graph` must be a `dependency_graph`.")
  mode <- .depgraph_normalize_constraint_mode(match.arg(mode))
  strategy <- match.arg(strategy)

  result <- switch(
    mode,
    subject = .derive_subject_constraints(graph, samples = samples),
    batch = .derive_batch_constraints(graph, samples = samples),
    study = .derive_study_constraints(graph, samples = samples),
    time = .derive_time_constraints(graph, samples = samples),
    composite = if (identical(strategy, "strict")) {
      .derive_composite_strict_constraints(graph, samples = samples, via = via)
    } else {
      .derive_rule_based_constraints(graph, samples = samples, priority = priority, via = via)
    }
  )

  if (!isTRUE(include_warnings)) {
    result$metadata$warnings <- character()
  }

  result
}

grouping_vector <- function(x) {
  .depgraph_assert(inherits(x, "split_constraint"), "`x` must be a `split_constraint`.")
  groups <- as.character(x$sample_map$group_id)
  names(groups) <- x$sample_map$sample_id
  groups
}
