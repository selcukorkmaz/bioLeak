# Validation helpers for dependency_graph objects.

.depgraph_empty_issue_data <- function() {
  data.frame(
    issue_id = character(),
    level = character(),
    severity = character(),
    code = character(),
    message = character(),
    node_ids = I(list()),
    edge_ids = I(list()),
    details = I(list()),
    stringsAsFactors = FALSE
  )
}

.new_validation_issue <- function(level, severity, code, message, node_ids = character(), edge_ids = character(), details = list()) {
  .depgraph_assert(level %in% .depgraph_validation_levels, paste0("Unsupported validation level: ", level))
  .depgraph_assert(severity %in% .depgraph_validation_severities, paste0("Unsupported validation severity: ", severity))

  data.frame(
    issue_id = NA_character_,
    level = level,
    severity = severity,
    code = code,
    message = as.character(message)[1L],
    node_ids = I(list(as.character(node_ids))),
    edge_ids = I(list(as.character(edge_ids))),
    details = I(list(.depgraph_normalize_attr_entry(details, context = "`details`"))),
    stringsAsFactors = FALSE
  )
}

.depgraph_bind_issues <- function(issues) {
  if (length(issues) == 0L) {
    return(.depgraph_empty_issue_data())
  }

  out <- do.call(rbind, issues)
  row.names(out) <- NULL
  if (nrow(out) == 0L) {
    return(.depgraph_empty_issue_data())
  }
  out$issue_id <- paste0("issue_", seq_len(nrow(out)))
  out
}

.depgraph_issue_count_table <- function(x, column) {
  if (nrow(x) == 0L) {
    return(data.frame(value = character(), n = integer(), stringsAsFactors = FALSE))
  }
  .depgraph_count_table(x[[column]])
}

.summarize_validation_issues <- function(issues) {
  list(
    n_issues = nrow(issues),
    by_level = .depgraph_issue_count_table(issues, "level"),
    by_severity = .depgraph_issue_count_table(issues, "severity"),
    by_code = .depgraph_issue_count_table(issues, "code")
  )
}

.depgraph_old_checks_to_levels <- function(checks) {
  checks <- unique(as.character(checks))
  levels <- character()

  if (any(checks %in% c("ids", "references", "cardinality"))) {
    levels <- c(levels, "structural")
  }
  if (any(checks %in% c("schema", "time", "cardinality"))) {
    levels <- c(levels, "semantic")
  }

  unique(levels)
}

.validate_structural <- function(graph) {
  node_data <- graph$nodes$data
  edge_data <- graph$edges$data
  issues <- list()

  if (anyNA(node_data$node_id) || any(!nzchar(node_data$node_id))) {
    issues[[length(issues) + 1L]] <- .new_validation_issue(
      level = "structural",
      severity = "error",
      code = "missing_node_id",
      message = "Node table contains missing or empty `node_id` values."
    )
  }

  if (anyDuplicated(node_data$node_id) > 0L) {
    dup_ids <- unique(node_data$node_id[duplicated(node_data$node_id)])
    issues[[length(issues) + 1L]] <- .new_validation_issue(
      level = "structural",
      severity = "error",
      code = "duplicate_node_id",
      message = paste0("Node table contains duplicated `node_id` values: ", paste(dup_ids, collapse = ", ")),
      node_ids = dup_ids
    )
  }

  if (anyNA(edge_data$edge_id) || any(!nzchar(edge_data$edge_id))) {
    issues[[length(issues) + 1L]] <- .new_validation_issue(
      level = "structural",
      severity = "error",
      code = "missing_edge_id",
      message = "Edge table contains missing or empty `edge_id` values."
    )
  }

  if (anyDuplicated(edge_data$edge_id) > 0L) {
    dup_ids <- unique(edge_data$edge_id[duplicated(edge_data$edge_id)])
    issues[[length(issues) + 1L]] <- .new_validation_issue(
      level = "structural",
      severity = "error",
      code = "duplicate_edge_id",
      message = paste0("Edge table contains duplicated `edge_id` values: ", paste(dup_ids, collapse = ", ")),
      edge_ids = dup_ids
    )
  }

  if (nrow(edge_data) > 0L) {
    triplets <- paste(edge_data$from, edge_data$to, edge_data$edge_type, sep = "||")
    if (anyDuplicated(triplets) > 0L) {
      dup_rows <- duplicated(triplets) | duplicated(triplets, fromLast = TRUE)
      dup_edges <- edge_data$edge_id[dup_rows]
      issues[[length(issues) + 1L]] <- .new_validation_issue(
        level = "structural",
        severity = "error",
        code = "duplicate_edge_triplet",
        message = "Edge table contains duplicate `(from, to, edge_type)` relations.",
        edge_ids = dup_edges
      )
    }
  }

  valid_nodes <- node_data$node_id
  bad_from <- setdiff(unique(edge_data$from), valid_nodes)
  bad_to <- setdiff(unique(edge_data$to), valid_nodes)
  if (length(bad_from) > 0L) {
    issues[[length(issues) + 1L]] <- .new_validation_issue(
      level = "structural",
      severity = "error",
      code = "missing_source_node",
      message = paste0("Edges reference missing source nodes: ", paste(bad_from, collapse = ", ")),
      node_ids = bad_from
    )
  }
  if (length(bad_to) > 0L) {
    issues[[length(issues) + 1L]] <- .new_validation_issue(
      level = "structural",
      severity = "error",
      code = "missing_target_node",
      message = paste0("Edges reference missing target nodes: ", paste(bad_to, collapse = ", ")),
      node_ids = bad_to
    )
  }

  bad_node_types <- setdiff(unique(node_data$node_type), .depgraph_node_types)
  if (length(bad_node_types) > 0L) {
    issues[[length(issues) + 1L]] <- .new_validation_issue(
      level = "structural",
      severity = "error",
      code = "unsupported_node_type",
      message = paste0("Unsupported node types found: ", paste(bad_node_types, collapse = ", "))
    )
  }

  type_map <- stats::setNames(node_data$node_type, node_data$node_id)
  if (nrow(edge_data) > 0L) {
    for (i in seq_len(nrow(edge_data))) {
      edge_type <- edge_data$edge_type[[i]]
      schema_row <- .depgraph_edge_type_schema(edge_type)
      if (nrow(schema_row) == 0L) {
        issues[[length(issues) + 1L]] <- .new_validation_issue(
          level = "structural",
          severity = "error",
          code = "unsupported_edge_type",
          message = paste0("Unknown edge type encountered: ", edge_type),
          edge_ids = edge_data$edge_id[[i]]
        )
        next
      }

      observed_from <- unname(type_map[[edge_data$from[[i]]]])
      observed_to <- unname(type_map[[edge_data$to[[i]]]])
      if (!identical(observed_from, schema_row$from_type[[1L]]) || !identical(observed_to, schema_row$to_type[[1L]])) {
        issues[[length(issues) + 1L]] <- .new_validation_issue(
          level = "structural",
          severity = "error",
          code = "invalid_edge_signature",
          message = paste0(
            "Edge `", edge_data$edge_id[[i]], "` violates schema for `", edge_type,
            "` (expected ", schema_row$from_type[[1L]], " -> ", schema_row$to_type[[1L]],
            ", found ", observed_from, " -> ", observed_to, ")."
          ),
          node_ids = c(edge_data$from[[i]], edge_data$to[[i]]),
          edge_ids = edge_data$edge_id[[i]]
        )
      }
    }
  }

  for (edge_type in .depgraph_single_target_edges) {
    subset_idx <- edge_data$edge_type == edge_type
    if (!any(subset_idx)) {
      next
    }
    edge_subset <- edge_data[subset_idx, c("edge_id", "from", "to"), drop = FALSE]
    counts <- stats::aggregate(
      edge_subset$to,
      by = list(from = edge_subset$from),
      FUN = function(x) length(unique(x))
    )
    bad <- counts[counts$x > 1L, , drop = FALSE]
    if (nrow(bad) > 0L) {
      offending_edges <- edge_subset$edge_id[edge_subset$from %in% bad$from]
      issues[[length(issues) + 1L]] <- .new_validation_issue(
        level = "structural",
        severity = "error",
        code = "single_target_violation",
        message = paste0(
          "Edge type `", edge_type, "` links a source node to multiple targets: ",
          paste(bad$from, collapse = ", ")
        ),
        node_ids = bad$from,
        edge_ids = offending_edges,
        details = list(edge_type = edge_type)
      )
    }
  }

  .depgraph_bind_issues(issues)
}

.validate_semantic <- function(graph) {
  node_data <- graph$nodes$data
  edge_data <- graph$edges$data
  issues <- list()

  for (i in seq_len(nrow(node_data))) {
    attr_check <- .depgraph_validate_node_attrs(
      node_type = node_data$node_type[[i]],
      attrs = node_data$attrs[[i]]
    )

    if (length(attr_check$missing_required) > 0L) {
      issues[[length(issues) + 1L]] <- .new_validation_issue(
        level = "semantic",
        severity = "error",
        code = "missing_required_node_attr",
        message = paste0(
          "Node `", node_data$node_id[[i]], "` is missing required attributes: ",
          paste(attr_check$missing_required, collapse = ", ")
        ),
        node_ids = node_data$node_id[[i]]
      )
    }

    if (length(attr_check$unknown_attrs) > 0L) {
      issues[[length(issues) + 1L]] <- .new_validation_issue(
        level = "semantic",
        severity = "warning",
        code = "non_canonical_node_attr",
        message = paste0(
          "Node `", node_data$node_id[[i]], "` contains non-canonical attributes: ",
          paste(attr_check$unknown_attrs, collapse = ", ")
        ),
        node_ids = node_data$node_id[[i]]
      )
    }
  }

  sample_ids <- node_data$node_id[node_data$node_type == "Sample"]
  for (relation in names(.depgraph_semantic_relation_rules)) {
    rule <- .depgraph_semantic_relation_rules[[relation]]
    relation_edges <- edge_data[edge_data$edge_type == relation, c("edge_id", "from", "to"), drop = FALSE]

    counts <- stats::setNames(integer(length(sample_ids)), sample_ids)
    targets <- vector("list", length(sample_ids))
    names(targets) <- sample_ids

    if (nrow(relation_edges) > 0L) {
      observed_counts <- stats::aggregate(
        relation_edges$to,
        by = list(from = relation_edges$from),
        FUN = function(x) length(unique(x))
      )
      counts[observed_counts$from] <- observed_counts$x
      for (sample_id in unique(relation_edges$from)) {
        targets[[sample_id]] <- unique(relation_edges$to[relation_edges$from == sample_id])
      }
    }

    if (!is.null(rule$missing_code) && !is.null(rule$missing_severity) && rule$min_targets > 0L) {
      missing_samples <- names(counts)[counts < rule$min_targets]
      if (length(missing_samples) > 0L) {
        issues[[length(issues) + 1L]] <- .new_validation_issue(
          level = "semantic",
          severity = rule$missing_severity,
          code = rule$missing_code,
          message = paste0(
            "Samples missing required relation `", relation, "`: ",
            paste(missing_samples, collapse = ", ")
          ),
          node_ids = missing_samples,
          details = list(edge_type = relation)
        )
      }
    }

    allow_multiple <- identical(relation, "sample_belongs_to_subject") &&
      .depgraph_validation_override(graph, "allow_multi_subject_samples", default = FALSE)

    if (!allow_multiple && !is.null(rule$multiple_code) && !is.null(rule$multiple_severity) && is.finite(rule$max_targets)) {
      multi_samples <- names(counts)[counts > rule$max_targets]
      if (length(multi_samples) > 0L) {
        edge_ids <- relation_edges$edge_id[relation_edges$from %in% multi_samples]
        related_nodes <- unique(c(
          multi_samples,
          unlist(targets[multi_samples], use.names = FALSE)
        ))
        issues[[length(issues) + 1L]] <- .new_validation_issue(
          level = "semantic",
          severity = rule$multiple_severity,
          code = rule$multiple_code,
          message = paste0(
            "Samples violate relation cardinality for `", relation, "`: ",
            paste(multi_samples, collapse = ", ")
          ),
          node_ids = related_nodes,
          edge_ids = edge_ids,
          details = list(edge_type = relation)
        )
      }
    }
  }

  feature_nodes <- node_data[node_data$node_type == "FeatureSet", , drop = FALSE]
  if (nrow(feature_nodes) > 0L) {
    for (i in seq_len(nrow(feature_nodes))) {
      derivation_scope <- .depgraph_extract_attr(feature_nodes$attrs[[i]], "derivation_scope")
      if (!is.na(derivation_scope) && !derivation_scope %in% .depgraph_allowed_derivation_scopes) {
        issues[[length(issues) + 1L]] <- .new_validation_issue(
          level = "semantic",
          severity = "warning",
          code = "invalid_featureset_derivation_scope",
          message = paste0(
            "FeatureSet node `", feature_nodes$node_id[[i]], "` has invalid `derivation_scope`: ",
            derivation_scope
          ),
          node_ids = feature_nodes$node_id[[i]]
        )
      }
    }
  }

  outcome_nodes <- node_data[node_data$node_type == "Outcome", , drop = FALSE]
  if (nrow(outcome_nodes) > 0L) {
    for (i in seq_len(nrow(outcome_nodes))) {
      observation_level <- .depgraph_extract_attr(outcome_nodes$attrs[[i]], "observation_level")
      if (!is.na(observation_level) && !observation_level %in% .depgraph_allowed_observation_levels) {
        issues[[length(issues) + 1L]] <- .new_validation_issue(
          level = "semantic",
          severity = "warning",
          code = "invalid_outcome_observation_level",
          message = paste0(
            "Outcome node `", outcome_nodes$node_id[[i]], "` has invalid `observation_level`: ",
            observation_level
          ),
          node_ids = outcome_nodes$node_id[[i]]
        )
      }
    }
  }

  time_nodes <- node_data[node_data$node_type == "Timepoint", , drop = FALSE]
  if (nrow(time_nodes) > 0L) {
    for (i in seq_len(nrow(time_nodes))) {
      value <- .depgraph_extract_attr(time_nodes$attrs[[i]], "time_index")
      converted <- suppressWarnings(as.numeric(value))
      if (!is.na(value) && is.na(converted)) {
        issues[[length(issues) + 1L]] <- .new_validation_issue(
          level = "semantic",
          severity = "warning",
          code = "invalid_time_index",
          message = paste0("Timepoint node `", time_nodes$node_id[[i]], "` has a non-numeric `time_index`."),
          node_ids = time_nodes$node_id[[i]]
        )
      }
    }
  }

  precedes <- edge_data[edge_data$edge_type == "timepoint_precedes", , drop = FALSE]
  if (nrow(precedes) > 0L && any(precedes$from == precedes$to)) {
    issues[[length(issues) + 1L]] <- .new_validation_issue(
      level = "semantic",
      severity = "error",
      code = "timepoint_self_loop",
      message = "`timepoint_precedes` contains self-loops.",
      edge_ids = precedes$edge_id[precedes$from == precedes$to]
    )
  }

  .depgraph_bind_issues(issues)
}

.validate_leakage <- function(graph) {
  node_data <- graph$nodes$data
  edge_data <- graph$edges$data
  issues <- list()

  subject_edges <- edge_data[edge_data$edge_type == "sample_belongs_to_subject", c("edge_id", "from", "to"), drop = FALSE]
  if (nrow(subject_edges) > 0L) {
    sample_counts <- stats::aggregate(
      subject_edges$from,
      by = list(subject = subject_edges$to),
      FUN = function(x) length(unique(x))
    )
    repeated <- sample_counts[sample_counts$x > 1L, , drop = FALSE]
    for (i in seq_len(nrow(repeated))) {
      subject_id <- repeated$subject[[i]]
      sample_ids <- unique(subject_edges$from[subject_edges$to == subject_id])
      edge_ids <- subject_edges$edge_id[subject_edges$to == subject_id]
      issues[[length(issues) + 1L]] <- .new_validation_issue(
        level = "leakage",
        severity = "advisory",
        code = "repeated_subject_samples",
        message = paste0("Subject `", subject_id, "` is linked to multiple samples."),
        node_ids = c(subject_id, sample_ids),
        edge_ids = edge_ids
      )
    }
  }

  study_edges <- edge_data[edge_data$edge_type == "sample_from_study", c("edge_id", "from", "to"), drop = FALSE]
  if (nrow(subject_edges) > 0L && nrow(study_edges) > 0L) {
    subject_study <- merge(
      subject_edges[, c("from", "to")],
      study_edges[, c("from", "to")],
      by = "from",
      suffixes = c("_subject", "_study")
    )
    if (nrow(subject_study) > 0L) {
      study_counts <- stats::aggregate(
        subject_study$to_study,
        by = list(subject = subject_study$to_subject),
        FUN = function(x) length(unique(x))
      )
      cross_study <- study_counts[study_counts$x > 1L, , drop = FALSE]
      for (i in seq_len(nrow(cross_study))) {
        subject_id <- cross_study$subject[[i]]
        linked_samples <- unique(subject_study$from[subject_study$to_subject == subject_id])
        issues[[length(issues) + 1L]] <- .new_validation_issue(
          level = "leakage",
          severity = "warning",
          code = "subject_cross_study_overlap",
          message = paste0("Subject `", subject_id, "` appears across multiple studies."),
          node_ids = c(subject_id, linked_samples)
        )
      }
    }
  }

  feature_nodes <- node_data[node_data$node_type == "FeatureSet", , drop = FALSE]
  if (nrow(feature_nodes) > 0L) {
    for (i in seq_len(nrow(feature_nodes))) {
      derivation_scope <- .depgraph_extract_attr(feature_nodes$attrs[[i]], "derivation_scope")
      if (!is.na(derivation_scope) && identical(derivation_scope, "per_dataset")) {
        issues[[length(issues) + 1L]] <- .new_validation_issue(
          level = "leakage",
          severity = "advisory",
          code = "per_dataset_featureset",
          message = paste0("FeatureSet `", feature_nodes$node_id[[i]], "` was derived at the full-dataset scope."),
          node_ids = feature_nodes$node_id[[i]]
        )
      }
    }
  }

  feature_edges <- edge_data[edge_data$edge_type == "sample_uses_featureset", c("edge_id", "from", "to"), drop = FALSE]
  if (nrow(feature_edges) > 0L) {
    feature_counts <- stats::aggregate(
      feature_edges$from,
      by = list(featureset = feature_edges$to),
      FUN = function(x) length(unique(x))
    )
    shared <- feature_counts[feature_counts$x > 1L, , drop = FALSE]
    for (i in seq_len(nrow(shared))) {
      featureset_id <- shared$featureset[[i]]
      linked_samples <- unique(feature_edges$from[feature_edges$to == featureset_id])
      edge_ids <- feature_edges$edge_id[feature_edges$to == featureset_id]
      issues[[length(issues) + 1L]] <- .new_validation_issue(
        level = "leakage",
        severity = "advisory",
        code = "shared_featureset_provenance",
        message = paste0("FeatureSet `", featureset_id, "` is shared across multiple samples."),
        node_ids = c(featureset_id, linked_samples),
        edge_ids = edge_ids
      )
    }
  }

  has_time_edges <- any(edge_data$edge_type == "sample_collected_at_timepoint")
  if (has_time_edges) {
    has_precedes <- any(edge_data$edge_type == "timepoint_precedes")
    time_nodes <- node_data[node_data$node_type == "Timepoint", , drop = FALSE]
    time_indices <- vapply(
      time_nodes$attrs,
      function(x) {
        value <- .depgraph_extract_attr(x, "time_index")
        if (length(value) == 0L || is.na(value)) {
          return(NA_character_)
        }
        as.character(value)
      },
      FUN.VALUE = character(1)
    )
    if (!has_precedes && (length(time_indices) == 0L || all(is.na(time_indices)))) {
      issues[[length(issues) + 1L]] <- .new_validation_issue(
        level = "leakage",
        severity = "warning",
        code = "missing_time_ordering",
        message = "Timepoint-linked samples exist but no explicit time ordering metadata were found."
      )
    }
  }

  batch_edges <- edge_data[edge_data$edge_type == "sample_processed_in_batch", c("edge_id", "from", "to"), drop = FALSE]
  n_samples <- sum(node_data$node_type == "Sample")
  if (nrow(batch_edges) > 0L && n_samples > 0L) {
    batch_counts <- stats::aggregate(
      batch_edges$from,
      by = list(batch = batch_edges$to),
      FUN = function(x) length(unique(x))
    )
    threshold <- max(3L, ceiling(n_samples * 0.5))
    heavy <- batch_counts[batch_counts$x >= threshold, , drop = FALSE]
    for (i in seq_len(nrow(heavy))) {
      batch_id <- heavy$batch[[i]]
      issues[[length(issues) + 1L]] <- .new_validation_issue(
        level = "leakage",
        severity = "advisory",
        code = "heavy_batch_reuse",
        message = paste0("Batch `", batch_id, "` is reused across many samples."),
        node_ids = batch_id,
        edge_ids = batch_edges$edge_id[batch_edges$to == batch_id]
      )
    }
  }

  .depgraph_bind_issues(issues)
}

validate_graph <- function(graph, checks = c("ids", "references", "cardinality", "schema", "time"), error_on_fail = FALSE, levels = NULL, severities = NULL) {
  .depgraph_assert(inherits(graph, "dependency_graph"), "`graph` must be a `dependency_graph`.")

  selected_levels <- levels
  if (is.null(selected_levels)) {
    if (missing(checks)) {
      selected_levels <- .depgraph_validation_levels
    } else {
      selected_levels <- .depgraph_old_checks_to_levels(checks)
      if (length(selected_levels) == 0L) {
        selected_levels <- .depgraph_validation_levels
      }
    }
  }
  selected_levels <- unique(as.character(selected_levels))
  .depgraph_assert(all(selected_levels %in% .depgraph_validation_levels), "`levels` contains unsupported values.")

  selected_severities <- if (is.null(severities)) .depgraph_validation_severities else unique(as.character(severities))
  .depgraph_assert(all(selected_severities %in% .depgraph_validation_severities), "`severities` contains unsupported values.")

  issue_sets <- list()
  if ("structural" %in% selected_levels) {
    issue_sets[[length(issue_sets) + 1L]] <- .validate_structural(graph)
  }
  if ("semantic" %in% selected_levels) {
    issue_sets[[length(issue_sets) + 1L]] <- .validate_semantic(graph)
  }
  if ("leakage" %in% selected_levels) {
    issue_sets[[length(issue_sets) + 1L]] <- .validate_leakage(graph)
  }

  issues <- .depgraph_bind_issues(issue_sets)
  if (nrow(issues) > 0L) {
    issues <- issues[issues$severity %in% selected_severities, , drop = FALSE]
    row.names(issues) <- NULL
  }

  metrics <- list(
    n_nodes = nrow(graph$nodes$data),
    n_edges = nrow(graph$edges$data),
    n_samples = sum(graph$nodes$data$node_type == "Sample")
  )

  report <- depgraph_validation_report(
    graph_name = graph$metadata$graph_name,
    issues = issues,
    metrics = metrics,
    metadata = list(
      levels = selected_levels,
      severities = selected_severities
    )
  )

  if (!isTRUE(report$valid) && isTRUE(error_on_fail)) {
    stop(paste(c("Graph validation failed.", report$errors), collapse = "\n"), call. = FALSE)
  }

  report
}

validate_depgraph <- function(graph, checks = c("ids", "references", "cardinality", "schema", "time"), error_on_fail = FALSE, levels = NULL, severities = NULL) {
  validate_graph(
    graph = graph,
    checks = checks,
    error_on_fail = error_on_fail,
    levels = levels,
    severities = severities
  )
}
