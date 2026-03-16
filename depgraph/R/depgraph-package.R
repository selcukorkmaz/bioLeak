# Core package constants and internal helpers for depgraph.

.depgraph_schema_version <- "0.1.0"

.depgraph_node_types <- c(
  "Sample", "Subject", "Batch", "Study",
  "Timepoint", "Assay", "FeatureSet", "Outcome"
)

.depgraph_prefix_map <- c(
  Sample = "sample",
  Subject = "subject",
  Batch = "batch",
  Study = "study",
  Timepoint = "timepoint",
  Assay = "assay",
  FeatureSet = "featureset",
  Outcome = "outcome"
)

.depgraph_node_schema <- list(
  Sample = list(
    key_field = "sample_id",
    required_attrs = character(),
    optional_attrs = c(
      "subject_id", "batch_id", "study_id", "timepoint_id",
      "assay_id", "featureset_id", "outcome_id", "sample_role",
      "collection_date", "source_file"
    )
  ),
  Subject = list(
    key_field = "subject_id",
    required_attrs = character(),
    optional_attrs = c(
      "subject_type", "species", "sex", "site", "enrollment_study"
    )
  ),
  Batch = list(
    key_field = "batch_id",
    required_attrs = character(),
    optional_attrs = c(
      "batch_type", "platform", "processing_center", "run_date",
      "operator"
    )
  ),
  Study = list(
    key_field = "study_id",
    required_attrs = character(),
    optional_attrs = c(
      "study_name", "study_type", "institution", "country",
      "study_phase"
    )
  ),
  Timepoint = list(
    key_field = "timepoint_id",
    required_attrs = character(),
    optional_attrs = c(
      "time_index", "time_value", "time_unit", "collection_datetime",
      "visit_label"
    )
  ),
  Assay = list(
    key_field = "assay_id",
    required_attrs = character(),
    optional_attrs = c(
      "modality", "platform", "vendor", "protocol_version"
    )
  ),
  FeatureSet = list(
    key_field = "featureset_id",
    required_attrs = character(),
    optional_attrs = c(
      "featureset_name", "modality", "feature_count",
      "derivation_scope", "derivation_version", "generated_from"
    )
  ),
  Outcome = list(
    key_field = "outcome_id",
    required_attrs = character(),
    optional_attrs = c(
      "outcome_name", "outcome_type", "outcome_value",
      "outcome_scale", "observation_level"
    )
  )
)

.depgraph_allowed_derivation_scopes <- c(
  "external",
  "per_training_split",
  "per_dataset"
)

.depgraph_allowed_observation_levels <- c(
  "sample",
  "subject"
)

.depgraph_validation_levels <- c(
  "structural",
  "semantic",
  "leakage"
)

.depgraph_validation_severities <- c(
  "error",
  "warning",
  "advisory"
)

.depgraph_edge_schema <- data.frame(
  edge_type = c(
    "sample_belongs_to_subject",
    "sample_processed_in_batch",
    "sample_from_study",
    "sample_collected_at_timepoint",
    "sample_measured_by_assay",
    "sample_uses_featureset",
    "sample_has_outcome",
    "subject_has_outcome",
    "timepoint_precedes",
    "featureset_generated_from_study",
    "featureset_generated_from_batch"
  ),
  from_type = c(
    "Sample", "Sample", "Sample", "Sample", "Sample",
    "Sample", "Sample", "Subject", "Timepoint", "FeatureSet",
    "FeatureSet"
  ),
  to_type = c(
    "Subject", "Batch", "Study", "Timepoint", "Assay",
    "FeatureSet", "Outcome", "Outcome", "Timepoint", "Study",
    "Batch"
  ),
  stringsAsFactors = FALSE
)

.depgraph_single_target_edges <- c(
  "sample_collected_at_timepoint",
  "sample_measured_by_assay"
)

.depgraph_semantic_relation_rules <- list(
  sample_belongs_to_subject = list(
    min_targets = 0L,
    max_targets = 1L,
    missing_code = NULL,
    multiple_code = "sample_multiple_subject_assignments",
    missing_severity = NULL,
    multiple_severity = "error"
  ),
  sample_from_study = list(
    min_targets = 1L,
    max_targets = 1L,
    missing_code = "sample_missing_study_assignment",
    multiple_code = "sample_multiple_study_assignments",
    missing_severity = "warning",
    multiple_severity = "error"
  ),
  sample_processed_in_batch = list(
    min_targets = 0L,
    max_targets = 1L,
    missing_code = NULL,
    multiple_code = "sample_multiple_batch_assignments",
    missing_severity = NULL,
    multiple_severity = "warning"
  )
)

.depgraph_empty_node_data <- function() {
  data.frame(
    node_id = character(),
    node_type = character(),
    node_key = character(),
    label = character(),
    attrs = I(list()),
    stringsAsFactors = FALSE
  )
}

.depgraph_empty_edge_data <- function() {
  data.frame(
    edge_id = character(),
    from = character(),
    to = character(),
    edge_type = character(),
    attrs = I(list()),
    stringsAsFactors = FALSE
  )
}

.depgraph_empty_attrs <- function(n) {
  rep(list(list()), n)
}

.depgraph_assert <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

.depgraph_match_node_type <- function(type) {
  .depgraph_assert(length(type) == 1L && !is.na(type), "`type` must be length 1.")
  match_idx <- match(tolower(type), tolower(.depgraph_node_types))
  .depgraph_assert(!is.na(match_idx), paste0("Unsupported node type: ", type))
  .depgraph_node_types[[match_idx]]
}

.depgraph_prefix_for_type <- function(type) {
  type <- .depgraph_match_node_type(type)
  .depgraph_prefix_map[[type]]
}

.depgraph_node_type_schema <- function(type) {
  type <- .depgraph_match_node_type(type)
  .depgraph_node_schema[[type]]
}

.depgraph_edge_type_schema <- function(edge_type) {
  edge_type <- as.character(edge_type)[1L]
  .depgraph_edge_schema[.depgraph_edge_schema$edge_type == edge_type, , drop = FALSE]
}

.depgraph_make_node_id <- function(type, key, prefix = TRUE) {
  key <- as.character(key)
  .depgraph_assert(all(!is.na(key) & nzchar(key)), "Node keys must be non-missing strings.")
  if (!isTRUE(prefix)) {
    return(key)
  }
  paste0(.depgraph_prefix_for_type(type), ":", key)
}

.depgraph_normalize_attr_entry <- function(x, context = "`attrs`") {
  if (is.null(x)) {
    return(list())
  }

  .depgraph_assert(is.list(x), paste0(context, " entries must be lists."))
  if (length(x) == 0L) {
    return(list())
  }

  attr_names <- names(x)
  .depgraph_assert(
    !is.null(attr_names) && all(!is.na(attr_names)) && all(nzchar(attr_names)),
    paste0(context, " entries must be named lists.")
  )
  .depgraph_assert(
    anyDuplicated(attr_names) == 0L,
    paste0(context, " entries cannot contain duplicate names.")
  )

  as.list(x)
}

.depgraph_make_attr_list <- function(data, attr_cols = NULL) {
  if (is.null(attr_cols) || length(attr_cols) == 0L) {
    return(I(.depgraph_empty_attrs(nrow(data))))
  }

  missing_cols <- setdiff(attr_cols, names(data))
  .depgraph_assert(
    length(missing_cols) == 0L,
    paste0("Missing attribute columns: ", paste(missing_cols, collapse = ", "))
  )

  attrs <- lapply(seq_len(nrow(data)), function(i) {
    as.list(data[i, attr_cols, drop = FALSE])
  })
  I(attrs)
}

.depgraph_validate_node_attrs <- function(node_type, attrs) {
  schema <- .depgraph_node_type_schema(node_type)
  attr_names <- names(attrs)

  missing_required <- setdiff(schema$required_attrs, attr_names)
  unknown_attrs <- setdiff(attr_names, c(schema$required_attrs, schema$optional_attrs))

  list(
    missing_required = missing_required,
    unknown_attrs = unknown_attrs
  )
}

.depgraph_build_igraph <- function(node_data, edge_data) {
  vertices <- node_data[, c("node_id", "node_type", "node_key", "label"), drop = FALSE]
  edges_for_graph <- edge_data[, c("from", "to", "edge_id", "edge_type"), drop = FALSE]

  igraph::graph_from_data_frame(
    d = edges_for_graph,
    vertices = vertices,
    directed = TRUE
  )
}

.depgraph_build_caches <- function(node_data, edge_data) {
  list(
    sample_ids = node_data$node_key[node_data$node_type == "Sample"],
    node_index = stats::setNames(seq_len(nrow(node_data)), node_data$node_id),
    edge_index = stats::setNames(seq_len(nrow(edge_data)), edge_data$edge_id)
  )
}

.depgraph_normalize_node_data <- function(data) {
  if (is.null(data)) {
    return(.depgraph_empty_node_data())
  }

  .depgraph_assert(is.data.frame(data), "`data` must be a data.frame.")

  required_cols <- c("node_id", "node_type", "node_key", "label")
  missing_cols <- setdiff(required_cols, names(data))
  .depgraph_assert(
    length(missing_cols) == 0L,
    paste0("Node data is missing columns: ", paste(missing_cols, collapse = ", "))
  )

  if (!"attrs" %in% names(data)) {
    data$attrs <- I(.depgraph_empty_attrs(nrow(data)))
  }

  .depgraph_assert(is.list(data$attrs), "`attrs` must be a list column.")

  data$node_id <- as.character(data$node_id)
  data$node_type <- as.character(data$node_type)
  data$node_key <- as.character(data$node_key)
  data$label <- as.character(data$label)
  data$attrs <- I(lapply(data$attrs, .depgraph_normalize_attr_entry, context = "`attrs`"))

  .depgraph_assert(
    all(!is.na(data$node_id) & nzchar(data$node_id)),
    "Node data contains missing or empty `node_id` values."
  )
  .depgraph_assert(
    all(!is.na(data$node_key) & nzchar(data$node_key)),
    "Node data contains missing or empty `node_key` values."
  )
  .depgraph_assert(
    all(!is.na(data$label) & nzchar(data$label)),
    "Node data contains missing or empty `label` values."
  )

  node_types <- vapply(data$node_type, .depgraph_match_node_type, character(1), USE.NAMES = FALSE)
  data$node_type <- node_types

  data[, c("node_id", "node_type", "node_key", "label", "attrs"), drop = FALSE]
}

.depgraph_normalize_edge_data <- function(data) {
  if (is.null(data)) {
    return(.depgraph_empty_edge_data())
  }

  .depgraph_assert(is.data.frame(data), "`data` must be a data.frame.")

  required_cols <- c("edge_id", "from", "to", "edge_type")
  missing_cols <- setdiff(required_cols, names(data))
  .depgraph_assert(
    length(missing_cols) == 0L,
    paste0("Edge data is missing columns: ", paste(missing_cols, collapse = ", "))
  )

  if (!"attrs" %in% names(data)) {
    data$attrs <- I(.depgraph_empty_attrs(nrow(data)))
  }

  .depgraph_assert(is.list(data$attrs), "`attrs` must be a list column.")

  data$edge_id <- as.character(data$edge_id)
  data$from <- as.character(data$from)
  data$to <- as.character(data$to)
  data$edge_type <- as.character(data$edge_type)
  data$attrs <- I(lapply(data$attrs, .depgraph_normalize_attr_entry, context = "`attrs`"))

  .depgraph_assert(
    all(!is.na(data$edge_id) & nzchar(data$edge_id)),
    "Edge data contains missing or empty `edge_id` values."
  )
  .depgraph_assert(
    all(!is.na(data$from) & nzchar(data$from)),
    "Edge data contains missing or empty `from` values."
  )
  .depgraph_assert(
    all(!is.na(data$to) & nzchar(data$to)),
    "Edge data contains missing or empty `to` values."
  )
  .depgraph_assert(
    all(!is.na(data$edge_type) & nzchar(data$edge_type)),
    "Edge data contains missing or empty `edge_type` values."
  )

  data[, c("edge_id", "from", "to", "edge_type", "attrs"), drop = FALSE]
}

.depgraph_count_table <- function(x) {
  counts <- sort(table(x), decreasing = TRUE)
  data.frame(
    value = names(counts),
    n = as.integer(counts),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

.depgraph_bind_data <- function(objects, field) {
  if (inherits(objects, c("graph_node_set", "graph_edge_set"))) {
    objects <- list(objects)
  }
  .depgraph_assert(is.list(objects) && length(objects) > 0L, "`objects` must be a non-empty list.")
  pieces <- lapply(objects, function(x) x[[field]])
  out <- do.call(rbind, pieces)
  row.names(out) <- NULL
  out
}

.depgraph_extract_attr <- function(attrs, key) {
  value <- attrs[[key]]
  if (is.null(value) || length(value) == 0L) {
    return(NA)
  }
  value[[1L]]
}

.depgraph_validation_override <- function(graph, key, default = FALSE) {
  overrides <- graph$metadata$validation_overrides
  if (is.null(overrides) || is.null(overrides[[key]])) {
    return(default)
  }
  isTRUE(overrides[[key]])
}
