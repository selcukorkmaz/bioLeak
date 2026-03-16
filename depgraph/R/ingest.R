# Metadata ingestion and canonical node/edge builders.

ingest_metadata <- function(data, col_map = NULL, dataset_name = NULL, strict = TRUE) {
  .depgraph_assert(is.data.frame(data), "`data` must be a data.frame.")

  out <- data
  if (!is.null(col_map)) {
    .depgraph_assert(is.character(col_map) && !is.null(names(col_map)), "`col_map` must be a named character vector.")
    missing_cols <- unname(col_map)[!unname(col_map) %in% names(out)]
    .depgraph_assert(
      length(missing_cols) == 0L,
      paste0("Mapped columns not found in `data`: ", paste(missing_cols, collapse = ", "))
    )
    rename_idx <- match(unname(col_map), names(out))
    names(out)[rename_idx] <- names(col_map)
  }

  required_cols <- "sample_id"
  if (isTRUE(strict)) {
    missing_required <- setdiff(required_cols, names(out))
    .depgraph_assert(
      length(missing_required) == 0L,
      paste0("Missing required metadata columns: ", paste(missing_required, collapse = ", "))
    )
  }

  id_like <- intersect(
    names(out),
    c(
      "sample_id", "subject_id", "batch_id", "study_id",
      "timepoint_id", "assay_id", "featureset_id", "outcome_id"
    )
  )

  for (col in id_like) {
    out[[col]] <- as.character(out[[col]])
  }

  attr(out, "dataset_name") <- dataset_name
  out
}

create_nodes <- function(data, type, id_col, label_col = NULL, attr_cols = NULL, prefix = TRUE, dedupe = TRUE) {
  .depgraph_assert(is.data.frame(data), "`data` must be a data.frame.")
  .depgraph_assert(id_col %in% names(data), paste0("Missing `id_col`: ", id_col))

  type <- .depgraph_match_node_type(type)
  schema <- .depgraph_node_type_schema(type)

  if (identical(id_col, schema$key_field) && is.null(attr_cols)) {
    attr_cols <- intersect(schema$optional_attrs, names(data))
  }

  keep_cols <- unique(c(id_col, label_col, attr_cols))
  keep_cols <- keep_cols[!is.na(keep_cols) & nzchar(keep_cols)]
  selected <- data[, keep_cols, drop = FALSE]
  selected <- selected[!is.na(selected[[id_col]]) & nzchar(selected[[id_col]]), , drop = FALSE]

  if (isTRUE(dedupe) && anyDuplicated(selected[[id_col]]) > 0L) {
    dup_ids <- unique(selected[[id_col]][duplicated(selected[[id_col]])])
    dup_rows <- selected[selected[[id_col]] %in% dup_ids, , drop = FALSE]
    if (nrow(unique(dup_rows)) > length(dup_ids)) {
      warning("Duplicate node identifiers with conflicting attributes were collapsed to the first row.", call. = FALSE)
    }
    selected <- selected[!duplicated(selected[[id_col]]), , drop = FALSE]
  }

  labels <- if (is.null(label_col)) selected[[id_col]] else selected[[label_col]]
  label_missing <- is.na(labels) | !nzchar(as.character(labels))
  if (any(label_missing)) {
    labels[label_missing] <- selected[[id_col]][label_missing]
  }

  node_data <- data.frame(
    node_id = .depgraph_make_node_id(type, selected[[id_col]], prefix = prefix),
    node_type = type,
    node_key = as.character(selected[[id_col]]),
    label = as.character(labels),
    attrs = .depgraph_make_attr_list(selected, attr_cols),
    stringsAsFactors = FALSE
  )

  graph_node_set(node_data, source = list(type = type, id_col = id_col))
}

create_edges <- function(data, from_col, to_col, from_type, to_type, relation, attr_cols = NULL, allow_missing = FALSE, dedupe = TRUE) {
  .depgraph_assert(is.data.frame(data), "`data` must be a data.frame.")
  .depgraph_assert(from_col %in% names(data), paste0("Missing `from_col`: ", from_col))
  .depgraph_assert(to_col %in% names(data), paste0("Missing `to_col`: ", to_col))

  from_type <- .depgraph_match_node_type(from_type)
  to_type <- .depgraph_match_node_type(to_type)
  relation <- as.character(relation)[1L]

  schema_row <- .depgraph_edge_type_schema(relation)
  if (nrow(schema_row) == 1L) {
    .depgraph_assert(
      identical(schema_row$from_type[[1L]], from_type) && identical(schema_row$to_type[[1L]], to_type),
      paste0("Relation `", relation, "` expects ", schema_row$from_type[[1L]], " -> ", schema_row$to_type[[1L]], ".")
    )
  }

  keep_cols <- unique(c(from_col, to_col, attr_cols))
  keep_cols <- keep_cols[!is.na(keep_cols) & nzchar(keep_cols)]
  selected <- data[, keep_cols, drop = FALSE]

  missing_endpoint <- is.na(selected[[from_col]]) | !nzchar(as.character(selected[[from_col]])) |
    is.na(selected[[to_col]]) | !nzchar(as.character(selected[[to_col]]))

  if (any(missing_endpoint)) {
    if (isTRUE(allow_missing)) {
      selected <- selected[!missing_endpoint, , drop = FALSE]
    } else {
      stop("Missing endpoint values found while creating edges.", call. = FALSE)
    }
  }

  edge_data <- data.frame(
    from = .depgraph_make_node_id(from_type, selected[[from_col]], prefix = TRUE),
    to = .depgraph_make_node_id(to_type, selected[[to_col]], prefix = TRUE),
    edge_type = relation,
    attrs = .depgraph_make_attr_list(selected, attr_cols),
    stringsAsFactors = FALSE
  )

  if (isTRUE(dedupe) && nrow(edge_data) > 0L) {
    keep <- !duplicated(edge_data[, c("from", "to", "edge_type"), drop = FALSE])
    edge_data <- edge_data[keep, , drop = FALSE]
  }

  edge_data$edge_id <- paste0(edge_data$edge_type, ":", seq_len(nrow(edge_data)))
  edge_data <- edge_data[, c("edge_id", "from", "to", "edge_type", "attrs"), drop = FALSE]

  graph_edge_set(edge_data, source = list(relation = relation, from_col = from_col, to_col = to_col))
}
