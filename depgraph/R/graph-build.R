# Graph assembly helpers.

build_dependency_graph <- function(nodes, edges, graph_name = NULL, dataset_name = NULL, validate = TRUE, validation_overrides = list()) {
  node_data <- .depgraph_bind_data(nodes, "data")
  edge_data <- .depgraph_bind_data(edges, "data")

  .depgraph_assert(any(node_data$node_type == "Sample"), "A dependency graph must contain at least one `Sample` node.")
  .depgraph_assert(anyDuplicated(node_data$node_id) == 0L, "Duplicate `node_id` values found in node sets.")
  .depgraph_assert(anyDuplicated(edge_data$edge_id) == 0L, "Duplicate `edge_id` values found in edge sets.")
  .depgraph_assert(
    all(edge_data$from %in% node_data$node_id),
    "Edge data contains `from` references that are not present in the node table."
  )
  .depgraph_assert(
    all(edge_data$to %in% node_data$node_id),
    "Edge data contains `to` references that are not present in the node table."
  )

  graph_obj <- dependency_graph(
    nodes = graph_node_set(node_data),
    edges = graph_edge_set(edge_data),
    graph = NULL,
    metadata = list(
      graph_name = graph_name,
      dataset_name = dataset_name,
      created_at = Sys.time(),
      schema_version = .depgraph_schema_version,
      validation_overrides = validation_overrides
    )
  )

  if (isTRUE(validate)) {
    validation <- validate_graph(graph_obj)
    if (!isTRUE(validation$valid)) {
      stop(
        paste(
          c("Graph validation failed.", validation$errors, validation$warnings),
          collapse = "\n"
        ),
        call. = FALSE
      )
    }
  }

  graph_obj
}

build_depgraph <- function(nodes, edges, graph_name = NULL, dataset_name = NULL, validate = TRUE, validation_overrides = list()) {
  build_dependency_graph(
    nodes = nodes,
    edges = edges,
    graph_name = graph_name,
    dataset_name = dataset_name,
    validate = validate,
    validation_overrides = validation_overrides
  )
}

as_igraph <- function(x) {
  .depgraph_assert(inherits(x, "dependency_graph"), "`x` must be a `dependency_graph`.")
  x$graph
}
