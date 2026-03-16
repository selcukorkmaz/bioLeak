# Constructors for core depgraph S3 classes.

graph_node_set <- function(data = NULL, schema_version = .depgraph_schema_version, source = list()) {
  data <- .depgraph_normalize_node_data(data)
  structure(
    list(
      data = data,
      schema_version = as.character(schema_version),
      source = source
    ),
    class = "graph_node_set"
  )
}

graph_edge_set <- function(data = NULL, schema_version = .depgraph_schema_version, source = list()) {
  data <- .depgraph_normalize_edge_data(data)
  structure(
    list(
      data = data,
      schema_version = as.character(schema_version),
      source = source
    ),
    class = "graph_edge_set"
  )
}

dependency_graph <- function(nodes, edges, graph, metadata = list(), caches = list()) {
  .depgraph_assert(inherits(nodes, "graph_node_set"), "`nodes` must be a `graph_node_set`.")
  .depgraph_assert(inherits(edges, "graph_edge_set"), "`edges` must be a `graph_edge_set`.")
  if (missing(graph) || is.null(graph)) {
    graph <- .depgraph_build_igraph(nodes$data, edges$data)
  }
  .depgraph_assert(inherits(graph, "igraph"), "`graph` must be an `igraph` object.")

  if (length(caches) == 0L) {
    caches <- .depgraph_build_caches(nodes$data, edges$data)
  }

  structure(
    list(
      nodes = nodes,
      edges = edges,
      graph = graph,
      metadata = utils::modifyList(
        list(
          graph_name = NULL,
          dataset_name = NULL,
          created_at = Sys.time(),
          schema_version = .depgraph_schema_version
        ),
        metadata
      ),
      caches = caches
    ),
    class = "dependency_graph"
  )
}

new_depgraph_nodes <- function(data = NULL, schema_version = .depgraph_schema_version, source = list()) {
  graph_node_set(data = data, schema_version = schema_version, source = source)
}

new_depgraph_edges <- function(data = NULL, schema_version = .depgraph_schema_version, source = list()) {
  graph_edge_set(data = data, schema_version = schema_version, source = source)
}

new_depgraph <- function(nodes, edges, graph = NULL, metadata = list(), caches = list()) {
  dependency_graph(nodes = nodes, edges = edges, graph = graph, metadata = metadata, caches = caches)
}

graph_query_result <- function(query = "", params = list(), nodes = NULL, edges = NULL, table = NULL, metadata = list()) {
  structure(
    list(
      query = as.character(query)[1L],
      params = params,
      nodes = if (is.null(nodes)) data.frame(stringsAsFactors = FALSE) else nodes,
      edges = if (is.null(edges)) data.frame(stringsAsFactors = FALSE) else edges,
      table = if (is.null(table)) data.frame(stringsAsFactors = FALSE) else table,
      metadata = metadata
    ),
    class = "graph_query_result"
  )
}

dependency_constraint <- function(constraint_id, relation_types, sample_map, transitive = TRUE, metadata = list()) {
  .depgraph_assert(is.data.frame(sample_map), "`sample_map` must be a data.frame.")
  structure(
    list(
      constraint_id = as.character(constraint_id)[1L],
      relation_types = as.character(relation_types),
      sample_map = sample_map,
      transitive = isTRUE(transitive),
      metadata = metadata
    ),
    class = "dependency_constraint"
  )
}

split_constraint <- function(strategy, sample_map, recommended_bioleak_args = list(), metadata = list()) {
  .depgraph_assert(is.data.frame(sample_map), "`sample_map` must be a data.frame.")
  structure(
    list(
      strategy = as.character(strategy)[1L],
      sample_map = sample_map,
      recommended_bioleak_args = recommended_bioleak_args,
      metadata = metadata
    ),
    class = "split_constraint"
  )
}

leakage_constraint <- function(issue_type, severity, affected_samples, evidence = NULL, recommendation = "", metadata = list()) {
  structure(
    list(
      issue_type = as.character(issue_type)[1L],
      severity = as.character(severity)[1L],
      affected_samples = as.character(affected_samples),
      evidence = if (is.null(evidence)) data.frame(stringsAsFactors = FALSE) else evidence,
      recommendation = as.character(recommendation)[1L],
      metadata = metadata
    ),
    class = "leakage_constraint"
  )
}

depgraph_validation_report <- function(graph_name = NULL, issues = NULL, metrics = list(), metadata = list()) {
  if (is.null(issues)) {
    issues <- data.frame(
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

  summary <- .summarize_validation_issues(issues)
  errors <- unique(issues$message[issues$severity == "error"])
  warnings <- unique(issues$message[issues$severity == "warning"])
  advisories <- unique(issues$message[issues$severity == "advisory"])

  structure(
    list(
      graph_name = graph_name,
      valid = length(errors) == 0L,
      issues = issues,
      summary = summary,
      metadata = utils::modifyList(
        list(
          created_at = Sys.time(),
          schema_version = .depgraph_schema_version
        ),
        metadata
      ),
      errors = errors,
      warnings = warnings,
      advisories = advisories,
      metrics = utils::modifyList(
        list(
          n_issues = nrow(issues),
          n_errors = length(errors),
          n_warnings = length(warnings),
          n_advisories = length(advisories)
        ),
        metrics
      )
    ),
    class = "depgraph_validation_report"
  )
}

bioleak_split_spec <- function(sample_data = NULL, group_var = "group_id", block_vars = character(), time_var = NULL, ordering_required = FALSE, constraint_mode = NULL, constraint_strategy = NULL, recommended_resampling = NULL, metadata = list()) {
  if (is.null(sample_data)) {
    sample_data <- data.frame(
      sample_id = character(),
      sample_node_id = character(),
      group_id = character(),
      primary_group = character(),
      batch_group = character(),
      study_group = character(),
      timepoint_id = character(),
      time_index = numeric(),
      order_rank = integer(),
      stringsAsFactors = FALSE
    )
  }

  .depgraph_assert(is.data.frame(sample_data), "`sample_data` must be a data.frame.")

  structure(
    list(
      sample_data = sample_data,
      group_var = as.character(group_var)[1L],
      block_vars = as.character(block_vars),
      time_var = if (is.null(time_var)) NULL else as.character(time_var)[1L],
      ordering_required = isTRUE(ordering_required),
      constraint_mode = if (is.null(constraint_mode)) NULL else as.character(constraint_mode)[1L],
      constraint_strategy = if (is.null(constraint_strategy)) NULL else as.character(constraint_strategy)[1L],
      recommended_resampling = if (is.null(recommended_resampling)) NULL else as.character(recommended_resampling)[1L],
      metadata = metadata
    ),
    class = "bioleak_split_spec"
  )
}

bioleak_split_spec_validation <- function(issues = NULL, metadata = list()) {
  if (is.null(issues)) {
    issues <- data.frame(
      issue_id = character(),
      severity = character(),
      code = character(),
      message = character(),
      n_affected = integer(),
      details = I(list()),
      stringsAsFactors = FALSE
    )
  }

  summary <- list(
    by_severity = .depgraph_count_table(if (nrow(issues) == 0L) character() else issues$severity),
    by_code = .depgraph_count_table(if (nrow(issues) == 0L) character() else issues$code)
  )
  valid <- !any(issues$severity %in% "error")

  structure(
    list(
      valid = valid,
      issues = issues,
      summary = summary,
      metadata = metadata
    ),
    class = "bioleak_split_spec_validation"
  )
}

leakage_risk_summary <- function(overview = character(), diagnostics = NULL, validation_summary = list(), constraint_summary = list(), split_spec_summary = list(), metadata = list()) {
  if (is.null(diagnostics)) {
    diagnostics <- data.frame(
      severity = character(),
      category = character(),
      message = character(),
      source = character(),
      n_affected = integer(),
      stringsAsFactors = FALSE
    )
  }

  structure(
    list(
      overview = as.character(overview),
      diagnostics = diagnostics,
      validation_summary = validation_summary,
      constraint_summary = constraint_summary,
      split_spec_summary = split_spec_summary,
      metadata = metadata
    ),
    class = "leakage_risk_summary"
  )
}
