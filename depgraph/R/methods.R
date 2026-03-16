# S3 methods for depgraph core classes.

print.graph_node_set <- function(x, ...) {
  cat(
    "<graph_node_set>",
    nrow(x$data), "nodes across",
    length(unique(x$data$node_type)), "types",
    "(schema", x$schema_version, ")\n"
  )
  invisible(x)
}

summary.graph_node_set <- function(object, ...) {
  list(
    schema_version = object$schema_version,
    n_nodes = nrow(object$data),
    node_types = .depgraph_count_table(object$data$node_type)
  )
}

as.data.frame.graph_node_set <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$data
}

print.graph_edge_set <- function(x, ...) {
  cat(
    "<graph_edge_set>",
    nrow(x$data), "edges across",
    length(unique(x$data$edge_type)), "relations",
    "(schema", x$schema_version, ")\n"
  )
  invisible(x)
}

summary.graph_edge_set <- function(object, ...) {
  list(
    schema_version = object$schema_version,
    n_edges = nrow(object$data),
    edge_types = .depgraph_count_table(object$data$edge_type)
  )
}

as.data.frame.graph_edge_set <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$data
}

print.dependency_graph <- function(x, ...) {
  graph_name <- x$metadata$graph_name %||% "<unnamed>"
  cat(
    "<dependency_graph>",
    graph_name,
    "\n  Nodes:", igraph::vcount(x$graph),
    "\n  Edges:", igraph::ecount(x$graph),
    "\n"
  )
  invisible(x)
}

summary.dependency_graph <- function(object, ...) {
  list(
    graph_name = object$metadata$graph_name,
    dataset_name = object$metadata$dataset_name,
    schema_version = object$metadata$schema_version,
    n_nodes = nrow(object$nodes$data),
    n_edges = nrow(object$edges$data),
    node_types = .depgraph_count_table(object$nodes$data$node_type),
    edge_types = .depgraph_count_table(object$edges$data$edge_type)
  )
}

plot.dependency_graph <- function(x, ...) {
  graphics::plot(x$graph, ...)
  invisible(x)
}

print.graph_query_result <- function(x, ...) {
  cat("<graph_query_result>", x$query, "\n")
  cat("  Rows:", nrow(x$table), "\n")
  invisible(x)
}

summary.graph_query_result <- function(object, ...) {
  list(
    query = object$query,
    n_nodes = nrow(object$nodes),
    n_edges = nrow(object$edges),
    n_rows = nrow(object$table)
  )
}

as.data.frame.graph_query_result <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$table
}

print.dependency_constraint <- function(x, ...) {
  cat("<dependency_constraint>", x$constraint_id, "\n")
  cat("  Samples:", nrow(x$sample_map), "\n")
  invisible(x)
}

summary.dependency_constraint <- function(object, ...) {
  group_col <- intersect(c("group_id", "split_unit", "block_id"), names(object$sample_map))
  n_groups <- if (length(group_col) == 0L) NA_integer_ else length(unique(object$sample_map[[group_col[[1L]]]]))
  list(
    constraint_id = object$constraint_id,
    relation_types = object$relation_types,
    n_samples = nrow(object$sample_map),
    n_groups = n_groups,
    transitive = object$transitive
  )
}

as.data.frame.dependency_constraint <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$sample_map
}

print.split_constraint <- function(x, ...) {
  cat("<split_constraint>", x$strategy, "\n")
  cat("  Samples:", nrow(x$sample_map), "\n")
  if ("group_id" %in% names(x$sample_map)) {
    cat("  Groups:", length(unique(x$sample_map$group_id)), "\n")
  }
  warnings <- x$metadata$warnings %||% character()
  if (length(warnings) > 0L) {
    cat("  Warnings:", length(warnings), "\n")
  }
  invisible(x)
}

summary.split_constraint <- function(object, ...) {
  warnings <- object$metadata$warnings %||% character()
  list(
    strategy = object$strategy,
    mode = object$metadata$mode %||% NA_character_,
    n_samples = nrow(object$sample_map),
    n_groups = if ("group_id" %in% names(object$sample_map)) length(unique(object$sample_map$group_id)) else NA_integer_,
    n_warnings = length(warnings),
    columns = names(object$sample_map)
  )
}

as.data.frame.split_constraint <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$sample_map
}

print.leakage_constraint <- function(x, ...) {
  cat("<leakage_constraint>", x$issue_type, "(", x$severity, ")\n")
  cat("  Affected samples:", length(x$affected_samples), "\n")
  invisible(x)
}

summary.leakage_constraint <- function(object, ...) {
  list(
    issue_type = object$issue_type,
    severity = object$severity,
    n_affected_samples = length(object$affected_samples),
    recommendation = object$recommendation
  )
}

as.data.frame.leakage_constraint <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$evidence
}

print.depgraph_validation_report <- function(x, ...) {
  cat("<depgraph_validation_report>")
  if (!is.null(x$graph_name) && nzchar(x$graph_name)) {
    cat(" ", x$graph_name)
  }
  cat(
    "\n  Valid:", x$valid,
    "\n  Issues:", nrow(x$issues),
    "\n"
  )

  severity_counts <- x$summary$by_severity
  if (nrow(severity_counts) > 0L) {
    cat("  By severity:\n")
    for (i in seq_len(nrow(severity_counts))) {
      cat("   -", severity_counts$value[[i]], ":", severity_counts$n[[i]], "\n")
    }
  }
  invisible(x)
}

summary.depgraph_validation_report <- function(object, ...) {
  object$summary
}

as.data.frame.depgraph_validation_report <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$issues
}

print.bioleak_split_spec <- function(x, ...) {
  cat("<bioleak_split_spec>", x$constraint_mode %||% "<unknown>", "\n")
  cat("  Samples:", nrow(x$sample_data), "\n")
  cat("  Groups:", length(unique(x$sample_data[[x$group_var]])), "\n")
  if (!is.null(x$recommended_resampling)) {
    cat("  Recommended resampling:", x$recommended_resampling, "\n")
  }
  invisible(x)
}

summary.bioleak_split_spec <- function(object, ...) {
  list(
    constraint_mode = object$constraint_mode,
    constraint_strategy = object$constraint_strategy,
    n_samples = nrow(object$sample_data),
    n_groups = length(unique(object$sample_data[[object$group_var]])),
    block_vars = object$block_vars,
    time_var = object$time_var,
    ordering_required = object$ordering_required,
    recommended_resampling = object$recommended_resampling
  )
}

as.data.frame.bioleak_split_spec <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$sample_data
}

print.bioleak_split_spec_validation <- function(x, ...) {
  cat("<bioleak_split_spec_validation>\n")
  cat("  Valid:", x$valid, "\n")
  cat("  Issues:", nrow(x$issues), "\n")
  invisible(x)
}

summary.bioleak_split_spec_validation <- function(object, ...) {
  object$summary
}

as.data.frame.bioleak_split_spec_validation <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$issues
}

print.leakage_risk_summary <- function(x, ...) {
  cat("<leakage_risk_summary>\n")
  cat("  Overview:", paste(x$overview, collapse = " "), "\n")
  cat("  Diagnostics:", nrow(x$diagnostics), "\n")
  invisible(x)
}

summary.leakage_risk_summary <- function(object, ...) {
  list(
    overview = object$overview,
    n_diagnostics = nrow(object$diagnostics),
    by_source = .depgraph_count_table(if (nrow(object$diagnostics) == 0L) character() else object$diagnostics$source),
    by_severity = .depgraph_count_table(if (nrow(object$diagnostics) == 0L) character() else object$diagnostics$severity)
  )
}

as.data.frame.leakage_risk_summary <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$diagnostics
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) {
    return(y)
  }
  x
}
