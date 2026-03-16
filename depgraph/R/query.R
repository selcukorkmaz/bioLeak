# Query helpers for dependency_graph objects.

.depgraph_resolve_node_ids <- function(graph, node_ids) {
  node_data <- graph$nodes$data
  node_ids <- unique(as.character(node_ids))
  .depgraph_assert(length(node_ids) > 0L, "`node_ids` must contain at least one value.")
  .depgraph_assert(
    all(node_ids %in% node_data$node_id),
    paste0("Unknown node IDs: ", paste(setdiff(node_ids, node_data$node_id), collapse = ", "))
  )
  node_ids
}

.depgraph_resolve_sample_node_ids <- function(graph, samples = NULL) {
  sample_nodes <- graph$nodes$data[graph$nodes$data$node_type == "Sample", , drop = FALSE]
  if (is.null(samples)) {
    return(sample_nodes$node_id)
  }

  samples <- unique(as.character(samples))
  matched <- sample_nodes$node_id[sample_nodes$node_id %in% samples]
  matched <- unique(c(
    matched,
    sample_nodes$node_id[match(samples, sample_nodes$node_key, nomatch = 0L)]
  ))
  matched <- matched[!is.na(matched) & nzchar(matched)]

  .depgraph_assert(length(matched) > 0L, "No matching sample nodes found.")
  matched
}

.depgraph_subset_edges <- function(graph, edge_types = NULL, node_ids = NULL) {
  edge_data <- graph$edges$data

  if (!is.null(edge_types)) {
    edge_types <- unique(as.character(edge_types))
    edge_data <- edge_data[edge_data$edge_type %in% edge_types, , drop = FALSE]
  }

  if (!is.null(node_ids)) {
    node_ids <- unique(as.character(node_ids))
    edge_data <- edge_data[edge_data$from %in% node_ids | edge_data$to %in% node_ids, , drop = FALSE]
  }

  edge_data
}

.depgraph_query_result <- function(graph, query, params, node_ids = NULL, edge_ids = NULL, table = NULL, metadata = list()) {
  node_data <- graph$nodes$data
  edge_data <- graph$edges$data

  nodes <- if (is.null(node_ids)) {
    node_data[0, , drop = FALSE]
  } else {
    node_data[node_data$node_id %in% unique(node_ids), , drop = FALSE]
  }

  edges <- if (is.null(edge_ids)) {
    edge_data[0, , drop = FALSE]
  } else {
    edge_data[edge_data$edge_id %in% unique(edge_ids), , drop = FALSE]
  }

  graph_query_result(
    query = query,
    params = params,
    nodes = nodes,
    edges = edges,
    table = if (is.null(table)) data.frame(stringsAsFactors = FALSE) else table,
    metadata = metadata
  )
}

.depgraph_filter_graph_by_edge_type <- function(graph, edge_types = NULL) {
  if (is.null(edge_types)) {
    return(graph$graph)
  }

  edge_types <- unique(as.character(edge_types))
  keep_ids <- graph$edges$data$edge_id[graph$edges$data$edge_type %in% edge_types]
  keep_edges <- igraph::E(graph$graph)[igraph::E(graph$graph)$edge_id %in% keep_ids]
  igraph::subgraph_from_edges(graph$graph, eids = keep_edges, delete.vertices = FALSE)
}

.depgraph_edge_between <- function(edge_data, from, to, undirected = FALSE) {
  matches <- edge_data$from == from & edge_data$to == to
  if (isTRUE(undirected)) {
    matches <- matches | (edge_data$from == to & edge_data$to == from)
  }
  edge_data[matches, , drop = FALSE]
}

.depgraph_format_paths <- function(graph, path_node_ids, edge_data, query_name, params) {
  node_data <- graph$nodes$data
  rows <- list()
  used_node_ids <- character()
  used_edge_ids <- character()

  if (length(path_node_ids) == 0L) {
    return(.depgraph_query_result(graph, query_name, params, table = data.frame(
      path_id = character(),
      step = integer(),
      node_id = character(),
      node_type = character(),
      edge_id = character(),
      edge_type = character(),
      stringsAsFactors = FALSE
    )))
  }

  row_idx <- 1L
  for (i in seq_along(path_node_ids)) {
    nodes_in_path <- path_node_ids[[i]]
    if (length(nodes_in_path) == 0L) {
      next
    }

    for (step_idx in seq_along(nodes_in_path)) {
      node_id <- nodes_in_path[[step_idx]]
      node_row <- node_data[node_data$node_id == node_id, , drop = FALSE]
      edge_id <- NA_character_
      edge_type <- NA_character_

      if (step_idx > 1L) {
        prev_node <- nodes_in_path[[step_idx - 1L]]
        edge_row <- .depgraph_edge_between(edge_data, prev_node, node_id, undirected = identical(params$mode, "all"))
        if (nrow(edge_row) > 0L) {
          edge_id <- edge_row$edge_id[[1L]]
          edge_type <- edge_row$edge_type[[1L]]
          used_edge_ids <- c(used_edge_ids, edge_id)
        }
      }

      rows[[row_idx]] <- data.frame(
        path_id = paste0("path_", i),
        step = step_idx,
        node_id = node_id,
        node_type = node_row$node_type[[1L]],
        edge_id = edge_id,
        edge_type = edge_type,
        stringsAsFactors = FALSE
      )
      used_node_ids <- c(used_node_ids, node_id)
      row_idx <- row_idx + 1L
    }
  }

  table <- do.call(rbind, rows)
  row.names(table) <- NULL

  .depgraph_query_result(
    graph = graph,
    query = query_name,
    params = params,
    node_ids = unique(used_node_ids),
    edge_ids = unique(used_edge_ids),
    table = table
  )
}

.depgraph_edge_types_for_via <- function(via, edge_types = NULL) {
  if (!is.null(edge_types)) {
    return(unique(as.character(edge_types)))
  }

  via <- unique(vapply(via, .depgraph_match_node_type, character(1), USE.NAMES = FALSE))
  .depgraph_edge_schema$edge_type[
    .depgraph_edge_schema$from_type == "Sample" &
      .depgraph_edge_schema$to_type %in% via
  ]
}

.depgraph_shared_dependency_table <- function(graph, via, samples = NULL, edge_types = NULL) {
  node_data <- graph$nodes$data
  edge_data <- graph$edges$data
  sample_nodes <- .depgraph_resolve_sample_node_ids(graph, samples)
  selected_edge_types <- .depgraph_edge_types_for_via(via, edge_types)

  edges <- edge_data[
    edge_data$edge_type %in% selected_edge_types &
      edge_data$from %in% sample_nodes,
    ,
    drop = FALSE
  ]

  rows <- list()
  row_idx <- 1L
  sample_key_map <- stats::setNames(node_data$node_key, node_data$node_id)
  node_type_map <- stats::setNames(node_data$node_type, node_data$node_id)

  if (nrow(edges) == 0L) {
    return(data.frame(
      sample_id_1 = character(),
      sample_id_2 = character(),
      sample_node_id_1 = character(),
      sample_node_id_2 = character(),
      shared_node_id = character(),
      shared_node_type = character(),
      edge_type = character(),
      stringsAsFactors = FALSE
    ))
  }

  split_targets <- split(edges, list(edges$to, edges$edge_type), drop = TRUE)
  for (group in split_targets) {
    sample_group <- sort(unique(group$from))
    if (length(sample_group) < 2L) {
      next
    }

    pairs <- utils::combn(sample_group, 2L, simplify = FALSE)
    for (pair in pairs) {
      rows[[row_idx]] <- data.frame(
        sample_id_1 = sample_key_map[[pair[[1L]]]],
        sample_id_2 = sample_key_map[[pair[[2L]]]],
        sample_node_id_1 = pair[[1L]],
        sample_node_id_2 = pair[[2L]],
        shared_node_id = group$to[[1L]],
        shared_node_type = node_type_map[[group$to[[1L]]]],
        edge_type = group$edge_type[[1L]],
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1L
    }
  }

  if (length(rows) == 0L) {
    return(data.frame(
      sample_id_1 = character(),
      sample_id_2 = character(),
      sample_node_id_1 = character(),
      sample_node_id_2 = character(),
      shared_node_id = character(),
      shared_node_type = character(),
      edge_type = character(),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, rows)
  row.names(out) <- NULL
  out
}

.depgraph_project_sample_dependencies <- function(graph, via, edge_types = NULL) {
  shared <- .depgraph_shared_dependency_table(graph, via = via, edge_types = edge_types)
  if (nrow(shared) == 0L) {
    return(shared)
  }

  projection <- unique(shared[, c("sample_node_id_1", "sample_node_id_2"), drop = FALSE])
  projection$projection_edge_id <- paste0("projection_", seq_len(nrow(projection)))
  projection
}

.depgraph_shared_dependency_edge_ids <- function(graph, shared_table) {
  if (nrow(shared_table) == 0L) {
    return(character())
  }

  edge_data <- graph$edges$data
  edge_ids <- character()
  for (i in seq_len(nrow(shared_table))) {
    relevant_edges <- edge_data[
      edge_data$edge_type == shared_table$edge_type[[i]] &
        edge_data$to == shared_table$shared_node_id[[i]] &
        edge_data$from %in% c(shared_table$sample_node_id_1[[i]], shared_table$sample_node_id_2[[i]]),
      "edge_id",
      drop = TRUE
    ]
    edge_ids <- c(edge_ids, relevant_edges)
  }

  unique(edge_ids)
}

query_node_type <- function(graph, node_types, ids = NULL) {
  .depgraph_assert(inherits(graph, "dependency_graph"), "`graph` must be a `dependency_graph`.")
  node_types <- unique(vapply(node_types, .depgraph_match_node_type, character(1), USE.NAMES = FALSE))

  node_data <- graph$nodes$data
  table <- node_data[node_data$node_type %in% node_types, , drop = FALSE]
  if (!is.null(ids)) {
    ids <- unique(as.character(ids))
    table <- table[table$node_id %in% ids, , drop = FALSE]
  }

  .depgraph_query_result(
    graph = graph,
    query = "query_node_type",
    params = list(node_types = node_types, ids = ids),
    node_ids = table$node_id,
    table = table
  )
}

query_edge_type <- function(graph, edge_types, node_ids = NULL) {
  .depgraph_assert(inherits(graph, "dependency_graph"), "`graph` must be a `dependency_graph`.")
  edge_types <- unique(as.character(edge_types))
  if (!is.null(node_ids)) {
    node_ids <- .depgraph_resolve_node_ids(graph, node_ids)
  }

  table <- .depgraph_subset_edges(graph, edge_types = edge_types, node_ids = node_ids)
  touched_nodes <- unique(c(table$from, table$to))

  .depgraph_query_result(
    graph = graph,
    query = "query_edge_type",
    params = list(edge_types = edge_types, node_ids = node_ids),
    node_ids = touched_nodes,
    edge_ids = table$edge_id,
    table = table
  )
}

query_neighbors <- function(graph, node_ids, edge_types = NULL, node_types = NULL, direction = c("out", "in", "all")) {
  .depgraph_assert(inherits(graph, "dependency_graph"), "`graph` must be a `dependency_graph`.")
  direction <- match.arg(direction)
  seed_ids <- .depgraph_resolve_node_ids(graph, node_ids)
  g <- .depgraph_filter_graph_by_edge_type(graph, edge_types)

  rows <- list()
  row_idx <- 1L
  used_nodes <- character()
  used_edges <- character()
  node_data <- graph$nodes$data
  edge_data <- .depgraph_subset_edges(graph, edge_types = edge_types)

  for (seed in seed_ids) {
    neigh <- igraph::neighbors(g, v = seed, mode = direction)
    neigh_ids <- igraph::as_ids(neigh)
    if (!is.null(node_types) && length(neigh_ids) > 0L) {
      node_types <- unique(vapply(node_types, .depgraph_match_node_type, character(1), USE.NAMES = FALSE))
      keep <- node_data$node_type[match(neigh_ids, node_data$node_id)] %in% node_types
      neigh_ids <- neigh_ids[keep]
    }

    for (neighbor in neigh_ids) {
      edge_rows <- if (identical(direction, "out")) {
        edge_data[edge_data$from == seed & edge_data$to == neighbor, , drop = FALSE]
      } else if (identical(direction, "in")) {
        edge_data[edge_data$from == neighbor & edge_data$to == seed, , drop = FALSE]
      } else {
        .depgraph_edge_between(edge_data, seed, neighbor, undirected = TRUE)
      }

      if (nrow(edge_rows) == 0L) {
        next
      }

      neighbor_row <- node_data[node_data$node_id == neighbor, , drop = FALSE]
      for (j in seq_len(nrow(edge_rows))) {
        rows[[row_idx]] <- data.frame(
          seed_node_id = seed,
          node_id = neighbor,
          node_type = neighbor_row$node_type[[1L]],
          edge_id = edge_rows$edge_id[[j]],
          edge_type = edge_rows$edge_type[[j]],
          direction = direction,
          stringsAsFactors = FALSE
        )
        used_nodes <- c(used_nodes, seed, neighbor)
        used_edges <- c(used_edges, edge_rows$edge_id[[j]])
        row_idx <- row_idx + 1L
      }
    }
  }

  table <- if (length(rows) == 0L) {
    data.frame(
      seed_node_id = character(),
      node_id = character(),
      node_type = character(),
      edge_id = character(),
      edge_type = character(),
      direction = character(),
      stringsAsFactors = FALSE
    )
  } else {
    out <- do.call(rbind, rows)
    row.names(out) <- NULL
    out
  }

  .depgraph_query_result(
    graph = graph,
    query = "query_neighbors",
    params = list(node_ids = seed_ids, edge_types = edge_types, node_types = node_types, direction = direction),
    node_ids = unique(used_nodes),
    edge_ids = unique(used_edges),
    table = table
  )
}

query_paths <- function(graph, from, to, edge_types = NULL, node_types = NULL, mode = c("out", "in", "all"), max_length = NULL) {
  .depgraph_assert(inherits(graph, "dependency_graph"), "`graph` must be a `dependency_graph`.")
  mode <- match.arg(mode)
  from_ids <- .depgraph_resolve_node_ids(graph, from)
  to_ids <- .depgraph_resolve_node_ids(graph, to)
  filtered_graph <- .depgraph_filter_graph_by_edge_type(graph, edge_types)
  edge_data <- .depgraph_subset_edges(graph, edge_types = edge_types)
  path_nodes <- list()
  path_idx <- 1L

  for (from_id in from_ids) {
    for (to_id in to_ids) {
      paths <- igraph::all_simple_paths(
        filtered_graph,
        from = from_id,
        to = to_id,
        mode = mode,
        cutoff = if (is.null(max_length)) -1 else max_length
      )
      if (length(paths) == 0L) {
        next
      }

      for (path in paths) {
        node_ids_in_path <- igraph::as_ids(path)
        if (!is.null(node_types)) {
          allowed_types <- unique(vapply(node_types, .depgraph_match_node_type, character(1), USE.NAMES = FALSE))
          path_types <- graph$nodes$data$node_type[match(node_ids_in_path, graph$nodes$data$node_id)]
          if (!all(path_types %in% allowed_types)) {
            next
          }
        }
        path_nodes[[path_idx]] <- node_ids_in_path
        path_idx <- path_idx + 1L
      }
    }
  }

  .depgraph_format_paths(
    graph = graph,
    path_node_ids = path_nodes,
    edge_data = edge_data,
    query_name = "query_paths",
    params = list(from = from_ids, to = to_ids, edge_types = edge_types, node_types = node_types, mode = mode, max_length = max_length)
  )
}

query_shortest_paths <- function(graph, from, to, edge_types = NULL, node_types = NULL, mode = c("out", "in", "all")) {
  .depgraph_assert(inherits(graph, "dependency_graph"), "`graph` must be a `dependency_graph`.")
  mode <- match.arg(mode)
  from_ids <- .depgraph_resolve_node_ids(graph, from)
  to_ids <- .depgraph_resolve_node_ids(graph, to)
  filtered_graph <- .depgraph_filter_graph_by_edge_type(graph, edge_types)
  edge_data <- .depgraph_subset_edges(graph, edge_types = edge_types)
  path_nodes <- list()
  path_idx <- 1L

  for (from_id in from_ids) {
    for (to_id in to_ids) {
      shortest <- igraph::shortest_paths(filtered_graph, from = from_id, to = to_id, mode = mode)$vpath
      if (length(shortest) == 0L || length(shortest[[1L]]) == 0L) {
        next
      }
      node_ids_in_path <- igraph::as_ids(shortest[[1L]])
      if (!is.null(node_types)) {
        allowed_types <- unique(vapply(node_types, .depgraph_match_node_type, character(1), USE.NAMES = FALSE))
        path_types <- graph$nodes$data$node_type[match(node_ids_in_path, graph$nodes$data$node_id)]
        if (!all(path_types %in% allowed_types)) {
          next
        }
      }
      path_nodes[[path_idx]] <- node_ids_in_path
      path_idx <- path_idx + 1L
    }
  }

  .depgraph_format_paths(
    graph = graph,
    path_node_ids = path_nodes,
    edge_data = edge_data,
    query_name = "query_shortest_paths",
    params = list(from = from_ids, to = to_ids, edge_types = edge_types, node_types = node_types, mode = mode)
  )
}

detect_dependency_components <- function(graph, via = c("Subject", "Batch", "Study", "Timepoint", "Assay", "FeatureSet", "Outcome"), edge_types = NULL, min_size = 1) {
  .depgraph_assert(inherits(graph, "dependency_graph"), "`graph` must be a `dependency_graph`.")
  shared <- .depgraph_shared_dependency_table(graph, via = via, edge_types = edge_types)
  projection <- .depgraph_project_sample_dependencies(graph, via = via, edge_types = edge_types)
  sample_nodes <- graph$nodes$data[graph$nodes$data$node_type == "Sample", , drop = FALSE]

  sample_graph <- if (nrow(projection) == 0L) {
    igraph::make_empty_graph(n = nrow(sample_nodes), directed = FALSE)
  } else {
    g <- igraph::graph_from_data_frame(
      d = data.frame(
        from = projection$sample_node_id_1,
        to = projection$sample_node_id_2,
        stringsAsFactors = FALSE
      ),
      vertices = data.frame(name = sample_nodes$node_id, stringsAsFactors = FALSE),
      directed = FALSE
    )
    g
  }
  igraph::V(sample_graph)$name <- sample_nodes$node_id

  components <- igraph::components(sample_graph)
  table <- data.frame(
    sample_id = sample_nodes$node_key,
    sample_node_id = sample_nodes$node_id,
    component_id = paste0("component_", components$membership[match(sample_nodes$node_id, names(components$membership))]),
    component_size = components$csize[components$membership[match(sample_nodes$node_id, names(components$membership))]],
    stringsAsFactors = FALSE
  )

  table <- table[table$component_size >= min_size, , drop = FALSE]
  keep_nodes <- table$sample_node_id
  projection <- projection[
    projection$sample_node_id_1 %in% keep_nodes &
      projection$sample_node_id_2 %in% keep_nodes,
    ,
    drop = FALSE
  ]
  shared <- shared[
    shared$sample_node_id_1 %in% keep_nodes &
      shared$sample_node_id_2 %in% keep_nodes,
    ,
    drop = FALSE
  ]
  edge_ids <- .depgraph_shared_dependency_edge_ids(graph, shared)

  .depgraph_query_result(
    graph = graph,
    query = "detect_dependency_components",
    params = list(via = via, edge_types = edge_types, min_size = min_size),
    node_ids = keep_nodes,
    table = table,
    metadata = list(
      n_components = length(unique(table$component_id)),
      projection_edges = projection
    ),
    edge_ids = edge_ids
  )
}

detect_shared_dependencies <- function(graph, via = c("Subject", "Batch", "Study", "Timepoint"), samples = NULL) {
  .depgraph_assert(inherits(graph, "dependency_graph"), "`graph` must be a `dependency_graph`.")
  table <- .depgraph_shared_dependency_table(graph, via = via, samples = samples)

  node_ids <- unique(c(
    table$sample_node_id_1,
    table$sample_node_id_2,
    table$shared_node_id
  ))

  edge_ids <- .depgraph_shared_dependency_edge_ids(graph, table)

  .depgraph_query_result(
    graph = graph,
    query = "detect_shared_dependencies",
    params = list(via = via, samples = samples),
    node_ids = node_ids,
    edge_ids = unique(edge_ids),
    table = table
  )
}
