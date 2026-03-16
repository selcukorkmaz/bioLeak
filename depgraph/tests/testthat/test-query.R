make_query_graph <- function() {
  meta <- data.frame(
    sample_id = c("S1", "S2", "S3", "S4"),
    subject_id = c("P1", "P1", "P2", "P3"),
    batch_id = c("B1", "B2", "B1", "B3"),
    study_id = c("ST1", "ST1", "ST2", "ST2"),
    timepoint_id = c("T0", "T1", "T0", "T0"),
    outcome_id = c("O1", "O1", "O2", "O3"),
    stringsAsFactors = FALSE
  )

  sample_nodes <- create_nodes(meta, type = "Sample", id_col = "sample_id")
  subject_nodes <- create_nodes(meta, type = "Subject", id_col = "subject_id")
  batch_nodes <- create_nodes(meta, type = "Batch", id_col = "batch_id")
  study_nodes <- create_nodes(meta, type = "Study", id_col = "study_id")
  time_nodes <- create_nodes(
    transform(meta, time_index = c(0L, 1L, 0L, 0L)),
    type = "Timepoint",
    id_col = "timepoint_id",
    attr_cols = c("time_index")
  )
  outcome_nodes <- create_nodes(
    data.frame(
      outcome_id = c("O1", "O2", "O3"),
      observation_level = c("subject", "subject", "subject"),
      stringsAsFactors = FALSE
    ),
    type = "Outcome",
    id_col = "outcome_id",
    attr_cols = c("observation_level")
  )

  subject_edges <- create_edges(meta, "sample_id", "subject_id", "Sample", "Subject", "sample_belongs_to_subject")
  batch_edges <- create_edges(meta, "sample_id", "batch_id", "Sample", "Batch", "sample_processed_in_batch")
  study_edges <- create_edges(meta, "sample_id", "study_id", "Sample", "Study", "sample_from_study")
  time_edges <- create_edges(meta, "sample_id", "timepoint_id", "Sample", "Timepoint", "sample_collected_at_timepoint")
  subject_outcome_edges <- create_edges(
    data.frame(
      subject_id = c("P1", "P2", "P3"),
      outcome_id = c("O1", "O2", "O3"),
      stringsAsFactors = FALSE
    ),
    "subject_id", "outcome_id", "Subject", "Outcome", "subject_has_outcome"
  )

  build_dependency_graph(
    nodes = list(sample_nodes, subject_nodes, batch_nodes, study_nodes, time_nodes, outcome_nodes),
    edges = list(subject_edges, batch_edges, study_edges, time_edges, subject_outcome_edges),
    graph_name = "query_graph",
    dataset_name = "QueryDemo"
  )
}

test_that("query_node_type returns graph_query_result with tidy node rows", {
  graph <- make_query_graph()
  result <- query_node_type(graph, node_types = "Subject")

  expect_s3_class(result, "graph_query_result")
  expect_true(all(result$table$node_type == "Subject"))
  expect_equal(nrow(result$table), 3L)
  expect_identical(as.data.frame(result), result$table)
})

test_that("query_edge_type returns graph_query_result with tidy edge rows", {
  graph <- make_query_graph()
  result <- query_edge_type(graph, edge_types = "sample_from_study")

  expect_s3_class(result, "graph_query_result")
  expect_true(all(result$table$edge_type == "sample_from_study"))
  expect_equal(nrow(result$table), 4L)
})

test_that("query_neighbors respects direction and edge type filtering", {
  graph <- make_query_graph()
  out_result <- query_neighbors(
    graph,
    node_ids = "sample:S1",
    edge_types = "sample_belongs_to_subject",
    direction = "out"
  )
  in_result <- query_neighbors(
    graph,
    node_ids = "subject:P1",
    edge_types = "sample_belongs_to_subject",
    direction = "in"
  )

  expect_s3_class(out_result, "graph_query_result")
  expect_equal(out_result$table$node_id, "subject:P1")
  expect_true(all(in_result$table$seed_node_id == "subject:P1"))
  expect_setequal(in_result$table$node_id, c("sample:S1", "sample:S2"))
})

test_that("query_paths returns stepwise path rows", {
  graph <- make_query_graph()
  result <- query_paths(
    graph,
    from = "sample:S1",
    to = "outcome:O1",
    edge_types = c("sample_belongs_to_subject", "subject_has_outcome"),
    max_length = 2
  )

  expect_s3_class(result, "graph_query_result")
  expect_true(all(c("path_id", "step", "node_id", "node_type", "edge_id", "edge_type") %in% names(result$table)))
  expect_true(any(result$table$node_id == "subject:P1"))
  expect_true(any(result$table$node_id == "outcome:O1"))
})

test_that("query_shortest_paths returns shortest path rows", {
  graph <- make_query_graph()
  result <- query_shortest_paths(
    graph,
    from = "sample:S1",
    to = "outcome:O1",
    edge_types = c("sample_belongs_to_subject", "subject_has_outcome")
  )

  expect_s3_class(result, "graph_query_result")
  expect_equal(length(unique(result$table$path_id)), 1L)
  expect_true(any(result$table$node_id == "subject:P1"))
})

test_that("detect_dependency_components uses projected sample dependencies", {
  graph <- make_query_graph()
  result <- detect_dependency_components(graph, via = c("Subject", "Batch"))

  expect_s3_class(result, "graph_query_result")
  expect_true(all(c("sample_id", "sample_node_id", "component_id", "component_size") %in% names(result$table)))

  comp_sizes <- stats::setNames(result$table$component_size, result$table$sample_id)
  expect_equal(comp_sizes[["S1"]], 3L)
  expect_equal(comp_sizes[["S2"]], 3L)
  expect_equal(comp_sizes[["S3"]], 3L)
  expect_equal(comp_sizes[["S4"]], 1L)
  expect_true(all(result$edges$edge_id %in% graph$edges$data$edge_id))
  expect_true(is.data.frame(result$metadata$projection_edges))
  expect_true(all(c("sample_node_id_1", "sample_node_id_2", "projection_edge_id") %in% names(result$metadata$projection_edges)))
})

test_that("detect_shared_dependencies finds direct shared subject and batch links", {
  graph <- make_query_graph()
  subject_result <- detect_shared_dependencies(graph, via = "Subject")
  batch_result <- detect_shared_dependencies(graph, via = "Batch")

  expect_s3_class(subject_result, "graph_query_result")
  expect_true(all(c(
    "sample_id_1", "sample_id_2", "sample_node_id_1", "sample_node_id_2",
    "shared_node_id", "shared_node_type", "edge_type"
  ) %in% names(subject_result$table)))

  expect_true(any(
    subject_result$table$sample_id_1 == "S1" &
      subject_result$table$sample_id_2 == "S2" &
      subject_result$table$shared_node_id == "subject:P1"
  ))

  expect_true(any(
    batch_result$table$sample_id_1 == "S1" &
      batch_result$table$sample_id_2 == "S3" &
      batch_result$table$shared_node_id == "batch:B1"
  ))
})
