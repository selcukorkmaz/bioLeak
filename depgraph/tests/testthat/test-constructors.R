test_that("create_nodes and create_edges build valid core objects", {
  meta <- data.frame(
    sample_id = c("S1", "S2"),
    subject_id = c("P1", "P2"),
    batch_id = c("B1", "B1"),
    study_id = c("ST1", "ST1"),
    stringsAsFactors = FALSE
  )

  samples <- create_nodes(meta, type = "Sample", id_col = "sample_id")
  subjects <- create_nodes(meta, type = "Subject", id_col = "subject_id")
  edges <- create_edges(meta, "sample_id", "subject_id", "Sample", "Subject", "sample_belongs_to_subject")

  graph <- build_dependency_graph(list(samples, subjects), list(edges))
  validation <- validate_graph(graph)

  expect_s3_class(samples, "graph_node_set")
  expect_s3_class(edges, "graph_edge_set")
  expect_s3_class(graph, "dependency_graph")
  expect_true(validation$valid)
  expect_named(samples$data$attrs[[1]], c("subject_id", "batch_id", "study_id"))
  expect_equal(igraph::vcount(as_igraph(graph)), nrow(graph$nodes$data))
  expect_equal(igraph::ecount(as_igraph(graph)), nrow(graph$edges$data))
})

test_that("compatibility aliases delegate to the primary constructors", {
  meta <- data.frame(
    sample_id = c("S1", "S2"),
    subject_id = c("P1", "P2"),
    stringsAsFactors = FALSE
  )

  samples <- create_nodes(meta, type = "Sample", id_col = "sample_id")
  subjects <- new_depgraph_nodes(create_nodes(meta, type = "Subject", id_col = "subject_id")$data)
  edges <- new_depgraph_edges(
    create_edges(meta, "sample_id", "subject_id", "Sample", "Subject", "sample_belongs_to_subject")$data
  )
  full_nodes <- new_depgraph_nodes(rbind(samples$data, subjects$data))

  graph <- new_depgraph(nodes = full_nodes, edges = edges)
  graph2 <- build_depgraph(nodes = list(samples, subjects), edges = list(edges))

  expect_s3_class(subjects, "graph_node_set")
  expect_s3_class(edges, "graph_edge_set")
  expect_s3_class(graph, "dependency_graph")
  expect_s3_class(graph2, "dependency_graph")
  expect_equal(igraph::vcount(as_igraph(graph2)), 4)
})

test_that("build_dependency_graph stores validation overrides in metadata", {
  meta <- data.frame(
    sample_id = c("S1"),
    subject_id = c("P1"),
    stringsAsFactors = FALSE
  )

  samples <- create_nodes(meta, type = "Sample", id_col = "sample_id")
  subjects <- create_nodes(meta, type = "Subject", id_col = "subject_id")
  edges <- create_edges(meta, "sample_id", "subject_id", "Sample", "Subject", "sample_belongs_to_subject")

  graph <- build_dependency_graph(
    list(samples, subjects),
    list(edges),
    validation_overrides = list(allow_multi_subject_samples = TRUE)
  )

  expect_true(isTRUE(graph$metadata$validation_overrides$allow_multi_subject_samples))
})

test_that("low-level constructors reject malformed attrs", {
  expect_error(
    graph_node_set(
      data.frame(
        node_id = "sample:S1",
        node_type = "Sample",
        node_key = "S1",
        label = "S1",
        attrs = I(list(list("bad"))),
        stringsAsFactors = FALSE
      )
    ),
    "named lists"
  )

  expect_error(
    graph_edge_set(
      data.frame(
        edge_id = "sample_belongs_to_subject:1",
        from = "sample:S1",
        to = "subject:P1",
        edge_type = "sample_belongs_to_subject",
        attrs = I(list(list("bad"))),
        stringsAsFactors = FALSE
      )
    ),
    "named lists"
  )
})
