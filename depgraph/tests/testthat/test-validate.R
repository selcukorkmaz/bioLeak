test_that("validate_graph catches unsupported edge types", {
  nodes <- graph_node_set(
    data.frame(
      node_id = c("sample:S1", "subject:P1"),
      node_type = c("Sample", "Subject"),
      node_key = c("S1", "P1"),
      label = c("S1", "P1"),
      attrs = I(list(list(), list())),
      stringsAsFactors = FALSE
    )
  )

  edges <- graph_edge_set(
    data.frame(
      edge_id = "mystery:1",
      from = "sample:S1",
      to = "subject:P1",
      edge_type = "mystery_relation",
      attrs = I(list(list())),
      stringsAsFactors = FALSE
    )
  )

  graph <- dependency_graph(nodes = nodes, edges = edges, graph = NULL)
  validation <- validate_graph(graph)

  expect_false(validation$valid)
  expect_match(validation$errors[[1]], "Unknown edge type", fixed = TRUE)
})

test_that("validate_graph returns a structured validation report", {
  meta <- data.frame(
    sample_id = c("S1", "S2"),
    subject_id = c("P1", "P2"),
    stringsAsFactors = FALSE
  )

  samples <- create_nodes(meta, type = "Sample", id_col = "sample_id")
  subjects <- create_nodes(meta, type = "Subject", id_col = "subject_id")
  edges <- create_edges(meta, "sample_id", "subject_id", "Sample", "Subject", "sample_belongs_to_subject")
  graph <- build_dependency_graph(list(samples, subjects), list(edges))

  validation <- validate_graph(graph)

  expect_s3_class(validation, "depgraph_validation_report")
  expect_true(is.data.frame(validation$issues))
  expect_true(all(c("issue_id", "level", "severity", "code", "message", "node_ids", "edge_ids", "details") %in% names(validation$issues)))
  expect_true(is.logical(validation$valid))
  expect_true(is.list(validation$summary))
  expect_true(is.list(validation$metrics))
  expect_equal(as.data.frame(validation), validation$issues)
})

test_that("validate_graph enforces structural single-target cardinality", {
  nodes <- graph_node_set(
    data.frame(
      node_id = c("sample:S1", "assay:A1", "assay:A2"),
      node_type = c("Sample", "Assay", "Assay"),
      node_key = c("S1", "A1", "A2"),
      label = c("S1", "A1", "A2"),
      attrs = I(list(list(), list(), list())),
      stringsAsFactors = FALSE
    )
  )

  edges <- graph_edge_set(
    data.frame(
      edge_id = c("sample_measured_by_assay:1", "sample_measured_by_assay:2"),
      from = c("sample:S1", "sample:S1"),
      to = c("assay:A1", "assay:A2"),
      edge_type = c("sample_measured_by_assay", "sample_measured_by_assay"),
      attrs = I(list(list(), list())),
      stringsAsFactors = FALSE
    )
  )

  graph <- dependency_graph(nodes = nodes, edges = edges, graph = NULL)
  validation <- validate_graph(graph)

  expect_false(validation$valid)
  expect_true(any(validation$issues$code == "single_target_violation"))
})

test_that("validate_graph warns on non-canonical node attributes", {
  nodes <- graph_node_set(
    data.frame(
      node_id = "sample:S1",
      node_type = "Sample",
      node_key = "S1",
      label = "S1",
      attrs = I(list(list(custom_field = "value"))),
      stringsAsFactors = FALSE
    )
  )

  edges <- graph_edge_set(
    data.frame(
      edge_id = character(),
      from = character(),
      to = character(),
      edge_type = character(),
      attrs = I(list()),
      stringsAsFactors = FALSE
    )
  )

  graph <- dependency_graph(nodes = nodes, edges = edges, graph = NULL)
  validation <- validate_graph(graph)

  expect_true(validation$valid)
  expect_true(any(grepl("non-canonical attributes", validation$warnings, fixed = TRUE)))
})

test_that("validate_depgraph is a compatibility alias", {
  meta <- data.frame(
    sample_id = c("S1"),
    subject_id = c("P1"),
    stringsAsFactors = FALSE
  )

  samples <- create_nodes(meta, type = "Sample", id_col = "sample_id")
  subjects <- create_nodes(meta, type = "Subject", id_col = "subject_id")
  edges <- create_edges(meta, "sample_id", "subject_id", "Sample", "Subject", "sample_belongs_to_subject")
  graph <- build_dependency_graph(list(samples, subjects), list(edges))

  validation <- validate_depgraph(graph)

  expect_s3_class(validation, "depgraph_validation_report")
  expect_identical(validation$valid, validate_graph(graph)$valid)
})

test_that("semantic validation flags invalid FeatureSet derivation scope", {
  nodes <- graph_node_set(
    data.frame(
      node_id = c("sample:S1", "featureset:FS1"),
      node_type = c("Sample", "FeatureSet"),
      node_key = c("S1", "FS1"),
      label = c("S1", "FS1"),
      attrs = I(list(list(), list(derivation_scope = "invalid_scope"))),
      stringsAsFactors = FALSE
    )
  )

  edges <- graph_edge_set(
    data.frame(
      edge_id = "sample_uses_featureset:1",
      from = "sample:S1",
      to = "featureset:FS1",
      edge_type = "sample_uses_featureset",
      attrs = I(list(list())),
      stringsAsFactors = FALSE
    )
  )

  graph <- dependency_graph(nodes = nodes, edges = edges, graph = NULL)
  validation <- validate_graph(graph)

  expect_true(any(validation$issues$code == "invalid_featureset_derivation_scope"))
  expect_true(any(validation$issues$severity == "warning"))
})

test_that("semantic validation flags missing study assignment when study structure exists", {
  nodes <- graph_node_set(
    data.frame(
      node_id = c("sample:S1", "sample:S2", "study:ST1"),
      node_type = c("Sample", "Sample", "Study"),
      node_key = c("S1", "S2", "ST1"),
      label = c("S1", "S2", "ST1"),
      attrs = I(list(list(), list(), list())),
      stringsAsFactors = FALSE
    )
  )

  edges <- graph_edge_set(
    data.frame(
      edge_id = "sample_from_study:1",
      from = "sample:S1",
      to = "study:ST1",
      edge_type = "sample_from_study",
      attrs = I(list(list())),
      stringsAsFactors = FALSE
    )
  )

  graph <- dependency_graph(nodes = nodes, edges = edges, graph = NULL)
  validation <- validate_graph(graph)

  expect_true(any(validation$issues$code == "sample_missing_study_assignment"))
  expect_true(any(validation$issues$severity == "warning"))
})

test_that("semantic validation flags multi-subject samples unless explicitly allowed", {
  nodes <- graph_node_set(
    data.frame(
      node_id = c("sample:S1", "subject:P1", "subject:P2"),
      node_type = c("Sample", "Subject", "Subject"),
      node_key = c("S1", "P1", "P2"),
      label = c("S1", "P1", "P2"),
      attrs = I(list(list(), list(), list())),
      stringsAsFactors = FALSE
    )
  )

  edges <- graph_edge_set(
    data.frame(
      edge_id = c("sample_belongs_to_subject:1", "sample_belongs_to_subject:2"),
      from = c("sample:S1", "sample:S1"),
      to = c("subject:P1", "subject:P2"),
      edge_type = c("sample_belongs_to_subject", "sample_belongs_to_subject"),
      attrs = I(list(list(), list())),
      stringsAsFactors = FALSE
    )
  )

  graph <- dependency_graph(nodes = nodes, edges = edges, graph = NULL)
  validation <- validate_graph(graph)

  expect_false(validation$valid)
  expect_true(any(validation$issues$code == "sample_multiple_subject_assignments"))

  allowed_graph <- dependency_graph(
    nodes = nodes,
    edges = edges,
    graph = NULL,
    metadata = list(validation_overrides = list(allow_multi_subject_samples = TRUE))
  )
  allowed_validation <- validate_graph(allowed_graph)

  expect_false(any(allowed_validation$issues$code == "sample_multiple_subject_assignments"))
})

test_that("leakage validation reports repeated-subject and per-dataset feature risks", {
  nodes <- graph_node_set(
    data.frame(
      node_id = c("sample:S1", "sample:S2", "subject:P1", "featureset:FS1"),
      node_type = c("Sample", "Sample", "Subject", "FeatureSet"),
      node_key = c("S1", "S2", "P1", "FS1"),
      label = c("S1", "S2", "P1", "FS1"),
      attrs = I(list(
        list(),
        list(),
        list(),
        list(derivation_scope = "per_dataset")
      )),
      stringsAsFactors = FALSE
    )
  )

  edges <- graph_edge_set(
    data.frame(
      edge_id = c(
        "sample_belongs_to_subject:1",
        "sample_belongs_to_subject:2",
        "sample_uses_featureset:1",
        "sample_uses_featureset:2"
      ),
      from = c("sample:S1", "sample:S2", "sample:S1", "sample:S2"),
      to = c("subject:P1", "subject:P1", "featureset:FS1", "featureset:FS1"),
      edge_type = c(
        "sample_belongs_to_subject",
        "sample_belongs_to_subject",
        "sample_uses_featureset",
        "sample_uses_featureset"
      ),
      attrs = I(list(list(), list(), list(), list())),
      stringsAsFactors = FALSE
    )
  )

  graph <- dependency_graph(nodes = nodes, edges = edges, graph = NULL)
  validation <- validate_graph(graph)

  expect_true(any(validation$issues$code == "repeated_subject_samples"))
  expect_true(any(validation$issues$code == "per_dataset_featureset"))
  expect_true(any(validation$issues$severity == "advisory"))
})
