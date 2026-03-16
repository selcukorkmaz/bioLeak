test_that("summary and print methods expose stable metadata", {
  meta <- data.frame(
    sample_id = c("S1", "S2"),
    subject_id = c("P1", "P2"),
    stringsAsFactors = FALSE
  )

  samples <- create_nodes(meta, type = "Sample", id_col = "sample_id")
  subjects <- create_nodes(meta, type = "Subject", id_col = "subject_id")
  edges <- create_edges(meta, "sample_id", "subject_id", "Sample", "Subject", "sample_belongs_to_subject")
  graph <- build_dependency_graph(list(samples, subjects), list(edges), graph_name = "demo")

  sample_summary <- summary(samples)
  edge_summary <- summary(edges)
  graph_summary <- summary(graph)

  expect_equal(sample_summary$schema_version, "0.1.0")
  expect_equal(edge_summary$schema_version, "0.1.0")
  expect_equal(graph_summary$schema_version, "0.1.0")
  expect_equal(graph_summary$graph_name, "demo")
  expect_true(length(capture.output(print(graph))) > 0L)
  expect_true(length(capture.output(print(samples))) > 0L)
  expect_true(length(capture.output(print(edges))) > 0L)
  expect_s3_class(as.data.frame(samples), "data.frame")
  expect_s3_class(as.data.frame(edges), "data.frame")
})
