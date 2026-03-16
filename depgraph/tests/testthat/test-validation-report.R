test_that("validation report methods summarize issues cleanly", {
  meta <- data.frame(
    sample_id = c("S1", "S2"),
    subject_id = c("P1", "P1"),
    stringsAsFactors = FALSE
  )

  samples <- create_nodes(meta, type = "Sample", id_col = "sample_id")
  subjects <- create_nodes(meta, type = "Subject", id_col = "subject_id")
  edges <- create_edges(meta, "sample_id", "subject_id", "Sample", "Subject", "sample_belongs_to_subject")
  graph <- build_dependency_graph(list(samples, subjects), list(edges))
  validation <- validate_graph(graph)

  validation_summary <- summary(validation)

  expect_true(length(capture.output(print(validation))) > 0L)
  expect_true(is.list(validation_summary))
  expect_true(all(c("by_level", "by_severity", "by_code") %in% names(validation_summary)))
  expect_s3_class(as.data.frame(validation), "data.frame")
})
