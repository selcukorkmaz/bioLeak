make_bioleak_graph <- function(time_source = c("index", "precedence")) {
  time_source <- match.arg(time_source)

  meta <- data.frame(
    sample_id = c("S1", "S2", "S3", "S4", "S5"),
    subject_id = c("P1", "P1", "P2", "P3", "P4"),
    batch_id = c("B1", "B2", "B1", "B3", NA),
    study_id = c("ST1", "ST1", "ST2", "ST2", "ST3"),
    timepoint_id = c("T0", "T1", "T0", "T2", NA),
    stringsAsFactors = FALSE
  )

  sample_nodes <- create_nodes(meta, type = "Sample", id_col = "sample_id")
  subject_nodes <- create_nodes(meta, type = "Subject", id_col = "subject_id")
  batch_nodes <- create_nodes(meta, type = "Batch", id_col = "batch_id")
  study_nodes <- create_nodes(meta, type = "Study", id_col = "study_id")

  time_data <- unique(meta[!is.na(meta$timepoint_id), "timepoint_id", drop = FALSE])
  time_nodes <- if (identical(time_source, "index")) {
    create_nodes(
      transform(time_data, time_index = c(0L, 1L, 2L)),
      type = "Timepoint",
      id_col = "timepoint_id",
      attr_cols = c("time_index")
    )
  } else {
    create_nodes(time_data, type = "Timepoint", id_col = "timepoint_id")
  }

  edge_sets <- list(
    create_edges(meta, "sample_id", "subject_id", "Sample", "Subject", "sample_belongs_to_subject"),
    create_edges(meta, "sample_id", "batch_id", "Sample", "Batch", "sample_processed_in_batch", allow_missing = TRUE),
    create_edges(meta, "sample_id", "study_id", "Sample", "Study", "sample_from_study"),
    create_edges(meta, "sample_id", "timepoint_id", "Sample", "Timepoint", "sample_collected_at_timepoint", allow_missing = TRUE)
  )

  if (identical(time_source, "precedence")) {
    edge_sets <- c(edge_sets, list(
      create_edges(
        data.frame(
          from_timepoint = c("T0", "T1"),
          to_timepoint = c("T1", "T2"),
          stringsAsFactors = FALSE
        ),
        "from_timepoint", "to_timepoint", "Timepoint", "Timepoint", "timepoint_precedes"
      )
    ))
  }

  build_dependency_graph(
    nodes = list(sample_nodes, subject_nodes, batch_nodes, study_nodes, time_nodes),
    edges = edge_sets,
    graph_name = paste0("bioleak_graph_", time_source),
    dataset_name = "BioLeakDemo"
  )
}

test_that("as_bioleak_split_spec translates constraints into canonical sample data", {
  graph <- make_bioleak_graph("index")
  constraint <- derive_split_constraints(graph, mode = "batch")
  spec <- as_bioleak_split_spec(constraint, graph = graph)

  expect_s3_class(spec, "bioleak_split_spec")
  expect_true(all(c(
    "sample_id", "sample_node_id", "group_id", "primary_group", "batch_group",
    "study_group", "timepoint_id", "time_index", "order_rank"
  ) %in% names(spec$sample_data)))
  expect_equal(spec$group_var, "group_id")
  expect_true("batch_group" %in% spec$block_vars)
  expect_equal(spec$recommended_resampling, "blocked_cv")
  expect_equal(spec$sample_data$batch_group[spec$sample_data$sample_id == "S1"], "B1")
  expect_equal(spec$sample_data$study_group[spec$sample_data$sample_id == "S1"], "ST1")
  expect_equal(spec$sample_data$order_rank[spec$sample_data$sample_id == "S4"], 3L)
  expect_equal(spec$metadata$source_mode, "batch")
  expect_equal(spec$metadata$source_strategy, "batch")
})

test_that("as_bioleak_split_spec preserves time-aware ordering semantics", {
  graph <- make_bioleak_graph("precedence")
  constraint <- derive_split_constraints(graph, mode = "time")
  spec <- as_bioleak_split_spec(constraint, graph = graph)

  expect_s3_class(spec, "bioleak_split_spec")
  expect_identical(spec$time_var, "order_rank")
  expect_true(isTRUE(spec$ordering_required))
  expect_equal(spec$sample_data$order_rank[spec$sample_data$sample_id == "S2"], 2L)
  expect_equal(spec$recommended_resampling, "ordered_split")
})

test_that("validate_bioleak_split_spec detects missing ordering and duplicate samples", {
  graph <- make_bioleak_graph("index")
  constraint <- derive_split_constraints(graph, mode = "time")
  spec <- as_bioleak_split_spec(constraint, graph = graph)
  broken <- spec
  broken$sample_data$order_rank[] <- NA_integer_
  broken$sample_data <- rbind(broken$sample_data, broken$sample_data[1, , drop = FALSE])

  validation <- validate_bioleak_split_spec(broken)

  expect_s3_class(validation, "bioleak_split_spec_validation")
  expect_false(validation$valid)
  expect_true(any(validation$issues$code == "duplicate_sample_id"))
  expect_true(any(validation$issues$code == "missing_required_ordering"))
  expect_true(length(capture.output(print(validation))) > 0L)
})

test_that("summarize_leakage_risks combines graph, constraint, and split-spec diagnostics", {
  graph <- make_bioleak_graph("index")
  constraint <- derive_split_constraints(graph, mode = "composite", strategy = "strict", via = c("Subject", "Batch"))
  spec <- as_bioleak_split_spec(constraint, graph = graph)
  summary_obj <- summarize_leakage_risks(graph, constraint = constraint, split_spec = spec)

  expect_s3_class(summary_obj, "leakage_risk_summary")
  expect_true(is.character(summary_obj$overview))
  expect_true(all(c("severity", "category", "message", "source", "n_affected") %in% names(summary_obj$diagnostics)))
  expect_true(any(summary_obj$diagnostics$source == "validation"))
  expect_true(any(summary_obj$diagnostics$source == "constraint"))
  expect_true(any(summary_obj$diagnostics$source == "split_spec"))
  expect_true(any(summary_obj$diagnostics$category == "blocking_available"))
  expect_true(any(summary_obj$diagnostics$category == "ordering_available"))
  expect_true("block_vars" %in% names(summary_obj$split_spec_summary))
  expect_true("ordering_required" %in% names(summary_obj$split_spec_summary))
  expect_true(length(capture.output(print(summary_obj))) > 0L)
})
