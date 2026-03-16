make_constraint_graph <- function(time_source = c("index", "precedence", "missing")) {
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

  subject_edges <- create_edges(
    meta, "sample_id", "subject_id", "Sample", "Subject", "sample_belongs_to_subject"
  )
  batch_edges <- create_edges(
    meta, "sample_id", "batch_id", "Sample", "Batch", "sample_processed_in_batch", allow_missing = TRUE
  )
  study_edges <- create_edges(
    meta, "sample_id", "study_id", "Sample", "Study", "sample_from_study"
  )
  time_edges <- create_edges(
    meta, "sample_id", "timepoint_id", "Sample", "Timepoint", "sample_collected_at_timepoint", allow_missing = TRUE
  )

  edge_sets <- list(subject_edges, batch_edges, study_edges, time_edges)
  if (identical(time_source, "precedence")) {
    precedence_edges <- create_edges(
      data.frame(
        from_timepoint = c("T0", "T1"),
        to_timepoint = c("T1", "T2"),
        stringsAsFactors = FALSE
      ),
      "from_timepoint", "to_timepoint", "Timepoint", "Timepoint", "timepoint_precedes"
    )
    edge_sets <- c(edge_sets, list(precedence_edges))
  }

  build_dependency_graph(
    nodes = list(sample_nodes, subject_nodes, batch_nodes, study_nodes, time_nodes),
    edges = edge_sets,
    graph_name = paste0("constraint_graph_", time_source),
    dataset_name = "ConstraintDemo"
  )
}

test_that("subject constraints group samples by subject and expose grouping vectors", {
  graph <- make_constraint_graph("index")
  result <- derive_split_constraints(graph, mode = "subject")

  expect_s3_class(result, "split_constraint")
  expect_true(all(c(
    "sample_id", "sample_node_id", "group_id", "constraint_type", "group_label", "explanation"
  ) %in% names(result$sample_map)))
  expect_true(all(result$sample_map$constraint_type == "subject"))
  expect_equal(result$sample_map$group_id[result$sample_map$sample_id == "S1"], result$sample_map$group_id[result$sample_map$sample_id == "S2"])
  expect_equal(result$recommended_bioleak_args$group_var, "group_id")

  groups <- grouping_vector(result)
  expect_true(is.character(groups))
  expect_identical(unname(groups[["S1"]]), result$sample_map$group_id[result$sample_map$sample_id == "S1"])
})

test_that("batch constraints warn on missing batches and keep assignments deterministic", {
  graph <- make_constraint_graph("index")
  result <- derive_split_constraints(graph, mode = "batch")

  expect_s3_class(result, "split_constraint")
  expect_equal(result$sample_map$group_id[result$sample_map$sample_id == "S1"], result$sample_map$group_id[result$sample_map$sample_id == "S3"])
  expect_false(result$sample_map$group_id[result$sample_map$sample_id == "S2"] == result$sample_map$group_id[result$sample_map$sample_id == "S1"])
  expect_true(any(grepl("missing batch", result$metadata$warnings, ignore.case = TRUE)))
  expect_true(grepl("unlinked", result$sample_map$explanation[result$sample_map$sample_id == "S5"], ignore.case = TRUE))
})

test_that("study constraints respect sample subsets", {
  graph <- make_constraint_graph("index")
  result <- derive_split_constraints(graph, mode = "study", samples = c("S1", "sample:S4"))

  expect_s3_class(result, "split_constraint")
  expect_setequal(result$sample_map$sample_id, c("S1", "S4"))
  expect_equal(nrow(result$sample_map), 2L)
  expect_true(all(result$sample_map$constraint_type == "study"))
})

test_that("time constraints use time_index ordering when available", {
  graph <- make_constraint_graph("index")
  result <- derive_split_constraints(graph, mode = "time")

  expect_s3_class(result, "split_constraint")
  expect_true(all(c("time_index", "timepoint_id", "order_rank") %in% names(result$sample_map)))
  expect_equal(result$sample_map$time_index[result$sample_map$sample_id == "S1"], 0)
  expect_equal(result$sample_map$order_rank[result$sample_map$sample_id == "S1"], 1)
  expect_equal(result$sample_map$order_rank[result$sample_map$sample_id == "S4"], 3)
  expect_true(any(grepl("missing timepoint", result$metadata$warnings, ignore.case = TRUE)))
})

test_that("time constraints can derive ordering from timepoint precedence edges", {
  graph <- make_constraint_graph("precedence")
  result <- derive_split_constraints(graph, mode = "time")

  expect_s3_class(result, "split_constraint")
  expect_equal(result$sample_map$order_rank[result$sample_map$sample_id == "S1"], 1)
  expect_equal(result$sample_map$order_rank[result$sample_map$sample_id == "S2"], 2)
  expect_equal(result$sample_map$order_rank[result$sample_map$sample_id == "S4"], 3)
  expect_true(all(is.na(result$sample_map$time_index[result$sample_map$sample_id %in% c("S1", "S2", "S4")])))
})

test_that("composite strict constraints use transitive dependency closure", {
  graph <- make_constraint_graph("index")
  result <- derive_split_constraints(graph, mode = "composite", strategy = "strict", via = c("Subject", "Batch"))

  expect_s3_class(result, "split_constraint")
  expect_equal(result$sample_map$group_id[result$sample_map$sample_id == "S1"], result$sample_map$group_id[result$sample_map$sample_id == "S2"])
  expect_equal(result$sample_map$group_id[result$sample_map$sample_id == "S1"], result$sample_map$group_id[result$sample_map$sample_id == "S3"])
  expect_false(result$sample_map$group_id[result$sample_map$sample_id == "S4"] == result$sample_map$group_id[result$sample_map$sample_id == "S1"])
  expect_equal(result$metadata$strategy, "strict")
})

test_that("composite rule-based constraints follow priority fallback rules", {
  graph <- make_constraint_graph("index")
  result <- derive_split_constraints(
    graph,
    mode = "composite",
    strategy = "rule_based",
    priority = c("batch", "study", "subject", "time")
  )

  expect_s3_class(result, "split_constraint")
  expect_equal(result$sample_map$group_id[result$sample_map$sample_id == "S1"], result$sample_map$group_id[result$sample_map$sample_id == "S3"])
  expect_false(result$sample_map$group_id[result$sample_map$sample_id == "S2"] == result$sample_map$group_id[result$sample_map$sample_id == "S1"])
  expect_equal(result$sample_map$constraint_type[result$sample_map$sample_id == "S5"], "study")
  expect_equal(result$sample_map$group_label[result$sample_map$sample_id == "S5"], "ST3")
  expect_equal(result$metadata$strategy, "rule_based")
  expect_true(isTRUE(result$recommended_bioleak_args$ordering_required))
  expect_equal(result$recommended_bioleak_args$time_var, "order_rank")
})
