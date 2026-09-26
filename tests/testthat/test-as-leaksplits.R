# Tests for as_leaksplits(), the splitGraph split_spec adapter.
# Before 0.3.8 the site/region/platform/assay/relatedness/spatial modes failed
# with "subscript out of bounds" and composite failed with "'primary_axis'
# must be a list".

fake_spec <- function(mode, n = 24, n_groups = 6, ordering_required = FALSE,
                      extra = list()) {
  sd <- data.frame(
    sample_id = paste0("S", seq_len(n)),
    group_id = paste0("G", rep(seq_len(n_groups), length.out = n)),
    stringsAsFactors = FALSE
  )
  for (nm in names(extra)) sd[[nm]] <- extra[[nm]]
  structure(list(sample_data = sd, group_var = "group_id",
                 constraint_mode = mode, time_var = NULL,
                 ordering_required = ordering_required),
            class = "split_spec")
}

obs_frame <- function(n = 24, seed = 1) {
  set.seed(seed)
  data.frame(sample_id = paste0("S", seq_len(n)),
             y = factor(rep(c("a", "b"), length.out = n)),
             x1 = rnorm(n), stringsAsFactors = FALSE)
}

groups_respected <- function(splits, groups) {
  all(vapply(splits@indices, function(f) {
    !length(intersect(groups[f$train], groups[f$test]))
  }, logical(1)))
}

test_that("grouping modes map to subject_grouped folds on group_id", {
  dat <- obs_frame()
  for (mode in c("subject", "site", "region", "platform", "assay",
                 "relatedness", "spatial", "composite")) {
    spec <- fake_spec(mode)
    sp <- suppressMessages(as_leaksplits(spec, dat, outcome = "y", v = 3, progress = FALSE))
    expect_s4_class(sp, "LeakSplits")
    expect_identical(sp@mode, "subject_grouped", info = mode)
    expect_identical(sp@info$group, "group_id", info = mode)
    grp <- spec$sample_data$group_id[match(dat$sample_id, spec$sample_data$sample_id)]
    expect_true(groups_respected(sp, grp), info = mode)
  }
})

test_that("batch and study specs use their blocking columns", {
  dat <- obs_frame()
  b <- rep(c("B1", "B2", "B3"), length.out = 24)
  spec_b <- fake_spec("batch", extra = list(batch_group = b))
  sp_b <- suppressMessages(as_leaksplits(spec_b, dat, outcome = "y", v = 3, progress = FALSE))
  expect_identical(sp_b@mode, "batch_blocked")
  expect_identical(sp_b@info$batch, "batch_group")

  st <- rep(c("ST1", "ST2", "ST3"), each = 8)
  spec_s <- fake_spec("study", extra = list(study_group = st))
  sp_s <- suppressMessages(as_leaksplits(spec_s, dat, outcome = "y", progress = FALSE))
  expect_identical(sp_s@mode, "study_loocv")
  expect_identical(sp_s@info$study, "study_group")
})

test_that("composite specs that require ordering give a clear error", {
  expect_error(
    as_leaksplits(fake_spec("composite", ordering_required = TRUE), obs_frame(),
                  outcome = "y", v = 3),
    "requires ordering"
  )
})

test_that("a spec with a single group gives a clear error", {
  expect_error(
    as_leaksplits(fake_spec("composite", n_groups = 1), obs_frame(),
                  outcome = "y", v = 3),
    "places all samples in 1 group"
  )
})

test_that("unknown future modes fall back to grouped folds with a warning", {
  expect_warning(
    sp <- suppressMessages(as_leaksplits(fake_spec("galaxy"), obs_frame(),
                                         outcome = "y", v = 3, progress = FALSE)),
    "not recognised"
  )
  expect_identical(sp@mode, "subject_grouped")
})

test_that("input validation errors are informative", {
  dat <- obs_frame()
  expect_error(as_leaksplits(list(), dat, outcome = "y"), "split_spec")
  expect_error(as_leaksplits(fake_spec("subject"), dat, outcome = "nope"), "outcome")
  expect_error(as_leaksplits(fake_spec("subject"), dat, outcome = "y",
                             sample_id_col = "id"), "not found")
  dat_short <- dat[1:20, ]
  expect_error(as_leaksplits(fake_spec("subject"), dat_short, outcome = "y"),
               "check ID alignment")
  dat_clash <- dat
  dat_clash$group_id <- "x"
  expect_error(as_leaksplits(fake_spec("subject"), dat_clash, outcome = "y"),
               "already has column")
})

test_that("real splitGraph specs are accepted for every mode", {
  skip_if_not_installed("splitGraph")
  meta <- data.frame(
    sample_id    = paste0("S", 1:12),
    subject_id   = rep(paste0("P", 1:6), each = 2),
    study_id     = rep(c("ST1", "ST2"), each = 6),
    batch_id     = rep(c("B1", "B2", "B3"), times = 4),
    site_id      = rep(c("NYC", "BOS", "SFO"), times = 4),
    region_id    = rep(c("cortex", "liver"), times = 6),
    platform_id  = rep(c("HiSeq", "NovaSeq"), each = 6),
    assay_id     = rep(c("rna", "atac", "rna"), times = 4),
    timepoint_id = rep(c("T1", "T2"), 6),
    time_index   = rep(c(1, 2), 6),
    stringsAsFactors = FALSE
  )
  pairs <- data.frame(id1 = c("P1", "P3"), id2 = c("P2", "P4"),
                      kinship = c(0.25, 0.5), stringsAsFactors = FALSE)
  coords <- data.frame(sample_id = meta$sample_id,
                       x = c(0, 0.1, 5, 5.1, 10, 10.1, 15, 15.1, 20, 20.1, 25, 25.1),
                       y = 0)
  g <- splitGraph::graph_from_metadata(meta)
  g <- splitGraph::add_edges(g, list(
    splitGraph::relatedness_edges_from_kinship(pairs, threshold = 0.1),
    splitGraph::spatial_edges_from_coords(coords, radius = 0.5)
  ))
  set.seed(1)
  dat <- data.frame(sample_id = paste0("S", 1:12), y = rbinom(12, 1, 0.5),
                    stringsAsFactors = FALSE)

  for (mode in c("subject", "batch", "study", "time", "site", "region",
                 "platform", "assay", "relatedness", "spatial")) {
    spec <- splitGraph::as_split_spec(
      splitGraph::derive_split_constraints(g, mode = mode), graph = g)
    sp <- suppressMessages(as_leaksplits(spec, data = dat, outcome = "y", v = 2,
                                         progress = FALSE))
    expect_s4_class(sp, "LeakSplits")
    if (!mode %in% c("batch", "study", "time")) {
      grp <- spec$sample_data$group_id[match(dat$sample_id, spec$sample_data$sample_id)]
      expect_true(groups_respected(sp, grp), info = mode)
    }
  }

  spec_rb <- splitGraph::as_split_spec(
    splitGraph::derive_split_constraints(g, mode = "composite", strategy = "rule_based",
                                         via = c("subject", "site")),
    graph = g)
  sp_rb <- suppressMessages(as_leaksplits(spec_rb, data = dat, outcome = "y", v = 2,
                                          progress = FALSE))
  expect_identical(sp_rb@mode, "subject_grouped")
  grp <- spec_rb$sample_data$group_id[match(dat$sample_id, spec_rb$sample_data$sample_id)]
  expect_true(groups_respected(sp_rb, grp))

  # strict composite over subject + site links every sample into one component
  spec_st <- splitGraph::as_split_spec(
    splitGraph::derive_split_constraints(g, mode = "composite", strategy = "strict",
                                         via = c("subject", "site")),
    graph = g)
  expect_error(as_leaksplits(spec_st, data = dat, outcome = "y", v = 2),
               "places all samples in 1 group")
})
