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

test_that("composite specs that require ordering give a clear error when time ordering is missing", {
  expect_error(
    as_leaksplits(fake_spec("composite", ordering_required = TRUE), obs_frame(),
                  outcome = "y", v = 3),
    "requires ordering"
  )
})

test_that("ordered composite and subject specs enforce both time ordering and group disjointness", {
  dat <- obs_frame(n = 24)
  # Balanced longitudinal panel: 6 groups, each observed at time 1, 2, 3, 4
  ord_rank <- rep(1:4, each = 6)
  spec_ord <- fake_spec(
    "composite", n = 24, n_groups = 6, ordering_required = TRUE,
    extra = list(order_rank = ord_rank, timepoint_id = paste0("T", ord_rank))
  )
  spec_ord$time_var <- "order_rank"

  sp <- suppressMessages(as_leaksplits(spec_ord, dat, outcome = "y", v = 3, progress = FALSE))
  expect_s4_class(sp, "LeakSplits")
  expect_identical(sp@mode, "time_series")
  expect_identical(sp@info$group, "group_id")
  expect_identical(sp@info$time, "order_rank")

  cd <- sp@info$coldata
  for (f in sp@indices) {
    expect_equal(length(intersect(cd$group_id[f$train], cd$group_id[f$test])), 0L)
    expect_true(all(cd$order_rank[f$train] < min(cd$order_rank[f$test])))
  }
})

test_that("as_leaksplits forwards custom stratum_var and supports use_block_constraints", {
  dat <- obs_frame(n = 24)
  b <- rep(paste0("B", 1:6), each = 4)
  s <- rep(c("S1", "S2"), length.out = 24)
  spec_blk <- fake_spec(
    "subject", n = 24, n_groups = 12,
    extra = list(batch_group = b, site_group = b, stratum = s)
  )
  spec_blk$block_vars <- c("batch_group", "site_group")
  spec_blk$stratum_var <- "stratum"

  sp_strat <- suppressMessages(as_leaksplits(spec_blk, dat, outcome = "y", v = 3, progress = FALSE))
  expect_identical(sp_strat@info$stratify_by, "stratum")
  expect_true(isTRUE(sp_strat@info$stratify))

  sp_comb <- suppressMessages(as_leaksplits(
    spec_blk, dat, outcome = "y", v = 3,
    use_block_constraints = TRUE, progress = FALSE
  ))
  expect_identical(sp_comb@mode, "combined")
  expect_equal(length(sp_comb@info$constraints), 3L)
  cd <- sp_comb@info$coldata
  for (f in sp_comb@indices) {
    expect_equal(length(intersect(cd$group_id[f$train], cd$group_id[f$test])), 0L)
    expect_equal(length(intersect(cd$batch_group[f$train], cd$batch_group[f$test])), 0L)
    expect_equal(length(intersect(cd$site_group[f$train], cd$site_group[f$test])), 0L)
  }
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

  # min_cut composite severs secondary bridges so as_leaksplits + fit_resample succeed cleanly
  spec_mc <- splitGraph::as_split_spec(
    splitGraph::derive_split_constraints(g, mode = "composite", strategy = "min_cut",
                                         via = c("subject", "site"), k = 3L),
    graph = g)
  set.seed(7)
  dat_shuf <- data.frame(
    sample_id = sample(paste0("S", 1:12)),
    y = factor(rep(c("ctrl", "case"), 6)),
    x1 = rnorm(12),
    x2 = rnorm(12),
    stringsAsFactors = FALSE
  )
  sp_mc <- suppressMessages(as_leaksplits(spec_mc, data = dat_shuf, outcome = "y", v = 3, progress = FALSE))
  expect_identical(as.character(sp_mc@info$coldata$sample_id), dat_shuf$sample_id)
  grp_shuf <- spec_mc$sample_data$group_id[match(dat_shuf$sample_id, spec_mc$sample_data$sample_id)]
  expect_true(groups_respected(sp_mc, grp_shuf))

  custom_glm <- list(
    glm = list(
      fit = function(x, y, task, weights, ...) {
        suppressWarnings(stats::glm(y ~ ., data = as.data.frame(x),
                                    family = stats::binomial(), weights = weights))
      },
      predict = function(object, newdata, task, ...) {
        as.numeric(suppressWarnings(stats::predict(object,
                                                   newdata = as.data.frame(newdata),
                                                   type = "response")))
      }
    )
  )

  fit_mc <- fit_resample(dat_shuf, outcome = "y", splits = sp_mc,
                         learner = "glm", custom_learners = custom_glm, metrics = "auc")
  expect_false("sample_id" %in% fit_mc@feature_names)
  expect_false(any(grepl("^sample_id", fit_mc@feature_names)))
})

test_that("as_leaksplits handles partial NA block_vars, duplicate IDs, and compact grouped time_series in fit_resample", {
  custom_glm <- list(
    glm = list(
      fit = function(x, y, task, weights, ...) {
        suppressWarnings(stats::glm(y ~ ., data = as.data.frame(x),
                                    family = stats::binomial(), weights = weights))
      },
      predict = function(object, newdata, task, ...) {
        as.numeric(suppressWarnings(stats::predict(object,
                                                   newdata = as.data.frame(newdata),
                                                   type = "response")))
      }
    )
  )

  dat <- obs_frame(n = 36)
  dat$x2 <- rnorm(36)
  # Partial NA in batch_group (unlinked samples) with non-chaining components
  grp_clean <- rep(paste0("G", 1:9), each = 4)
  b_partial <- c(rep(c("B1", "B2", "B3"), each = 8), rep(NA_character_, 12))
  spec_na <- fake_spec("subject", n = 36, n_groups = 9,
                       extra = list(group_id = grp_clean, batch_group = b_partial))
  spec_na$block_vars <- "batch_group"
  sp_na <- suppressMessages(as_leaksplits(
    spec_na, dat, outcome = "y", v = 3, use_block_constraints = TRUE, progress = FALSE
  ))
  expect_identical(sp_na@mode, "combined")
  fit_comb <- fit_resample(dat, outcome = "y", splits = sp_na,
                           learner = "glm", custom_learners = custom_glm, metrics = "auc")
  expect_false(any(grepl("sample_id|group_id|batch_group", fit_comb@feature_names)))

  # Duplicate sample_id in data is rejected even when nrow matches
  dat_dup <- dat
  dat_dup$sample_id[36] <- dat_dup$sample_id[1]
  expect_error(as_leaksplits(fake_spec("subject", n = 36), dat_dup, outcome = "y"), "check ID alignment")

  # Compact grouped time_series split resolves group exclusion inside fit_resample
  ord_rank <- rep(1:6, each = 6)
  grp_ts <- rep(paste0("G", 1:6), each = 6)
  spec_ts <- fake_spec("composite", n = 36, n_groups = 6, ordering_required = TRUE,
                       extra = list(group_id = grp_ts, order_rank = ord_rank))
  spec_ts$time_var <- "order_rank"
  sp_ts_compact <- suppressMessages(as_leaksplits(
    spec_ts, dat, outcome = "y", v = 3, compact = TRUE, progress = FALSE
  ))
  expect_true(isTRUE(sp_ts_compact@info$compact))
  fit_ts <- fit_resample(dat, outcome = "y", splits = sp_ts_compact,
                         learner = "glm", custom_learners = custom_glm, metrics = "auc")
  expect_s4_class(fit_ts, "LeakFit")
})
