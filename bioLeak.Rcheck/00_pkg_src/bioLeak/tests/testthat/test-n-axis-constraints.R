test_that("constraints with 2 axes produces identical results to primary_axis/secondary_axis", {
  set.seed(42)
  df <- data.frame(
    subject = rep(1:20, each = 3),
    batch = rep(rep(c("A", "B", "C", "D"), each = 5), 3),
    outcome = rbinom(60, 1, 0.5),
    x1 = rnorm(60),
    x2 = rnorm(60)
  )

  splits_legacy <- suppressWarnings(make_split_plan(
    df, outcome = "outcome",
    mode = "combined",
    primary_axis = list(type = "subject", col = "subject"),
    secondary_axis = list(type = "batch", col = "batch"),
    v = 3, seed = 1, progress = FALSE
  ))

  splits_new <- suppressWarnings(make_split_plan(
    df, outcome = "outcome",
    mode = "combined",
    constraints = list(
      list(type = "subject", col = "subject"),
      list(type = "batch", col = "batch")
    ),
    v = 3, seed = 1, progress = FALSE
  ))

  expect_equal(length(splits_legacy@indices), length(splits_new@indices))
  for (i in seq_along(splits_legacy@indices)) {
    expect_equal(splits_legacy@indices[[i]]$train, splits_new@indices[[i]]$train)
    expect_equal(splits_legacy@indices[[i]]$test, splits_new@indices[[i]]$test)
  }
})

test_that("3-axis constraint (subject + batch + site) enforces all three", {
  set.seed(42)
  # Design: 9 groups of 5 subjects each. Each group has a unique batch-site pair.
  # 3 batches x 3 sites, nested so batch i maps to sites {i} only.
  # This makes the 3 axes perfectly aligned in 3 blocks, enabling v=3 splits.
  n_groups <- 9
  subjects_per_group <- 5
  group_batch <- rep(paste0("B", 1:3), each = 3)
  group_site <- rep(paste0("S", 1:3), times = 3)
  df <- data.frame(
    subject = rep(seq_len(n_groups), each = subjects_per_group),
    batch = rep(group_batch, each = subjects_per_group),
    site = rep(group_site, each = subjects_per_group),
    outcome = rbinom(n_groups * subjects_per_group, 1, 0.5),
    x1 = rnorm(n_groups * subjects_per_group),
    stringsAsFactors = FALSE
  )
  splits <- suppressWarnings(make_split_plan(
    df, outcome = "outcome",
    mode = "combined",
    constraints = list(
      list(type = "subject", col = "subject"),
      list(type = "batch", col = "batch"),
      list(type = "study", col = "site")
    ),
    v = 3, progress = FALSE
  ))
  expect_s4_class(splits, "LeakSplits")

  for (fold in splits@indices) {
    if (is.null(fold$train) || is.null(fold$test)) next
    # No subject overlap
    expect_length(intersect(unique(df$subject[fold$train]),
                            unique(df$subject[fold$test])), 0)
    # No batch overlap
    expect_length(intersect(unique(df$batch[fold$train]),
                            unique(df$batch[fold$test])), 0)
    # No site overlap
    expect_length(intersect(unique(df$site[fold$train]),
                            unique(df$site[fold$test])), 0)
  }
})

test_that("check_split_overlap validates all N axes", {
  set.seed(42)
  n_groups <- 9
  subjects_per_group <- 5
  group_batch <- rep(paste0("B", 1:3), each = 3)
  group_site <- rep(paste0("S", 1:3), times = 3)
  df <- data.frame(
    subject = rep(seq_len(n_groups), each = subjects_per_group),
    batch = rep(group_batch, each = subjects_per_group),
    site = rep(group_site, each = subjects_per_group),
    outcome = rbinom(n_groups * subjects_per_group, 1, 0.5),
    x1 = rnorm(n_groups * subjects_per_group),
    stringsAsFactors = FALSE
  )
  splits <- suppressWarnings(make_split_plan(
    df, outcome = "outcome",
    mode = "combined",
    constraints = list(
      list(type = "subject", col = "subject"),
      list(type = "batch", col = "batch"),
      list(type = "study", col = "site")
    ),
    v = 3, progress = FALSE
  ))
  result <- check_split_overlap(splits)
  expect_true(all(result$pass))
  expect_true(all(c("subject", "batch", "site") %in% result$col))
})

test_that("legacy primary_axis/secondary_axis still works", {
  set.seed(42)
  df <- data.frame(
    subject = rep(1:20, each = 3),
    batch = rep(rep(c("A", "B", "C", "D"), each = 5), 3),
    outcome = rbinom(60, 1, 0.5),
    x1 = rnorm(60)
  )
  splits <- suppressWarnings(make_split_plan(
    df, outcome = "outcome",
    mode = "combined",
    primary_axis = list(type = "subject", col = "subject"),
    secondary_axis = list(type = "batch", col = "batch"),
    v = 3, progress = FALSE
  ))
  expect_s4_class(splits, "LeakSplits")
  # constraints should be populated in info
  expect_equal(length(splits@info$constraints), 2)
})

test_that("error when both constraints and primary_axis provided", {
  df <- data.frame(subject = 1:10, batch = rep("A", 10),
                   outcome = rbinom(10, 1, 0.5), x1 = rnorm(10))
  expect_error(
    make_split_plan(df, outcome = "outcome", mode = "combined",
                    primary_axis = list(type = "subject", col = "subject"),
                    constraints = list(
                      list(type = "subject", col = "subject"),
                      list(type = "batch", col = "batch")
                    ),
                    v = 3, progress = FALSE),
    "Cannot specify both"
  )
})

test_that("error when constraints has fewer than 2 elements", {
  df <- data.frame(subject = 1:10, batch = rep("A", 10),
                   outcome = rbinom(10, 1, 0.5), x1 = rnorm(10))
  expect_error(
    make_split_plan(df, outcome = "outcome", mode = "combined",
                    constraints = list(
                      list(type = "subject", col = "subject")
                    ),
                    v = 3, progress = FALSE),
    "at least 2 elements"
  )
})

test_that("error when constraint col doesn't exist", {
  df <- data.frame(subject = 1:10, batch = rep("A", 10),
                   outcome = rbinom(10, 1, 0.5), x1 = rnorm(10))
  expect_error(
    make_split_plan(df, outcome = "outcome", mode = "combined",
                    constraints = list(
                      list(type = "subject", col = "subject"),
                      list(type = "batch", col = "nonexistent")
                    ),
                    v = 3, progress = FALSE),
    "not found"
  )
})

test_that("constraints info slot stores all axes", {
  set.seed(42)
  n_groups <- 9
  subjects_per_group <- 5
  group_batch <- rep(paste0("B", 1:3), each = 3)
  group_site <- rep(paste0("S", 1:3), times = 3)
  df <- data.frame(
    subject = rep(seq_len(n_groups), each = subjects_per_group),
    batch = rep(group_batch, each = subjects_per_group),
    site = rep(group_site, each = subjects_per_group),
    outcome = rbinom(n_groups * subjects_per_group, 1, 0.5),
    x1 = rnorm(n_groups * subjects_per_group),
    stringsAsFactors = FALSE
  )
  splits <- suppressWarnings(make_split_plan(
    df, outcome = "outcome",
    mode = "combined",
    constraints = list(
      list(type = "subject", col = "subject"),
      list(type = "batch", col = "batch"),
      list(type = "study", col = "site")
    ),
    v = 3, progress = FALSE
  ))
  expect_equal(length(splits@info$constraints), 3)
  # backward compat: primary_axis and secondary_axis still populated
  expect_equal(splits@info$primary_axis$col, "subject")
  expect_equal(splits@info$secondary_axis$col, "batch")
})
