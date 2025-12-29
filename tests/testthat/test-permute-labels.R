test_that("quantile break cache and majority level helpers behave", {
  vals <- c(1, 2, 3, 4, NA)
  br1 <- bioLeak:::.get_cached_quantile_breaks(vals, probs = c(0, 0.5, 1))
  br2 <- bioLeak:::.get_cached_quantile_breaks(vals, probs = c(0, 0.5, 1))
  expect_equal(br1, br2)

  expect_equal(bioLeak:::.majority_level(c("a", "b", "a", NA)), "a")
  expect_true(is.na(bioLeak:::.majority_level(c(NA, NA))))
})

test_that("grouped and within-group permutation functions preserve values", {
  set.seed(1)
  y <- c(1, 2, 3, 4, 5, 6)
  subj <- c("s1", "s1", "s2", "s2", "s3", "s3")
  perm <- bioLeak:::.permute_subject_grouped(y, subj)
  expect_equal(sort(perm), sort(y))

  group <- c("g1", "g1", "g2", "g2", "g2", "g1")
  perm2 <- bioLeak:::.permute_within_group(y, group)
  expect_equal(sort(perm2), sort(y))

  perm3 <- bioLeak:::.permute_within_batch(y, group)
  perm4 <- bioLeak:::.permute_within_study(y, group)
  expect_equal(sort(perm3), sort(y))
  expect_equal(sort(perm4), sort(y))
})

test_that("permute_labels_factory returns per-fold permutations", {
  set.seed(1)
  df <- make_class_df(20)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "subject_grouped", group = "subject",
                              v = 4, stratify = TRUE, seed = 1)
  perm_fun <- bioLeak:::.permute_labels_factory(
    cd = df, outcome = "outcome", mode = "subject_grouped",
    folds = splits@indices, perm_stratify = TRUE,
    time_block = "circular", block_len = 2, seed = 1,
    group_col = "subject", batch_col = "batch", study_col = "study"
  )
  out <- perm_fun(1)
  expect_equal(length(out), length(splits@indices))
  expect_equal(length(out[[1]]), length(splits@indices[[1]]$test))
})

test_that("permute_labels_factory warns for small numeric stratification", {
  df <- make_class_df(12)
  df$outcome <- rnorm(12)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "batch_blocked", batch = "batch",
                              v = 3, stratify = FALSE, seed = 1)
  expect_warning({
    perm_fun <- bioLeak:::.permute_labels_factory(
      cd = df, outcome = "outcome", mode = "batch_blocked",
      folds = splits@indices, perm_stratify = TRUE,
      time_block = "circular", block_len = 2, seed = 1,
      batch_col = "batch"
    )
  }, "requires at least 20")
  out <- perm_fun(1)
  expect_equal(length(out), length(splits@indices))
})

test_that("permute_labels_factory handles time-series permutations", {
  df <- make_class_df(30)
  splits <- make_splits_quiet(df, outcome = "outcome",
                              mode = "time_series", time = "time",
                              v = 4, seed = 1)
  perm_fun <- bioLeak:::.permute_labels_factory(
    cd = df, outcome = "outcome", mode = "time_series",
    folds = splits@indices, perm_stratify = FALSE,
    time_block = "stationary", block_len = 3, seed = 1
  )
  out <- perm_fun(1)
  expect_equal(length(out), length(splits@indices))
  expect_equal(length(out[[1]]), length(splits@indices[[1]]$test))
})

test_that("permute_labels_factory errors on invalid metadata", {
  df <- data.frame(outcome = c(NA, NA), x = 1:2)
  splits <- list(list(test = 1:2, train = integer(0)))
  expect_error(bioLeak:::.permute_labels_factory(
    cd = df, outcome = "outcome", mode = "subject_grouped",
    folds = splits, perm_stratify = FALSE,
    time_block = "circular", block_len = 2, seed = 1
  ), "only NA")
})
