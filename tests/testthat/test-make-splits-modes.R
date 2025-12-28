test_that("make_splits respects grouping constraints", {
  set.seed(1)
  df <- data.frame(
    outcome = rep(c(0, 1), each = 10),
    subject = rep(1:10, each = 2),
    batch = rep(letters[1:4], length.out = 20),
    study = rep(LETTERS[1:5], length.out = 20),
    time = seq_len(20),
    x1 = rnorm(20),
    x2 = rnorm(20)
  )

  splits_subject <- make_splits(df, outcome = "outcome",
                                mode = "subject_grouped", group = "subject",
                                v = 5, repeats = 1, stratify = TRUE,
                                seed = 1, progress = FALSE)
  for (fold in splits_subject@indices) {
    tr <- df$subject[fold$train]
    te <- df$subject[fold$test]
    expect_equal(length(intersect(unique(tr), unique(te))), 0)
  }

  splits_batch <- make_splits(df, outcome = "outcome",
                              mode = "batch_blocked", batch = "batch",
                              v = 4, repeats = 1, stratify = FALSE,
                              seed = 1, progress = FALSE)
  for (fold in splits_batch@indices) {
    tr <- df$batch[fold$train]
    te <- df$batch[fold$test]
    expect_equal(length(intersect(unique(tr), unique(te))), 0)
  }

  splits_study <- make_splits(df, outcome = "outcome",
                              mode = "study_loocv", study = "study",
                              seed = 1, progress = FALSE)
  for (fold in splits_study@indices) {
    te_study <- unique(df$study[fold$test])
    expect_equal(length(te_study), 1)
    expect_false(te_study %in% unique(df$study[fold$train]))
  }

  splits_time <- make_splits(df, outcome = "outcome",
                             mode = "time_series", time = "time",
                             v = 4, horizon = 1, seed = 1, progress = FALSE)
  expect_true(length(splits_time@indices) > 0)
  for (fold in splits_time@indices) {
    tmin <- min(df$time[fold$test])
    expect_true(all(df$time[fold$train] <= (tmin - 1)))
  }
})
