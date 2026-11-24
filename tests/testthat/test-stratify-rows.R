# J8 (JSS review): stratify_method = "rows" balances class row counts across
# folds for grouped designs with unequal group sizes and mixed-outcome groups.

.j8_data <- function(seed, n_groups = 40) {
  set.seed(seed)
  sizes <- sample(c(1:3, 10:20), n_groups, replace = TRUE)
  p1 <- sample(c(0.05, 0.5, 0.95), n_groups, replace = TRUE)
  do.call(rbind, lapply(seq_len(n_groups), function(g) {
    data.frame(subject = paste0("s", g),
               y = factor(ifelse(stats::runif(sizes[g]) < p1[g], "case", "ctrl"),
                          levels = c("case", "ctrl")),
               x = stats::rnorm(sizes[g]))
  }))
}

# Largest deviation of a fold's case share (in rows) from the overall share,
# and fold-size spread relative to the mean fold size.
.j8_imbalance <- function(df, plan) {
  p_all <- mean(df$y == "case")
  dev <- vapply(plan@indices, function(f) abs(mean(df$y[f$test] == "case") - p_all), numeric(1))
  sz <- vapply(plan@indices, function(f) length(f$test), numeric(1))
  c(class = max(dev), size = stats::sd(sz) / mean(sz))
}

test_that("stratify_method = 'rows' improves row-level class balance over 'groups' and random", {
  res <- t(vapply(1:20, function(s) {
    df <- .j8_data(s)
    pr <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject",
                          v = 5, stratify = TRUE, stratify_method = "rows", seed = s, progress = FALSE)
    pg <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject",
                          v = 5, stratify = TRUE, seed = s, progress = FALSE)
    pn <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject",
                          v = 5, stratify = FALSE, seed = s, progress = FALSE)
    c(rows = .j8_imbalance(df, pr), groups = .j8_imbalance(df, pg), none = .j8_imbalance(df, pn))
  }, numeric(6)))
  m <- colMeans(res)
  expect_lt(m["rows.class"], 0.5 * m["groups.class"])
  expect_lt(m["rows.class"], 0.5 * m["none.class"])
  expect_lt(m["rows.size"], 0.5 * m["groups.size"])
})

test_that("stratify_method = 'rows' keeps groups intact, uses all folds and is reproducible", {
  df <- .j8_data(3)
  p1 <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject",
                        v = 5, repeats = 2, stratify = TRUE, stratify_method = "rows",
                        seed = 9, progress = FALSE)
  p2 <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject",
                        v = 5, repeats = 2, stratify = TRUE, stratify_method = "rows",
                        seed = 9, progress = FALSE)
  expect_identical(p1@indices, p2@indices)
  expect_equal(length(p1@indices), 10L)
  for (f in p1@indices) {
    expect_length(intersect(df$subject[f$test], df$subject[f$train]), 0L)
  }
  for (r in 1:2) {
    test_rows <- unlist(lapply(Filter(function(f) f$repeat_id == r, p1@indices), `[[`, "test"))
    expect_setequal(test_rows, seq_len(nrow(df)))
  }
  expect_identical(p1@info$stratify_method, "rows")
  # repeats differ
  r1 <- lapply(Filter(function(f) f$repeat_id == 1, p1@indices), `[[`, "test")
  r2 <- lapply(Filter(function(f) f$repeat_id == 2, p1@indices), `[[`, "test")
  expect_false(identical(r1, r2))
})

test_that("default stratify_method leaves existing plans unchanged", {
  df <- .j8_data(4)
  p_def <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject",
                           v = 5, stratify = TRUE, seed = 2, progress = FALSE)
  p_grp <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject",
                           v = 5, stratify = TRUE, stratify_method = "groups", seed = 2,
                           progress = FALSE)
  expect_identical(p_def@indices, p_grp@indices)
  expect_identical(p_def@info$stratify_method, "groups")
})

test_that("stratify_method = 'rows' works for batch_blocked, combined and compact plans", {
  df <- .j8_data(5)
  df$batch <- paste0("b", as.integer(factor(df$subject)) %% 12)
  pb <- make_split_plan(df, outcome = "y", mode = "batch_blocked", batch = "batch", v = 4,
                        stratify = TRUE, stratify_method = "rows", seed = 1, progress = FALSE)
  expect_identical(pb@info$stratify_method, "rows")
  for (f in pb@indices) expect_length(intersect(df$batch[f$test], df$batch[f$train]), 0L)

  df$site <- paste0("site", as.integer(factor(df$subject)) %% 8)
  pc <- make_split_plan(df, outcome = "y", mode = "combined",
                        constraints = list(list(type = "subject", col = "subject"),
                                           list(type = "batch", col = "site")),
                        v = 4, stratify = TRUE, stratify_method = "rows", seed = 1,
                        progress = FALSE)
  expect_identical(pc@info$stratify_method, "rows")
  for (f in pc@indices) {
    expect_length(intersect(df$subject[f$test], df$subject[f$train]), 0L)
    expect_length(intersect(df$site[f$test], df$site[f$train]), 0L)
  }

  pk <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject",
                        v = 5, stratify = TRUE, stratify_method = "rows", seed = 1,
                        compact = TRUE, progress = FALSE)
  pe <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject",
                        v = 5, stratify = TRUE, stratify_method = "rows", seed = 1,
                        progress = FALSE)
  fa <- pk@info$fold_assignments[[1]]
  expect_identical(lapply(pe@indices, `[[`, "test"),
                   lapply(seq_len(5), function(k) which(fa == k)))
})

test_that("stratify_method = 'rows' applies when majority labels collapse to one class", {
  # every subject is mostly 'ctrl', but cases exist at row level
  set.seed(1)
  df <- data.frame(subject = rep(paste0("s", 1:20), each = 5),
                   y = factor(rep(c("case", "ctrl", "ctrl", "ctrl", "ctrl"), 20)))
  df$y[1:5] <- "ctrl"
  expect_warning(
    pg <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject",
                          v = 5, stratify = TRUE, seed = 1, progress = FALSE),
    "only one class"
  )
  pr <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject",
                        v = 5, stratify = TRUE, stratify_method = "rows", seed = 1,
                        progress = FALSE)
  expect_true(isTRUE(pr@info$stratify))
})

test_that("stratify_method = 'rows' warns when stratify is FALSE and is validated", {
  df <- .j8_data(6)
  expect_warning(
    make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject",
                    v = 5, stratify_method = "rows", progress = FALSE),
    "stratify_method"
  )
  expect_error(
    make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject",
                    v = 5, stratify = TRUE, stratify_method = "bogus", progress = FALSE)
  )
})

test_that("nested plans accept stratify_method = rows", {
  df <- .j8_data(7)
  p <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject",
                       v = 3, stratify = TRUE, stratify_method = "rows", nested = TRUE,
                       seed = 1, progress = FALSE)
  expect_identical(p@info$stratify_method, "rows")
  expect_false(is.null(p@info$inner))
})

test_that("refit permutations re-draw plans with the recorded stratify_method", {
  set.seed(1)
  df <- .j8_data(8)
  df$x2 <- stats::rnorm(nrow(df))
  sp <- make_split_plan(df, outcome = "y", mode = "subject_grouped", group = "subject",
                        v = 4, stratify = TRUE, stratify_method = "rows", seed = 1,
                        progress = FALSE)
  fit <- NULL
  utils::capture.output(
    fit <- fit_resample(df[, c("subject", "y", "x", "x2")], outcome = "y", splits = sp,
                        learner = "glm", custom_learners = make_custom_learners(),
                        metrics = "auc", refit = FALSE, seed = 1)
  )
  orig <- make_split_plan
  seen <- character(0)
  local_mocked_bindings(make_split_plan = function(...) {
    seen <<- c(seen, list(...)$stratify_method %||% "default")
    orig(...)
  })
  utils::capture.output(suppressWarnings(
    audit_leakage(fit, metric = "auc", B = 3, perm_refit = TRUE, perm_folds = "redraw",
                  target_scan = FALSE, target_scan_multivariate = FALSE)
  ))
  expect_gte(length(seen), 3L)
  expect_true(all(seen == "rows"))
})
