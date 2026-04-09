test_that(".bio_capture_provenance returns list with all required keys", {
  prov <- bioLeak:::.bio_capture_provenance()
  expect_true(is.list(prov))
  required_keys <- c("r_version", "packages", "platform", "timestamp",
                      "git_sha", "hardware")
  expect_true(all(required_keys %in% names(prov)))
})

test_that("packages is a data.frame with package and version columns", {
  prov <- bioLeak:::.bio_capture_provenance()
  expect_true(is.data.frame(prov$packages))
  expect_true("package" %in% names(prov$packages))
  expect_true("version" %in% names(prov$packages))
  expect_gt(nrow(prov$packages), 0)
})

test_that("git_sha is NA when capture_git = FALSE", {
  prov <- bioLeak:::.bio_capture_provenance(capture_git = FALSE)
  expect_true(is.na(prov$git_sha))
})

test_that("hardware is NULL when capture_hardware = FALSE", {
  prov <- bioLeak:::.bio_capture_provenance(capture_hardware = FALSE)
  expect_null(prov$hardware)
})

test_that("hardware is a list with expected keys when capture_hardware = TRUE", {
  prov <- bioLeak:::.bio_capture_provenance(capture_hardware = TRUE)
  expect_true(is.list(prov$hardware))
  expect_true(all(c("sysname", "nodename", "machine") %in% names(prov$hardware)))
})

test_that("fit_resample output has provenance in info slot", {
  set.seed(42)
  df <- data.frame(
    subject = rep(1:30, each = 3),
    outcome = rep(c(0, 1), length.out = 90),
    x1 = rnorm(90),
    x2 = rnorm(90)
  )
  splits <- make_split_plan(df, outcome = "outcome",
                            mode = "subject_grouped", group = "subject",
                            v = 3, progress = FALSE)
  fit <- fit_resample(df, outcome = "outcome", splits = splits,
                      learner = "glmnet", preprocess = list(
                        impute = list(method = "median"),
                        normalize = list(method = "zscore")
                      ))
  expect_true("provenance" %in% names(fit@info))
  prov <- fit@info$provenance
  expect_true(is.list(prov))
  expect_true(nchar(prov$r_version) > 0)
  expect_true(is.data.frame(prov$packages))
})
