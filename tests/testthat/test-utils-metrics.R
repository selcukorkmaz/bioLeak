test_that("cindex_pairwise computes concordance correctly", {
  pred <- c(0.1, 0.2, 0.3)
  truth <- c(1, 2, 3)
  expect_equal(bioLeak:::.cindex_pairwise(pred, truth), 1)

  pred2 <- c(0.3, 0.2, 0.1)
  expect_equal(bioLeak:::.cindex_pairwise(pred2, truth), 0)

  expect_true(is.na(bioLeak:::.cindex_pairwise(1, 1)))
  expect_true(is.na(bioLeak:::.cindex_pairwise(c(0.1, 0.2), c(1, 1))))
})

test_that("survival C-index treats predictions as risk scores (higher = shorter survival)", {
  skip_if_not_installed("survival")
  set.seed(1)
  x  <- rnorm(300)
  y  <- survival::Surv(rexp(300, rate = exp(1.5 * x)), rbinom(300, 1, 0.8))
  expect_gt(bioLeak:::.cindex_survival(1.5 * x, y), 0.7)    # true risk score
  expect_lt(bioLeak:::.cindex_survival(-1.5 * x, y), 0.3)   # inverted score
})

test_that("fit_resample survival C-index is above 0.5 for an informative feature", {
  skip_if_not_installed("survival")
  skip_if_not_installed("glmnet")
  set.seed(2)
  n <- 200
  df <- data.frame(id = seq_len(n), x1 = rnorm(n), x2 = rnorm(n))
  df$time <- rexp(n, rate = exp(1.2 * df$x1)); df$status <- rbinom(n, 1, 0.8)
  sp <- make_split_plan(df, outcome = c("time", "status"), mode = "subject_grouped",
                        group = "id", v = 4, progress = FALSE)
  fit <- suppressWarnings(fit_resample(df, outcome = c("time", "status"), splits = sp,
                                       learner = "glmnet", metrics = "cindex"))
  expect_gt(mean(fit@metrics$cindex, na.rm = TRUE), 0.6)
  if (requireNamespace("censored", quietly = TRUE)) {
    suppressWarnings(suppressPackageStartupMessages(library(censored)))
    fit_ph <- suppressWarnings(fit_resample(df, outcome = c("time", "status"), splits = sp,
                                            learner = parsnip::proportional_hazards(),
                                            metrics = "cindex"))
    expect_gt(mean(fit_ph@metrics$cindex, na.rm = TRUE), 0.6)
  }
})

test_that("parsnip survival predictions are oriented as risk for every engine", {
  skip_if_not_installed("survival")
  skip_if_not_installed("censored")
  skip_if_not_installed("rpart")
  suppressWarnings(suppressPackageStartupMessages(library(censored)))
  set.seed(11)
  n <- 240
  df <- data.frame(id = seq_len(n), x1 = rnorm(n), x2 = rnorm(n))
  df$time <- rexp(n, rate = 0.1 * exp(1.2 * df$x1)); df$status <- rbinom(n, 1, 0.8)
  sp <- make_split_plan(df, outcome = c("time", "status"), mode = "subject_grouped",
                        group = "id", v = 4, progress = FALSE)
  # survreg "time" rises with survival; censored rpart "time" rises with hazard.
  specs <- list(
    survreg = parsnip::survival_reg() |> parsnip::set_engine("survival"),
    rpart = parsnip::decision_tree() |> parsnip::set_mode("censored regression") |>
      parsnip::set_engine("rpart"))
  for (nm in names(specs)) {
    f <- suppressWarnings(fit_resample(df, outcome = c("time", "status"), splits = sp,
                                       learner = specs[[nm]], metrics = "cindex"))
    expect_gt(mean(f@metrics$cindex, na.rm = TRUE), 0.6, label = nm)
  }
  # under the null, choosing the orientation on the training fold does not bias C
  set.seed(12)
  null_c <- vapply(1:6, function(r) {
    d <- data.frame(id = seq_len(n), x1 = rnorm(n), x2 = rnorm(n),
                    time = rexp(n, 0.1), status = rbinom(n, 1, 0.8))
    s <- make_split_plan(d, outcome = c("time", "status"), mode = "subject_grouped",
                         group = "id", v = 4, progress = FALSE, seed = r)
    f <- suppressWarnings(fit_resample(d, outcome = c("time", "status"), splits = s,
                                       learner = specs$survreg, metrics = "cindex"))
    mean(f@metrics$cindex, na.rm = TRUE)
  }, numeric(1))
  expect_lt(abs(mean(null_c) - 0.5), 0.05)
})

test_that("survival direction that cannot be decided gives NA predictions with a reason", {
  skip_if_not_installed("survival")
  skip_if_not_installed("censored")
  skip_if_not_installed("rpart")
  suppressWarnings(suppressPackageStartupMessages(library(censored)))
  set.seed(3)
  n <- 120
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
  y <- survival::Surv(rexp(n, 0.1 * exp(1.2 * X$x1)), rbinom(n, 1, 0.8))
  spec <- parsnip::decision_tree() |> parsnip::set_mode("censored regression") |>
    parsnip::set_engine("rpart")
  fit <- parsnip::fit_xy(spec, x = X[1:90, ], y = y[1:90])
  ok <- bioLeak:::.bio_surv_risk_pred(fit, X[91:120, ], X[1:90, ], y[1:90])
  expect_gt(survival::concordance(y[91:120] ~ ok, reverse = TRUE)$concordance, 0.6)
  all_cens <- survival::Surv(as.numeric(y[1:90, 1]), rep(0, 90))
  na <- bioLeak:::.bio_surv_risk_pred(fit, X[91:120, ], X[1:90, ], all_cens)
  expect_true(all(is.na(na)))
  expect_match(attr(na, "surv_undecided"), "no comparable pairs")
  expect_null(attr(ok, "surv_undecided"))
})
