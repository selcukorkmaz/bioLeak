test_that("recipe leakage validator flags trained recipes", {
  skip_if_not_installed("recipes")

  df <- data.frame(
    outcome = factor(rep(c(0, 1), each = 5), levels = c(0, 1)),
    x1 = rnorm(10),
    x2 = rnorm(10)
  )

  rec <- recipes::recipe(outcome ~ x1 + x2, data = df)
  rec_trained <- recipes::prep(rec, training = df, retain = TRUE)

  expect_warning(
    bioLeak:::.bio_validate_recipe_graph(rec_trained, context = "test", mode = "warn"),
    "trained recipe"
  )
  expect_error(
    bioLeak:::.bio_validate_recipe_graph(rec_trained, context = "test", mode = "error"),
    "trained recipe"
  )
})

test_that("recipe leakage validator flags negative lag steps", {
  skip_if_not_installed("recipes")

  df <- data.frame(
    outcome = factor(rep(c(0, 1), each = 5), levels = c(0, 1)),
    x1 = rnorm(10),
    x2 = rnorm(10)
  )
  rec_lag <- recipes::recipe(outcome ~ x1 + x2, data = df) |>
    recipes::step_lag(x1, lag = -1)

  expect_warning(
    bioLeak:::.bio_validate_recipe_graph(rec_lag, context = "test", mode = "warn"),
    "negative lag"
  )
})

test_that("workflow leakage validator flags trained workflows", {
  skip_if_not_installed("workflows")
  skip_if_not_installed("parsnip")

  df <- data.frame(
    outcome = factor(rep(c(0, 1), each = 6), levels = c(0, 1)),
    x1 = rnorm(12),
    x2 = rnorm(12)
  )
  wf <- workflows::workflow() |>
    workflows::add_model(
      parsnip::logistic_reg(mode = "classification") |>
        parsnip::set_engine("glm")
    ) |>
    workflows::add_formula(outcome ~ x1 + x2)
  wf_trained <- workflows::fit(wf, data = df)

  expect_warning(
    bioLeak:::.bio_validate_workflow_graph(wf_trained, context = "test", mode = "warn"),
    "trained workflow"
  )
})
