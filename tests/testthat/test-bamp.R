fast_mcmc <- list(number_of_iterations = 2000, burn_in = 1000, step = 20, tuning = 400)

test_that("bamp AP model returns apc object with expected structure", {
  skip_on_cran()
  data(apc, package = "bamp", envir = environment())
  model <- bamp(cases, population,
                age = "rw1", period = "rw1", cohort = " ",
                periods_per_agegroup = 5,
                mcmc.options = fast_mcmc,
                parallel = FALSE, verbose = FALSE)

  expect_s3_class(model, "apc")
  expect_false(is.null(model$samples))
  expect_false(is.null(model$model))
  expect_equal(model$model$age,    "rw1")
  expect_equal(model$model$period, "rw1")
  expect_equal(model$model$cohort, " ")
})

test_that("effects.apc returns list with age and period components", {
  skip_on_cran()
  data(apc, package = "bamp", envir = environment())
  model <- bamp(cases, population,
                age = "rw1", period = "rw1", cohort = " ",
                periods_per_agegroup = 5,
                mcmc.options = fast_mcmc,
                parallel = FALSE, verbose = FALSE)

  eff <- effects(model)
  expect_type(eff, "list")
  expect_false(is.null(eff$age))
  expect_false(is.null(eff$period))
})
