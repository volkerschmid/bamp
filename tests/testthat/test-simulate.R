test_that("apcSimulate returns correct structure and dimensions", {
  set.seed(42)
  noa <- 5; nop <- 6; ppag <- 2
  noc <- ppag * (noa - 1) + nop
  age    <- rep(0, noa)
  period <- rep(0, nop)
  cohort <- rep(0, noc)

  result <- apcSimulate(-3, age, period, cohort, ppag, 1000)

  expect_named(result, c("cases", "population"))
  expect_equal(dim(result$cases),      c(nop, noa))
  expect_equal(dim(result$population), c(nop, noa))
})

test_that("apcSimulate cases are within [0, population]", {
  set.seed(42)
  noa <- 5; nop <- 6; ppag <- 2
  noc <- ppag * (noa - 1) + nop
  result <- apcSimulate(-3, rep(0, noa), rep(0, nop), rep(0, noc), ppag, 1000)

  expect_true(all(result$cases >= 0))
  expect_true(all(result$cases <= result$population))
})

test_that("apcSimulate stops when cohort has wrong length", {
  expect_error(
    apcSimulate(-3, rep(0, 5), rep(0, 6), rep(0, 10), 2, 1000),
    regexp = "cohort"
  )
})

test_that("apcSimulate works with scalar population", {
  set.seed(1)
  noa <- 4; nop <- 5; ppag <- 1
  noc <- ppag * (noa - 1) + nop
  result <- apcSimulate(0, rep(0, noa), rep(0, nop), rep(0, noc), ppag, 1e5)

  expect_equal(dim(result$cases), c(nop, noa))
  # with intercept=0 approx 50% probability; rates should be near 0.5
  rates <- result$cases / result$population
  expect_true(all(rates > 0.3 & rates < 0.7))
})
