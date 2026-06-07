test_that("coh returns correct cohort indices", {
  # oldest age group, first period = first cohort
  expect_equal(coh(10, 1, 10, 5), 1)
  # youngest age group, last period = last cohort
  expect_equal(coh(1, 8, 10, 5), 53)
  # general: (noa - agegroup) * ppag + period
  expect_equal(coh(3, 2, 5, 2), (5 - 3) * 2 + 2)
})

test_that("coh is vectorized over age groups", {
  result <- coh(1:3, 1, 5, 2)
  expect_equal(result, c((5-1)*2+1, (5-2)*2+1, (5-3)*2+1))
})

test_that("coh is vectorized over periods", {
  result <- coh(5, 1:3, 5, 2)
  expect_equal(result, 1:3)
})
