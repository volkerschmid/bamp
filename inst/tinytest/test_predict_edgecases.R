## Regression tests for predict_apc / bamp edge cases (Bugs 1-4)
##  Bug 1: predict_apc used exp(x)/(1+exp(x)) -> NaN for large logits (now plogis)
##  Bug 2: age-period model (no cohort) could not be predicted (abind / undefined c2)
##  Bug 3: cohort=NULL errored "argument is of length zero" (het check before default)
##  Bug 4: non-integer population -> rbinom(size=non-integer) -> NA -> crash
suppressWarnings(suppressMessages(library(bamp)))

set.seed(1)
A <- 10L; P <- 12L
pop <- matrix(1000, P, A)
lograte <- outer(seq(-6, -3, length.out = P), seq(0, 2, length.out = A), `+`)
cases <- matrix(stats::rpois(P * A, pop * plogis(lograte)), P, A)
opts <- list(number_of_iterations = 800, burn_in = "auto", step = "auto", tuning = 100)

## Bug 1: plogis is a finite, stable sigmoid where exp/(1+exp) overflowed to NaN
big <- c(-800, -50, 0, 50, 800)
expect_true(all(is.finite(plogis(big))))
expect_true(anyNA(exp(big) / (1 + exp(big))))          # the old formula did NaN

## Bug 3: cohort = NULL must be accepted (normalised to " ")
fit0 <- tryCatch(suppressWarnings(bamp(cases, pop, age = "rw1", period = "rw1",
          cohort = NULL, periods_per_agegroup = 1, mcmc.options = opts)),
          error = function(e) e)
expect_false(inherits(fit0, "error"))

## Bug 2: age-period (no cohort) fit AND predict, with finite rates
fitap <- suppressWarnings(bamp(cases, pop, age = "rw1", period = "rw1", cohort = " ",
          periods_per_agegroup = 1, mcmc.options = opts))
pp <- tryCatch(predict_apc(fitap, periods = 4,
          population = rbind(pop, pop[rep(P, 4), ]), quantiles = 0.5),
          error = function(e) e)
expect_false(inherits(pp, "error"))
if (!inherits(pp, "error")) {
  rt <- if (!is.null(pp$pr)) pp$pr else pp$rate
  expect_true(all(is.finite(rt)))
}

## Bug 4: non-integer population (e.g. PCLM-ungrouped exposure) must predict
fit <- suppressWarnings(bamp(cases, pop, age = "rw1", period = "rw1", cohort = "rw1",
          periods_per_agegroup = 1, mcmc.options = opts))
popNI <- rbind(pop, pop[rep(P, 4), ]) + 0.37           # non-integer exposure
pp4 <- tryCatch(predict_apc(fit, periods = 4, population = popNI, quantiles = 0.5),
          error = function(e) e)
expect_false(inherits(pp4, "error"))
if (!inherits(pp4, "error")) {
  rt4 <- if (!is.null(pp4$pr)) pp4$pr else pp4$rate
  expect_true(all(is.finite(rt4)))
}
