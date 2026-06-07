# Changelog

## bamp 2.1.7

- Effects (age, period, cohort) are now computed automatically inside
  [`bamp()`](https://volkerschmid.github.io/bamp/reference/bamp.md) and
  stored in the returned object (`model$effects`), so a separate call to
  [`effects.apc()`](https://volkerschmid.github.io/bamp/reference/effects.apc.md)
  is no longer needed for the default median summary.
- Fixed bug in
  [`effects.apc()`](https://volkerschmid.github.io/bamp/reference/effects.apc.md):
  cache check was looking for `x$effect` (singular) instead of
  `x$effects` (plural), and referenced an undefined variable in the
  cache condition.
- MCMC chains now use warm starts: on a restart the previous sample
  values are used as initial values, reducing burn-in time for automatic
  convergence checking. First-run behaviour is unchanged.

## bamp 2.1.6

- Fixed bug in cohort heterogeneity hyperparameter check
  (`cohort="rw2+het"` was not recognized correctly).
- Fixed typo in `period_covariate` handling that silently prevented
  vector coercion.
- Fixed
  [`checkConvergence()`](https://volkerschmid.github.io/bamp/reference/checkConvergence.md):
  cohort convergence check used age hyperparameter instead of cohort
  hyperparameter.
- Fixed hyperparameter display in
  [`checkConvergence()`](https://volkerschmid.github.io/bamp/reference/checkConvergence.md)
  info output.
- Fixed memory leaks in C++ MCMC code (raw pointer allocations per
  iteration).
- Fixed `GetRNGstate`/`PutRNGstate` pairing in random number generation.
- Added Cholesky decomposition error check with informative message.
- Replaced fork-based `mclapply` with socket-based
  `makeCluster`/`parLapply` for reliable parallel execution on all
  platforms, including macOS GUI environments (RStudio). Falls back to
  sequential if cluster setup fails.
- MCMC chains that fail due to numerical errors are now discarded
  cleanly instead of continuing with corrupted values.
- Cholesky decomposition now recovers from non-positive-definite
  precision matrices (common with RW2 priors during tuning) by adding
  increasing diagonal regularization and retrying, instead of aborting
  the chain.
- Added unit tests.
- Updated CITATION to use
  [`bibentry()`](https://rdrr.io/r/utils/bibentry.html).

## bamp 2.1.5

CRAN release: 2026-06-06

- bugfix RW2 by Chris Kypridemos
- Fixed population dimension handling in predict_apc.R by pietrosa

## bamp 2.1.3

CRAN release: 2022-10-18

- invalid UTF-8 in comment removed

## bamp 2.1.2

- Adapted to R 4.2

## bamp 2.1.1

CRAN release: 2022-05-05

- USE_FC_LEN_T for Fortran code

## bamp 2.1.0

CRAN release: 2021-06-10

- Better default settings (burn in times ten, more informative prior for
  age).
- Add warnings for failed convergence checks and removed chains in
  bamp(), including suggestions.
- Add warnings for failed convergence checks in print.apc().
- Fixed unwanted doubling of MCMC when verbose=2.

## bamp 2.0.8

CRAN release: 2020-01-23

- Minor bug fixes, fix “additional issues”.

## bamp 2.0.7

CRAN release: 2019-05-16

- Better initial setting for restarting iterations; helps with RW2
  priors.

## bamp 2.0.6

CRAN release: 2019-01-08

- Introductory vignette renamed (double vignette name warning from
  CRAN).

## bamp 2.0.5

CRAN release: 2018-12-12

- Removed ambiguities (mail Brian Ripley) and clean up in C code.

## bamp 2.0.4

CRAN release: 2018-11-30

- Added examples to all functions.

## bamp 2.0.3

- Added more details to help pages.

## bamp 2.0.2

- Reference in description changed.

## bamp 2.0.1

- Smaller vignettes.

## bamp 2.0.0

- R package.
