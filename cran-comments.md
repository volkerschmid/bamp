## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* macOS Tahoe 26.5.1, R 4.6.0 (aarch64-apple-darwin25.4.0, local)
* Windows (win-builder, R-devel)

## Summary of changes in 2.2.0

This release fixes several bugs and improves robustness:

* Age/period/cohort effects are now computed automatically inside `bamp()`
  and stored in the returned object, so a separate call to `effects.apc()`
  is no longer needed.
* MCMC chains now use warm starts on restart, reducing burn-in time for
  automatic convergence checking.
* Fixed cache check in `effects.apc()` (wrong field name `x$effect` vs
  `x$effects`).
* Fixed cohort heterogeneity hyperparameter check (`cohort="rw2+het"`
  was not recognized correctly).
* Fixed typo in `period_covariate` handling that silently prevented
  vector coercion.
* Fixed `checkConvergence()`: cohort convergence check used age
  hyperparameter instead of cohort hyperparameter.
* Fixed memory leaks in C++ MCMC code.
* Fixed `GetRNGstate`/`PutRNGstate` pairing in random number generation.
* Added Cholesky decomposition error check with informative message;
  non-positive-definite matrices are now handled by diagonal
  regularization and retry instead of aborting the chain.
* Replaced fork-based `mclapply` with socket-based `makeCluster`/
  `parLapply` for reliable parallel execution on all platforms,
  including macOS GUI environments (RStudio).
* MCMC chains that fail due to numerical errors are now discarded
  cleanly.
* Added unit tests.
* Updated CITATION to use `bibentry()`.
