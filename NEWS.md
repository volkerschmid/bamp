# bamp 3.0.0

Major release. The Polya-Gamma sampler and its native C engine, DIC-based `selectModel()`, and the identifiability-aware convergence diagnostics in this release were contributed by Chris Kypridemos.

* New Polya-Gamma Gibbs sampler (`method = "pg"`, now the default): a joint
  Polya-Gamma data-augmentation sampler with exact full conditionals and no
  Metropolis tuning. Supports overdispersion, age/period/cohort heterogeneity
  and period/cohort covariates natively. The legacy Taylor sampler remains
  available via `method = "taylor"`.
* Native C implementation of the Polya-Gamma sampler (`pg_engine = "C"`,
  the default) with an equivalent pure-R reference (`pg_engine = "R"`);
  the two agree to numerical tolerance.
* `mcmc.options` values `number_of_iterations`, `burn_in` and `step` may now be
  set to `"auto"` (the default), which chooses the MCMC length from the rarity
  of the data. Any value given as a number is used exactly as before.
* New `selectModel()`: automatic APC model selection by DIC.
* New `prior_scale` argument for `method = "pg"`.
* `checkConvergence()` now assesses the identified quantities (smoothing
  precisions, intercept and the fitted linear predictor per Lexis cell), which
  are invariant to the age-period-cohort trend aliasing, rather than the raw
  effect chains that drift along the non-identified trend.
* `effects.apc()` and `plot.apc()` gain a `convention` argument that fixes the
  non-identified linear trend to a chosen display gauge, making the effect
  curves reproducible between runs; `plot.apc()` also handles any number of
  quantiles and zero-covariate models.
* `print.apc()` now also reports the intercept (5%/50%/95% quantiles).
* Fixed `predict_apc(periods = 0)` crash (downward-sequence off-by-one).
* Fixed `predict_apc()` logit overflow: use `plogis()` instead of
  `exp(x)/(1+exp(x))`, which returned `NaN` for large forecast logits and
  crashed downstream.
* Fixed `bamp()`: `age`/`period`/`cohort = NULL` is now accepted (previously
  errored with "argument is of length zero").
* Fixed `predict_apc()`: age-period models without a cohort effect can now be
  predicted, and non-integer population/exposure no longer produces `NA`.
* `bamp(..., method = "pg")`: a chain that fails under forked parallelism
  (`parallel::mclapply`) now reports its actual error instead of the opaque
  "subscript out of bounds" that resulted from silently indexing into the
  failed chain's result.
* Fixed `bamp()`'s convergence warning being tied to `verbose` backwards: it
  now prints the "did not converge" message when `verbose = TRUE` and stays
  silent when `verbose = FALSE`, as intended.

# bamp 2.2.0

* Effects (age, period, cohort) are now computed automatically inside `bamp()` and stored in the returned object (`model$effects`), so a separate call to `effects.apc()` is no longer needed for the default median summary.
* Fixed bug in `effects.apc()`: cache check was looking for `x$effect` (singular) instead of `x$effects` (plural), and referenced an undefined variable in the cache condition.
* MCMC chains now use warm starts: on a restart the previous sample values are used as initial values, reducing burn-in time for automatic convergence checking. First-run behaviour is unchanged.
* Fixed bug in cohort heterogeneity hyperparameter check (`cohort="rw2+het"` was not recognized correctly).
* Fixed typo in `period_covariate` handling that silently prevented vector coercion.
* Fixed `checkConvergence()`: cohort convergence check used age hyperparameter instead of cohort hyperparameter.
* Fixed hyperparameter display in `checkConvergence()` info output.
* Fixed memory leaks in C++ MCMC code (raw pointer allocations per iteration).
* Fixed `GetRNGstate`/`PutRNGstate` pairing in random number generation.
* Added Cholesky decomposition error check with informative message.
* Replaced fork-based `mclapply` with socket-based `makeCluster`/`parLapply` for reliable parallel execution on all platforms, including macOS GUI environments (RStudio). Falls back to sequential if cluster setup fails.
* MCMC chains that fail due to numerical errors are now discarded cleanly instead of continuing with corrupted values.
* Cholesky decomposition now recovers from non-positive-definite precision matrices (common with RW2 priors during tuning) by adding increasing diagonal regularization and retrying, instead of aborting the chain.
* Added unit tests.
* Updated CITATION to use `bibentry()`.

# bamp 2.1.5
* bugfix RW2 by Chris Kypridemos
* Fixed population dimension handling in predict_apc.R by pietrosa

# bamp 2.1.3
* invalid UTF-8 in comment removed

# bamp 2.1.2
* Adapted to R 4.2

# bamp 2.1.1
* USE_FC_LEN_T for Fortran code

# bamp 2.1.0
* Better default settings (burn in times ten, more informative prior for age).
* Add warnings for failed convergence checks and removed chains in bamp(), including suggestions.
* Add warnings for failed convergence checks in print.apc().
* Fixed unwanted doubling of MCMC when verbose=2.

# bamp 2.0.8
* Minor bug fixes, fix "additional issues".

# bamp 2.0.7
* Better initial setting for restarting iterations; helps with RW2 priors.

# bamp 2.0.6
* Introductory vignette renamed (double vignette name warning from CRAN).

# bamp 2.0.5
* Removed ambiguities (mail Brian Ripley) and clean up in C code.

# bamp 2.0.4
* Added examples to all functions.

# bamp 2.0.3
* Added more details to help pages.

# bamp 2.0.2
* Reference in description changed.

# bamp 2.0.1
* Smaller vignettes.

# bamp 2.0.0
* R package.
