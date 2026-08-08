## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* macOS Tahoe 26.5.2, R 4.6.0 (aarch64-apple-darwin23, local)
* Windows (win-builder, R-devel)

## Summary of changes in 3.0.0

This is a major release, adding a new default MCMC engine and several bug fixes:

* New Polya-Gamma Gibbs sampler (`method = "pg"`, now the default): a joint
  data-augmentation sampler with exact full conditionals and no Metropolis
  tuning. Supports overdispersion, age/period/cohort heterogeneity and
  period/cohort covariates natively. The legacy Taylor expansion proposal sampler remains
  available via `method = "taylor"`.
* Native C implementation of the Polya-Gamma sampler (`pg_engine = "C"`,
  the default) with an equivalent pure-R reference (`pg_engine = "R"`).
* `mcmc.options` values `number_of_iterations`, `burn_in` and `step` may now
  be set to `"auto"` (the default), choosing the MCMC length from the
  rarity of the data.
* New `selectModel()`: automatic APC model selection by DIC.
* New `prior_scale` argument for `method = "pg"`.
* `checkConvergence()` now assesses identified quantities (smoothing
  precisions, intercept, fitted linear predictor per Lexis cell) instead of
  the raw effect chains, which drift along the non-identified
  age-period-cohort trend.
* `effects.apc()` and `plot.apc()` gain a `convention` argument fixing the
  non-identified linear trend to a chosen display gauge, making the effect
  curves reproducible between runs; `plot.apc()` also handles any number of
  quantiles and zero-covariate models.
* Fixed `predict_apc(periods = 0)` crash (downward-sequence off-by-one).
* Fixed `predict_apc()` logit overflow (`plogis()` instead of
  `exp(x)/(1+exp(x))`, which returned `NaN` for large forecast logits).
* Fixed `bamp()` to accept `age`/`period`/`cohort = NULL`.
* Fixed `predict_apc()` for age-period (no-cohort) models and non-integer
  population/exposure.
* `bamp(..., method = "pg")` now surfaces the real error message from a
  failed MCMC chain instead of an opaque "subscript out of bounds".

## Downstream dependencies

There are no reverse dependencies on CRAN.
