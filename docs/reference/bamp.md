# Bayesian Age-Period-Cohort Modeling and Prediction (bamp)

Bayesian Age-Period-Cohort Modeling for the analyze of incidence or
mortality data on the Lexis diagram. For each pixel in the Lexis diagram
(that is for a specific age group and specific period) data must be
available on the number of persons under risk (population number) and
the number of disease cases (typically cancer incidence or mortality). A
hierarchical model is assumed with a binomial model in the first-stage.
As smoothing priors for the age, period and cohort parameters random
walks of first and second order (RW1 or RW2) available. Deviance
information criterion and effective number of parameters is computed for
model comparison. Note that there is a non-identifiability in the
likelihood of the APC-model, see e.g. Clayton and Schifflers (1987,
DOI:10.1002/sim.4780060406), which indices some problems in interpreting
the latent effects. Only for RW1 model, the parameters are (weakly)
identifiable. Period and age groups do not need to be on the same grid,
for example periods can be in one year intervals and age groups in five
year intervals.  
Additionally to the model described in Knorr-Held and Rainer (2001,
DOI:10.1093/biostatistics/2.1.109), `bamp` can handle

- models with and without global heterogeneity parameter
  (overdispersion),

- models with additional age, period and/or cohort heterogeneity,

- additional covariates.

## Usage

``` r
bamp(
  cases,
  population,
  age,
  period,
  cohort,
  overdisp = FALSE,
  period_covariate = NULL,
  cohort_covariate = NULL,
  periods_per_agegroup,
  mcmc.options = list(number_of_iterations = 1e+05, burn_in = 50000, step = 50, tuning =
    500),
  hyperpar = list(age = c(1, 0.5), period = c(1, 5e-04), cohort = c(1, 5e-04), overdisp =
    c(1, 0.05)),
  dic = TRUE,
  parallel = TRUE,
  verbose = FALSE
)
```

## Arguments

- cases:

  number of cases

- population:

  population number

- age:

  prior for age groups ("rw1", "rw2", "rw1+het", "rw2+het", " ")

- period:

  prior for periods ("rw1", "rw2", "rw1+het", "rw2+het", " ")

- cohort:

  prior for cohorts ("rw1", "rw2", "rw1+het", "rw2+het", " ")

- overdisp:

  logical, add overdispersion to model

- period_covariate:

  covariate for period

- cohort_covariate:

  covariate for cohort

- periods_per_agegroup:

  periods per age group

- mcmc.options:

  list of options for MCMC.

  - burn_in: number of iterations used as burnin at the beginning of the
    algorithm, these iterations will be removed.

  - step: Step size, for example default is 50, so only every 50th
    iterations will be stored.

  - tuning: number of iterations for automatic tuning. Depending on the
    model, the MCMC algorithm will tune certain parameters for more
    efficient MCMC chains. After tuning, the algorithm is restarted.

- hyperpar:

  list of hyper parameters. The hyper prior for the precision (inverse
  variance) in the random walk priors is a Gamma distribution with
  parameters \\a\\ and \\b\\; expected value is \\a/b\\, variance is
  \\a/b^2\\. Weak hyper parameters are suggested, defaults are \\a=1,
  b=0.5\\ for age, \\a=1, b=0.0005\\ for period and cohort effects and
  \\a=1, b=0.05\\ for overdispersion (if added). It is recommended to
  choose the hyper priors depending on the model, in particular on the
  order of the random walk.

- dic:

  logical. If true. DIC will be computed

- parallel:

  logical, should computation be done in parallel. This uses the
  parallel package, which does not allow parallel computing under
  Windows.

- verbose:

  verbose mode

## Details

This functions returns an
[`apc`](https://volkerschmid.github.io/bamp/reference/apc.md) object.
Only samples from the posterior are computed, point estimates and
credible intervals will be computed in
[`effects.apc`](https://volkerschmid.github.io/bamp/reference/effects.apc.md),
[`print.apc`](https://volkerschmid.github.io/bamp/reference/print.apc.md)
and
[`plot.apc`](https://volkerschmid.github.io/bamp/reference/plot.apc.md).
[`predict_apc`](https://volkerschmid.github.io/bamp/reference/predict_apc.md)
can be used for for prediction of the future rates and number of cases
and for a retrospective prediction for model checking.

## See also

`vignette("modeling", package = "bamp")`

## Examples

``` r
if (FALSE) { # \dontrun{
data(apc)
model <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", periods_per_agegroup = 5)
} # }
```
