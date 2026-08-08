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
year intervals.\
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
  mcmc.options = list(number_of_iterations = "auto", burn_in = "auto", step = "auto",
    tuning = 500),
  hyperpar = list(age = c(1, 0.5), period = c(1, 5e-04), cohort = c(1, 5e-04), overdisp =
    c(1, 0.05)),
  dic = TRUE,
  parallel = TRUE,
  verbose = FALSE,
  method = c("pg", "taylor"),
  prior_scale = FALSE,
  pg_engine = c("C", "R")
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

  - step: Step size, so only every step-th iteration is stored.

  - tuning: number of iterations for automatic tuning (used by
    `method="taylor"`). Depending on the model, the MCMC algorithm will
    tune certain parameters for more efficient MCMC chains. After
    tuning, the algorithm is restarted.

  Each of `number_of_iterations`, `burn_in` and `step` may be a number
  or the string `"auto"` (the default). `"auto"` chooses the value from
  the data: rare or zero-heavy counts (whose rare-event cells mix more
  slowly) get more iterations (from 40000 for well-populated data up to
  120000 when almost every cell is empty or has very few events),
  `burn_in` defaults to half the iterations, and `step` is set to keep
  about 1000 stored samples per chain. Any value given as a number is
  used exactly as supplied, so explicit settings reproduce the previous
  behaviour.

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

  should the chains be run in parallel. `TRUE`/`FALSE`, or a number
  giving the requested number of cores (capped at the number of chains).
  Uses the `parallel` package: forked workers
  ([`mclapply`](https://rdrr.io/r/parallel/mclapply.html)) on Unix and
  macOS, and – for `method = "pg"` – a PSOCK cluster on Windows (where
  forking is unavailable), so the default engine now runs in parallel on
  all platforms. (The legacy `method = "taylor"` engine still runs
  serially on Windows.) Parallel runs are reproducible: the per-chain
  seeds are drawn in the main process, so a given
  [`set.seed()`](https://rdrr.io/r/base/Random.html) yields the same
  result serially or in parallel.

- verbose:

  verbose mode

- method:

  MCMC engine. `"pg"` (default) is a joint sampler that combines
  Polya-Gamma data augmentation (Polson, Scott & Windle 2013) with a
  Laplace (Newton) Metropolis-Hastings refinement: each sweep draws the
  intercept and the age, period and cohort effects jointly in one exact
  Gibbs step and then refines them with a joint Newton proposal against
  the true binomial likelihood. It has no Metropolis tuning, never
  restarts on low acceptance and does not prune chains; it is markedly
  more robust for RW2 priors and converges the high-population,
  rare-event cells of incidence/mortality data that the Gibbs step alone
  mixes only slowly. It natively supports all of the package's models –
  RW1/RW2 priors, heterogeneity (`"rw1+het"`/`"rw2+het"`),
  overdispersion and period/cohort covariates. The Polya-Gamma weights
  use a normal approximation that is essentially exact for the large
  population counts of incidence/mortality data, so it typically needs
  far fewer iterations than the legacy sampler. `"taylor"` is the
  original block Metropolis-Hastings sampler with taylor expansion
  proposals (the default in versions 2.x); it remains available and can
  be faster on well-behaved (non rare-event) data, but it can fail to
  converge or prune all chains on sparse/zero-cell data.

- prior_scale:

  logical; only used by `method="pg"`. If `TRUE`, the intrinsic
  random-walk structure matrices are scaled to unit generalised variance
  (Sorbye & Rue 2014) so that a single hyper-prior is comparable across
  random-walk orders, grid sizes and data sets. The default is `FALSE`,
  which keeps the same prior parameterisation (and the same default
  hyper-parameters) as `method="taylor"`; if you set it to `TRUE` you
  should choose hyper-parameters appropriate for the scaled prior. See
  ‘Scaling the random-walk priors’ below for the rationale and benefits,
  and the examples for a short demonstration.

- pg_engine:

  implementation of the `method="pg"` sampler, one of `"C"` (default) or
  `"R"`. Both run the identical algorithm and, for a given seed, produce
  the same draws to floating-point tolerance; the `"C"` engine is a
  compiled port of the inner loop (no extra package dependency) and is
  roughly twice as fast. `"R"` is the readable reference implementation,
  kept for verification. Ignored for `method="taylor"`.

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

## Scaling the random-walk priors (`prior_scale`)

Each age, period and cohort effect has an intrinsic Gaussian
(random-walk) prior with precision (smoothing) parameter \\\kappa\\: the
effect vector \\x\\ has density proportional to
\\\exp(-\tfrac{1}{2}\kappa\\ x'Kx)\\, where \\K=D'D\\ is built from the
first- or second-order difference operator \\D\\. A
\\\mathrm{Gamma}(a,b)\\ hyper-prior is placed on \\\kappa\\. The
difficulty is that the smoothness implied by a given \\\kappa\\ is
governed not by \\\kappa\\ alone but by the marginal variance of the
effect, the generalised inverse of \\\kappa K\\; and the eigenvalues of
\\K\\ grow with the number of time points and with the random-walk
order. The *same* hyper-prior on \\\kappa\\ therefore implies very
different prior smoothness for, say, an RW1 over 10 periods and an RW2
over 50 cohorts. A hyper-prior tuned on one model silently means
something different on another, which is one reason a fixed default can
behave inconsistently across data sets.

With `prior_scale = TRUE` the structure matrix \\K\\ is rescaled so that
the geometric mean of the (generalised) marginal variances equals one
(Sorbye and Rue, 2014, DOI:10.1080/01621459.2013.866549). After scaling,
\\1/\sqrt{\kappa}\\ is, to a good approximation, the marginal standard
deviation of a typical effect element *on the log-odds (logit) scale*,
independently of the random-walk order, the number of age/period/cohort
points and the grid spacing.

Benefits: (i) **portable hyper-priors** – one \\\mathrm{Gamma}(a,b)\\
encodes the same smoothness belief across RW1/RW2 and across data sets
of different size; (ii) an **interpretable prior** – you can set
\\(a,b)\\ to express a belief about \\1/\sqrt{\kappa}\\ as a prior
effect standard deviation on the logit scale; (iii) **fairer model
comparison** (e.g. RW1 vs RW2 by DIC), because the prior is not
implicitly penalising one model far more than another. Scaling affects
only the smooth random-walk blocks; the i.i.d. heterogeneity components
and overdispersion already have an interpretable scale and are
unchanged.

The default is `prior_scale = FALSE` so that `method = "pg"` reproduces
the prior parameterisation (and default `hyperpar`) of the legacy
`method = "taylor"` engine. If you turn scaling on you should set
`hyperpar` for the scaled prior, where \\\kappa \approx
1/\mathrm{variance}\\; using the unscaled defaults with
`prior_scale = TRUE` would impose a different (and probably unintended)
amount of smoothing. Scaling is most worthwhile when fitting many models
or data sets and you want one coherent, interpretable prior across all
of them. The example below shows the effect concretely.

## See also

`vignette("modeling", package = "bamp")`

## Examples

``` r
if (FALSE) { # \dontrun{
data(apc)
model <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", periods_per_agegroup = 5)
} # }

## Demonstration of prior_scale (no MCMC, runs instantly): for a fixed
## precision kappa, report the geometric-mean prior marginal standard
## deviation of a random-walk effect on the logit scale, with and without
## Sorbye-Rue scaling, across random-walk orders and grid sizes.
prior_sd <- function(L, order, kappa = 1, scale = FALSE) {
  K <- crossprod(diff(diag(L), differences = order))   # structure matrix D'D
  if (scale) {                                          # Sorbye-Rue unit-variance scaling
    e <- eigen(K, symmetric = TRUE); keep <- e$values > max(e$values) * 1e-9
    V <- e$vectors[, keep, drop = FALSE]
    Sigma <- V %*% diag(1 / e$values[keep], sum(keep)) %*% t(V)
    K <- K * exp(mean(log(diag(Sigma))))
  }
  e <- eigen(K, symmetric = TRUE); keep <- e$values > max(e$values) * 1e-9
  V <- e$vectors[, keep, drop = FALSE]
  Sig <- V %*% diag(1 / e$values[keep], sum(keep)) %*% t(V) / kappa
  sqrt(exp(mean(log(diag(Sig)))))                       # geometric-mean marginal SD
}
grid <- expand.grid(order = 1:2, L = c(10, 25, 50))
data.frame(grid,
           unscaled = round(mapply(prior_sd, grid$L, grid$order, scale = FALSE), 3),
           scaled   = round(mapply(prior_sd, grid$L, grid$order, scale = TRUE), 3))
#>   order  L unscaled scaled
#> 1     1 10    1.224      1
#> 2     2 10    1.349      1
#> 3     1 25    1.943      1
#> 4     2 25    5.195      1
#> 5     1 50    2.749      1
#> 6     2 50   14.646      1
## With prior_scale = FALSE the same kappa = 1 implies an effect SD ranging
## from ~1.2 to ~14.6 across these models; with prior_scale = TRUE it is 1.0
## throughout, so a single hyper-prior on kappa means the same smoothness for
## every random-walk order and grid size.
```
