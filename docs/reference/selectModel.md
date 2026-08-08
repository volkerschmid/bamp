# Automatic model selection for age-period-cohort models

Searches over age-period-cohort model specifications and returns the one
best supported by the data, by Deviance Information Criterion (DIC). It
answers the practical questions "is a first- or second-order random walk
more appropriate for each effect?", "are the data overdispersed?" and
(optionally) "is extra heterogeneity warranted?" without the user
fitting every model by hand.

## Usage

``` r
selectModel(
  cases,
  population,
  periods_per_agegroup,
  age = NULL,
  period = NULL,
  cohort = NULL,
  overdispersion = NULL,
  try_heterogeneity = FALSE,
  dic_margin = 4,
  psrf_tol = 1.1,
  screen = list(number_of_iterations = 10000, burn_in = 5000, step = 5, tuning = 200),
  final = "auto",
  refit = TRUE,
  hyperpar = list(age = c(1, 0.5), period = c(1, 5e-04), cohort = c(1, 5e-04), overdisp =
    c(1, 0.05), age_het = c(1, 0.05), period_het = c(1, 0.05), cohort_het = c(1, 0.05)),
  parallel = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- cases:

  number of cases (matrix, periods x age groups), as in
  [`bamp`](https://volkerschmid.github.io/bamp/reference/bamp.md).

- population:

  population number, as in
  [`bamp`](https://volkerschmid.github.io/bamp/reference/bamp.md).

- periods_per_agegroup:

  periods per age group.

- age, period, cohort:

  optional fixed value for an effect (`"rw1"`, `"rw2"`, `"rw1+het"`,
  `"rw2+het"` or `" "` for absent). If `NULL` (default) the effect is
  present and its random-walk order is searched over `"rw1"`/`"rw2"`.

- overdispersion:

  optional fixed logical. If `NULL` (default), whether to include
  overdispersion is part of the search.

- try_heterogeneity:

  logical; if `TRUE` the search may also add heterogeneity (`"+het"`) to
  an effect. Default `FALSE`.

- dic_margin:

  minimum DIC improvement required to adopt a more complex model
  (parsimony threshold). Default 4 (a conventional "clearly better" DIC
  difference).

- psrf_tol:

  convergence tolerance: a fit counts as converged if its maximum
  fitted-value Gelman-Rubin statistic is at or below this. Default 1.1.

- screen:

  list of MCMC settings used for the comparison fits, kept moderate for
  speed (default 10000 iterations, 5000 burn-in, step 5). Passed as
  `mcmc.options`. If no candidate converges at this length the search
  warns and selects by DIC only; increase these settings and re-run.

- final:

  list of MCMC settings used to refit the selected model, or `"auto"` to
  use the data-adaptive default of
  [`bamp`](https://volkerschmid.github.io/bamp/reference/bamp.md).

- refit:

  logical; if `TRUE` (default) refit the selected model with the `final`
  settings and return it; if `FALSE` return the screening fit of the
  selected model.

- hyperpar:

  hyper-parameter list passed to
  [`bamp`](https://volkerschmid.github.io/bamp/reference/bamp.md);
  defaults include heterogeneity hyper-parameters so `"+het"` models can
  be fitted.

- parallel:

  passed to
  [`bamp`](https://volkerschmid.github.io/bamp/reference/bamp.md)
  (chains run in parallel).

- verbose:

  logical; if `TRUE` (default) report progress and the running best
  model.

- ...:

  further arguments passed to
  [`bamp`](https://volkerschmid.github.io/bamp/reference/bamp.md) (e.g.
  `prior_scale`, `pg_engine`).

## Value

A list (class `"apcselect"`) with elements

- `table`: a data frame of every specification fitted, with its DIC,
  effective number of parameters `pD`, mean deviance, convergence flag,
  maximum fitted-value PSRF and fit time, ordered by DIC.

- `best`: the selected specification (a named list).

- `model`: the fitted
  [`apc`](https://volkerschmid.github.io/bamp/reference/apc.md) object
  for the selected specification (refitted with `final` settings if
  `refit = TRUE`).

- `path`: the sequence of specifications adopted by the greedy search.

## Details

The search is a **greedy forward selection by complexity**. It starts
from the simplest model – a first-order random walk (`"rw1"`) for every
effect that is present, no overdispersion – and at each round considers
every candidate that is exactly one step more complex than the current
best: an effect upgraded from `"rw1"` to `"rw2"`, overdispersion
switched on, or (if `try_heterogeneity = TRUE`) heterogeneity added to
an effect. All candidates in a round are fitted and the one with the
lowest DIC is adopted, but **only if** it improves DIC by at least
`dic_margin` *and* its chains converged; otherwise the search stops.
This costs a handful of fits rather than the full grid, follows an
interpretable path, and – through the margin – prefers the simpler model
unless the data clearly favour the more complex one. Each distinct
specification is fitted at most once (results are cached).

Model comparison uses DIC (lower is better), which rewards fit and
penalises effective complexity (`pD`); see
[`bamp`](https://volkerschmid.github.io/bamp/reference/bamp.md). A
specification that did not converge is never selected, because a low DIC
from a chain that has not mixed is not trustworthy – convergence is
judged via the same criterion as
[`checkConvergence`](https://volkerschmid.github.io/bamp/reference/checkConvergence.md):
the maximum Gelman-Rubin statistic over the smoothing precisions of the
effects present and the fitted log-odds across Lexis cells must be at or
below `psrf_tol`.

For speed and fairness all candidates are fitted with the same short
`screen` MCMC settings; the selected model is then optionally refitted
(`refit = TRUE`) with the longer `final` settings before being returned.
Fitting uses `method = "pg"`, which is robust on the sparse, rare-event
data where the legacy Taylor sampler can fail to converge.

Pin an axis to exclude it from the search by passing a fixed value: e.g.
`age = "rw2"` fixes the age effect (it is not searched), `age = " "`
removes the age effect entirely, and `overdispersion = FALSE` forbids
overdispersion. Any axis left `NULL` is searched.

## See also

[`bamp`](https://volkerschmid.github.io/bamp/reference/bamp.md),
[`checkConvergence`](https://volkerschmid.github.io/bamp/reference/checkConvergence.md)

## Examples

``` r
if (FALSE) { # \dontrun{
data(apc)
sel <- selectModel(cases, population, periods_per_agegroup = 5)
sel$table          # ranked comparison of the models tried
sel$best           # the chosen specification
plot(sel$model)    # the refitted best model
} # }
```
