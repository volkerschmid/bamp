# Plot apc object

Plot apc object

## Usage

``` r
# S3 method for class 'apc'
plot(
  x,
  quantiles = c(0.05, 0.5, 0.95),
  convention = c("age", "period", "cohort", "none"),
  combined = FALSE,
  ...
)
```

## Arguments

- x:

  apc object

- quantiles:

  quantiles to plot. Default: `c(0.05,0.5,0.95)` is median and 90%
  credible interval.

- convention:

  display-layer gauge convention for the linear trend (drift) in a full
  age-period-cohort model; one of `"age"` (default), `"period"`,
  `"cohort"` or `"none"`. In a full APC model the three effects are
  identifiable only up to a shared linear trend; this fixes that single
  degree of freedom, removing the run-to-run drift that is the dominant
  source of non-reproducibility in the plotted curves (residual
  curvature and Monte-Carlo noise remain). `"age"` shows the age effect
  as curvature about a zero trend and puts the drift in the period and
  cohort effects; `"none"` plots the raw, un-gauged effects. It is
  display-only (the fitted rates and predictions are unchanged) and is
  ignored for models that are not full APC. See
  [`effects.apc`](https://volkerschmid.github.io/bamp/reference/effects.apc.md).

- combined:

  logical. For heterogeneity models, if `TRUE` plot the full effect
  (smooth + iid heterogeneity); default `FALSE` plots the smooth
  component only. Ignored for models without heterogeneity.

- ...:

  Additional arguments will be ignored

## Value

plot

## Details

Plot of age, period and cohort effects from apc objects. If covariates
have been used for period/cohort, a second plot with covariate, absolute
effect and relative effect is created. Absolute effect is relative
effect times covariate.

## Examples

``` r
if (FALSE) { # \dontrun{
data(apc)
model <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", periods_per_agegroup = 5)
plot(model)
} # }
```
