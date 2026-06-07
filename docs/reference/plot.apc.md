# Plot apc object

Plot apc object

## Usage

``` r
# S3 method for class 'apc'
plot(x, quantiles = c(0.05, 0.5, 0.95), ...)
```

## Arguments

- x:

  apc object

- quantiles:

  quantiles to plot. Default: `c(0.05,0.5,0.95)` is median and 90%
  credible interval.

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
