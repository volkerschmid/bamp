# Effects from Fitted APC Model

Effects from Fitted APC Model

## Usage

``` r
# S3 method for class 'apc'
effects(object, mean = FALSE, quantiles = 0.5, update = FALSE, ...)
```

## Arguments

- object:

  an apc object

- mean:

  logical. If TRUE, mean effects are computed

- quantiles:

  Scalar or vector of quantiles to compute (only if mean=FALSE)

- update:

  logical. If TRUE, the apc object including the effects is returned

- ...:

  Additional arguments will be ignored

## Value

List of age, period, cohort effects or apc object including effects (if
update=TRUE)

## Examples

``` r
if (FALSE) { # \dontrun{
data(apc)
model <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", periods_per_agegroup = 5)
effects(model)
} # }
```
