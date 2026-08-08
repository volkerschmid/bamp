# Check apc object, whether MCMC has converged

This function uses Gelman and Rubin's R (potential scale reduction
factor) to check convergence. All checked quantities should have R\<1.1.
[`bamp`](https://volkerschmid.github.io/bamp/reference/bamp.md) runs at
least four MCMC chains by default (more if `parallel` is more than
four).

## Usage

``` r
checkConvergence(x, info = FALSE, level = 2, auto = FALSE)
```

## Arguments

- x:

  An apc object

- info:

  logical; print more information (including the raw per-effect
  diagnostic, which is affected by the age-period-cohort identifiability
  and should not be used on its own, see Details)

- level:

  level of check; 1 uses point estimate, 2 uses upper C.I.

- auto:

  logical; should be TRUE if called automatically from
  [`bamp`](https://volkerschmid.github.io/bamp/reference/bamp.md)

## Value

logical; TRUE if check is fine.

## Details

In an age-period-cohort model the age, period and cohort effects are
linearly dependent (Clayton and Schifflers, 1987): a linear trend can be
moved between the three effects without changing the likelihood. The
individual effect chains can therefore drift along this non-identified
direction even when the model has fully converged, which makes a naive
Gelman-R on the raw effects report spurious non-convergence.

`checkConvergence` therefore assesses the quantities that are actually
identified: the smoothing precisions and the fitted linear predictor
(log-odds) in every cell of the Lexis diagram, which is invariant to the
trend re-allocation. With `info=TRUE` the raw per-effect diagnostic is
also printed for reference.

## Examples

``` r
if (FALSE) { # \dontrun{
data(apc)
model <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", periods_per_agegroup = 5)
checkConvergence(model)
} # }
```
