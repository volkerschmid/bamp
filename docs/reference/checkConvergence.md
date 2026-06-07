# Check apc object, whether MCMC has converged

This functions uses Gelman and Rubin's R to check convergence for all
main parameters. All parameters should have R\<1.1.
[`bamp`](https://volkerschmid.github.io/bamp/reference/bamp.md) runs at
least four MCMC chains by default (more if parallel is more than four).

## Usage

``` r
checkConvergence(x, info = FALSE, level = 2, auto = FALSE)
```

## Arguments

- x:

  An apc object

- info:

  logical; print more information

- level:

  level of check; 1 uses point point estimation, 2 uses upper C.I.

- auto:

  logical; should be TRUE if called automatically from
  [`bamp`](https://volkerschmid.github.io/bamp/reference/bamp.md) \#'

## Value

logical; TRUE if check is fine.

## Examples

``` r
if (FALSE) { # \dontrun{
data(apc)
model <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", periods_per_agegroup = 5)
checkConvergence(model)
} # }
```
