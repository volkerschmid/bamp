# Effects from Fitted APC Model

Effects from Fitted APC Model

## Usage

``` r
# S3 method for class 'apc'
effects(
  object,
  mean = FALSE,
  quantiles = 0.5,
  update = FALSE,
  convention = c("age", "period", "cohort", "none"),
  combined = FALSE,
  ...
)
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

- convention:

  display-layer gauge convention for the linear trend (drift) in a full
  age-period-cohort model; one of `"age"` (default), `"period"`,
  `"cohort"` or `"none"`. See Details.

- combined:

  logical. For heterogeneity models (`"rw1+het"` / `"rw2+het"`), if
  `TRUE` the returned effect is the full effect (smooth + iid
  heterogeneity component); if `FALSE` (default) only the smooth
  component is returned. Ignored for models without heterogeneity.

- ...:

  Additional arguments will be ignored

## Value

List of age, period, cohort effects or apc object including effects (if
update=TRUE)

## Details

In a full age-period-cohort model the age, period and cohort effects are
identifiable only up to one shared linear trend (drift), because a
linear trend can be moved between the three effects without changing the
fitted rates (Clayton and Schifflers, 1987). When the raw posterior
samples are summarised directly, this non-identified trend makes the
effect curves drift between MCMC iterations and between runs, so they
are not reproducible.

`convention` fixes that single degree of freedom for display: `"age"`
removes the linear slope of the age effect (age is shown as curvature
about a zero trend, the drift is shown in the period and cohort
effects); `"period"` and `"cohort"` pin the corresponding effect's slope
to zero instead; `"none"` returns the raw, un-gauged effects (previous
behaviour, curves may differ between runs). The gauge is applied per
MCMC draw before the quantiles are computed, and only for full APC
models – for models without all three effects there is no trend aliasing
and the argument is ignored.

Fixing the gauge removes the run-to-run linear-trend (drift) component,
which is typically the dominant source of non-reproducibility in the
effect curves; the residual curvature and Monte-Carlo sampling noise are
unaffected, so two independent runs are made much closer but need not
agree exactly. The zero-slope property holds for every individual MCMC
draw; because quantiles are non-linear, the summarised median/quantile
curve has only approximately (not exactly) zero slope.

The gauge is display-only: it never modifies the stored samples
(`object$samples`), and the fitted rates, the predictions from
[`predict_apc`](https://volkerschmid.github.io/bamp/reference/predict_apc.md)
and the DIC are invariant to it – only the way the common linear trend
is split among the three curves changes. The convention actually used is
recorded in `attr(result, "gauge_convention")`.

## Examples

``` r
if (FALSE) { # \dontrun{
data(apc)
model <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", periods_per_agegroup = 5)
effects(model)
} # }
```
