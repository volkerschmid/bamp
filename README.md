# Bayesian Age-period-cohort Modeling and Prediction

BAMP is a software package to analyze incidence or mortality data on the Lexis diagram, using a Bayesian version of an age-period-cohort model. Such models have been described in, e.g., [Berzuini and Clayton (1994)](https://doi.org/10.1002/sim.4780130804),  [Besag, J.E., P.J. Green, D.M. Higdon and K.L. Mengersen (1995)](https://doi.org/10.1214/ss/1177010123) and [Knorr-Held and Rainer (2001)](https://doi.org/10.1093/biostatistics/2.1.109). For each pixel in the Lexis diagram (that is for a specific age group and specific period) data must be available on the number of persons under risk (population number) and the number of disease cases (typically cancer incidence or mortality). A hierarchical model is assumed with a binomial model in the first-stage.

As smoothing priors for the age, period and cohort parameters random walks of first and second order (RW1 or RW2) available. BAMP also allows to drop one or more of the latent components, for example to drop the cohort effect and to analyze a age-period model. Additional unstructured prior distributions are assumed for each pixel in the Lexis diagram. Note that there is a nonidentifiability in the likelihood of the APC-model, see [Clayton and Schifflers (1987)](https://doi.org/10.1002/sim.4780060406), which induces some problems in interpreting the latent effects. The RW1 model is (weakly) identifiable; for the full RW2 model, `effects.apc()`/`plot.apc()` provide a `convention` argument that fixes the non-identified linear trend to a chosen display gauge, making the effect curves interpretable and reproducible between runs.

BAMP covers a range of models:

- AP and AC models,
- models with and without global heterogeneity parameter (overdispersion),
- models with additional age, period and/or cohort heterogeneity,
- models including covariates in the period or cohort effect,

The package includes features like

- data does not need to be on the same grid for age groups and periods, for example period can be in one year intervals and age group in five year intervals,
- prediction of future rates and number of cases,
- retrospective prediction for model checking,
- automatic model selection based on DIC (`selectModel()`).

Since version 3.0.0, BAMP uses a Polya-Gamma Gibbs sampler (`method = "pg"`, the default), a joint data-augmentation sampler with exact full conditionals and no Metropolis tuning; the legacy Taylor-expansion sampler remains available via `method = "taylor"`.

There are some graphical routines available in order to

- plot estimated age, period and cohort effects
- compare observed and fitted rates
- predict rates
- assess the pointwise credible intervals of the unstructured parameters. This helps to identify variation in the data, which is not supported by the age, period and cohort parameters.

## BAMP R package

[The bamp R package is available on CRAN.](https://CRAN.R-project.org/package=bamp)

<!-- badges: start -->
[![R-CMD-check](https://github.com/volkerschmid/bamp/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/volkerschmid/bamp/actions/workflows/check-standard.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/bamp)](https://CRAN.R-project.org/package=bamp)
[![CRAN_Download_Badge](https://cranlogs.r-pkg.org/badges/bamp)](https://cran.r-project.org/package=bamp)
[![CRAN_Download_Badge](https://cranlogs.r-pkg.org/badges/grand-total/bamp)](https://cran.r-project.org/package=bamp)
<!-- badges: end -->
