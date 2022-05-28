<img src="https://volkerschmid.github.io/bamp/figures/bamplogo8.png" align="right" />

# Bayesian Age-period-cohort Modeling and Prediction

BAMP is a software package to analyze incidence or mortality data on the Lexis diagram, using a Bayesian version of an age-period-cohort model. Such models have been described in, e.g., [Berzuini and Clayton (1994)](https://doi.org/10.1002/sim.4780130804),  [Besag, J.E., P.J. Green, D.M. Higdon and K.L. Mengersen (1995)](https://doi.org/10.1214/ss/1177010123) and [Knorr-Held and Rainer (2001)](https://doi.org/10.1093/biostatistics/2.1.109). For each pixel in the Lexis diagram (that  is for a specific age group and specific period) data must be available on the number of persons under risk (population number) and the number of disease cases (typically cancer incidence or mortality). A hierarchical model is assumed with a binomial model in the first-stage.

As smoothing priors for the age, period and cohort parameters random walks of first and second order (RW1 or RW2) available. BAMP also allows to drop one or more of the latent components, for example to drop the cohort effect and to analyze a age-period model. Additional unstructured prior distributions are assumed for each pixel in the Lexis diagram. Note that there is a nonidentifiability in the likelihood of the APC-model, see [Clayton and Schifflers (1987)](https://doi.org/10.1002/sim.4780060406), which indices some problems in interpreting the latent effects. Only for RW1 model, the parameters are (weakly) identifiable.

BAMP has several features which are described more detailed in [Knorr-Held and Rainer (2001)](https://doi.org/10.1093/biostatistics/2.1.109):

- The data does not need to be on the same grid, for example period can be in one year intervals and age group in five year intervals.
- BAMP allows for prediction of the future number of cases
- BAMP allows for a retrospective prediction for model checking

Additionally to the model described in [Knorr-Held and Rainer (2001)](https://doi.org/10.1093/biostatistics/2.1.109), BAMP can handle
- AP and AC models
- models with and without global heterogenity parameter (overdispersion)
- models with additional age, period and/or cohort heterogenity
- including covariates (still in development)
Detail about this feature can be found in [Schmid (2004 - in German)](https://edoc.ub.uni-muenchen.de/3000/)

There are some graphical routines available in order to

- plot estimated age, period and cohort effects (only for RW1 model)
- compare observed and fitted rates
- predict rates
- assess the "significance" of the unstructured parameters. This helps  to identify variation in the data, which is not supported by the age, period and cohort parameters. 

## BAMP old version (1.3)

[Find the older standalone version here.](https://volkerschmid.github.io/bamp/articles/standaloneversion/)

## BAMP R package (2.1)

[The bamp R package is available on CRAN.](https://CRAN.R-project.org/package=bamp)

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/bamp)](https://CRAN.R-project.org/package=bamp)
[![CRAN_Download_Badge](https://cranlogs.r-pkg.org/badges/bamp)](https://cran.r-project.org/package=bamp)
[![CRAN_Download_Badge](https://cranlogs.r-pkg.org/badges/grand-total/bamp)](https://cran.r-project.org/package=bamp)
<!-- badges: end -->
