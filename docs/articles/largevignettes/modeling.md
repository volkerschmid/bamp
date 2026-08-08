# Bayesian Age-Period-Cohort Modeling

This vignette walks through the range of age-period-cohort (APC) models
that [`bamp()`](https://volkerschmid.github.io/bamp/reference/bamp.md)
can fit: different smoothness priors for the age, period and cohort
effects, a model that omits an effect entirely, and models where the
period or cohort effect is scaled by a known covariate. For a first
introduction to
[`bamp()`](https://volkerschmid.github.io/bamp/reference/bamp.md), see
*vignette(“bamp”, package=“bamp”)*; for forecasting, see
*vignette(“prediction”, package=“bamp”)*.

## Data example

BAMP includes a data example: `cases` and `population`, matrices of
observed case counts and population at risk, cross-classified by
calendar year (rows) and age group (columns). It also bundles `cov_p`
and `cov_c`, example period and cohort covariates used in the last two
models below.

[`data`](https://rdrr.io/r/utils/data.html)`(``apc``)`` `[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(``cases``[``,``1``]``,type``=``"l"``,ylim``=`[`range`](https://rdrr.io/r/base/range.html)`(``cases``)``, ylab``=``"cases"``, xlab``=``"year"``, main``=``"cases per age group"``)`` ``for`` ``(``i`` ``in`` ``2``:``8``)`[`lines`](https://rdrr.io/r/graphics/lines.html)`(``cases``[``,``i``]``, col``=``i``)`

![](modeling_files/figure-html/loadplot-1.png)

## APC model with random-walk first-order priors

A first-order random walk (`"rw1"`) prior assumes each effect changes
smoothly, without a systematic trend, between neighbouring age groups,
periods or cohorts. This is the most common default.

`model1`` ``<-`` `[`bamp`](https://volkerschmid.github.io/bamp/reference/bamp.md)`(``cases``, ``population``, age``=``"rw1"``, period``=``"rw1"``, cohort``=``"rw1"``,`` `` periods_per_agegroup ``=`` ``5``)`

[`bamp()`](https://volkerschmid.github.io/bamp/reference/bamp.md)
automatically checks MCMC convergence using the Gelman-Rubin diagnostic.
We can also run this check manually:

[`checkConvergence`](https://volkerschmid.github.io/bamp/reference/checkConvergence.md)`(``model1``)`

    ## [1] TRUE

[`print()`](https://volkerschmid.github.io/bamp/reference/print.apc.md)
summarises the fitted model: posterior estimates of the smoothing
(precision) parameters, and the deviance and DIC, which are useful for
comparing the models fitted in this vignette:

[`print`](https://volkerschmid.github.io/bamp/reference/print.apc.md)`(``model1``)`

    ## 
    ##  Model:
    ## age (rw1)  - period (rw1)  - cohort (rw1) model
    ## Deviance:     230.12
    ## pD:            36.60
    ## DIC:          266.72
    ## 
    ## 
    ##  Hyper parameters:                 5%           50%          95%         
    ## age                              0.436        1.108        2.295
    ## period                          48.842      128.964      320.787
    ## cohort                          35.190       61.080      103.190
    ## 
    ## 
    ## Markov Chains convergence checked succesfully using Gelman's R (potential scale reduction factor).

The fitted age, period and cohort effects can be plotted with point-wise
posterior quantiles:

[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(``model1``)`

![](modeling_files/figure-html/plot_model-1.png)![](modeling_files/figure-html/plot_model-2.png)![](modeling_files/figure-html/plot_model-3.png)

By default the plot shows the median and a 90% interval; other quantiles
can be requested explicitly:

[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(``model1``, quantiles ``=`` `[`c`](https://rdrr.io/r/base/c.html)`(``0.025``,``0.1``,``0.5``,``0.9``,``0.975``)``)`

![](modeling_files/figure-html/plot_model_with_more_quantiles-1.png)![](modeling_files/figure-html/plot_model_with_more_quantiles-2.png)![](modeling_files/figure-html/plot_model_with_more_quantiles-3.png)

## APC model with random-walk second-order priors

A second-order random walk (`"rw2"`) prior penalises curvature rather
than level, so it produces effects that follow a smooth, roughly linear
trend instead of reverting towards a constant. RW2 priors are less
informative than RW1, so the chains mix more slowly and need
substantially more iterations to converge; here we also set the
smoothing hyperparameters explicitly instead of using the defaults.

`model2`` ``<-`` `[`bamp`](https://volkerschmid.github.io/bamp/reference/bamp.md)`(``cases``, ``population``, age``=``"rw2"``, period``=``"rw2"``, cohort``=``"rw2"``,`` `` periods_per_agegroup ``=`` ``5``,`` `` mcmc.options``=`[`list`](https://rdrr.io/r/base/list.html)`(``"number_of_iterations"``=``200000``, ``"burn_in"``=``100000``, ``"step"``=``50``, ``"tuning"``=``500``)``,`` `` hyperpar``=`[`list`](https://rdrr.io/r/base/list.html)`(``"age"``=`[`c`](https://rdrr.io/r/base/c.html)`(``1``,``.5``)``, ``"period"``=`[`c`](https://rdrr.io/r/base/c.html)`(``1``,``0.05``)``, ``"cohort"``=`[`c`](https://rdrr.io/r/base/c.html)`(``1``,``0.05``)``)``)`

[`checkConvergence`](https://volkerschmid.github.io/bamp/reference/checkConvergence.md)`(``model2``)`

    ## [1] TRUE

[`print`](https://volkerschmid.github.io/bamp/reference/print.apc.md)`(``model2``)`

    ## 
    ##  Model:
    ## age (rw2)  - period (rw2)  - cohort (rw2) model
    ## Deviance:     233.58
    ## pD:            35.72
    ## DIC:          269.30
    ## 
    ## 
    ##  Hyper parameters:                 5%           50%          95%         
    ## age                              1.090        3.024        6.847
    ## period                          16.052       41.275       91.027
    ## cohort                          24.825       46.878       83.398
    ## 
    ## 
    ## Markov Chains convergence checked succesfully using Gelman's R (potential scale reduction factor).

[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(``model2``)`

![](modeling_files/figure-html/Model2_results-1.png)![](modeling_files/figure-html/Model2_results-2.png)![](modeling_files/figure-html/Model2_results-3.png)

## APC model without a period effect

Age, period and cohort are linearly dependent – any one is determined by
the other two – so a model with all three as unrestricted smooth trends
is not identifiable. One common way to sidestep this is to drop the
effect of least interest – here the period effect (`period=" "`) – and
keep only age and cohort:

`model3``<-`[`bamp`](https://volkerschmid.github.io/bamp/reference/bamp.md)`(``cases``, ``population``, age``=``"rw1"``, period``=``" "``, cohort``=``"rw2"``,`` `` periods_per_agegroup ``=`` ``5``)`

[`checkConvergence`](https://volkerschmid.github.io/bamp/reference/checkConvergence.md)`(``model3``)`

    ## [1] TRUE

[`print`](https://volkerschmid.github.io/bamp/reference/print.apc.md)`(``model3``)`

    ## 
    ##  Model:
    ## age (rw1) cohort (rw2) model
    ## Deviance:     275.82
    ## pD:            29.69
    ## DIC:          305.52
    ## 
    ## 
    ##  Hyper parameters:                 5%           50%          95%         
    ## age                              0.299        0.724        1.497
    ## cohort                          36.892       73.124      139.543
    ## 
    ## 
    ## Markov Chains convergence checked succesfully using Gelman's R (potential scale reduction factor).

[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(``model3``)`

![](modeling_files/figure-html/model3-1.png)![](modeling_files/figure-html/model3-2.png)

## APC model with a cohort covariate

Instead of dropping an effect, a known covariate can scale it:
`cohort_covariate` multiplies the cohort effect by a given value at each
cohort (here the bundled example covariate `cov_c`), which is useful
when an external, cohort-specific risk factor is available.

`(``model4``<-`[`bamp`](https://volkerschmid.github.io/bamp/reference/bamp.md)`(``cases``, ``population``, age``=``"rw1"``, period``=``"rw1"``, cohort``=``"rw1"``,`` `` cohort_covariate ``=`` ``cov_c``, periods_per_agegroup ``=`` ``5``)``)`

    ## 
    ##  Model:
    ## age (rw1)  - period (rw1)  - cohort (rw1) model
    ## Deviance:     230.43
    ## pD:            37.39
    ## DIC:          267.81
    ## 
    ## 
    ##  Hyper parameters:                 5%           50%          95%         
    ## age                              0.408        1.021        2.101
    ## period                          55.546      142.037      349.491
    ## cohort                          31.623       54.483       92.487
    ## 
    ## 
    ## Markov Chains convergence checked succesfully using Gelman's R (potential scale reduction factor).

[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(``model4``)`

![](modeling_files/figure-html/model4-1.png)![](modeling_files/figure-html/model4-2.png)![](modeling_files/figure-html/model4-3.png)![](modeling_files/figure-html/model4-4.png)![](modeling_files/figure-html/model4-5.png)![](modeling_files/figure-html/model4-6.png)

## APC model with a period covariate

Analogously, `period_covariate` scales the period effect by a known,
period-specific covariate (here `cov_p`):

`(``model5``<-`[`bamp`](https://volkerschmid.github.io/bamp/reference/bamp.md)`(``cases``, ``population``, age``=``"rw1"``, period``=``"rw1"``, cohort``=``"rw1"``,`` `` period_covariate ``=`` ``cov_p``, periods_per_agegroup ``=`` ``5``)``)`

    ## 
    ##  Model:
    ## age (rw1)  - period (rw1)  - cohort (rw1) model
    ## Deviance:     229.64
    ## pD:            36.71
    ## DIC:          266.35
    ## 
    ## 
    ##  Hyper parameters:                 5%           50%          95%         
    ## age                              0.442        1.098        2.226
    ## period                          54.273      147.975      384.812
    ## cohort                          34.560       60.147       98.309
    ## 
    ## 
    ## Markov Chains convergence checked succesfully using Gelman's R (potential scale reduction factor).

[`plot`](https://rdrr.io/r/graphics/plot.default.html)`(``model5``)`

![](modeling_files/figure-html/model5-1.png)![](modeling_files/figure-html/model5-2.png)![](modeling_files/figure-html/model5-3.png)![](modeling_files/figure-html/model5-4.png)![](modeling_files/figure-html/model5-5.png)![](modeling_files/figure-html/model5-6.png)
