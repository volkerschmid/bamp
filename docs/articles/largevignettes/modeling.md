# Bayesian Age-Period-Cohort Modeling

## Data example

BAMP includes a data example.

``` r

data(apc)
plot(cases[,1],type="l",ylim=range(cases), ylab="cases", xlab="year", main="cases per age group")
for (i in 2:8)lines(cases[,i], col=i)
```

![](modeling_files/figure-html/loadplot-1.png)

## APC model with random walk first order prior

``` r

model1 <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1",
              periods_per_agegroup = 5)
```

bamp() automatically performs a check for MCMC convergence using Gelman
and Rubin’s convergence diagnostic. We can manually check the
convergence again:

``` r

checkConvergence(model1)
```

    ## [1] TRUE

Now we have a look at the model results. This includes estimates of
smoothing parameters and deviance and DIC:

``` r

print(model1)
```

    ## 
    ##  Model:
    ## age (rw1)  - period (rw1)  - cohort (rw1) model
    ## Deviance:     231.44
    ## pD:            36.91
    ## DIC:          268.35
    ## 
    ## 
    ##  Hyper parameters:                 5%           50%          95%         
    ## age                              0.353        0.939        2.008
    ## period                          68.286      197.863      618.072
    ## cohort                          34.381       59.508       96.741
    ## 
    ## 
    ## Markov Chains convergence checked succesfully using Gelman's R (potential scale reduction factor).

We can plot the main APC effects using point-wise quantiles:

``` r

plot(model1)
```

![](modeling_files/figure-html/plot_model-1.png)![](modeling_files/figure-html/plot_model-2.png)![](modeling_files/figure-html/plot_model-3.png)

More quantiles are possible:

``` r

plot(model1, quantiles = c(0.025,0.1,0.5,0.9,0.975))
```

![](modeling_files/figure-html/plot_model_with_more_quantiles-1.png)![](modeling_files/figure-html/plot_model_with_more_quantiles-2.png)![](modeling_files/figure-html/plot_model_with_more_quantiles-3.png)

## APC model with random walk second order prior

``` r

model2 <- bamp(cases, population, age="rw2", period="rw2", cohort="rw2",
              periods_per_agegroup = 5,
              mcmc.options=list("number_of_iterations"=200000, "burn_in"=100000, "step"=50, "tuning"=500),
              hyperpar=list("age"=c(1,.5), "period"=c(1,0.05), "cohort"=c(1,0.05)))
```

``` r

checkConvergence(model2)
```

    ## [1] TRUE

``` r

print(model2)
```

    ## 
    ##  Model:
    ## age (rw2)  - period (rw2)  - cohort (rw2) model
    ## Deviance:     234.27
    ## pD:            36.77
    ## DIC:          271.04
    ## 
    ## 
    ##  Hyper parameters:                 5%           50%          95%         
    ## age                              1.054        2.863        6.478
    ## period                          16.626       41.725       89.427
    ## cohort                          23.759       45.271       82.147
    ## 
    ## 
    ## Markov Chains convergence checked succesfully using Gelman's R (potential scale reduction factor).

``` r

plot(model2)
```

![](modeling_files/figure-html/Model2_results-1.png)![](modeling_files/figure-html/Model2_results-2.png)![](modeling_files/figure-html/Model2_results-3.png)

``` r

model3<-bamp(cases, population, age="rw1", period=" ", cohort="rw2",
              periods_per_agegroup = 5)
```

    ## hain. Please check for convergence using checkConvergence() and maybe change your model settings (maybe add overdispersion).

``` r

checkConvergence(model3)
```

    ## [1] TRUE

``` r

print(model3)
```

    ## 
    ##  Model:
    ## age (rw1) cohort (rw2) model
    ## Deviance:     276.65
    ## pD:            30.23
    ## DIC:          306.88
    ## 
    ## 
    ##  Hyper parameters:                 5%           50%          95%         
    ## age                              0.278        0.737        1.469
    ## cohort                          38.698       74.630      140.661
    ## 
    ## 
    ## Markov Chains convergence checked succesfully using Gelman's R (potential scale reduction factor).

``` r

plot(model3)
```

![](modeling_files/figure-html/model3-1.png)![](modeling_files/figure-html/model3-2.png)

``` r

(model4<-bamp(cases, population, age="rw1", period="rw1", cohort="rw1",
             cohort_covariate = cov_c, periods_per_agegroup = 5))
plot(model4)
```

![](modeling_files/figure-html/model4-1.png)![](modeling_files/figure-html/model4-2.png)![](modeling_files/figure-html/model4-3.png)![](modeling_files/figure-html/model4-4.png)![](modeling_files/figure-html/model4-5.png)![](modeling_files/figure-html/model4-6.png)

``` r

(model5<-bamp(cases, population, age="rw1", period="rw1", cohort="rw1",
             period_covariate = cov_p, periods_per_agegroup = 5))
plot(model5)
```

![](modeling_files/figure-html/model5-1.png)![](modeling_files/figure-html/model5-2.png)![](modeling_files/figure-html/model5-3.png)![](modeling_files/figure-html/model5-4.png)![](modeling_files/figure-html/model5-5.png)![](modeling_files/figure-html/model5-6.png)
