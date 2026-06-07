# Short Introduction to BAMP

## Data example

BAMP includes a data example.

``` r

data(apc)
plot(cases[,1],type="l",ylim=range(cases), ylab="cases", xlab="year", main="cases per age group")
for (i in 2:8)lines(cases[,i], col=i)
```

![](bamp_files/figure-html/loadplot-1.png)

For simulating APC data, see *vignette(“simulation”, package=“bamp”)*.

## APC model with random walk first order prior

``` r

model1 <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1",
              periods_per_agegroup = 5)
```

    ## Warning: MCMC chains did not converge!

bamp() automatically performs a check for MCMC convergence using Gelman
and Rubin’s convergence diagnostic. We can manually check the
convergence again:

``` r

checkConvergence(model1)
```

    ## Warning: MCMC chains did not converge!

    ## [1] FALSE

Now we have a look at the model results. This includes estimates of
smoothing parameters and deviance and DIC:

``` r

print(model1)
```

    ## 
    ## WARNING! Markov Chains have apparently not converged! DO NOT TRUST THIS MODEL!
    ## 
    ##  Model:
    ## age (rw1)  - period (rw1)  - cohort (rw1) model
    ## Deviance:     231.52
    ## pD:            36.86
    ## DIC:          268.38
    ## 
    ## 
    ##  Hyper parameters:                 5%           50%          95%         
    ## age                              0.358        0.890        1.887
    ## period                          72.847      204.735      663.295
    ## cohort                          34.838       59.817       97.881

We can plot the main APC effects using point-wise quantiles:

``` r

plot(model1)
```

![](bamp_files/figure-html/plot_model-1.png)![](bamp_files/figure-html/plot_model-2.png)![](bamp_files/figure-html/plot_model-3.png)

More quantiles are possible:

``` r

plot(model1, quantiles = c(0.025,0.1,0.5,0.9,0.975))
```

![](bamp_files/figure-html/plot_model_with_more_quantiles-1.png)![](bamp_files/figure-html/plot_model_with_more_quantiles-2.png)![](bamp_files/figure-html/plot_model_with_more_quantiles-3.png)

For other models see *vignette(“modeling”,package=“bamp”)*.

## Prediction

Using the prior assumption of a random walk for the period and cohort
effect, one can predict cases for upcoming years.

``` r

pred <- predict_apc(object=model1, periods=3)
```

``` r

m<-max(pred$pr[2,,])
plot(pred$pr[2,,8],type="l", ylab="probability", xlab="year", ylim=c(0,m))
for (i in 7:1)
  lines(pred$pr[2,,i],col=8-i)
legend(1,m,col=8:1,legend=paste("Age group",1:8),lwd=2,cex=0.6)
lines(c(10.5,10.5),c(0,1),lty=2)
```

![](bamp_files/figure-html/unnamed-chunk-1-1.png)

More details see *vignette(“prediction”,package=“bamp”)*.
