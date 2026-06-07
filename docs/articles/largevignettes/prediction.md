# Bayesian Age-Period-Cohort Prediction

## Prediction

Using the prior assumption of a random walk for the period and cohort
effect, one can predict cases for upcoming years.

Here, we use the included data example.

``` r

data(apc)
plot(cases[,1],type="l",ylim=range(cases), ylab="cases", xlab="year", main="cases per age group")
for (i in 2:8)lines(cases[,i], col=i)
```

![](prediction_files/figure-html/loadplot-1.png)

We us only nine years and predict the last year.

``` r

model0 <- bamp(cases[-10,], population[-10,], age="rw1", period="rw1", cohort="rw1",
              periods_per_agegroup = 5)
```

``` r

model0<-predict_apc(object=model0, periods=1, population=population, update = TRUE)
```

Plot of predicted cases with credible intervals and true data

``` r

ts.plot(t(model0$predicted$cases_period), lty=c(2,1,2))
points(apply(cases,1,sum), pch=19)
```

![](prediction_files/figure-html/unnamed-chunk-1-1.png)

Plot period and cohort effects including prediction of year 10.

``` r

ts.plot(t(model0$predicted$period), lty=c(2,1,2))
```

![](prediction_files/figure-html/unnamed-chunk-2-1.png)

``` r

ts.plot(t(model0$predicted$cohort), lty=c(2,1,2))
```

![](prediction_files/figure-html/unnamed-chunk-2-2.png)
