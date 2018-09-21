## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache=TRUE)
library(bamp)

## ----loadplot------------------------------------------------------------
data(apc)
plot(cases[,1],type="l",ylim=range(cases), ylab="cases", xlab="year", main="cases per age group")
for (i in 2:8)lines(cases[,i], col=i)

## ----Model1_with_RW1_priors----------------------------------------------
model1 <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1",
              periods_per_agegroup = 5)

## ----check_model---------------------------------------------------------
checkConvergence(model1)

## ----print_model---------------------------------------------------------
print(model1)

## ----plot_model----------------------------------------------------------
plot(model1)

## ----plot_model_with_more_quantiles--------------------------------------
plot(model1, quantiles = c(0.025,0.1,0.5,0.9,0.975))

## ----predmodel-----------------------------------------------------------
model2 <- bamp(cases[-10,], population[-10,], age="rw1", period="rw1", cohort="rw1",
              periods_per_agegroup = 5)

## ------------------------------------------------------------------------
#ts.plot(t(prediction$cases_period), lty=c(2,1,2))
#points(apply(cases,1,sum), pch=19)

## ----prdiction-----------------------------------------------------------
predion<-predict_apc(object=model1, periods=1, quantiles=c(0.05,0.5,0.95), update=TRUE)

## ----prediction----------------------------------------------------------
prection<-predict_apc(model2)

## ------------------------------------------------------------------------
#ts.plot(t(prediction$period), lty=c(2,1,2))
#ts.plot(t(prediction$cohort), lty=c(2,1,2))

