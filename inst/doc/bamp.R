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

## ----Model2_with_RW2_priors----------------------------------------------
model2 <- bamp(cases, population, age="rw2", period="rw2", cohort="rw2",
              periods_per_agegroup = 5)

## ----Model2_results------------------------------------------------------
checkConvergence(model2)
print(model2)
plot(model2)

## ----model3--------------------------------------------------------------
model3<-bamp(cases, population, age="rw1", period=" ", cohort="rw2",
              periods_per_agegroup = 5)
checkConvergence(model3)
print(model3)
plot(model3)

## ----model4--------------------------------------------------------------
(model4<-bamp(cases, population, age="rw1", period="rw1", cohort="rw1",
             cohort_covariate = cov_c, periods_per_agegroup = 5, 
             parallel = FALSE, 
             mcmc.options=list("number_of_iterations"=50000)
             ))
plot(model4)

## ----model5--------------------------------------------------------------
(model5<-bamp(cases, population, age="rw1", period="rw1", cohort="rw1",
             period_covariate = cov_p, periods_per_agegroup = 5, 
             parallel = FALSE))
plot(model5)

## ----predmodel-----------------------------------------------------------
model0 <- bamp(cases[-10,], population[-10,], age="rw1", period="rw1", cohort="rw1",
              periods_per_agegroup = 5)

## ----prediction----------------------------------------------------------
prediction<-predict_apc(object=model0, periods=1, population=population)

## ------------------------------------------------------------------------
ts.plot(t(prediction$cases_period), lty=c(2,1,2))
points(apply(cases,1,sum), pch=19)

## ------------------------------------------------------------------------
ts.plot(t(prediction$period), lty=c(2,1,2))
ts.plot(t(prediction$cohort), lty=c(2,1,2))

