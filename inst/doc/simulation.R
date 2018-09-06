## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache=TRUE, out.width = "90%")
library(bamp)

## ----prediction----------------------------------------------------------
prediction<-predict_apc(simmod, periods=5, population=array(1e6,c(20,10)))
plot(prediction$cases_period[2,], ylim=range(prediction$cases_period),ylab="",pch=19)
points(prediction$cases_period[1,],pch="–",cex=2)
points(prediction$cases_period[3,],pch="–",cex=2)
for (i in 1:20)lines(rep(i,3),prediction$cases_period[,i])
plot(prediction$period[2,])

## ----bamp_rw2------------------------------------------------------------
simmodrw2 <- bamp(cases = simdata$cases, population = simdata$population, age = "rw2", 
period = "rw2", cohort = "rw2", periods_per_agegroup =periods_per_agegroup)

## ----check_print_plot_rw2------------------------------------------------
print(simmodrw2)
checkConvergence(simmodrw2)
plot(simmodrw2)

## ----simulate_covariate_for_period---------------------------------------
cov_p<-rnorm(15,period,.1)

## ----bamp_with_covariate-------------------------------------------------
simmod2 <- bamp(cases = simdata$cases, population = simdata$population, age = "rw1", 
period = "rw1", cohort = "rw1", periods_per_agegroup =periods_per_agegroup,
period_covariate = cov_p)

