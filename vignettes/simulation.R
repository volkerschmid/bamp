## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache=TRUE, out.width = "90%")
library(bamp)

## ----simulate_age--------------------------------------------------------
age=2*sqrt(seq(1,20,length=10))
age<- age-mean(age)
plot(age, type="l")

## ----simulate_period-----------------------------------------------------
period=15:1
period[8:15]<-8:15
period<-period/5
period<-period-mean(period)
plot(period, type="l")

## ----simulate_cohort-----------------------------------------------------
periods_per_agegroup=5
number_of_cohorts <- periods_per_agegroup*(10-1)+15
cohort<-rep(0,60)
cohort[1:15]<-(14:0)
cohort[16:30]<- (1:15)/2
cohort[31:60]<- 8
cohort<-cohort/10
cohort<-cohort-mean(cohort)
plot(cohort, type="l")

## ----simulate_data-------------------------------------------------------
simdata<-apcSimulate(-10, age, period, cohort, periods_per_agegroup, 1e6)
print(simdata$cases)

## ----bamp----------------------------------------------------------------
simmod <- bamp(cases = simdata$cases, population = simdata$population, age = "rw1", 
period = "rw1", cohort = "rw1", periods_per_agegroup =periods_per_agegroup)

## ----check_print_plot----------------------------------------------------
print(simmod)
checkConvergence(simmod)
plot(simmod)

## ----plot_comparison-----------------------------------------------------
effects<-effects(simmod)
effects2<-effects(simmod, mean=TRUE)
#par(mfrow=c(3,1))
plot(age, type="l")
lines(effects$age, col="blue")
lines(effects2$age, col="green")
plot(period, type="l")
lines(effects$period, col="blue")
lines(effects2$period, col="green")
plot(cohort, type="l")
lines(effects$cohort, col="blue")
lines(effects2$cohort, col="green")

## ----prediction----------------------------------------------------------
prediction<-predict_apc(simmod, periods=5, population=array(1e6,c(20,10)))

## ----prediction1---------------------------------------------------------
plot(prediction$cases_period[2,], ylim=range(prediction$cases_period),ylab="",pch=19)
points(prediction$cases_period[1,],pch="–",cex=2)
points(prediction$cases_period[3,],pch="–",cex=2)
for (i in 1:20)lines(rep(i,3),prediction$cases_period[,i])

## ----prediction4---------------------------------------------------------
plot(prediction$period[2,])

