## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache=TRUE)
library(bamp)

## ----loadplot------------------------------------------------------------
data(apc)
plot(cases[,1],type="l",ylim=range(cases), ylab="cases", xlab="year", main="cases per age group")
for (i in 2:8)lines(cases[,i], col=i)

## ----predmodel-----------------------------------------------------------
model0 <- bamp(cases[-10,], population[-10,], age="rw1", period="rw1", cohort="rw1",
              periods_per_agegroup = 5)

## ----prediction----------------------------------------------------------
model0<-predict_apc(object=model0, periods=1, population=population, update = TRUE)

## ------------------------------------------------------------------------
ts.plot(t(model0$predicted$cases_period), lty=c(2,1,2))
points(apply(cases,1,sum), pch=19)

## ------------------------------------------------------------------------
ts.plot(t(model0$predicted$period), lty=c(2,1,2))
ts.plot(t(model0$predicted$cohort), lty=c(2,1,2))

