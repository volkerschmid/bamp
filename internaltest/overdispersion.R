
library(bamp)
data(apc)

model <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1",
               periods_per_agegroup = 5, overdisp=TRUE,
               mcmc.options=list(number_of_iterations = 6000, burn_in = 3000, step = 50,
                                tuning = 500), parallel=2)
#


m<-predict(model, periods=2)
