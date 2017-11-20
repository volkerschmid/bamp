setwd("~/Dropbox/temp/bamp/")
library(bamp)
data(apc)
model1 <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1",
               periods_per_agegroup = 5,  verbose=TRUE)

print(model1)
samples<-model1$samples

model2 <- bamp(cases, population, age="rw2", period="rw2", cohort="rw2",
               periods_per_agegroup = 5,  verbose=TRUE)

