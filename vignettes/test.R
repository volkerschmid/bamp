
library(bamp)
data(apc)

model <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", overdisp = TRUE,
               periods_per_agegroup = 5,  verbose=TRUE)

