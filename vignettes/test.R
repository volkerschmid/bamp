
library(bamp)
data(apc)

system.time(model1 <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", overdisp = FALSE,
               periods_per_agegroup = 5,  verbose=2, parallel=4))

system.time(model2 <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", overdisp = TRUE,
              periods_per_agegroup = 5,  verbose=2))

print(model1)
print(model2)
