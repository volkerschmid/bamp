
library(bamp)
data(apc)

model1 <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", overdisp = FALSE,
               periods_per_agegroup = 5,  verbose=2, parallel=6)

model2 <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", overdisp = TRUE,
              periods_per_agegroup = 5,  verbose=2, parallel=6)

print(model1)
print(model2)
