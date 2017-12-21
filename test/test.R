
library(bamp)
data(apc)

system.time(model1 <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", overdisp = FALSE,
               periods_per_agegroup = 5,  verbose=2, mcmc.options=list("number_of_iterations"=20000)))

system.time(model2 <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", overdisp = TRUE,
              periods_per_agegroup = 5,  verbose=2, mcmc.options=list("number_of_iterations"=20000)))

print(model1)
print(model2)

system.time(model3 <- bamp(cases, population, age="rw2", period="rw2", cohort="rw2", overdisp = FALSE,
                           periods_per_agegroup = 5,  verbose=2, mcmc.options=list("number_of_iterations"=40000)))

system.time(model4 <- bamp(cases, population, age="rw2", period="rw2", cohort="rw2", overdisp = TRUE,
                           periods_per_agegroup = 5,  verbose=2, mcmc.options=list("number_of_iterations"=105000)))

print(model3)
print(model4)
