library(bamp)

# Input data (US)
counts<-as.matrix(read.table("us_data_2014.txt",
                   row.names=1, header=F))
pop<-as.matrix(read.table("us_pop_2014.txt",
                row.names=1, header=F))
agegroups <- c("25-29", "30-34", "35-39", "40-44", "45-49", "50-54", 
               "55-59", "60-64", "65-69", "70-74", "75-79", "80-84" )
age <- seq(27,82,by=5)

#APC model with random walk first order prior
#model1 <- bamp(counts, pop, age="rw1", period="rw1", cohort="rw1",
#               periods_per_agegroup = 5,  
#               mcmc.options = list(number_of_iterations = 750,
#                                   burn_in = 600, step = 50, tuning = 500),
#               verbose=TRUE, parallel = FALSE)
#bamp() automatically performs a check for MCMC convergence
checkConvergence(model1)

#Ergebnis
print(model1)
plot(model1)

#More quantiles are possible:
plot(model1, quantiles = c(0.025,0.1,0.5,0.9,0.975))

## rw2
model2 <- bamp(counts, pop, age="rw2", period="rw2", cohort="rw2",
  periods_per_agegroup = 5, mcmc.options = list(number_of_iterations = 5000,
                                                burn_in = 2000, step = 50, tuning = 1000),
  verbose=TRUE, parallel = FALSE)
  
  