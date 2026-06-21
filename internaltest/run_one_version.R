library(bamp)

data(apc)

t1 <- proc.time()
model1 <- bamp(cases, population, age = "rw1", period = "rw1", cohort = "rw1",
               periods_per_agegroup = 5)
elapsed_fit <- (proc.time() - t1)[["elapsed"]]

cat(sprintf("TIME fit: %.3f s\n", elapsed_fit))
checkConvergence(model1)
print(model1)

t2 <- proc.time()
pred <- predict_apc(object = model1, periods = 3)
elapsed_pred <- (proc.time() - t2)[["elapsed"]]

cat(sprintf("TIME predict: %.3f s\n", elapsed_pred))
