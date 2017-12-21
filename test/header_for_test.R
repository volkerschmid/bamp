data(apc)
        age= period= cohort="rw1"
overdisp=FALSE
         period_covariate=NULL
cohort_covariate=NULL
periods_per_agegroup = 5
verbose=2
parallel=FALSE
mcmc.options=list("number_of_iterations"=20000, "burn_in"=5000, "step"=50, "tuning"=500)
hyperpar=list("age"=c(1,0.0005), "period"=c(1,0.0005), "cohort"=c(1,0.0005), "overdisp"=c(1,0.05))
         dic=TRUE
