#' Bayesian Age-Period-Cohort Modelling and Prediction (bamp)
#'
#' @param cases number of cases
#' @param population population number
#' @param age prior for age groups ("rw1", "rw2", "rw1+het", "rw2+het", " ")
#' @param period prior for periods ("rw1", "rw2", "rw1+het", "rw2+het", " ")
#' @param cohort prior for cohorts ("rw1", "rw2", "rw1+het", "rw2+het", " ")
#' @param overdisp logical, add overdispersion to model
#' @param periods_per_agegroup periods per agegroup
#' @param period_covariate covariate for period
#' @param cohort_covariate covariate for cohort
#' @param mcmc.options Options for MCMC, see description
#' @param hyperpar Hyper parameters, see description
#' @param parallel logical, should computation be done in parallel
#' @param verbose verbose mode
#'
#' @param mcmc.options 
#' @param hyperpar 
#' 
#' @description 
#'
#' @useDynLib bamp
#' @export
#' @import coda
#'
bamp <-
function(cases, population,
        age, period, cohort, overdisp=FALSE,
        period_covariate=NULL, cohort_covariate=NULL,
        periods_per_agegroup,
        mcmc.options=list("number_of_iterations"=105000, "burn_in"=5000, "step"=1, "tuning"=500),
        hyperpar=list("age"=c(1,0.0005), "period"=c(1,0.0005), "cohort"=c(1,0.0005), "overdisp"=c(1,0.05)),
        parallel=TRUE, verbose=FALSE){
  output=apc()

  age_hyperpar_a=hyperpar$age[1]
  age_hyperpar_b=hyperpar$age[2]
  period_hyperpar_a=hyperpar$period[1]
  period_hyperpar_b=hyperpar$period[2]
  cohort_hyperpar_a=hyperpar$cohort[1]
  cohort_hyperpar_b=hyperpar$cohort[2]
  z_hyperpar_a=hyperpar$overdisp[1]
  z_hyperpar_b=hyperpar$overdisp[2]
  if (age=="rw1+het"|age=="rw2+het")
    {
    age_hyperpar_a2=hyperpar$age_het[1]
    age_hyperpar_b2=hyperpar$age_het[2]
  }
  else
  {
    age_hyperpar_a2=1
    age_hyperpar_b2=1
  }
  if (period=="rw1+het"|period=="rw2+het")
  {
    period_hyperpar_a2=hyperpar$period_het[1]
  period_hyperpar_b2=hyperpar$period_het[2]
  }
  else
  {
    period_hyperpar_a2=1
    period_hyperpar_b2=1
  }
  if (cohort=="rw1+het"|period=="rw2+het")
  {
    cohort_hyperpar_a2=hyperpar$cohort_het[1]
  cohort_hyperpar_b2=hyperpar$cohort_het[2]
  }
  else
  {
    cohort_hyperpar_a2=1
    cohort_hyperpar_b2=1
  }
  
  chains=4
  if (parallel>1)chains=parallel
  if(!parallel)chain=1
  
  # have been options before
  if (is.null(mcmc.options$burn_in))
  {
    burn_in=5000
  }
  else
  {
      burn_in=mcmc.options$burn_in
  }
  if(is.null(mcmc.options$number_of_iterations))
    {
    number_of_iterations=100000+burn_in
  }
  else
  {
        number_of_iterations=mcmc.options$number_of_iterations
  }
  if (is.null(mcmc.options$tuning))
    {
    tuning=500
    }
  else
    {
      tuning=mcmc.options$tuning
    }
  if (is.null(mcmc.options$step))
  {
    step=1
    }
  else
  {
      step=mcmc.options$step
      }
  dataorder = 0
  number_of_agegroups=dim(cases)[2]
  number_of_periods=dim(cases)[1]
  cohort_start = 1
  period_start = 1

  z_mode=ifelse(overdisp,1,0)

  model<-list()
  if (is.null(age))age=" "
  age_block = switch(age,
    "rw1" = 1,
    "rw2" = 2,
    "rw1+het" = 3,
    "rw2+het" = 4,
    " " = 0
  )
  if (is.null(period))period=" "
  period_block = switch(period,
                     "rw1" = 1,
                     "rw2" = 2,
                     "rw1+het" = 3,
                     "rw2+het" = 4,
                     " " = 0
  )
  if (is.null(cohort))cohort=" "
  cohort_block = switch(cohort,
                     "rw1" = 1,
                     "rw2" = 2,
                     "rw1+het" = 3,
                     "rw2+het" = 4,
                     " " = 0
  )

  model$age=age
  model$period=period
  model$cohort=cohort

  #if (!is.null(age_covariate))age_block=age_block=age_block+7
  if (!is.null(period_covariate))period_block=period_block=period_block+7
  if (!is.null(cohort_covariate))cohort_block=cohort_block=cohort_block+7


  zentrieren <-
    function(mat, my){
      summe <- rowSums(mat)

      for(i in 1:dim(mat)[1]){
        mat[i,] <- mat[i,] - summe[i]/dim(mat)[2]
        my[i] <- my[i] + summe[i]/dim(mat)[2]
      }


      return(list(mat, my))
    }

##################################################################################################################################
# check of the input and preparation of the variables (bamp.cc S.1-S.11)

  stopifnot(is.data.frame(cases) || is.matrix(cases))              # cases must be a data.frame or a matrix
  if(is.data.frame(cases)){
    cases <- as.matrix(cases)
  }

  stopifnot(is.data.frame(population) || is.matrix(population))    # population must be a data.frame or a matrix
  if(is.data.frame(population)){
    population <- as.matrix(population)
  }

  if(dataorder == 0 || dataorder == 1){                            # dataorder must be in (0,1)
  }else{
    stop("ERROR: Dataorder must be 0 or 1!")
  }

  if(dataorder == 0){
    cases <- t(cases)
    population <- t(population)
  }

  stopifnot(number_of_agegroups%%1 == 0 && number_of_agegroups > 0)# number_of_agegroups must be a whole number
  stopifnot(number_of_periods%%1 == 0 && number_of_periods > 0)    # number_of_periods must be a whole number
  stopifnot(is.numeric(periods_per_agegroup))                      # periods_per_agegroup must be numeric


  stopifnot(number_of_iterations%%1 == 0)                          # number_of_iterations must be a whole number
  stopifnot(burn_in%%1 == 0)                                       # burn_in must be a whole number
  if(number_of_iterations <= burn_in){                             # number_of_iterations must be bigger than burn_in
    stop("ERROR: Number of iterations must be bigger than burn in!")
  }

  stopifnot(step%%1 == 0)                                          # step must be a whole number
  stopifnot(number_of_iterations-burn_in >= step)                  # number_ofiterations-burn_in musst be bigger than step
  stopifnot(tuning%%1 == 0)                                        # tuning must be a whole number
  if(burn_in <= tuning){                                           # burn_in must be bigger than tuning
    stop("ERROR: Burn in must be bigger than tuning constant!")
  }

  if(age_block == 0 | age_block == 1 | age_block == 2 |            # age_block must be in (0,1,2,3,4,8,9)
     age_block == 3 | age_block == 4){
  }else{
     stop("ERROR: Age_Block must be 0, 1, 2, 3 or 4!")
  }
  if(age_block > 0){
    stopifnot(is.numeric(age_hyperpar_a) & age_hyperpar_a > 0)     # needed for gammadistribution
    stopifnot(is.numeric(age_hyperpar_b) & age_hyperpar_b > 0)     # needed for gammadistribution
  }
  if(age_block == 3 | age_block == 4){
    stopifnot(is.numeric(age_hyperpar_a2) & age_hyperpar_a2 > 0)   # needed for gammadistribution
    stopifnot(is.numeric(age_hyperpar_b2) & age_hyperpar_b2 > 0)   # needed for gammadistribution
  }

  if(period_block == 0 | period_block == 1 |                       # period_block must be in (0,1,2,3,4,8,9)
     period_block == 2 | period_block == 3 |
     period_block == 4 | period_block == 8 |
     period_block == 9){
  }else{
    stop("ERROR: Period Block must be 0, 1, 2, 3, 4, 8 or 9!")
  }

  if(period_block > 0){
    stopifnot(is.numeric(period_hyperpar_a) & period_hyperpar_a > 0)     # needed for gammadistribution
    stopifnot(is.numeric(period_hyperpar_b) & period_hyperpar_b > 0)     # needed for gammadistribution
  }
  if(period_block == 3 | period_block == 4){
    stopifnot(is.numeric(period_hyperpar_a2) & period_hyperpar_a2 > 0)   # needed for gammadistribution
    stopifnot(is.numeric(period_hyperpar_b2) & period_hyperpar_b2 > 0)   # needed for gammadistribution
  }

  if(cohort_block == 0 | cohort_block == 1 |                      # cohort_block must be in (0,1,2,3,4,8,9)
     cohort_block == 2 | cohort_block == 3 |
     cohort_block == 4 | cohort_block == 8 |
     cohort_block == 9){
  }else{
    stop("ERROR: Cohort Block must be 0, 1, 2, 3, 4, 8 or 9!")
  }
  if(cohort_block > 0){
    stopifnot(is.numeric(cohort_hyperpar_a) & cohort_hyperpar_a > 0)     # needed for gammadistribution
    stopifnot(is.numeric(cohort_hyperpar_b) & cohort_hyperpar_b > 0)     # needed for gammadistribution
  }
  if(cohort_block == 3 | cohort_block == 4){
    stopifnot(is.numeric(cohort_hyperpar_a2) & cohort_hyperpar_a2 > 0)   # needed for gammadistribution
    stopifnot(is.numeric(cohort_hyperpar_b2) & cohort_hyperpar_b2 > 0)   # needed for gammadistribution
  }

  if(z_mode == 0 | z_mode == 1){                                   # z_mode must be in (0,1)
  }else{
    stop("ERROR: Z Mode must be 0 or 1!")
  }
  if(z_mode == 1){
    stopifnot(is.numeric(z_hyperpar_a) & z_hyperpar_a > 0)         # needed for gammadistribution
    stopifnot(is.numeric(z_hyperpar_b) & z_hyperpar_b > 0)         # needed for gammadistribution
  }

# additional conditions and new variables
  if(number_of_periods <= periods_per_agegroup){               # number_of_periods must be bigger than periods_per_agegroup
    stop("ERROR: Need more periods!")
  }


  number_of_cohorts <- periods_per_agegroup*(number_of_agegroups-1)+number_of_periods # number_of_cohorts


  max_number_of_ap_combinations <- floor((number_of_periods/periods_per_agegroup)+0.999) # max_combinations


  number_of_extractions <- (number_of_iterations - burn_in)/step    # number_of_extractions

  period_plus <- 0
  if(period_block == 8 || period_block == 9){                       # conditions for period_block 8 and 9
    period_block <- period_block - 7
    period_data <- c()
    period_plus <- 1
  }

  cohort_plus <- 0
  if(cohort_block == 8 || cohort_block == 9){                       # conditions for cohort_block 8 and 9
    cohort_block <- cohort_block - 7
    cohort_data <- c()
    cohort_plus <- 1
  }

  cphi <- 1                                                         # set phi = 1
  cpsi <- 1                                                         # set psi = 1

  if(period_plus == 1){                                             # if period_plus = 1 (with period_covariate)
    if(is.vector(period_covariate)){
      cphi <- 0                                                     # set phi = 0
    }else{
      stop("ERROR:The period covariate must be in a vector form!")# period covariate must be in a vector form
    }
    if(period_start%%1 != 0 || period_start <= 0){
      stop("ERROR: Period start must be a positive integer!")     # period start must be a positive integer
    }else{
      if(period_start > length(period_covariate)){
        stop("ERROR: Period start is bigger than the count of the covariables!")# period start can not be bigger than the vector
      }else{
        period_covariate <- period_covariate[period_start:length(period_covariate)]
      }
    }
    if(length(period_covariate) < number_of_periods){
      stop("ERROR: Not enough observations of the period covariabe!")    # need more observations
    }else{
      period_data <- period_covariate[1:number_of_periods]  # collects the period data
    }
    for(i in 1:number_of_periods){
      cphi <- cphi + period_data[i]                                  # calculates a cphi
    }
    cphi <- cphi/number_of_periods
    for(i in 1:(number_of_periods)){
      period_data[i] <- period_data[i]/cphi                          # calculates the period data
    }
  }

  if(cohort_plus == 1){                                              # if cohort_plus = 1 (with cohort_covariate)
    if(is.vector(cohort_covariate)){
      cpsi <- 0                                                      # set psi = 0
    }else{
      stop("ERROR: The cohort covariate must be in a vector form!")# cohort covariate must be in a vector form
    }
  if(cohort_start%%1 != 0 || cohort_start <= 0){
    stop("ERROR: Cohort start must be a positive integer!")         # cohort start must be a positive integer
  }else{
    if(cohort_start > length(cohort_covariate)){
      stop("ERROR: Cohort start is bigger than the count of the covariables!")# cohort start can not be bigger than the vector
    }else{
      cohort_covariate <- cohort_covariate[cohort_start:length(cohort_covariate)]
    }
  }
    if(length(cohort_covariate) < number_of_cohorts){
      stop("ERROR: Not enough observations of the cohort covariabe!")  # need more observations
    }else{
      cohort_data <- cohort_covariate[1:number_of_cohorts]  # collects the cohort data
    }
    for(i in 1:number_of_cohorts){
      cpsi <- cpsi + cohort_data[i]                                        # calculates a cpsi
    }
    cpsi <- cpsi/number_of_cohorts
    for(i in 1:number_of_cohorts){
      cohort_data[i] <- cohort_data[i]/cpsi                                # calculates the cohort data
    }
  }

if (verbose)
  {
  ##################################################################################################################################
## output of settings (bamp.cc S.12-S-14)

  max_block <- max(age_block, (max(period_block, cohort_block)))
  if(max_block == 0){
    stop(cat("ERROR: No Block!\n"))
  }
  settings <- character()                                                 # AGE-Model ?
  if(age_block > 0){
    settings <- "AGE"
  }
  if(age_block*period_block > 0){
    settings <- paste(settings, "-", sep = "")
  }
  if(period_block > 0){
    settings <- paste(settings, "PERIOD", sep = "")                       # PERIOD-Model ?
  }
  if(age_block*cohort_block > 0 || period_block*cohort_block > 0){
    settings <- paste(settings, "-", sep = "")
  }
  if(cohort_block > 0){
    settings <- paste(settings, "COHORT", sep = "")                       # COHORT-Model ?
  }
  settings <- paste(settings, "Model", sep = " ")
  if(z_mode == 1){
    settings <- paste(settings, "with overdispersion", sep = " ")         # with overdispersion when z_mode = 1
  }
  if(period_block > 2){                                                   # with unstructured period effects for period_block > 2
    if(z_mode == 1){
      settings <- paste(settings, "and unstructured period effects", sep = " ")
    }else{
      settings <- paste(settings, "with unstructured period effects", sep = " ")
    }
  }
  if(age_block > 2){
    if(z_mode == 1 || period_block > 2){                                  # with unstructured age effects for age_block > 2
      settings <- paste(settings, "and unstructured age effects", sep = " ")
    }else{
      settings <- paste(settings, "with unstructured age effects", sep = " ")
    }
  }
  if(cohort_block > 2){
    if(z_mode == 1 || period_block > 2 || age_block > 2){                 # with unstructured cohort effects for cohort_block > 2
      settings <- paste(settings, "and unstructured cohort effects", sep = " ")
    }else{
      settings <- paste(settings, "with unstructured cohort effects", sep = " ")
    }
  }

  settings <- paste(settings, ".", sep = "")
  settings <- paste(settings, "\nPrioris: ", sep = "")                    # Prioris
  if(age_block > 0){
    settings <- paste(settings, "age effect", sep = " ")                  # age
  }
  if(age_block == 1 | age_block == 3){
    settings <- paste(settings, "RW 1", sep = " ")                        # RW 1
  }
  if(age_block == 2 | age_block == 4){
    settings <- paste(settings, "RW 2", sep = " ")                        # RW 2
  }
  if(age_block*period_block > 0){
    settings <- paste(settings, ",", sep = "")
  }
  if(period_block > 0){
    settings <- paste(settings, "period effect", sep = " ")               # period
  }
  if(period_block == 1 | period_block == 3){
    settings <- paste(settings, "RW 1", sep = " ")                        # RW 1
  }
  if(period_block == 2 | period_block == 4){
    settings <- paste(settings, "RW 2", sep = " ")                        # RW 2
  }
  if(cohort_block > 0){
    settings <- paste(settings, ", cohort effect", sep = "")              # cohort
  }
  if(cohort_block == 1 | cohort_block == 3){
    settings <- paste(settings, "RW 1", sep = " ")                        # RW 1
  }
  if(cohort_block == 2 | cohort_block == 4){
    settings <- paste(settings, "RW 2", sep = " ")                        # RW 2
  }
  settings <- paste(settings, "\n", sep = "")
  if(period_plus == 1){
    settings <- paste(settings, "Period effect with covariates ", names(period_covariate),   # period effect with covariates
                      " starting at position ", period_start, "\n", sep = "")
  }
  if(cohort_plus == 1){
    settings <- paste(settings, "Cohort effect with covariates ", names(cohort_covariate),   # cohort effect with covariates
                      " starting at position ", cohort_start, "\n", sep = "")
  }
  settings <- paste(settings, number_of_agegroups, " age groups, ", number_of_periods, " periods, ",
                    number_of_cohorts, " cohorts. ", sep = "")            # counts of agegroups, periods and cohorts
  settings <- paste(settings, ".", sep = "")
  settings <- paste(settings, "\n", sep = "")
  settings <- paste(settings, number_of_iterations, " iterations with ", burn_in, " burn in, using every ",
                    step, "th sample.", sep = "")                        # iterations, burn_in, step and tuning
  if(tuning > 0){
    settings <- paste(settings, " Tuning at iteration ", tuning, sep = "")
  }
  settings <- paste(settings, "\n", sep = "")

  done <- 0

  if(age_block > 0 || period_block > 0 || cohort_block > 0 || z_mode > 0){   # Hyper parameters
    settings <- paste(settings, "Hyper parameters: ", sep = "")
  }
  if(age_block > 0){                                                         # for age
    settings <- paste(settings, "age eff. (", age_hyperpar_a, ", ", age_hyperpar_b, ")", sep = "")
    done <- 1
    if(age_block == 3 | age_block == 4){
      settings <- paste(settings, " - unstr. age eff. (", age_hyperpar_a2, ", ", age_hyperpar_b2, ")", sep = "")
    }
  }
  if(period_block > 0){                                                      # for period
    if(done == 1){
      settings <- paste(settings, " - ", sep = "")
    }
    settings <- paste(settings, "period eff. (", period_hyperpar_a, ", ", period_hyperpar_b, ")", sep = "")
    done <- 1
    if(period_block == 3 | period_block == 4){
      settings <- paste(settings, " - unstr. period eff. (", period_hyperpar_a2, ", ", period_hyperpar_b2, ")", sep = "")
    }
  }
  if(cohort_block > 0){                                                      # for cohort
    if(done == 1){
      settings <- paste(settings, " - ", sep = "")
    }
    settings <- paste(settings, "cohort eff. (", cohort_hyperpar_a, ", ", cohort_hyperpar_b, ")", sep = "")
    done <-1
    if(cohort_block == 3 | cohort_block == 4){
      settings <- paste(settings, " - unstr. cohort eff. (", cohort_hyperpar_a2, ", ", cohort_hyperpar_b2, ")", sep = "")
    }
  }
  if(z_mode == 1){                                                           # for heterogeneity
    if(done == 1){
    settings <- paste(settings, " - ", sep = "")
    }
    settings <- paste(settings, "overdispersion (", z_hyperpar_a, ", ", z_hyperpar_b, ")", sep = "")
  }
  settings <- paste(settings, ".\n\n", sep = "")
}

##################################################################################################################################
# additional settings and variables



  # additional conditions for the data (TxJ-matrix for dataorder = 0, JxT-matrix for dataorder =1)
  if(dim(cases)[2]  != number_of_periods || dim(cases)[1] != number_of_agegroups){
    stop("ERROR: Cases data does not fit to settings!")
  }



  if(period_block == 3 || period_block == 4){
    delta <- 100                                                              # set delta = 100 for period_block 3 and 4
  }

  if(cohort_block == 3 || cohort_block == 4){
    delta <- 100                                                              # set delta = 100 for cohort_block 3 and 4
  }



##################################################################################################################################
  if (verbose)
    {
    ausgabe <- paste(settings, "Starting Iterations!\n\n")
    cat(ausgabe)
  }

 nr.samples<-floor(number_of_iterations/step)-floor(burn_in/step)
 number_of_cohorts = ceiling(periods_per_agegroup*(number_of_agegroups - 1)+number_of_periods)

 delta.sample<-kappa.sample<-kappa2.sample<-lambda.sample<-lambda2.sample<-ny.sample<-ny2.sample<-my.sample<-dev.sample<-rep(0,nr.samples)
 theta.sample<-theta2.sample<-rep(0,nr.samples*number_of_agegroups)
 phi.sample<-phi2.sample<-rep(0,nr.samples*number_of_periods)
 psi.sample<-psi2.sample<-rep(0,nr.samples*number_of_cohorts)
 #ksi.sample<-rep(0,nr.samples*number_of_agegroups*number_of_periods)

 blocks=c(age_block, period_block, cohort_block)
 numbers=c(number_of_agegroups, number_of_periods)
 numbersmcmc=c(number_of_iterations, burn_in, step, tuning)
 modelsettings=c(0, 0, z_mode)
 allhyper=c(age_hyperpar_a, age_hyperpar_b, period_hyperpar_a, period_hyperpar_b, cohort_hyperpar_a, cohort_hyperpar_b,
            age_hyperpar_a2, age_hyperpar_b2, period_hyperpar_a2, period_hyperpar_b2, cohort_hyperpar_a2, cohort_hyperpar_b2,
            z_hyperpar_a, z_hyperpar_b)

 singlerun<-function(i,cases,population,blocks,numbers,periods_per_agegroup,
                     numbersmcmc,modelsettings,allhyper,theta.sample,phi.sample,psi.sample,
                     theta2.sample,phi2.sample,psi2.sample,#ksi.sample,
                     delta.sample,kappa.sample,kappa2.sample,
                     lambda.sample,lambda2.sample,ny.sample,ny2.sample,my.sample,
                     dev.sample,verbose){
return(.C("bamp",
                          as.integer(cases),
                          as.integer(population),
                          as.integer(blocks),
                          as.integer(numbers),
                          as.double(periods_per_agegroup),
                          as.integer(numbersmcmc),
                          as.integer(modelsettings),
                          as.double(allhyper),
                          as.double(theta.sample), as.double(phi.sample), as.double(psi.sample),
                          as.double(theta2.sample), as.double(phi2.sample), as.double(psi2.sample),
                          #as.double(ksi.sample),
                          as.double(delta.sample), as.double(kappa.sample), as.double(kappa2.sample), as.double(lambda.sample),
                          as.double(lambda2.sample), as.double(ny.sample), as.double(ny2.sample), as.double(my.sample),
                          as.double(dev.sample),
                          as.integer(verbose)
                          )
)
}

if(parallel)results_list<-parallel::mclapply(1:chains,singlerun,cases,population,blocks,numbers,periods_per_agegroup,
                                                 numbersmcmc,modelsettings,allhyper,theta.sample,phi.sample,psi.sample,
                                                 theta2.sample,phi2.sample,psi2.sample,delta.sample,kappa.sample,kappa2.sample,
                                                 lambda.sample,lambda2.sample,ny.sample,ny2.sample,my.sample,dev.sample,verbose, mc.cores=chains)

 if(!parallel)results_list<-lapply(1:chains,singlerun,cases,population,blocks,numbers,periods_per_agegroup,
                                              numbersmcmc,modelsettings,allhyper,theta.sample,phi.sample,psi.sample,
                                              theta2.sample,phi2.sample,psi2.sample,delta.sample,kappa.sample,kappa2.sample,
                                              lambda.sample,lambda2.sample,ny.sample,ny2.sample,my.sample,dev.sample,verbose)

##################################################################################################################################


 theta<-phi<-psi<-theta2<-phi2<-psi2<-delta<-kappa<-kappa2<-lambda<-lambda2<-ny<-ny2<-my<-deviance<-vector("list",chains)

 for (i in 1:chains){
   deviance[[i]]=coda::mcmc(results_list[[i]][[23]])
 }
 deviance<-coda::as.mcmc.list(deviance)
 kick<-rep(TRUE, chains)
 deviance.mean<-unlist(lapply(deviance,mean))
 deviance.sd<-unlist(lapply(deviance,sd))
 dm<-median(deviance.mean)
 sd<-median(deviance.sd)
 kick<-((deviance.mean>=(dm-sd))&(deviance.mean<=(dm+sd)))
 

  for (i in (1:chains)[kick]){
   results=results_list[[i]]
   theta[[i]]=coda::mcmc(matrix(results[[9]],ncol=number_of_agegroups,byrow=TRUE))
   phi[[i]]=coda::mcmc(matrix(results[[10]],ncol=number_of_periods,byrow=TRUE))
   psi[[i]]=coda::mcmc(matrix(results[[11]],ncol=number_of_cohorts,byrow=TRUE))
   theta2[[i]]=coda::mcmc(matrix(results[[12]],ncol=number_of_agegroups,byrow=TRUE))
   phi2[[i]]=coda::mcmc(matrix(results[[13]],ncol=number_of_periods,byrow=TRUE))
   psi2[[i]]=coda::mcmc(matrix(results[[14]],ncol=number_of_cohorts,byrow=TRUE))
   kappa[[i]]=coda::mcmc(results[[16]])
   kappa2[[i]]=coda::mcmc(results[[17]])
   lambda[[i]]=coda::mcmc(results[[18]])
   lambda2[[i]]=coda::mcmc(results[[19]])
   ny[[i]]=coda::mcmc(results[[20]])
   ny2[[i]]=coda::mcmc(results[[21]])
   my[[i]]=coda::mcmc(results[[22]])
   deviance[[i]]=coda::mcmc(results[[23]])
   delta[[i]]=coda::mcmc(results[[15]])
   }


theta<-coda::as.mcmc.list(theta)
phi<-coda::as.mcmc.list(phi)
psi<-coda::as.mcmc.list(psi)
theta2<-coda::as.mcmc.list(theta2)
phi2<-coda::as.mcmc.list(phi2)
psi2<-coda::as.mcmc.list(psi2)
kappa<-coda::as.mcmc.list(kappa)
kappa2<-coda::as.mcmc.list(kappa2)
lambda<-coda::as.mcmc.list(lambda)
lambda2<-coda::as.mcmc.list(lambda2)
ny<-coda::as.mcmc.list(ny)
ny2<-coda::as.mcmc.list(ny2)
my<-coda::as.mcmc.list(my)
deviance<-coda::as.mcmc.list(deviance)


 samples=list("intercept"=my, "age"=theta, "period"=phi, "cohort"=psi)
 if (age_block==3)samples=c(samples,list("age2"=theta2))
 if (period_block==3)samples=c(samples,list("period2"=phi2))
 if (cohort_block==3)samples=c(samples,list("cohort2"=psi2))
 samples=c(samples, list("age_parameter"=kappa, "period_parameter"=lambda, "cohort_parameter"=ny))
 if (age_block==3)samples=c(samples,list("age2_parameter"=kappa2))
 if (period_block==3)samples=c(samples,list("period2_parameter"=lambda2))
 if (cohort_block==3)samples=c(samples,list("cohort2_parameter"=ny2))
 if (z_mode==1) samples=c(samples, list("delta"=delta))

 data=list("cases"=cases,"population"=population)

 output$model=model
 output$data=data
 output$samples=samples

 
 if (verbose)cat("\nComputing deviance and DIC.")
   
 theta.med<-summary(theta,0.5)$quantiles
 phi.med<-summary(phi,0.5)$quantiles
 psi.med<-summary(psi,0.5)$quantiles
 my.med<-summary(my,0.5)$quantiles
 if (z_mode==1)delta.med<-summary(delta,0.5)$quantiles
 
 devtemp=0.0;
  for(i in 1:number_of_agegroups){
       for(j in 1:number_of_periods){
         ksi <- my.med + theta.med[i] + phi.med[j] + psi.med[coh(i,j,number_of_agegroups,periods_per_agegroup)]
         pr <- exp(ksi)/(1+exp(ksi))
         ydach <- population[i,j]*pr
         if(cases[i,j]==0.0)
         {   
           devtemp=devtemp+2*((population[i,j]-cases[i,j])*log((population[i,j]-cases[i,j])/(population[i,j]-ydach)));
         }
         else
         {
           devtemp=devtemp+2*(cases[i,j]*log(cases[i,j]/ydach)+(population[i,j]-cases[i,j])*log((population[i,j]-cases[i,j])/(population[i,j]-ydach)));
         }
      }
  }
 
 med.deviance<-summary(deviance,0.5)$quantiles
 deviance <- list()
 deviance$med.deviance <- med.deviance
 deviance$deviance.med <- devtemp
 deviance$pD <- med.deviance-devtemp
 deviance$DIC <- 2*med.deviance-devtemp
 
 output$deviance=deviance
 
 # ksi_berechnen <-
 #   function(ksi, psi, vdb, noa, nop){
 #     for(i in 1:noa){
 #       for(j in 1:nop){
 #         ksi[j,i] <- psi[coh(i, j, noa, vdb)]
 #       }
 #     }
 #     return(ksi)
 #   }
 
# check(output)

 return(output)
}
