#' Prediction for age-period-cohort models
#'
#' @param object apc object
#' @param periods number of periods to predict
#' @param population matrix of (predicted) population, if NULL, population data from original bamp call will be used
#' @param quantiles vector of quantiles to compute
#' @param update boolean. If TRUE, object will be returned with results added to the object
#'
#' @return list with quantiles of predicted probabilities (pr), predicted cases (cases) and predicted cases per period (cases_period)
#' and a list samples with MCMC samples of pr, cases and cases_period. 
#' If update=TRUE, the apc object will be returned with this list (predicted) added.
#' 
#' @import parallel
#' @importFrom abind abind
#' @export
#'
predict_apc<-function(object, periods=0, population=NULL, quantiles=c(0.05,0.5,0.95), update=FALSE){
  ksi_prognose <-
    function(prepi, vdb, noa, nop, nop2, noc, zmode){
      my<-prepi[1]
      theta<-prepi[2:(noa+1)]
      phi<-prepi[(noa+2):(noa+nop2+1)]
      psi<-prepi[(noa+nop2+2):(noa+nop2+noc+1)]
      delta<-prepi[noa+nop2+noc+2]
      
      ksi<-array(0, c(nop2,noa))
      for(i in 1:noa){
        for(j in 1:nop2){
          ksi[j,i] <- my + theta[i] + phi[j] + psi[bamp::coh(i,j,noa,vdb)]
          
          if(zmode){
            ksi[j,i] <- ksi[j,i] + (rnorm(1, mean = 0, sd = 1)/sqrt(delta))
          }
        }
      }
      
      return(ksi)
    }
  
  predict_rw <-
    function(prepi, rw, n1, n2){
      lambda<-prepi[1]
      phi<-prepi[-1]
      if(rw == 1){
        for(i in (n1+1):n2){
          phi[i] <- (rnorm(1, mean = 0, sd = 1)/sqrt(lambda)) + phi[i-1]
        }
      }
      
      if(rw == 2){
        for(i in (n1+1):n2){
          phi[i] <- (rnorm(1, mean = 0, sd = 1)/sqrt(lambda)) + (2*phi[i-1] - phi[i-2])
        }
      }
      
      if(rw == 0){
        for(i in (n1+1):n2){
          phi[i] <- phi[i-1]
        }
      }
      
      
      return(phi)
    }
  
  phi<-psi<-NA
  
  a1<-dim(object$data$cases)[1]
  n1<-dim(object$data$cases)[2]
  n2<-n1+periods
  rwp<-rwc<-0
  
  if (!object$model$period=="")
  {
    rwp<-switch(object$model$period,
                rw1 = 1,
                rw2 = 2
    )
    ch<-length(object$samples$period)
      prep<-parallel::mclapply(1:ch, function(i,samples)cbind(object$samples$period_parameter[[i]],object$samples$period[[i]]), samples)
      phi<-parallel::mclapply(prep, function(prepi, rw, n1, n2){t(apply(prepi, 1, predict_rw, rw, n1, n2))}, rwp, n1, n2)
  }
  
  if (!object$model$cohort=="")
  {
    rwc<-switch(object$model$cohort,
                rw1 = 1,
                rw2 = 2
    )
    c1<-dim(object$samples$cohort[[1]])[2]
    c2<-bamp::coh(1,n2,a1,object$data$periods_per_agegroup)
    ch<-length(object$samples$cohort)
      prep<-parallel::mclapply(1:ch, function(i,samples)cbind(object$samples$cohort_parameter[[i]],object$samples$cohort[[i]]), samples)
      psi<-parallel::mclapply(prep, function(prepi, rw, n1, n2){t(apply(prepi, 1, predict_rw, rw, n1, n2))}, rwc, c1, c2)
  }
    
  nr.samples<-length(object$samples$intercept[[1]])
    theta<-if(object$model$age!=""){object$samples$age}else{NA}
    
      delta<-if(object$model$overdispersion){object$samples$overdispersion}else{NA}

      prepfx<- function(i, theta, phi, psi, my, delta, a1, n1, c1, nr){
      theta=if(any(is.na(theta))){array(0, c(nr, a1))}else{theta[[i]]}
      phi=if(any(is.na(phi))){array(0, c(nr, a1))}else{phi[[i]]}
      psi=if(any(is.na(psi))){array(0, c(nr, a1))}else{psi[[i]]}
      my=my[[i]]
      delta=if(any(is.na(delta))){rep(0, nr)}else{delta[[i]]}
      return(cbind(my,theta,phi,psi, delta))
    }
    
    prep<- parallel::mclapply(1:ch, prepfx, theta, phi, psi, object$samples$intercept, delta, a1, n1, c1, nr.samples)

  ksi<-parallel::mclapply(prep, function(prepi, vdb, noa, nop, nop2, noc, zmode){
    temp<-apply(prepi, 1, ksi_prognose, vdb, noa, nop, nop2, noc, zmode); return(array(temp,c(nop2,noa,dim(temp)[2])))}, object$data$periods_per_agegroup,
    a1, n1, n2, c2, object$model$overdispersion)

    
  ksi0<-ksi[[1]]
  if (ch>1)
  for (i in 2:ch)
    ksi0<-abind::abind(ksi0,ksi[[i]], along=3)
  
  pr <- exp(ksi0)/(1+exp(ksi0))
  
  if (is.null(population))population<-object$data$population
  
  n0<-min(dim(population)[1],n2)
  
  predictedcases <- array(apply(pr,3,function(pr1,n,n0){
    return(rbinom(n0*dim(n)[2],n[1:n0,],pr1[1:n0,]))},population,n0),c(n0,dim(pr)[2:3]))
  
  predictperiod<-apply(predictedcases,c(1,3),sum)

  qu_predicted_cases<-apply(predictedcases,1:2,quantile,quantiles)
  qu_predicted_pr<-apply(pr,1:2,quantile,quantiles)
  qu_predictperiod<-apply(predictperiod,1,quantile,quantiles)
  
  
  period<-phi[[1]]
  if (ch>1)
  for (i in 2:ch)
    period<-abind::abind(period,phi[[i]], along=1)
  cohort<-psi[[1]]
  if (ch>1)
  for (i in 2:ch)
    cohort<-abind::abind(cohort,psi[[i]], along=1)

  samples<-list(
    "pr"=pr,
    "cases"=predictedcases,
    "cases_period"=predictperiod,
    "period"=period,
    "cohort"=cohort
  )
  
  predicted<-list(
    "pr"=qu_predicted_pr,
    "cases"=qu_predicted_cases,
    "cases_period"=qu_predictperiod,
    "period"=apply(period,2,quantile,quantiles),
    "cohort"=apply(cohort,2,quantile,quantiles),
    "samples"=samples
  )
  
  if (!update){
    return(predicted)}
  else{
    object$predicted=predicted
    return(object)
  }
}
