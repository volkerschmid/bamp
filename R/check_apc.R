#' Check apc object, whether MCMC has converged
#'
#' @param x An apc object
#' @param info logical; print more information
#' @param level level of check; 1 uses point point estimation, 2 uses upper C.I.
#' @param auto logical; should be TRUE if called automatically from \code{\link{bamp}}
#' #'
#' @description 
#' This functions uses Gelman and Rubin's R to check convergence for all main parameters. 
#' All parameters should have R<1.1. 
#' \code{\link{bamp}} runs at least four MCMC chains by default (more if parallel is more than four).
#'
#' @import coda 
#' @return logical; TRUE if check is fine.
#' @export
#' @examples 
#' \dontrun{
#' data(apc)
#' model <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", periods_per_agegroup = 5)
#' checkConvergence(model)
#' }


checkConvergence<-function(x, info=FALSE, level=2, auto=FALSE)
{
  j<-level
  if(x$model$age!=" ")theta<-coda::gelman.diag(x$samples$age, multivariate = FALSE)
  if(x$model$period!=" ")phi<-coda::gelman.diag(x$samples$period, multivariate = FALSE)
  if(x$model$cohort!=" ")psi<-coda::gelman.diag(x$samples$cohort, multivariate = FALSE)
  if(x$model$age!=" ")kappa<-coda::gelman.diag(x$samples$age_parameter, multivariate = FALSE)
  if(x$model$period!=" ")lambda<-coda::gelman.diag(x$samples$period_parameter, multivariate = FALSE)
  if(x$model$cohort!=" ")ny<-coda::gelman.diag(x$samples$cohort_parameter, multivariate = FALSE)
  my<-coda::gelman.diag(x$samples$intercept, multivariate = FALSE)
  test<-FALSE
  if(x$model$age!=" ")test<-any(theta$psrf[,j]>1.1)|kappa$psrf[j]>1.1
  if(x$model$period!=" ")test<-test|any(phi$psrf[,j]>1.1)|lambda$psrf[j]>1.1
  if(x$model$cohort!=" ")test<-test|any(psi$psrf[,j]>1.1)|kappa$psrf[j]>1.1
  if (test|my$psrf[j]>1.1)
  {
    if(!auto){
      cat("Warning: MCMC chains did not converge!\n")
      if (info)
      {
        cat("(mean) R       Point est. Upper C.I. \n")
        cat("age effect    ",apply(theta$psrf,2,mean),"\n")
        cat("period effect ",apply(phi$psrf,2,mean),"\n")
        cat("cohort effect ",apply(psi$psrf,2,mean),"\n")
        cat("age hyperp.   ",apply(theta$psrf,2,mean),"\n")
        cat("period hyperp.",apply(phi$psrf,2,mean),"\n")
        cat("cohort hyperp.",apply(psi$psrf,2,mean),"\n")
      }
    }
    return=FALSE
  }
  else{
    return=TRUE
  }
  return(return)
}
