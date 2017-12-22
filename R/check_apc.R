#' Check apc object, whether MCMC has converged
#'
#' @param x An apc object
#' @param level level of check; 1 uses point point estimation, 2 uses upper C.I.
#' @param auto logical; should be TRUE if called automatically from bamp()
#'
#' @description 
#' bamp uses Gelman and Rubins R to check convergence for all main parameters. 
#' All parameters should have R<1.1. 
#' bamp() runs at least four MCM chains (more if parallel is more than four).
#' level=2 convergence check is not yet implemented.  
#'
#' @import coda 
#' @return logical; TRUE if check is fine.
#' @export
#'

checkConvergence<-function(x, level=2, auto=FALSE)
{
  j<-level
  theta<-coda::gelman.diag(x$samples$age, multivariate = FALSE)
  phi<-coda::gelman.diag(x$samples$period, multivariate = FALSE)
  psi<-coda::gelman.diag(x$samples$cohort, multivariate = FALSE)
  kappa<-coda::gelman.diag(x$samples$age_parameter, multivariate = FALSE)
  lambda<-coda::gelman.diag(x$samples$period_parameter, multivariate = FALSE)
  ny<-coda::gelman.diag(x$samples$cohort_parameter, multivariate = FALSE)
  my<-coda::gelman.diag(x$samples$intercept, multivariate = FALSE)
  if (any(theta$psrf[,j]>1.1)|any(phi$psrf[,j]>1.1)|any(psi$psrf[,j]>1.1)|kappa$psrf[j]>1.1|lambda$psrf[j]>1.1|ny$psrf[j]>1.1|my$psrf[j]>1.1)
  {
    if(!auto)cat("Warning: MCMC chains did not converge!\n")
    return=FALSE
  }
  else{
    if (!auto)
    {
      return=TRUE
    }
  }
  return(return)
}
