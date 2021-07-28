#' Print apc objects
#' @aliases print
#' @param x apc object
#' @param ... additional arguments will be ignored
#'
#' @return print
#' @export
#'
#' @examples 
#' \dontrun{
#' data(apc)
#' model <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", periods_per_agegroup = 5)
#' print(model)
#' }
print.apc<-function(x, ...)
{
  
  ## Check convergence
  
  cc <- checkConvergence(x, auto=TRUE)
  
  if (!cc)cat("\nWARNING! Markov Chains have apparently not converged! DO NOT TRUST THIS MODEL!\n")
  
  ## Model
  
  settings <- ""
  if (is.null(x$model$age))x$model$age=" "
  settings <- paste0(settings, switch(x$model$age,
                     "rw1" = "age (rw1) ",
                     "rw2" = "age (rw2) ",
                     "rw1+het" = "age (rw1 with heterogeneity) ",
                     "rw2+het" = "age (rw2 with heterogeneity) ",
                     " " = ""
  ))
  if (x$model$age !=" "&x$model$period!=" ")settings=paste0(settings," - ")
  if (is.null(x$model$period))x$model$period=" "
  settings <- paste0(settings, switch(x$model$period,
                      "rw1" = "period (rw1) ",
                      "rw2" = "period (rw2) ",
                      "rw1+het" = "period (rw1 with heterogeneity) ",
                      "rw2+het" = "period (rw2 with heterogeneity) ",
                      " " = ""
  ))
  if (x$model$cohort !=" "&x$model$period!=" ")settings=paste0(settings," - ")
  if (is.null(x$model$cohort))x$model$cohort=" "
  settings <- paste0(settings, switch(x$model$cohort,
                             "rw1" = "cohort (rw1)",
                             "rw2" = "cohort (rw2)",
                             "rw1+het" = "cohort (rw1 with heterogeneity)",
                             "rw2+het" = "cohort (rw2 with heterogeneity)",
                             " " = 0
  ))
  settings <- paste0(settings," model")
  if (!is.null(x$model$overdisp))if(x$model$overdisp)settings=paste0(settings," with overdispersion.")  
  cat("\n Model:\n");
  cat(paste0(settings,"\n"))
  
  ## Deviance
  
  cat(paste0("Deviance: ", format(x$deviance$mean.deviance,digits = 2,nsmall = 2, width=10,trim=FALSE),"\n"))
  cat(paste0("pD:       ", format(x$deviance$pD,digits = 2,nsmall = 2, width=10,trim=FALSE),"\n"))
  cat(paste0("DIC:      ", format(x$deviance$DIC,digits = 2,nsmall = 2, width=10,trim=FALSE),"\n\n"))
  
  ## Hyper parameters
  
  agepar <- quantile(unlist(x$samples$age_parameter), c(.05,.5,.95))
  perpar <- quantile(unlist(x$samples$period_parameter), c(.05,.5,.95))
  cohpar <- quantile(unlist(x$samples$cohort_parameter), c(.05,.5,.95))
  if(x$model$overdispersion)overdisp <- quantile(unlist(x$samples$overdispersion), c(.05,.5,.95))
  
  cat("\n","Hyper parameters:")
  cat(paste0(rep(" ",9)))
  cat(format(names(agepar),width=12))
  cat("\n")
  if (x$model$age!=" ")
  {cat("age   ", paste0(rep(" ",10)))
  cat(format(agepar,digits = 3, width=12,trim=FALSE,nsmall=3))
  cat("\n")}
  if (x$model$period!=" "){
    cat("period", paste0(rep(" ",10)))
  cat(format(perpar,digits = 3, width=12,trim=FALSE,nsmall=3))
  cat("\n")}
  if (x$model$cohort!=" "){
    cat("cohort", paste0(rep(" ",10)))
  cat(format(cohpar,digits = 3, width=12,trim=FALSE,nsmall=3))
  cat("\n")}
  if(x$model$overdispersion)
  {
    cat("overdispersion", paste0(rep(" ",6)))
    cat(format(overdisp,digits = 3, width=12,trim=FALSE,nsmall=3))
    cat("\n")
  }
  
  ## Convergence
  
  if (cc)cat("\n\nMarkov Chains convergence checked succesfully using Gelman's R (potential scale reduction factor).")
  
}
