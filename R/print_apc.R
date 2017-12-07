#' Print apc objects
#'
#' @param x apc object
#'
#' @return print
#' @export
#'
print.apc<-function(x, ...)
{
  ## Model
  
  settings <- ""
  if (is.null(x$model$age))x$model$age=" "
  settings <- paste0(settings, switch(x$model$age,
                     "rw1" = "age (rw1)",
                     "rw2" = "age (rw2)",
                     "rw1+het" = "age (rw1 with heterogeneity)",
                     "rw2+het" = "age (rw2 with heterogeneity)",
                     " " = ""
  ))
  if (x$model$age !=" "&x$model$period!=" ")settings=paste0(settings," - ")
  if (is.null(x$model$period))x$model$period=" "
  settings <- paste0(settings, switch(x$model$period,
                      "rw1" = "period (rw1)",
                      "rw2" = "period (rw2)",
                      "rw1+het" = "period (rw1 with heterogeneity)",
                      "rw2+het" = "period (rw2 with heterogeneity)",
                      " " = 0
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
  
  cat(paste0("Deviance: ", format(x$deviance$med.deviance,nsmall = 4, width=10),"\n"))
  cat(paste0("pD:       ", format(x$deviance$pD,nsmall = 4, width=10),"\n"))
  cat(paste0("DIC:      ", format(x$deviance$DIC,nsmall = 4, width=10),"\n\n"))
  
  ## Hyper parameters
  
  agepar <- quantile(unlist(x$samples$age_parameter), c(.05,.5,.95))
  perpar <- quantile(unlist(x$samples$period_parameter), c(.05,.5,.95))
  cohpar <- quantile(unlist(x$samples$cohort_parameter), c(.05,.5,.95))
  if(x$model$overdispersion)overdisp <- quantile(unlist(x$samples$overdispersion), c(.05,.5,.95))
  
  cat("\n","Hyper parameters:")
  cat(paste0(rep(" ",9)))
  cat(format(names(agepar),width=12))
  cat("\n")
  cat("age   ", paste0(rep(" ",10)))
  cat(format(agepar,digits = 4, width=12,trim=FALSE,nsmall=4))
  cat("\n")
  cat("period", paste0(rep(" ",10)))
  cat(format(perpar,digits = 4, width=12,trim=FALSE,nsmall=4))
  cat("\n")
  cat("cohort", paste0(rep(" ",10)))
  cat(format(cohpar,digits = 4, width=12,trim=FALSE,nsmall=4))
  cat("\n")
  if(x$model$overdispersion)
  {
    cat("overdispersion", paste0(rep(" ",6)))
    cat(format(overdisp,digits = 4, width=12,trim=FALSE,nsmall=4))
    cat("\n")
  }
}
