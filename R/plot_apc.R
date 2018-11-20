#' Plot apc object
#'
#' @param x apc object
#' @param quantiles quantiles to plot. Default: \code{c(0.05,0.5,0.95)} is median and 90\% credible interval.
#' @param ... Additional arguments will be ignored
#' 
#' @details Plot of age, period and cohort effects from apc objects. If covariates have been used for period/cohort, a second plot with covariate, absolute effect and relative effect is created. Absolute effect is relative effect times covariate. 
#' @import stats graphics
#' @return plot
#' @export
#' @examples
#' \dontrun{
#' data(apc)
#' model <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", periods_per_agegroup = 5)
#' plot(model)
#' }
plot.apc<-function(x, quantiles=c(0.05,0.5,0.95), ...)
{
  age<-as.array(x$samples$age)
  age<-apply(age,2,quantile,quantiles)
  period<-as.array(x$samples$period)
  period<-apply(period,2,quantile,quantiles)
  cohort<-as.array(x$samples$cohort)
  cohort<-apply(cohort,2,quantile,quantiles)

  q<-length(quantiles)

  lty=1
  if (length(quantiles)==3)lty=c(2,1,2)
  if (length(quantiles)==5)lty=c(3,2,1,2,3)
  #par(mfrow=c(3,1))
    if (x$model$age!=" ")
      {
      plot(age[1,],type="l",lty=lty[1],ylim=range(age),
         axes=FALSE,main="age",xlab="",ylab="")
    if (is.null(x$data$agegroups))axis(1,lwd=0)
    if (!is.null(x$data$agegroups))axis(1,lwd=0,at=1:length(x$data$agegroups),labels=x$data$agegroups)
    axis(2,lwd=0)
    for (i in 1:q)
      lines(age[i,],lty=lty[i])}
  if (x$model$period!=" ")
  {
    plot(period[1,],type="l",lty=lty[1],ylim=range(period),
         axes=FALSE,ylab="", xlab="", main="period")
    if (is.null(x$data$periods))axis(1,lwd=0)
    if (!is.null(x$data$periods))
      {
        axis(1,lwd=0,at=1:length(x$data$periods),labels=x$data$periods)
      }
    axis(2,lwd=0)
    for (i in 1:q)
      lines(period[i,],lty=lty[i])
  }
  if (x$model$cohort!=" ")
  {
    plot(cohort[1,],type="l",lty=lty[1],ylim=range(cohort),
         axes=FALSE,ylab="", xlab="", main="cohort")
    if (is.null(x$data$cohorts))axis(1,lwd=0)
    if (!is.null(x$data$cohorts))
    {
      axis(1,lwd=0,at=1:length(x$data$cohorts),labels=x$data$cohorts)
    }
    axis(2,lwd=0)
    for (i in 1:q)
      lines(cohort[i,],lty=lty[i])
  }
  
  if (!is.null(x$covariate))
  {
    if (!is.null(x$covariate$period))
    {
      c<-dim(period)[2]
      plot(x$covariate$period,type="l", main="period covariate", ylab="",
           xlim=c(1,c), axes=FALSE)
      if (is.null(x$data$periods))axis(1,lwd=0)
      if (!is.null(x$data$periods))
      {
        axis(1,lwd=0,at=1:length(x$data$periods),labels=x$data$periods)
      }
      axis(2,lwd=0)
      
      plot(period[1,],type="l",lty=lty[1],ylim=range(period),
           axes=FALSE,ylab="", xlab="", main="period effect")
      if (is.null(x$data$periods))axis(1,lwd=0)
      if (!is.null(x$data$periods))
      {
        axis(1,lwd=0,at=1:length(x$data$periods),labels=x$data$periods)
      }
      axis(2,lwd=0)
      for (i in 1:q)
        lines(period[i,],lty=lty[i])
      
      periodcov<-as.array(x$samples$period)
      for (i in 1:dim(periodcov)[1])
        for (j in 1:dim(periodcov)[3])
          periodcov[i,,j]<-periodcov[i,,j]/x$covariate$period[1:c]
      periodcov<-apply(periodcov,2,quantile,quantiles)
      
      plot(periodcov[1,],type="l",lty=lty[1],ylim=range(periodcov),
           axes=FALSE,ylab="", xlab="", main="raw period covariate effect")
      if (is.null(x$data$periods))axis(1,lwd=0)
      if (!is.null(x$data$periods))
      {
        axis(1,lwd=0,at=1:length(x$data$periods),labels=x$data$periods)
      }
      axis(2,lwd=0)
      for (i in 1:q)
        lines(periodcov[i,],lty=lty[i])
    }

    if (!is.null(x$covariate$cohort))
    {
      c<-dim(cohort)[2]
      plot(x$covariate$cohort,type="l", main="cohort covariate", ylab="",
           xlim=c(1,c), axes=FALSE)
      if (is.null(x$data$cohorts))axis(1,lwd=0)
      if (!is.null(x$data$cohorts))
      {
        axis(1,lwd=0,at=1:length(x$data$cohorts),labels=x$data$cohorts)
      }
      axis(2,lwd=0)
      
      plot(cohort[1,],type="l",lty=lty[1],ylim=range(cohort),
           axes=FALSE,ylab="", xlab="", main="cohort effect")
      if (is.null(x$data$cohorts))axis(1,lwd=0)
      if (!is.null(x$data$cohorts))
      {
        axis(1,lwd=0,at=1:length(x$data$cohorts),labels=x$data$cohorts)
      }
      axis(2,lwd=0)
      for (i in 1:q)
        lines(cohort[i,],lty=lty[i])
      
      cohortcov<-as.array(x$samples$cohort)
      for (i in 1:dim(cohortcov)[1])
        for (j in 1:dim(cohortcov)[3])
          cohortcov[i,,j]<-cohortcov[i,,j]/x$covariate$cohort[1:c]
      cohortcov<-apply(cohortcov,2,quantile,quantiles)
      
      plot(cohortcov[1,],type="l",lty=lty[1],ylim=range(cohortcov),
           axes=FALSE,ylab="", xlab="", main="raw effect of cohort covariate")
      if (is.null(x$data$cohorts))axis(1,lwd=0)
      if (!is.null(x$data$cohorts))
      {
        axis(1,lwd=0,at=1:length(x$data$cohorts),labels=x$data$cohorts)
      }
      axis(2,lwd=0)
      for (i in 1:q)
        lines(cohortcov[i,],lty=lty[i])
    }
  }

}