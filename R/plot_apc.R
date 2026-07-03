#' Plot apc object
#'
#' @param x apc object
#' @param quantiles quantiles to plot. Default: \code{c(0.05,0.5,0.95)} is median and 90\% credible interval.
#' @param convention display-layer gauge convention for the linear trend (drift)
#'   in a full age-period-cohort model; one of \code{"age"} (default),
#'   \code{"period"}, \code{"cohort"} or \code{"none"}. In a full APC model the
#'   three effects are identifiable only up to a shared linear trend; this fixes
#'   that single degree of freedom, removing the run-to-run drift that is the
#'   dominant source of non-reproducibility in the plotted curves (residual
#'   curvature and Monte-Carlo noise remain). \code{"age"} shows the age effect
#'   as curvature about a zero trend and puts the drift in the period and cohort
#'   effects; \code{"none"} plots the raw, un-gauged effects. It is display-only
#'   (the fitted rates and predictions are unchanged) and is ignored for models
#'   that are not full APC. See \code{\link{effects.apc}}.
#' @param combined logical. For heterogeneity models, if \code{TRUE} plot the
#'   full effect (smooth + iid heterogeneity); default \code{FALSE} plots the
#'   smooth component only. Ignored for models without heterogeneity.
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
plot.apc<-function(x, quantiles=c(0.05,0.5,0.95),
                   convention=c("age","period","cohort","none"), combined=FALSE, ...)
{
  ## guard preserved from volkerschmid/bamp 2.2.0: nothing to plot if every
  ## chain failed / no samples were kept.
  if (length(x$samples) == 0 || is.null(x$samples$age)) {
    message("No MCMC samples available.")
    return(invisible(NULL))
  }
  convention <- match.arg(convention)
  g <- .apc_regauge(x, convention)
  if (isTRUE(combined)) g <- .apc_add_het(x, g)  # add iid het to the smooth curve
  q<-length(quantiles)
  # summarise to a q x (n groups) matrix; force the matrix shape so a single
  # quantile (q==1) keeps a row dimension -- apply() would otherwise drop it to
  # a plain vector and the age[i,] indexing below would fail.
  qmat<-function(samp) matrix(apply(as.array(samp),2,quantile,quantiles),nrow=q)
  age<-qmat(g$age)
  period<-qmat(g$period)
  cohort<-qmat(g$cohort)

  # symmetric line types: solid (median) in the middle, increasingly dashed
  # toward the tails. Reproduces the historical c(2,1,2)/c(3,2,1,2,3) for q==3/5
  # and stays valid for any q (previously only q in {1,3,5} worked; for q==2/4
  # lty stayed scalar 1 and lty[i] became NA -> "invalid line type").
  center<-(q+1)/2
  lty<-pmin(1L+round(abs(seq_len(q)-center)),6L)
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
      
      periodcov<-as.array(g$period)
      for (i in 1:dim(periodcov)[1])
        for (j in 1:dim(periodcov)[3])
          periodcov[i,,j]<-periodcov[i,,j]/x$covariate$period[1:c]
      periodcov<-matrix(apply(periodcov,2,quantile,quantiles),nrow=q)

      plot(periodcov[1,],type="l",lty=lty[1],ylim=range(periodcov[is.finite(periodcov)]),
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
      
      cohortcov<-as.array(g$cohort)
      for (i in 1:dim(cohortcov)[1])
        for (j in 1:dim(cohortcov)[3])
          cohortcov[i,,j]<-cohortcov[i,,j]/x$covariate$cohort[1:c]
      cohortcov<-matrix(apply(cohortcov,2,quantile,quantiles),nrow=q)

      plot(cohortcov[1,],type="l",lty=lty[1],ylim=range(cohortcov[is.finite(cohortcov)]),
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