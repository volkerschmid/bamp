#' Plot apc object
#'
#' @param x apc object
#' @param quantiles quantiles to plot
#' @param ... Additional arguments will be ignored
#'
#' @import stats graphics
#' @return plot
#' @export
#'
plot.apc<-function(x, quantiles=c(.05,.5,.95), ...)
{
  age<-as.array(x$samples$age)
  if(x$model$age=="rw2"){
    n<-dim(age)[2]
    d<-(age[,n,]-age[,1,])/(n-1)
        for (i in n:1)
      age[,i,]<-age[,i,]-(i-1)*d
    m<-apply(age,c(1,3),mean)
    for (i in n:1)
      age[,i,]<-age[,i,]-m
    }
  age<-apply(age,2,quantile)
  period<-as.array(x$samples$period)
  if(x$model$period=="rw2"){
    n<-dim(period)[2]
    d<-(period[,n,]-period[,1,])/(n-1)
    for (i in n:1)
      period[,i,]<-period[,i,]-(i-1)*d
    m<-apply(period,c(1,3),mean)
    for (i in n:1)
      period[,i,]<-period[,i,]-m

  }
  period<-apply(period,2,quantile)
  cohort<-as.array(x$samples$cohort)
  if(x$model$cohort=="rw2"){
    n<-dim(cohort)[2]
    d<-(cohort[,n,]-cohort[,1,])/(n-1)
    for (i in n:1)
      cohort[,i,]<-cohort[,i,]-(i-1)*d
    m<-apply(cohort,c(1,3),mean)
    for (i in n:1)
      cohort[,i,]<-cohort[,i,]-m
  }
  cohort<-apply(cohort,2,quantile)

  q<-length(quantiles)

  lty=1
  if (length(quantiles)==3)lty=c(2,1,2)
  if (length(quantiles)==5)lty=c(3,2,1,2,3)
  par(mfrow=c(3,1))
    plot(age[1,],type="l",lty=lty[1],ylim=range(age),
         axes=FALSE,main="age",xlab="",ylab="")
    if (is.null(x$data$agegroups))axis(1,lwd=0)
    if (!is.null(x$data$agegroups))axis(1,lwd=0,at=1:length(x$data$agegroups),labels=x$data$agegroups)
    axis(2,lwd=0)
    for (i in 1:q)
      lines(age[i,],lty=lty[i])
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
