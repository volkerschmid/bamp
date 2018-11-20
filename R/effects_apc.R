#' Effects from Fitted APC Model
#'
#' @param object an apc object
#' @param mean logical. If TRUE, mean effects are computed
#' @param quantiles Scalar or vector of quantiles to compute (only if mean=FALSE)
#' @param update logical. If TRUE, the apc object including the effects is returned
#' @param ... Additional arguments will be ignored

#' @return List of age, period, cohort effects or apc object including effects (if update=TRUE)
#' @export
#' @examples
#' \dontrun{
#' data(apc)
#' model <- bamp(cases, population, age="rw1", period="rw1", cohort="rw1", periods_per_agegroup = 5)
#' effects(model)
#' }

effects.apc<-function(object, mean=FALSE, quantiles=0.5, update=FALSE, ...)
{
  x<-object
  #check, if we have done this before
  if (!is.null(x$effect))
  {
    if (all(attr(effects,"settings")==c(mean,quantiles)))return(x$effect) 
  }
  
  age=x$samples$age
  age=summary(age, quantiles=quantiles)
  if (mean)
  {
    age=age$statistics[,1]
  }
  else
  {
    age=age$quantiles
  }
  period=x$samples$period
  period=summary(period, quantiles=quantiles)
  if (mean)
  {
    period=period$statistics[,1]
  }
  else
  {
    period=period$quantiles
  }
  cohort=x$samples$cohort
  cohort=summary(cohort, quantiles=quantiles)
  if (mean)
  {
    cohort=cohort$statistics[,1]
  }
  else
  {
    cohort=cohort$quantiles
  }
  effects<-list("age"=age, "period"=period, "cohort"=cohort)
  attr(effects,"settings") <- c(mean,quantiles)
  if (update){x$effects<-effects; return(x)}
  return(effects)
}