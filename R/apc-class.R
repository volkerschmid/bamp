#' apc S3 class
#'
#' @return apc class
#' @export
#' @import coda
#' @description Class for (Bayesian) age-period-cohort objects
#' @details \code{\link{bamp}} will return an object of class apc. Available functions are
#' \itemize{
#' \item \code{\link{plot.apc}} plots main effects
#' \item \code{\link{print.apc}} print summary of model and effects 
#' \item \code{\link{effects.apc}} extract effects (mean, median and quantiles)}
#'
apc <- function()
{

  me <- list(
    samples=coda::mcmc.list(),
    model=list(),
    data=list()
  )

  ## Set the name for the class
  class(me) <- "apc"
  return(me)
}