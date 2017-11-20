#' apc S3 class
#'
#' @return apc class
#' @export
#' @import coda
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