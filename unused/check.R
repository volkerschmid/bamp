#' Check an object
#'
#' @param x An object
#' @return logical; TRUE if check is fine.
#' @export
#'
check <- function(x, ...)
{
  UseMethod("check",x)
}
