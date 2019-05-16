.onAttach <- function(lib, pkg)
{
    # startup message
    msg <- paste("bamp version", utils::packageVersion("bamp"))
    packageStartupMessage(msg)
    invisible()
}