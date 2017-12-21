#' Compute mean from mcmc.list object
#'
#' @param mcmclist 
#'
#' @return vector
#' @export
meanmcmc<-function(mcmclist){
  if(!is.null(dim(mcmclist)))
    return(apply(matrix(unlist(parallel::mclapply(mcmclist,function(x)apply(x,2,mean))),ncol=length(mcmclist)),1,mean))
  if(is.null(dim(mcmclist)))
   return(mean(unlist(mcmclist)))
}