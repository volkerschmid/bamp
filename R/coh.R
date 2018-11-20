#' Compute cohort index from age and period index
#'
#' @param agegroup age group index
#' @param period period index
#' @param noa number of age groups in total
#' @param periods_per_agegroup periods per age group
#'
#' @return cohort index
#' @export
#' @examples 
#' # last agegroup in first period equals first cohort
#' coh(10, 1, 10, 5)  
#' 
#' # first agegroup in last period equals last cohort 
#' coh(1, 8, 10, 5) 

coh <-
function(agegroup, period, noa, periods_per_agegroup){
  cohort <- (noa - agegroup)*periods_per_agegroup + period
  return(cohort)
}
