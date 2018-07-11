#' Compute cohort index from age and period index
#'
#' @param agegroup age group index
#' @param period period index
#' @param noa number of age groups in total
#' @param periods_per_agegroup periods per age group
#'
#' @return cohort index
#' @export

coh <-
function(agegroup, period, noa, periods_per_agegroup){
  cohort <- (noa - agegroup)*periods_per_agegroup + period
  return(cohort)
}
