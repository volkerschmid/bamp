#' Simulate from age-period-cohort model
#'
#' @description This functions simulates a data set of cases on the Lexis diagram from given age, period and cohort effects. Population numbers have to be given; can be one number for all age group/period combinations.
#' @param intercept Intercept
#' @param age Vector of effect for age groups
#' @param period Vector of effects for periods
#' @param cohort Vector of effect for cohorts
#' @param periods_per_agegroup Periods per age group
#' @param population Population number. Either a matrix or a scalar.
#' @seealso \code{vignette("simulation", package = "bamp")}
#' @return List with number of cases (matrix) and population numbers (matrix).
#' @export
#'
#' @examples age=sqrt(seq(5,0,length=10)); age<-1-age-mean(age)
#' period=15:1; period[8:15]<-8:15; period<-period/6; period<-period-mean(period)
#' periods_per_agegroup=5; number_of_cohorts <- periods_per_agegroup*(10-1)+15
#' cohort<-rep(0,60); cohort[1:10]<-10:1; cohort[41:60]<- -(1:20)/2; cohort<-cohort/10;
#' cohort<-cohort-mean(cohort)
#' simdata<-apcSimulate(-5, age, period, cohort, periods_per_agegroup, 1e6)
#' par(mfrow=c(3,1))
#' plot(age, type="l")
#' plot(period, type="l")
#' plot(cohort, type="l")
#' \dontrun{
#' simmod <- bamp(cases = simdata$cases, population = simdata$population, age = "rw1", 
#' period = "rw1", cohort = "rw1", periods_per_agegroup =periods_per_agegroup)
#' plot(simmod)
#' }

apcSimulate<-function(intercept, age, period, cohort, periods_per_agegroup, population)
{
  noa <- length(age)
  nop <- length(period)
  noc <- periods_per_agegroup*(noa-1)+nop
  if (noc!=length(cohort)){
    stop(paste0("Error: cohort needs to be of length ",noc))
  }
  n<-population
  if(is.null(dim(n)))n<-array(n,c(noa,nop))
  pr<- array(NA,c(noa, nop))
  for  (i in 1:noa)
    for (j in 1:nop)
    {
      ksi <- intercept + age[i] + period[j] + cohort[coh(i, j, noa, periods_per_agegroup)]
      pr[i,j] <- exp(ksi) / (1+exp(ksi))
    }
  y <-array(rbinom(n,n,pr), c(noa,nop))
  return(list("cases"=t(y), "population"=t(n)))
}



