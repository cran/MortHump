
#' Extrapolates non-extinct cohorts
#'
#' @description This is an internal function of package \pkg{MortHump} which extrapolates cohort mortality rates for non-extinct cohorts. It is meant to be called from the \link{HMD2MH} function.
#'
#' @param path path to the HMD folder for local access
#' @param max age at which the data should be right-censored
#' @param country HMD population letter code
#' @param year birth cohort for which the extrapolation should be performed
#' @param sex sex for which the extrapolation should be performed ("females" or "males")
#' @param cmx observed cohort mortality rates as given by a call to the \link{readHMD} function from the \pkg{HMDHFDplus} package
#'
#' @details The extrapolation is currently only available for females and males, not for the total population.
#' The function uses the model of Hyndman & Ullah (see reference), a variant of the Lee-Carter model which, among others, takes into account more than one dimension of the single value decomposition.
#'
#' The function first fits the model to period data, then extrapolates the age-specific death rates in the future, and reconstruct the cohort asdr by extracting the diagonal of the period asdr.
#'
#' @return Returns a vector of age-specific cohort mortality rates whose length depends on the \code{max} argument
#'
#'
#' @seealso \link{HMD2MH}, \link{fdm} from the \pkg{demography} package
#'
#' @references
#' Hyndman, R. J., & Ullah, M. S. (2007). Robust forecasting of mortality and fertility rates: a functional data approach. Computational Statistics & Data Analysis, 51(10), 4942-4956.
#'
#' @importFrom demography fdm
#' @importFrom demography forecast.fdm
#' @importFrom demography read.demogdata
#'
#' @export


xtp <- function(path, max, country, year, sex, cmx){

   if(sex %in% c("males","females") == F){warning("Extrapolation of non-extinct cohorts is currently only available for males and females.")}

    per <- read.demogdata(file = paste(path,country,"STATS","Mx_1x1.txt",sep="/"), popfile = paste(path,country,"STATS","Exposures_1x1.txt",sep="/"), type = "mortality", label = country)

    per$age <- 0:max
    per$rate$male <- per$rate$male[1:(max+1),]
    per$rate$female <- per$rate$female[1:(max+1),]
    per$pop$male <- per$pop$male[1:(max+1),]
    per$pop$female <- per$pop$female[1:(max+1),]

  if(any(colSums(is.na(per$rate$male)) > 10)){

    problem <- which(colSums(is.na(per$rate$male)) > 10)
    per$rate$male[,problem] <- per$rate$male[,min(problem)-1]
    per$rate$female[,problem] <- per$rate$female[,min(problem)-1]
    per$pop$male[,problem] <- per$pop$male[,min(problem)-1]
    per$pop$female[,problem] <- per$pop$female[,min(problem)-1]

    }

  end <- rev(per$year)[1]
  nonext <- (end - max):(end - 30)

  fun <- function(x){sum(is.na(x)) != 0 & is.na(x[1]) == F}
  av <- by(data = cmx$Female, INDICES = cmx$Year, FUN = fun)
  gens <- unique(cmx$Year)[av == T]

  cMx <- as.data.frame(matrix(nrow = max+1, ncol = length(gens)), row.names = 0:max)
  names(cMx) <- gens

  # projection
  s <- ifelse(sex == "males", "male", "female")
  fdm <- fdm(per, series = s)
  fc <- demography::forecast.fdm(fdm, h = max-30+1)

  if(sex == "males"){cmx.fc <- diag(fc$rate$male[(end-year+1):(max+1),])}
  if(sex == "females"){cmx.fc <- diag(fc$rate$female[(end-year+1):(max+1),])}

  cmx <- cmx[cmx$Year == year,]
  cmx[(end-year+1):(max+1),ifelse(sex == "males",4,3)] <- cmx.fc

  return(cmx)

}


