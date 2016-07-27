#' Format data from the Human Mortality Database for use in the MortHump package
#'
#' @param country \href{http://www.mortality.org/cgi-bin/hmd/hmd_download.php}{HMD population letter code}
#' @param year year (or birth cohort) for which the data is required
#' @param dim dimension (period or cohort)
#' @param xtra if TRUE extrapolates non-extinct cohorts
#' @param sex sex for which the data is required ("females", "males" or "total")
#' @param min age at which the data should be left-censored
#' @param max age at which the data should be right-censored
#' @param username personal HMD username for web access
#' @param password personal HMD password for web access
#' @param path path to the HMD folder for local access (containing the country folders)
#'
#' @details
#' Two methods are available to access the data, either online or locally. Web access is advised
#' to guarantee the most recent data, but can be slow. For a fast access, download first the \href{http://www.mortality.org/cgi-bin/hmd/hmd_download.php}{complete
#' zipped data files} (choose "All countries for the HMD"), and then indicate the path to the general folder.
#' The data collection method will be adapted depending if you indicate a path or a username and a password.
#'
#' If you choose cohort data, you have the option of either work with the observed data, or to extrapolate the non-extinct cohorts.
#' In the latter case, a variant of the Lee-Carter model is used to extrapolate the period mortality rates, which are used to identify cohort rates in the diagonal of the Lexis matrix.
#' This option uses a variant of the Lee-Carter model, namely the functional model proposed by Hyndman and Ullah (2007) and available in the \pkg{demography} package.
#'
#' The list of available countries and period coverage can be found on the \href{http://www.mortality.org/cgi-bin/hmd/DataAvailability.php}{HMD website}.
#' This list is also accessible with the function \code{getHMDcountries()} from the \pkg{HMDHFDplus} package.
#'
#' The \code{max} argument is designed to deal with either the presence of a mortality plateau among centenarians, which can be diffcult to capture with parametric models, or
#' a high level of stochasticity at old ages due to a small number of survivors. It is recommended to keep the value of \code{max} above 80 or 90.
#'
#' @return A data frame containing the following variables.
#' \describe{
#'   \item{\code{x}}{vector of ages.}
#'   \item{\code{d}}{vector of death counts.}
#'   \item{\code{n}}{vector of population exposures.}
#'   \item{\code{m}}{vector of death rates.}
#' }
#'
#' @seealso This function makes use of the functions \code{readHMD} and \code{readHMDweb} from the \pkg{HMDHFDplus} package.
#'
#' @references
#' Human Mortality Database. University of California, Berkeley (USA), and Max Planck Institute for Demographic Research (Germany). Available at \url{www.mortality.org} or \url{www.humanmortality.de}.
#'
#' @import HMDHFDplus
#'
#' @export


HMD2MH <- function(country, year, dim = "period", xtra = FALSE, sex, min = 0, max = NULL, username = NULL, password = NULL, path = NULL){

  if(all(is.null(c(username,password,path)))){warning("You must provide either valid username and password (web) or path to the folder (local).")}

  if(sex %in% c("males","females","total") == F){warning("sex must be one of 'females', 'males' or 'total'")}

  if(dim %in% c("period","cohort") == F){warning("dim must be one of 'period' or 'cohort'")}

  if(dim == "period"){

  if(is.null(path) == F){

    dx <- readHMD(filepath = file.path(path,country,"STATS","Deaths_1x1.txt"), fixup = T)
    nx <- readHMD(filepath = file.path(path,country,"STATS","Exposures_1x1.txt"), fixup = T)

    }

  if(all(is.null(c(username,password)) == F)){

    dx <- readHMDweb(username = username, password = password, CNTRY = country, item = "Deaths_1x1", fixup = T)
    nx <- readHMDweb(username = username, password = password, CNTRY = country, item = "Exposures_1x1", fixup = T)

    }

  dx <- dx[dx$Year == year,]
  nx <- nx[nx$Year == year,]

  if(is.null(max)){max <- max(dx$Age)}

  dx <- dx[dx$Age <= max,]
  nx <- nx[nx$Age <= max,]

  x <- dx$Age

  if(sex == "females"){dx <- dx[,3] ; nx <- nx[,3]}
  if(sex == "males"){dx <- dx[,4]; nx <- nx[,4]}
  if(sex == "total"){dx <- dx[,5]; nx <- nx[,5]}

  data <- data.frame(x = x, d = dx, n = nx)

  data$m <- dx / nx

  }

  if(dim == "cohort"){

    if(xtra == TRUE & is.null(path)){warning("The extrapolation option is only compatible with locally stored HMD data.")}

    cnx <- readHMD(filepath = paste(path,country,"STATS","cExposures_1x1.txt",sep = "/"), fixup = T)
    cmx <- readHMD(filepath = paste(path,country,"STATS","cMx_1x1.txt",sep = "/"), fixup = T)

    cnx <- cnx[cnx$Age <= max,]
    cmx <- cmx[cmx$Age <= max,]

    if(xtra == TRUE){cmx <- xtp(path = path, max = max, country = country, year = year, sex = sex, cmx = cmx)}else{
      cmx <- cmx[cmx$Year == year,]
    }
    cnx <- cnx[cnx$Year == year,]

    x <- cmx$Age

    if(sex == "females"){cmx <- cmx[,3] ; cnx <- cnx[,3]}
    if(sex == "males"){cmx <- cmx[,4]; cnx <- cnx[,4]}
    if(sex == "total"){cmx <- cmx[,5]; cnx <- cnx[,5]}

    if(xtra == TRUE){ for(i in which(is.na(cnx))){cnx[i] <- cnx[i-1] - cnx[i-1] * cmx[i]}  }

    cdx <- round(cmx * cnx)
    cmx <- cdx / cnx

    data <- data.frame(x = x, d = cdx, n = cnx, m = cmx)

  }

  return(data)

}
