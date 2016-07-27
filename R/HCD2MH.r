
#' @title Format data from the Human Cause-of-Death Database for use in the MortHump package
#'
#' @description
#' This is a data grabber for the Human Cause-of-death database that select and format age- and cause-specific death counts, exposures, rates and cause-of-death labels.
#' It formats the data so that they can be properly used by \code{codhump} function.
#'
#' @param country HCD population letter code
#' @param year year for which the data is required
#' @param sex sex for which the data is required ("females", "males" or "total")
#' @param list cause-of-death classification ("short", "interm" or "full")
#' @param unabr should the data be unabridged into single-year values
#' @param path path to the HCD folder for local access
#'
#' @details
#' This function assumes that all the HCD data are stored into a single folder, containing itself country-specific folders.
#' To constitute this folder, download the zipped \href{http://www.causesofdeath.org/cgi-bin/datazip.php}{datasets by country}, and unzip the country folders.
#' The list of available countries and period coverage can be found on the \href{http://www.causesofdeath.org/cgi-bin/data.php}{HCD website}.
#'
#' The \href{http://www.causesofdeath.org/docs/interm_list.pdf}{Intermediate} and \href{http://www.causesofdeath.org/docs/short_list.pdf}{short} cause-of-death classifications can be found on the HCD website.
#' The full list depends on each country and can be found on the country pages.
#'
#' @return HCD2MH returns a list of six elements:
#'  \describe{
#'    \item{mxc}{a matrix of cause- and age-specific death rates (age in rows, causes in columns)}
#'    \item{dxc}{a matrix of cause- and age-specific death counts (age in rows, causes in columns)}
#'    \item{nx}{a vector of age-specific exposures}
#'    \item{x}{a vector of mid-point for each age interval}
#'    \item{age}{a character vector of age labels}
#'    \item{inter}{a matrix containing the endpoints of the age intervals}
#'    \item{lab}{a character vector of short names for the causes of death}
#'  }
#'
#' @references
#'
#' Human Cause-of-Death Database. French Institute for Demographic Studies (France) and Max Planck Institute for Demographic Research (Germany). Available at \url{www.causeofdeath.org} (data downloaded in May 2016).
#'
#' @seealso
#'
#' \link{codhump}
#'
#' @export
#'

HCD2MH <- function(country, year, sex, list, unabr = FALSE, path){

  if(sex %in% c("males","females","total") == F){warning("sex must be one of 'females', 'males' or 'total'")}

  mxc <- read.csv(file = file.path(path,country,paste(country,"_m_",list,"_idr.csv",sep = "")))
  dxc <- read.csv(file = file.path(path,country,paste(country,"_d_",list,"_idr.csv",sep = "")))
  nx  <- read.csv(file = file.path(path,country,paste(country,"_e.csv",sep = "")))

  if(sex == "females"){sex <- 2}
  if(sex == "males"){sex <- 1}
  if(sex == "total"){sex <- 3}

  mxc <- mxc[mxc$year == year & mxc$sex == sex,]
  dxc <- dxc[dxc$year == year & dxc$sex == sex,]
  nx  <- nx[nx$year == year & nx$sex == sex,]

  af <- unique(c(mxc$agf,dxc$agf))
  if(length(af) > 1){warning("Two different age formats are present in the selected dataset.")}
  #data("agf")
  age <- as.character(agf[,af])

  if(0 %in% mxc$cause){mxc <- mxc[mxc$cause != 0,] ; dxc <- dxc[dxc$cause != 0,]}
  causes <- mxc$cause

  mxc <- as.data.frame(t(mxc[,7:ncol(mxc)]))
  suppressWarnings(mxc <- as.data.frame(apply(mxc, 2, as.numeric)))
  mxc <- mxc / 1e6
  mxc <- mxc[rowSums(is.na(mxc)) < ncol(mxc),]
  if(af == 2){mxc <- mxc[-19,]}
  if(af == 3){mxc <- mxc[-c(19,21),]}
  if(af == 4){mxc <- mxc[-c(19,21,23),]}
  names(mxc) <- causes
  mxc <- mxc[,colSums(is.na(mxc)) < nrow(mxc)]


  dxc <- as.data.frame(t(dxc[,8:ncol(dxc)]))
  suppressWarnings(dxc <- as.data.frame(apply(dxc, 2, as.numeric)))
  dxc <- dxc[rowSums(is.na(dxc)) < ncol(dxc),]
  if(af == 2){dxc <- dxc[-19,]}
  if(af == 3){dxc <- dxc[-c(19,21),]}
  if(af == 4){dxc <- dxc[-c(19,21,23),]}
  names(dxc) <- causes
  dxc <- dxc[,colSums(is.na(dxc)) < nrow(dxc)]


  nx <- t(nx[1,6:30])
  if(af == 1){nx <- nx[-c(20:length(nx))]}
  if(af == 2){nx <- nx[-c(19,22:length(nx))]}
  if(af == 3){nx <- nx[-c(19,21,24:length(nx))]}
  if(af == 4){nx <- nx[-c(19,21,23)]}

  if(af == 1){x <- age[-c(20:length(age))]}
  if(af == 2){x <- age[-c(19,22:length(age))]}
  if(af == 3){x <- age[-c(19,21,24:length(age))]}
  if(af == 4){x <- age[-c(19,21,23)]}
  age <- x
  plus <- grep(x = as.character(x), pattern = "+", fixed = TRUE)
  x[plus] <- sub(x = x[plus], pattern = "+", replacement = "", fixed = TRUE)
  inter <- do.call(rbind,strsplit(as.character(x), "-"))
  inter <- cbind(apply(inter,2,as.numeric))
  inter[,2] <- inter[,2] + 1
  inter[nrow(inter),2] <- inter[nrow(inter),1] + 10
  x <- rowMeans(inter)

  if(unabr == TRUE){

    ua <- unabridge(dxc = dxc, nx = nx, inter = inter)
    x <- ua$x + 0.5
    dxc <- ua$dxc
    nx <- ua$nx
    mxc <- dxc / nx
    inter <- cbind(x[-length(x)],x[-1])

    }

  #data("lists")
  if(list == "short"){lab <- as.data.frame(short)}
  if(list == "interm"){lab <- as.data.frame(interm)}
  if(list == "full"){message("The labels of the full list must be downloaded from the HCD website as they differ by country.") ; lab <- NA}
  lab <- lab[match(names(mxc),lab$no),]

  return(list(mxc = mxc, dxc = dxc, nx = nx, x = x, age = age, inter = inter, lab = lab))

}
