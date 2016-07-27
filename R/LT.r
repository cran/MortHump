#' @name LT
#' @title A function to build lifetables and extract certain lifetable statistics from data.
#'
#' @description Accepts either single age or 5-year abridged data. Accepts many optional arguments, such as differing methods for a(x) estimation, optional smoothing for M(x) or a(x) values, a changeable radix.
#'
#' @param Nx numeric vector of population exposures by age.
#' @param Dx numeric vector of death counts by age.
#' @param Mx numeric vector of central death rates (assumed in the function to be the lifetable m(x)), calculated as deaths/exposure.
#' @param ages optional, default = \code{"auto"}. If \code{"auto"}, the function tried to guess whether you have single-age data or 5-year abridged age groups. These are used in the lifetable labels, and do not enter into calculations. Otherwise, the user can specify a vector (character or numeric).
#' @param axmethod either \code{"keyfitz"}, \code{"schoen"}, \code{"preston"} or \code{"midpoint"}. Default = \code{"keyfitz"}, although this is not recommended for abridged ages. See comparisons in \code{\link{axEstimate}} examples. The user can also supply a numeric vector of a(x) values here (e.g. from a different estimation procedure or from a different population).
#' @param sex either \code{"male"} or \code{"female"} (default). It is only necessary to specify this if \code{"preston"} is the \code{axmethod}. It does not affect any other lifetable calculations.
#' @param mxsmooth logical, default = \code{TRUE}. Should the mx vector be smoothed? If \code{TRUE} and both \code{Nx} and \code{Dx} vectors are supplied (the ideal case), smoothing is done using the function \code{Mort1Dsmooth()} from Giancarlo Camarda's \code{MortalitySmooth} package. In this case, \code{Dx} values are smoothed using \code{log(Nx)} as an offset, and all other items are the function defaults. If \code{Mx} is provided instead of \code{Nx} and \code{Dx} a loess smoother is used, \code{loess}, with span set to .15 for single age data and .4 for 5-year abridged data. If these smoothing procedures are not satisfactory, the user may wish to pre-process the Mx estimate and specify \code{mxsmooth = FALSE}, or else leave it rough.
#' @param axsmooth logical, default = \code{TRUE}. Ignored if \code{mxsmooth = TRUE}. Should the a(x) values be calculated from a smoothed M(x) series? In this case, the M(x) series is smoothed within the \code{axEstimate()} function for a(x) estimation, but the smoothed M(x) function that was used is not returned. In general, it is better to smooth the M(x) function prior to putting it in this function, because the loess smoother used here has no weights or offset. If this is not possible, loess M(x) smoothing still produces more consistent and less erratic a(x) estimates. If \code{mxsmooth = FALSE} and \code{axsmooth = TRUE}, the Mx series is only smoothed for use in a(x) estimation, and does not affect any other lifetable calculations that are dependent on Mx.
#' @param radix The lifetable starting population at age 0, l(0). default = 1. Other common values are 1000 and 100000, although any value may be given.
#' @param verbose logical, default = \code{TRUE}. Should informative but possibly annoying messages be returned when the function does something that you might want to know about?
#' @param ... further arguments passed to or from other methods
#'
#' @details Either \code{Nx} must be specified along with \code{Dx}, *or* \code{Mx} must be specified directly. If smoothing is used, it is better to specify both \code{Nx} and \code{Dx}, since the \code{Nx} vector can be used as an offset in the \code{MortalitySmooth} smoother.
#'
#' @return a list is invisibly returned
#' \itemize{
#'   \item \code{LT} A \code{data.frame} with 11 columns and as many rows as you have ages. The columns are "Age" (character age labels, i.e. with "+"), "ages" (numeric), \code{"mx", "ax", "qx", "px", "lx", "dx", "Lx", "Tx", "ex"}.
#' All the individual components of the lifetable can be called and are unrounded when individually called. Some additional values are also available:
#'   \item \code{Age} character vector of ages. Denotes intervals in the case of an abridged table.
#'   \item \code{ages} numeric vector of ages. Left side of interval.
#'   \item \code{mx} the lifetable mx (may differ from a given Mx).
#'   \item \code{ax} Chiang's ax, either given by the user or estimated in \code{axEstimate}.
#'   \item \code{qx} typical lifetable qx. Death probability for interval x, x + n.
#'   \item \code{lx} typical lifetable lx. Number of radix individuals entering age x. l(0) = the radix population.
#'   \item \code{dx} typical lifetable dx. When \code{radix = 1} (default), this is the probability density function of deaths.
#'   \item \code{Lx} typical lifetable Lx. Lifetable exposure for the interval x,x+n.
#'   \item \code{Tx} typical lifetable Tx. Total number of years remaining to be lived by the cohort entering age x.
#'   \item \code{ex} typical lifetable ex. Life remaining life expectancy at age x. e(0) = life expectancy at birth. Two other estimates of e(0) are given below.
#'   \item \code{Sx} probability of surviving from age x until age x + n (L_{x+n}/L_{x}).
#'   \item \code{Widths} vector of age intervals (n).
#' }
#'
#' @references
#' The main reference for this function has been:
#'
#' Preston et al (2001). Demography: Measuring and Modelling Population Processes. Blackwell Publishing
#'
#' ax estimation also received input from:
#'
#' Chiang C.L.(1968) Introduction to Stochastic Processes in Biostatistics. New York: Wiley.
#'
#' Coale Anseley and Paul Demeny, with B Vaughan (1983). Regional Model Life Tables and Stable Populations. New York Academic Press.\\
#' Keyfitz, Nathan (1966) A Life Table that Agrees with the Data. Journal of the American Statistical Association, 61 (314):305-12. As described on page 44-45 of: Schoen R. (1978) Calculating lifetables by estimating Chiang's a from observed rates. Demography 15: 625-35.
#'
#' function calls \code{MortalitySmooth}: Carlo G Camarda (2009) MortalitySmooth: Smoothing Poisson counts with P-splines. (version 2.3 at the time of this writing) \url{http://CRAN.R-project.org/package=MortalitySmooth}.
#'
#' @seealso \code{\link{MortalitySmooth}}, \code{\link{axEstimate}}.
#'
#' @examples
#'
#' data(CHE2010m)
#' attach(CHE2010m)
#' lt <- LT(Nx = n, Dx = d, ages = x, mxsmooth = FALSE)
#' lt$ex[1] # life expectancy at birth
#'
#'
#' @author Tim Riffe
#'
#' @importFrom MortalitySmooth Mort1Dsmooth
#' @importFrom stats loess
#'
#' @export
#'

LT <- function(Nx=NULL, Dx=NULL, Mx = Dx/Nx, ages = 0:(length(Mx)-1), axmethod = "midpoint", sex = "female",
           mxsmooth = TRUE, axsmooth = TRUE, radix = 1, verbose = TRUE, ...){
    # the verbose function:
    Verb <- function(v, x){
      if (v) {
        cat(paste0(x, "\n"))
      }
    }
    # first a series of checks, messages and imputations to make sure the given Dx&Nx *or* Mx values are viable
    otherArgs <- list(...)
    if (length(otherArgs)>1){
      Verb(verbose,paste(unlist(otherArgs),collapse = ", "),"is/are not arguments to this function,\nalthough they may have been in the past. In that case, consider them deprecated")
    }
    if (length(Mx) == 0){
      # two checks that will stop the function
      # 1) in absence of Mx, need both Dx and Nx
      if (is.null(Nx) | is.null(Dx)){
        Verb(verbose,"you're missing the key arguments")
        stop("either specify Mx directly, or else specify both Nx and Dx.")
      }
      # 2) both Nx and Dx must be of equal length
      stopifnot(length(Nx) != length(Dx))

      # safety to avoid zeros in the denominator. will make Mx of 1 in those cases
      if (any(Nx == 0)) {
        Verb(verbose, "there was at least one 0 in your Nx vector, imputed with corresponding Dx.")
        Nx[Nx == 0] <- Dx[Nx == 0]
      }
      Mx <- Dx / Nx
    }

    # by this point we should have an Mx vector with no NAs: no more need for Nx,Dx
    # we want to be able to accept 0s...


    # N is just used for counting, to save space
    N                   <- length(Mx) #
    Widths              <- diff(ages)
    Widths              <- c(Widths, Widths[N - 1])
    # define character for Age in formatted lifetable
    if (all(Widths == 1)){
      Age <- as.character(ages)
      Age[N] <- paste0(Age[N],"+")
    } else {
      Age <- c(paste(ages[1:(N - 1)], ages[1:(N - 1)] + Widths[1:(N-1)] - 1,sep="-"),paste0(ages[N],"+"))
    }

    ages.mids.pre 		<- ages + Widths / 2
    ages.mids.pre[1] 	<- .1

    if(!is.null(Nx) & !is.null(Dx) & mxsmooth){
      # Giancarlo's package. I recommend either supplying an already-smoothed Mx vector (for complete control)
      # or else supplying Dx and Nx and setting mxsmooth to TRUE.
      # define some reasonable weights to block out certain cells if necessary
      w <- ifelse(Dx/Nx == 0 | is.null(Dx/Nx) | is.na(Dx/Nx) | log(Nx) < 0, 0, 1)
      fitBIC 			<- Mort1Dsmooth(x = ages.mids.pre, y = Dx, offset = log(Nx), w = w)
      Mx[2:N] 		<- (fitted(fitBIC) / Nx)[2:N]
    }

    if(is.null(Nx) & is.null(Dx) & mxsmooth){
      Verb(verbose,"mxsmooth was specified as TRUE, but since Mx was supplied directly, \nthere are no implicit weights (Nx). Function used a loess smoother \nto smooth out the Mx, but please be wary.")
      span 			<- ifelse(N > 30, .15, .4)
      logMx 			<- log(Mx)
      logMx[is.infinite(logMx)] <- NA
      logMx[2:N] 		<- predict(loess(logMx ~ ages.mids.pre, span = span, control = loess.control(surface = "interpolate")), newdata = ages[2:N])
      Mx 				<- exp(logMx)
    }

    # Assume that the lifetable mx is equal to the central death rate.
    mx                  <- Mx
    # these later 2 imputations should not be needed, but are there just in case.
    mx[is.na(mx)] 		<- 0
    mx[is.infinite(mx)] <- 0

    # we don't want 0s in m(x) for calculating a(x), because a(x) (except midpoint) is iterative
    # and we'd erroneously bring down neighboring age's ax values with zeros.
    # for later calulations, the zeros are taken 'as-is'
    Ind0 <- NULL
    if (length(axmethod) == 1 & min(mx) == 0){
      Ind0        <- mx == 0

      Verb(verbose, paste("\n\n*there were some ages (", ages[Ind0],
                          ") with no mortality.\nValues are imputed for calculating a(x), but the zeros are kept for the rest of the lifetable.\n"))
      span        <- ifelse(N > 30,.15,.4)
      logMx       <- log(mx)
      logMx[is.infinite(logMx)] <- NA
      Imp         <- exp(predict(loess(logMx ~ ages.mids.pre, span = span, control = loess.control(surface = "interpolate")), newdata = ages))[Ind0]
      if (any(is.na(Imp))){
        Imp     <- exp(spline(ages.mids.pre, logMx, xout=ages.mids.pre)$y[Ind0])
      }

      mx[Ind0]    <- Imp
    }
    ax <- NULL
    if (length(axmethod) == 1){
      if (axmethod %in% c("keyfitz", "schoen", "midpoint", "preston")){
        # here, the ax iteration is externalized to axEstimate()
        if (mxsmooth){
          Verb(verbose,"\nif mxsmooth = TRUE, then we don't smooth ax\n")
          axsmooth <- FALSE
        }
        ax <- axEstimate(Mx = mx, n = Widths, axsmooth = axsmooth, method = axmethod, sex = sex, verbose = verbose)
      }
      if (axmethod == "keyfitz" & !all(Widths == 1)){
        Verb(verbose, "\nIt appears you have an abridged lifetable, but have specified the keyfitz method of ax estimation.\nBe aware that this method presumes equal age intervals, as Preston et. al. (2001)\n warn on page 45. Consider using a different method or else specifying your own ax vector.\n Function continued nonetheless.")
      }
    }

    if (is.numeric(axmethod) & length(axmethod) == N) {
      ax <- axmethod
    }
    # last default # sex = "male"
    if (is.null(ax)){
      ax <- axEstimate(Mx = mx, n = Widths, axsmooth = axsmooth, method = "midpoint", sex = sex, verbose = verbose)
      Verb(verbose, "axmethod must be specified either as 'schoen','keyfitz','midpoint'\nor as a numeric vector the same length as Nx.\nThis wasn't the case, so the function defaulted to the midpoint method.")
    }

    # if zeros were imputed for ax estimation, then we put them back for the rest of the calculations
    if (is.null(Ind0)){
      mx[Ind0]        <- 0
    }

    # now start the lifetable columns:
    qx 			        <- (Widths * mx) / (1 + (Widths - ax) * mx)
    qx[N] 		        <- 1

    # can't have qx > 1, so we impute 1s where necessary: hopefully only at the penultimate, as the case may be
    qx[qx > 1] 	        <- 1
    px 			        <- 1 - qx

    lx                  <- c(radix, radix * cumprod(px))[1:N]
    dx 			        <- c(-diff(lx),lx[N])

    Lx 	                <- c(Widths[1:(N - 1)] * lx[2:N] + ax[1:(N - 1)] * dx[1:(N - 1)], lx[N] * ax[N])
    Lx[is.infinite(Lx)] <- 1
    Lx[is.na(Lx)] 	    <- 0

    Tx 			        <- rev(cumsum(rev(Lx)))
    ex 			        <- Tx / lx
    ex[N]               <- ifelse(mx[N] == 0, ax[N], {1 / mx[N]})

    # Sx is the pertinent output for projection matrices
    Sx   	            <- c((Lx[2:N] / Widths[2:N]) / (Lx[1:(N-1)] / Widths[1:(N - 1)]), Tx[N] / Tx[(N - 1)])

    # two not-very satisfying, and possibly redundant, substitutions:
    Sx[Lx == 0]   	    <- 0
    Sx[is.na(Sx)] 	    <- 0

    LT <- data.frame(cbind("Age" = Age, "ages" = ages, "mx" = round(mx, 5), "ax" = round(ax, digits = 2),
                           "qx" = round(qx, 5), "px" = round(px, 5), "lx" = round(lx, 5),
                           "dx" = round(dx, 5), "Lx" = round(Lx, 5), "Tx" = round(Tx, 5), "ex" = round(ex, 3)))
    # both LT as well as the individual pieces (not rounded) can be called
    out <- list(LT = LT,
                Age = Age,
                ages = ages,
                mx = mx,
                ax = ax,
                qx = qx,
                lx = lx,
                dx = dx,
                Lx = Lx,
                Tx = Tx,
                ex = ex,
                Sx = Sx,
                Widths = Widths)
    invisible(out)
  }



