
#' @title Split abridged death counts into single-age data
#'
#' @description
#' This function applies spline smoothing and interpolation in order to transform abridged age- and cause-specific death counts and exposures into single-age data.
#'
#' @param dxc age- and cause-specific death counts by age group
#' @param nx age-specific exposures by age group
#' @param inter matrix containing the begining and end of the age intervals
#' @param plot should a plot comparison of abridged and unabridged data be produced (logical)
#'
#' @details
#'
#' The method used is inspired by the \href{http://www.mortality.org/Public/Docs/MethodsProtocol.pdf}{method protocol} of the Human Mortality Database, albeit using a built-in spline smoother with fewer constraints.
#' The age- and cause-specific death counts are first transformed into age- and cause-specific cumulative death counts.
#' Smoothing cubic splines are then fitted to these cause-specific cumulative distributions, from which singe-year estimates are interpolated.
#' Single-age death counts are then computed by taking the first derivatives (slopes) of the interpolated cumulative values.
#'
#' The first age group (age 0) is treated separately in order to facilitate the estimation of the spline.
#' Practically, it is removed from the series before the interpolation and subsequently added to the unabridged death counts.
#'
#' Exposures are treated in a similar fashion, keeping the first age group together with the rest of the series.
#'
#' If \code{plot == TRUE}, two plots are produced, from the cumulative and original age-specific death counts respectively.
#' Both plots contain abridged (black) and unabridged (red) figures. In the second plot, the unabrided figures are inflated by the size of the original age interval to allow for a visual inspection of the quality of fit.
#' These graphs may include a few small negative death counts, as the spline fitted on the cumulative death counts is not constrained to be monotonically increasing. These are later corrected in the output.
#'
#' @return
#' The function returns a list of three objects.
#' \describe{
#'  \item{dxc}{a matrix of cause-specific single-age death counts (rows = ages, columns = causes)}
#'  \item{nx}{a vector of single-age exposures}
#'  \item{x}{a vector of the age (begining of the single-age intervals)}
#' }
#'
#' @examples
#'
#' data(USA2000m)
#'
#' unabr <- unabridge(dxc = USA2000m$dxc, nx = USA2000m$nx, inter = USA2000m$inter, plot = TRUE)
#'
#' @seealso
#'
#' \link{HCD2MH}, \link{smooth.spline}
#'
#' @export
#'

unabridge <- function(dxc, nx, inter, plot = FALSE){

    fun <- function(v, .inter){

    x <- c(.inter[1,1],.inter[,2])
    y <- c(0,cumsum(v))

    if(v[1] > 0.5 * sum(v)){

      xx <- min(x):(max(x)-1)
      yy <- rep(v / apply(.inter,1,diff), times = apply(.inter,1,diff))
      ycum <- c(0,cumsum(yy))

    }else{

    spl <- smooth.spline(x = x, y = y, all.knots = TRUE, spar = -1)
    ycum <- predict(spl, x = min(x):max(x))$y

    yy <- diff(ycum)
    xx <- min(x):(max(x)-1)

    }

    if(plot == TRUE){
    par(mfrow = c(1,2))
    plot(x, y, type = "s", las = 1, xlab = "age", ylab = "cdf")
    lines(min(x):max(x), ycum, lty = 3, col = 2, type = "s")

    plot(.inter[,1], diff(y), type = "s", xlab = "age", ylab = "pdf", las = 1)
    lines(xx, yy * rep(apply(.inter,1,diff), times = apply(.inter,1,diff)),
          type = "s", col = 2, lty = 3)
    par(mfrow = c(1,1))
    }

    return(yy)

    }

    dxc.unabr <- apply(dxc[-1,],2,fun,inter[-1,])

    dxc.unabr[dxc.unabr < 0] <- 0

    dxc.unabr <- rbind(dxc[1,],dxc.unabr)

    dxc.unabr <- as.data.frame(dxc.unabr)

    names(dxc.unabr) <- names(dxc)

    nx.unabr <- fun(nx, inter)
    nx.unabr[nx.unabr < 0] <- 0

    x <- inter[1,1]:(max(inter[,2])-1)

    return(list(dxc = dxc.unabr, x = x, nx = nx.unabr))

}
