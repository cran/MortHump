#' @title Confidence Intervals for morthump fit
#'
#' @description Computes confidence intervals for life expectancy, which can be used to assess the statistical significance of the young adult mortality hump.
#'
#' @param object the result of a \code{morthump} fit.
#' @param parm the parameter for which the confidence intevals are estimated, in this case e0.
#' @param method either \code{"chiang"} or \code{"MC"} (see details)
#' @param level the confidence level required.
#' @param iter number of iterations for the \code{"MC"} method.
#' @param ... further arguments passed to or from other methods.
#'
#' @details
#' This summary function applies to a \code{morthump} object, resulting from a call to the \code{morthump} function.
#' Using the fitted force of mortality and observed exposures, it computes confidence intervals for life expectancy at birth.
#' This value can then be compared with the number of years of life expectancy lost to the hump, as estimated by the \code{morthump} model, in order to assess whether the latter departs statistically from the stochastic noise.
#'
#' Two methods are available. When \code{method = "chiang"}, an algebraic solution is used (Chiang 1978).
#' When \code{method = "MC"}, a Monte-Carlo simulation approach is used (Andreev & Shkolnikov 2010). They usually reach very close results.
#'
#' @return
#'
#' \describe{
#'  \item{e0}{life expectancy estimated on the fitted force of mortality.}
#'  \item{lower}{lower bound for e0.}
#'  \item{upper}{upper bound for e0.}
#'  }
#'
#' @references
#' Chiang, C. L. (1978). Life table and mortality analysis. World Health Organization.
#' Andreev, E. M., & Shkolnikov, V. M. (2010). Spreadsheet for calculation of confidence limits for any life table or healthy-life table quantity. Rostock: Max Planck Institute for Demographic Research (MPIDR Technical Report, 5.
#' Brown, L. D., Cai, T. T., & DasGupta, A. (2001). Interval Estimation for a Binomial Proportion. Statistical Science, 16(2), 101-117.
#' Eayres, D., & Williams, E. (2004). Evaluation of methodologies for small area life expectancy estimation. Journal of epidemiology and community health, 58(3), 243-249.
#' Toson, B., & Baker, A. (2003). Life expectancy at birth: methodological options for small populations. National statistics methodological series, 33.
#'
#' @importFrom stats rbinom qnorm quantile
#' @export
#' @method confint morthump


confint.morthump <- function(object, parm = "e0", level = 0.95, method, iter = 1000,...){

    fit <- object
    loss <- summary(fit)$loss

    # function to shift a variable
    shift <- function(x,shift_by){
      #http://www.r-bloggers.com/generating-a-laglead-variables/
      stopifnot(is.numeric(shift_by))
      stopifnot(is.numeric(x))

      if (length(shift_by)>1)
        return(sapply(shift_by,shift, x=x))

      out<-NULL
      abs_shift_by=abs(shift_by)
      if (shift_by > 0 )
        out <- c(tail(x, -abs_shift_by), rep(NA, abs_shift_by))
      else if (shift_by < 0 )
        out <- c(rep(NA, abs_shift_by), head(x, -abs_shift_by))
      else
        out <- x
      out
    }

    # quick function to compute life expectancy (quicker than LT)
    LE <- function(mx){
      if (is.null(dim(mx))){
        mx <- as.matrix(mx)
      }
      N             <- nrow(mx)
      ax            <- mx * 0 + .5
      ax[1, ]       <- 0.1
      qx            <- mx / (1 + (1 - ax) * mx)
      qx[N, ]       <- 1
      qx[qx > 1]    <- 1
      px            <- 1 - qx
      lx            <- apply(px, 2,function(.px,.N){
        c(1,cumprod(.px)[-.N])
      }, .N = N)
      dx            <- lx * qx
      dx[N, ]       <- 1 - colSums(dx[-N, , drop = FALSE])
      Lx            <- rbind(
        lx[2:N, , drop = FALSE] + ax[1:(N -  1), , drop = FALSE] * dx[1:(N - 1), , drop = FALSE],
        mx[N,] / ax[N,]
      )
      Lx[is.infinite(Lx)] <- 1
      Lx[is.na(Lx)] <- 0
      Tx            <- apply(Lx, 2,function(.Lx){
        rev(cumsum(rev(.Lx)))
      })
      ex            <- Tx / lx
      return(ex[1, ])
    }

  data <- fit$data

  if(parm != "e0"){warning("Confidence intervals are only available for e0.")}
  if(sum(data$n) < 5000){message("Population under exposure lower than 5000. The variance estimates will not be reliable (see Eayres & William 2004 and Toson & Baker 2003:10).")}
  if(any(data$d == 0)){message("At least one age group has zero deaths. This may harm the computation of the variance, but the effect for populations over 5000 is minimal (see Toson & Baker 2003:17).")}

  stopifnot(method %in% c("MC", "chiang"))

  Mx    <- summary(fit, print = FALSE)$mhat
  n     <- length(Mx)
  Nx    <- round(data$n)
  Dx    <- Nx * Mx
  e0    <- LE(Mx)

  if(method == "MC"){
    rmx         <- t(mapply(function(.Mx, .N, iter){
                      rbinom(n = iter, prob = .Mx, size = .N)
                      },Mx,Nx,iter=iter)) / Nx
    re0         <- LE(rmx)
    lower       <- quantile(re0, probs = (1 - level) / 2)
    upper       <- quantile(re0, probs = 1 - (1 - level) / 2)
  }

  if(method == "chiang"){
    if(any(Dx < 5)){
      message("Number of deaths is too low in at least one age group to use Chiang's method. Wald's approximation for the variance of Mx needs at least 5 deaths for each age group (see Brown, Cai & DasGupta 2001).")
    }
    ax      <- rep(0.5, length(Mx))
    ax[1]   <- 0.1
    qx      <- Mx/(1 + (1 - ax) * Mx)
    qx[n]   <- 1
    qx[qx > 1] <- 1
    px      <- 1 - qx
    lx      <- c(1,cumprod(px)[-n])
    dx      <- -diff(lx)
    dx[n]   <- lx[n]

    Lx      <- c(lx[2:n] + ax[1:(n -  1)] * dx[1:(n - 1)], Mx[n] / ax[n])
    Lx[is.infinite(Lx)] <- 1
    Lx[is.na(Lx)] <- 0
    Tx      <- rev(cumsum(rev(Lx)))
    ex      <- Tx/lx
    varqx   <- (Mx * (1 - ax * Mx)) / (Nx * (1 + (1 - ax) * Mx) ^3) # Chiang 2
    yx      <- lx^2 * ((1 - ax) + c(shift(ex,1)[-n],0)) ^2 * varqx
    cumyx   <- rev(cumsum(rev(yx)))
    varex   <- cumyx / (lx^2)
    seex    <- sqrt(varex)
    lower   <- ex[1] + qnorm((1 - level) / 2) * seex[1]
    upper   <- ex[1] - qnorm((1 - level) / 2) * seex[1]
  }

  cat("Loss of life expectancy due to the hump (years):",loss)
  cat("\n")
  cat("Half-confidence interval of fitted life expectancy at birth:",round((e0-lower)/2,3))
  cat("\n")
  if(loss > (e0-lower)/2){cat("The hump is statistically significant")}else{
    cat("The hump is not statistically significant")
  }
  cat("\n")
  cat("Significance level:",level)

  return(list(e0 = e0, upper = upper, lower = lower))

}







