#' @title Summary method for morthump models
#'
#' @description Produces a summary of a fitted morthump model
#'
#' @param object the result of a \code{morthump} fit
#' @param x an object of class "summary.morthump", usually, a result of a call to summary.morthump
#' @param ... further arguments passed to or from other methods
#'
#' @details
#' This summary function applies to a \code{morthump} object, resulting from a call to the \code{morthump} function.
#' Using the fitted age-specific contributions to the hump, it computes summary measures about the intensity, centrality and dispersion of the hump.
#'
#' \bold{intensity}
#'
#' The indensity of the hump is measured in terms of the loss of life expectancy implied by the presence of the hump.
#' It is defined as the difference between life expectancy computed on the fitted force of mortality (sum of the three components), and life expectancy computed on the hump-free force of mortality.
#' This difference can be interpreted as the mean years of life lost in the population due to the presence of the hump, or alternatively as the potential gain in life expectancy that could be reached in the absence of the hump.
#'
#' Other measures of intensity are computed to offer a public health perspective to the presence of the hump.
#' The years of life lost (YLL), which represents the total number of years that could have been lived by those who died because of the hump.
#' The number of deaths due to the hump is simply the number of deaths that could have been averted in the absence of the hump.
#'
#' \bold{centrality}
#'
#' Three measures of centrality are currently proposed: the modal, median and mean age at death from the hump for each cause. By considering the hump as a density, standard summary measures can be applied to it.
#' The mode is defined as the age at which is observed the highest number of deaths due to the hump, and is therefore not sensitive to the proportion of deaths recorded at other ages.
#' The median is defined as the age at which half of the deaths due to the hump have occured. It is more sensitive to the overall shape of the hump, but not as much as the mean.
#' The mean is defined as the mean age at death for the people who die because of the hump, and is therefore sensitive to the general shape of the hump.
#' The difference between the mode and the mean can inform on the symetry of the hump, e.g. if the mode is located before the mean this suggests the presence of a flatter hump after the peak.
#'
#' \bold{dispersion}
#'
#' Dispersion is measured by the standard deviation of the age at which people die because of the hump.
#'
#' @return
#'
#' \describe{
#'  \item{mhat}{vector containing the fitted age-specific death rates.}
#'  \item{mhump}{vector containing the fitted age-specific contributions to the hump.}
#'  \item{pval}{when using a parametric model, this is the p-value associated to the F-test against a model with no hump (\code{hps}). It is formally the probability to reject the superiority of the model with a hump given that there is one, and can be intepreted as the probability that the observed hump only appears as the result of stochastic noise.}
#'  \item{pdf}{function that returns the value of the probability density function of the hump at a given age. Used to compute the mode.}
#'  \item{cdf}{function that returns the value of the cumulative density function of the hump at a given age. Used to compute any quantile of the hump, including the median.}
#'  \item{qtl}{function that returns quantiles of the hump (i.e. ages at which the hump has claimed a given percentage of its deaths).}
#'  \item{loss}{loss of life expectancy due to the hump.}
#'  \item{dh}{number of lives lost due to the hump.}
#'  \item{yll}{years of life lost due to the hump.}
#'  \item{mode}{modal age at death from the hump}
#'  \item{median}{median age at death from the hump}
#'  \item{mean}{mean age at death from the hump}
#'  \item{sd}{standard deviation for the age at death from the hump}
#'  \item{par}{parameters of the model (in case of a \code{sse} model, these are the coefficients of the base spline)}
#'  }
#'
#'
#' @export
#' @method summary morthump


summary.morthump <- function(object,...){

   fit <- object

   model <- attr(fit, "model")
   data <-  fit$data
   x <- data$x
   if(x[1] == 0){
     if(model %in% c("hp","hpk","hps")){x[1] <- 1e-5}else{x[1] <- 0.9}
   }


   # mx fitted, hump, and w/o hump
   mhat <- predict(object = fit, x = x, which = "fitted")
   mhump <- predict(object = fit, x = x, which = "hump")
   mhyp <- mhat - mhump

   par1 <- coef(object = fit)

  if(model %in% c("hp","hpk")){

     # Fisher test compared to a model w/o hump
     w <- fit$w
     m   <- data$m
     fit2 <- hps.fit(data = data, w = w)
     par2 <- coef(object = fit2)
     mhat2 <- predict(object = fit2)
     df1 <- length(par1) - length(par2)
     df2 <- length(x) - length(par)
     res1 <- sum(w * (m - mhat)^2, na.rm = TRUE)
     res2 <- sum(w * (m - mhat2)^2, na.rm = TRUE)
     test <- ((res2 - res1) / df1) / (res2 / df2)
     suppressWarnings(pval <- pf(test, df1 = df1, df2 = df2, lower.tail = FALSE))
     #pval <- stats:::anovalist.nls(fit,fit2)$"Pr(>F)"[2]

   }else{pval <- NA}

   # loss of life expectancy
   e0hat <- LT(Mx = mhat, type = "single-age", mxsmooth = FALSE, verbose = FALSE)$ex[1]
   e0hyp <- LT(Mx = mhyp, type = "single-age", mxsmooth = FALSE, verbose = FALSE)$ex[1]
   e0hump  <- round(e0hyp - e0hat,2)

   # loss of lives
   dh <- round(sum(data$n * mhump))

   # years of life lost
   yll <- round((data$n * mhump) %*% LT(Mx = mhat, type = "single-age", mxsmooth = FALSE, verbose = FALSE)$ex)

   if(model != "hps"){

   # functions
   pdf <- function(t){predict(object = fit, x = t, which = "hump", tol = 2) / sum(mhump)}
   pdf <- Vectorize(pdf)
   cdf <- Vectorize(FUN = function(t){integrate(f = function(a){pdf(a)}, lower = min(x), upper = t)$value})
   if(model == "sse"){cdf <- function(t){sum(pdf(seq(min(x),t,0.1))) / 10}}
   cdf <- Vectorize(cdf)
   qtl <- function(p){optimize(f = function(t){abs(cdf(t) - p)}, interval = range(x))$minimum}
   qtl <- Vectorize(qtl)

   # descriptive statistics
   mode.st <- x[which.max(predict(object = fit, which = "hump"))]
   mode <- optimize(f = function(x){-pdf(x)}, lower = mode.st * 0.8, upper = mode.st * 1.2)$minimum
   median <- qtl(0.5)
   mean <- as.numeric(pdf(x) %*% x)
   sd <- as.numeric(sqrt(pdf(x) %*% (data$x - mean)^2))

   }



   if(model != "hps"){
     out <- list(mhat = mhat, mhump = mhump, pval = pval, pdf = pdf, cdf = cdf, qtl = qtl, loss = e0hump, dh = dh, yll = yll, mode = mode, median = median, mean = mean, sd = sd, par = par1)
      }else{
     out <- list(mhat = mhat, mhump = mhump, pval = pval, loss = e0hump, dh = dh, yll = yll, par = par1, it = fit$it)
      }

   class(out) <- "summary.morthump"
   attr(out, "extra") <- list(model = model, type = fit$type, method = fit$method, it = fit$it, maxit = fit$maxit)

   invisible(out)
}

#' @rdname summary.morthump
#' @method print summary.morthump
#' @export

print.summary.morthump <- function(x, ...){

  model <- attr(x, "extra")$model
  type <- attr(x, "extra")$type

  cat("[[>>]] Information about the model used to fit the mortality schedule")


  cat("\nType of the model:",type)
  cat("\nName of the model:",model)
  cat("\nAlgorithm used for optimization:",attr(x,"extra")$method)
  if(type == "parametric"){
    cat("\nEstimated coeffients: ", paste0(names(x$par),"=",signif(x$par,3)))
    cat("\nF-test against a model w/o hump: ", signif(x$pval,3),"(H0: this model does not fit better than a model without a hump)")
  }else{
    cat("\nNumber of iterations:",attr(x, "extra")$it,"( max =",attr(x, "extra")$maxit,")")
  }

  cat("\n\n")

  if(model != "hps"){

    cat("[[>>]] Descriptive statistics about the hump")

    cat("\nIntensity")
    cat("\nLife expectancy at birth lost due to the hump: ", x$loss, "year")
    cat("\nDeaths due to the hump: ", x$dh)
    cat("\nYears of life lost due to the hump: ", x$yll)
    cat("\n")

    cat("\nCentrality")
    cat("\nMode:", round(x$mode,2))
    cat("\nMedian:", round(x$median,2))
    cat("\nMean:", round(x$mean,2))
    cat("\n")

    cat("\nDispersion")
    cat("\nStandard deviation:", round(x$sd,2))

  }

  invisible(x)
}

