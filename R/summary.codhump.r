#' @title Summary method for codhump models
#'
#' @description Produces a summary of a fitted codhump model
#'
#' @param object the result of a \code{codhump} fit
#' @param x an object of class "summary.codhump", usually, a result of a call to summary.codhump
#' @param ... further arguments passed to or from other methods
#'
#' @details
#' This summary function applies to a \code{codhump} object, resulting from a call to the \code{codhump} function.
#' Using the fitted age- and cause-specific contributions to the hump, it computes summary measures about the intensity, centrality and dispersion of the hump and its components.
#'
#' \bold{intensity}
#'
#' The indensity of the hump is measured in terms of the loss of life expectancy implied by the presence of the all-cause hump, as well as each cause-specific contributions.
#' For the all-cause mortality, it is defined as the difference between life expectancy computed on the fitted force of mortality (sum of both components), and life expectancy computed on the hump-free force of mortality (only the senescence component).
#' This difference can be interpreted as the mean years of life lost in the population due to the presence of the hump, or alternatively as the potential gain in life expectancy that could be reached in the absence of the hump.
#' This measure is extended to cause-specific contributions to the hump by substracting separately each of the hump components from the overall force of mortality and recomputing the life expectancy.
#' This counterfactual approach assumes implicitely the independance of the cause-specific contributions to the hump. In other words, it assumes that removing one cause's hump component will not affect the other causes' components.
#' In reality, one might argue that people who are under higher risk of death from one cause during early adulthood are also under a higher risk from other causes, and thus that the assumption of independence does not hold.
#' The loss of life expectancy should thus be interpreted with care and be considered as a purely statistical and theoretical exercice, and not be interpreted as a prediction of what could be gained if we could shelter people from the young adult excess mortality.
#'
#' \bold{centrality}
#'
#' Two measures of centrality are currently proposed: the modal and mean age at death from the hump for each cause. By considering the hump as a density, standard summary measures can be applied to it.
#' The mode is defined as the age at which is observed the highest number of deaths due to the hump, and is therefore not sensitive to the proportion of deaths recorded at other ages.
#' The mean is defined as the mean age at death for the people who die because of the hump, and is therefore sensitive to the general shape of the hump.
#' The difference between the mode and the mean can inform on the symetry of the hump, e.g. if the mode is located before the mean this suggests the presence of a flatter hump after the peak.
#'
#' \bold{dispersion}
#'
#' Dispersion is measured by the standard deviation of the age at which people die because of the hump.
#' The all-cause and the cause-specific dispersions are independant of each other and should be interpreted together with measures of centrality.
#' One may for instance observe a high all-cause dispersion combined with small cause-specific dispersions if the cause-specific humps are not centered on the same age.
#' Alternatively, if all the cause-specific humps are centered on the same age, the overall dispersion may be smaller than each of the cause-specific dispersions.
#' Missing values may be caused by remaining negative contributions to the hump, which can be solved by changing parameters in the fit of the \code{\link{codhump}}.
#'
#' @return
#'
#' \describe{
#'  \item{loss}{vector containing the losses in life expectancy due to the hump, first all causes taken together, then for each cause of death.}
#'  \item{mode}{vector containing the modal age at death from the hump, first all causes taken together, then for each cause of death}
#'  \item{mean}{vector containing the mean age at death from the hump, first all causes taken together, then for each cause of death}
#'  \item{sd}{vector containing the standard deviation of the age at death from the hump, first all causes taken together, then for each cause of death}
#'  }
#'
#'
#' @export
#' @method summary codhump


summary.codhump <- function(object, ...){

  fit <- object

  # loss of life expectancy : all-cause
  mhat <- fit$all$gamma1 + fit$all$gamma2
  mhyp <- mhat - fit$all$gamma1
  widths <- diff(fit$mxc$x) ; widths <- c(widths[1],widths)
  if(max(widths) == 1){ages <- fit$mxc$x}else{ages <- floor(fit$mxc$x/widths)*widths}
  e0hat <- LT(Mx = mhat, ages = ages, mxsmooth = FALSE, verbose = FALSE)$ex[1]
  e0hyp <- LT(Mx = mhyp, ages = ages, mxsmooth = FALSE, verbose = FALSE)$ex[1]
  e0hump  <- round(e0hyp - e0hat,2)

  # loss of life expectancy : cause-specific
  losses <- rep(NA, length(fit$par$typ))
  for(i in 1:length(fit$par$typ)){
    mhump <- fit$all$gamma1 - fit$fits[[i]]$gamma1
    mhyp <- mhat - mhump
    e0hyp <- LT(Mx = mhyp, ages = ages, mxsmooth = FALSE, verbose = FALSE)$ex[1]
    losses[i] <- round(e0hyp - e0hat, 2)
  }
  names(losses) <- names(fit$par$typ)

  # descriptive statistics
  x <- fit$mxc$x
  mh <- as.numeric(fit$all$gamma1)
  mode <- x[which.max(mh)]
  mean <- mh %*% x / sum(mh)
  sd <- sqrt(mh %*% (x - mean)^2 / sum(mh))
  modes <- means <- sds <- rep(NA, length(fit$par$typ))
  for(i in 1:length(modes)){
    modes[i] <- x[which.max(fit$decomp[,i])]
    means[i] <- fit$decomp[,i] %*% x / sum(fit$decomp[,i])
    suppressWarnings(sds[i]   <- sqrt(fit$decomp[,i] %*% (x - means[i])^2 / sum(fit$decomp[,i])))
  }
  names(sds) <- names(means) <- names(modes) <- names(fit$par$typ)

  out <- list(loss = c(e0hump,losses), mode = c(mode,modes), mean = c(mean,means), sd = c(sd,sds))

  class(out) <- "summary.codhump"
  attr(out, "extra") <- fit

  invisible(out)

}

#' @rdname summary.codhump
#' @method print summary.codhump
#' @export

print.summary.codhump <- function(x, ...){

fit <- attr(x, "extra")
data <- fit$par$data
share.abs <- function(c){round(sum(data$dxc[data$x >= 10 & data$x < 35,c]) / sum(data$dxc[data$x >= 10 & data$x < 35,]) * 100, 1)}

cat("[[>>]] Information about the model used to fit the cause- and age-specific death rates")
cat("\nNumber of causes contributing to the hump:",length(fit$par$typ))
cat("\nList of causes contributing to the hump:",names(fit$par$typ))
cat("\nShare of observed deaths between 10 and 34 (%):")
cat("\n",paste0(names(fit$par$typ)," = ",unlist(lapply(fit$par$typ, share.abs)),collapse = ", "))
cat("\nAge range used for the analysis:",min(fit$par$x.range),"-",max(fit$par$x.range))
cat("\nNumber of iterations:",length(fit$dif),"( max =",fit$par$maxit,")")

if(length(fit$dif) == fit$par$maxit){lab <- "Reached maximum number of iterations"}
if(diff(rev(fit$dif))[1] < 1e-3){lab <- "Reached a stable solution"}
if(fit$neg[length(fit$neg)] < fit$neg[length(fit$neg)-1]){lab <- "Could not approach the constraints any closer"}
cat("\nReason for stopping the optimization:",lab)
if(any(x$loss[-1] / x$loss[1] < 5e-2) | max(fit$neg) < -0.01){
  cat("\nConsider changing parameters to improve fit (see documentation from codhump)")
  cat(paste("\n>> Small contribution of :",names(fit$par$typ)[x$loss[-1] / x$loss[1] < 5e-2]))
  cat(paste("\n>> Negative contributions amount to",-round(max(fit$neg*100),1),"% of total hump"))
}


cat("\n\n")

cat("[[>>]] Descriptive statistics about the hump")

cat("\nIntensity")
cat("\nLoss of life expectancy due to the hump: ", x$loss[1], "year")
cat("\nLoss of life expectancy by cause (years):",paste0(paste0(names(fit$par$typ)," = "),x$loss[-1], collapse = ", "))
cat("\nLoss of life expectancy by cause (%):",paste0(paste0(names(fit$par$typ)," = "), round(x$loss[-1] / x$loss[1] * 100, 1), collapse = ", "))

cat("\nInaccuracy compared to overall hump:",round((sum(x$loss[-1]) - x$loss[1]) / x$loss[1] * 100, 1),"%")
cat("\n")

cat("\nCentrality")
cat("\nMode of overall hump:", round(x$mode[1],2))
cat("\nMode by cause:",paste0(paste0(names(fit$par$typ)," = "),x$mode[-1], collapse = ", "))
cat("\nMean of overall hump:", round(x$mean[1],2))
cat("\nMean by cause:",paste0(paste0(names(fit$par$typ)," = "),round(x$mean[-1],2), collapse = ", "))
cat("\n")

cat("\nDispersion")
cat("\nStandard deviation of overall hump:", round(x$sd[1],2))
cat("\nStandard deviation by cause:",paste0(paste0(names(fit$par$typ)," = "),round(x$sd[-1],2), collapse = ", "))



invisible(x)

}
