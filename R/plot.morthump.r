#' @title Plot diagnostic for a morthump object
#'
#' @description Two plots (selectable by \code{which}) are currently available. They represent the age-specific death rates, and the density of the hump.
#'
#' @param x morthump object resulting from a call to the morthump function
#' @param which type of plot to be produced (see details)
#' @param ... other parameters to be passed through to plotting functions
#'
#' @details
#'
#' \code{Plot 1} displays the observed age-specific death rates, along with the fitted force of mortality and a estimation of what it look like if there was no young adult mortality hump.
#'
#' \code{Plot 2} displays the age-specific contributions to the hump. These contributions are rescaled in order for the hump to sum up to 1.
#' Consequently, the surface can be interpreted as a probability density function (pdf), from which summary indicators can be computed such as the mean, mode and median age at death due to the hump.
#' Other summary indicators can be obtained by calling \link{summary.morthump}.
#'
#' @examples
#'
#' data("CHE2010m")
#' fit <- morthump(data = CHE2010m, model = "sse")
#' plot(fit, which = 1)
#'
#' @seealso
#'
#' \link{morthump}, \link{summary.morthump}
#'
#' @export
#' @method plot morthump


plot.morthump <- function(x, which,...){

  fit <- x

  model <- attr(fit, "model")
  data <- fit$data
  if(data$x[1] == 0 & model != "sse"){data$x[1] <- 1e-5}

  hump <- summary(fit, print = FALSE)

  if(which == 1){

  plot(x = data$x, y = log(data$m), type = "p", ylab = "mx (log scale)", xlab = "age", pch = 16, col = "grey", axes = FALSE, ...)
  axis(side = 1)
  yl <- 10^(floor(log10(min(data$m))):0)
  axis(side = 2, at = log(yl), labels = yl, las = 1)
  box()
  title("Age-specific death rates")

  lines(x = data$x, y = log(hump$mhat), lwd = 2, col = 1)

  if(model != "hps"){

  lines(x = data$x, y = log(hump$mhump), lwd = 1, lty = 2)

  polygon(x = c(data$x, rev(data$x)), y = log(c(hump$mhat - hump$mhump, rev(hump$mhat))), col = "#00000020")

  # text(labels = paste(round(hump$loss,2),"yr"), x = which.max(hump$mhump),y = log(predict(fit, x = which.max(hump$mhump)) / 2))

  legend(x = "bottomright", legend = c("observed","fitted","hump"), col = c(8,1,1), pch = c(16,NA,NA), lwd = c(NA,2,1), lty = c(NA,1,2))

  }else{legend(x = "bottomright", legend = c("observed","fitted"), col = c(8,1), pch = c(16,NA), lwd = c(NA,2), lty = c(NA,1))}

  }

  if(which == 2){

  if(model == "hps"){stop("The HPS model does not assume any hump")}

  mode <- hump$mode
  median <- hump$median
  mean <- hump$mean
  pdf <- hump$pdf
  cdf <- hump$cdf
  qtl <- hump$qtl
  cols <- brewer.pal(n = 3, name = "Set1")

  plot(x = data$x, y = pdf(data$x), type = "n", las = 1, xlab = "age", ylab = "pdf", xlim = range(which(pdf(data$x) > 1e-5)), ...)

  xx <- seq(min(data$x), max(data$x), 0.1)
  polygon(x = c(xx,max(xx),min(xx)),  y = c(pdf(xx),0,0), col = "grey", border = NA)

  # deal with overlapping labels
  xydist <- function(x1, x2){sqrt((x1 - x2)^2 + (pdf(x1) - pdf(x2))^2)}
  d1 <- xydist(mode, median)
  d2 <- xydist(median, mean)

  points(x = mode, y = pdf(mode), pch = 3, lwd = 2, col = cols[1])
  text(x = mode,  y = pdf(mode), labels = "mode", pos = ifelse(d1 < 1 & mode < median,2,4), col = cols[1])
  segments(x0 = mode, x1 = mode, y0 = 0, y1 = pdf(mode), lty = 3, col = cols[1])

  points(median, pdf(median), pch = 3, lwd = 2, col = cols[2])
  text(x = median,  y = ifelse(d2 < 1 & pdf(median) > pdf(mean), pdf(median)*1.02, pdf(median)), labels = "median", pos = ifelse(median > mode,4,2), col = cols[2])
  segments(x0 = median, x1 = median, y0 = 0, y1 = pdf(median), lty = 3, col = cols[2])

  points(mean, pdf(mean), pch = 3, lwd = 2, col = cols[3])
  text(x = mean,  y = ifelse(d2 < 1 & pdf(median) < pdf(mean), pdf(mean)*1.02, pdf(mean)), labels = "mean", pos = ifelse(mean > mode,4,2), col = cols[3])
  segments(x0 = mean, x1 = mean, y0 = 0, y1 = pdf(mean), lty = 3, col = cols[3])

  # legend("topright", legend = c(expression(paste(e[0]^-H - e[0],"=",round(hump$loss,2),"yr")),expression(paste(d^H,"=", hump$dh))), box.lty = 0)

  title("Density of the hump")

  }

}
