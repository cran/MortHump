#' @title Plot diagnostic for a codhump object
#'
#' @description Five plots (selectable by \code{which}) are currently available. Some are designed to estimate the pertinence of the choice of the parameters (plots 1 and 2),
#' others can serve as diagnostic to judge the quality of the fit (plot 3 and 4), and plot 5 illustrates the fitted cause- and age-specific contributions to the hump.
#'
#' @param x codhump object resulting from a call to the codhump function
#' @param which type of plot to be produced (see details)
#' @param horiz should horizontal display be favored over vertical in panel plots
#' @param ... other parameters to be passed through to plotting functions
#'
#' @details
#'
#' \code{Plot 1} displays the age-specific death rates for all-cause mortality (black), as well as the remaining group of causes that were not selected among the ones susceptible to contribute to the hump (grey).
#' If the remaining causes still display a hump, this is a sign that the typology chosen for the fit (\code{typ} argument) omits some other causes that contribute to the hump.
#' The colored lines represent the cause-deleted forces of mortality, i.e. the resulting age-specific death rates after the deletion of a given cause of death. The closer these lines are to the all-cause mortality, the smaller are their respective death rates.
#' If some of these lines follow very closely the all-cause death rates, it is worth considering eliminating them from the list of causes that are assumed to contribute to the hump, as their contribution is likely to be non-significant.
#' Keeping them probably adds more complexity to the model than it brings information.
#'
#' \code{Plot 2} displays the independent fits of SSE models on each cause-deleted forces of mortality. These serve in turn as starting values for the constrained simultaneous fit (step 4 of the codhump algorithm).
#' This plot can be used as an indicator of poor choices in certain parameters, especially \code{x.hump} and \code{lambda}. For instance, if the value of \code{x.hump} is too large, the hump component may interfer with the senesence component by "competing" with it in the middle ages.
#' Alternatively, if a given cause accounts for all of the hump, the respective cause-deleted hump component will be hard to fit with an SSE model, which can translate visualy by unreasonable components.
#'
#' \code{Plot 3} displays two diagnostic plots concerning the evolution of the optimization process.
#' The left panel indicates the evolution of the speed of convergence, which corresponds to the highest relative change in the coefficients of the splines of each cause and component.
#' It is expected that, after a certain number of iterations, the speed of convergence decreases as the current solutions approaches the optimal solution. Ultimately, the algorithm is stopped if this value falls below 1e-3.
#' The right panel indicates the amount of negative values within the cause- and age-specific contributions to the hump.
#' These negative values appear when the constraints are not respected (i.e. that the sum of each cause-specific contributions do not sum up to the overall hump, or that the cause-deleted forces of mortality are poorly approximated by the two components of the SSE model).
#' They are expected to disapear as the current solution approaches the optimal solution, but in some cases the solution may diverge from this consrtaint and the algorithm is stopped if the share of negative values starts increasing.
#'
#' \code{Plot 4} displays the constrained SSE models on each of the cause-deleted forces of mortality. It also suggests a comparison with the all-cause SSE model, and the resulting age-specific contributions to the hump.
#' On the plot, these contributions appear in the form of the areas between the all-cause and the cause-deleted humps, which, in a counterfactual reasoning, can be interpreted as the contributions of each cause to the hump.
#'
#' \code{Plot 5} displays the age-specific contributions to the hump by cause of death. These contributions are rescaled in order for the overall hump (black line) to sum up to 1.
#' Consequently, the overall curve can be interpreted as a probability density function (pdf).
#' This plot has two main functions: to control that the cause-specific contributions to the hump sum up to the overall hump, and to visualize the shape of these contributions. For summary indicators of these contributions, use \link{summary.codhump}.
#'
#' @examples
#'
#' data("USA2000m")
#' typ <- list()
#' typ$tac <- 93
#' typ$sui <- 100
#' typ$hom <- 101
#' typ$poi <- c(97,98)
#' typ$oac <- c(94:96,99,102)
#' fit <- codhump(data = USA2000m, typ = typ, start.correct = TRUE)
#' plot(fit, which = 5)
#'
#' @seealso
#'
#' \link{codhump}, \link{summary.codhump}
#'
#' @import RColorBrewer
#'
#' @export
#' @method plot codhump

plot.codhump <- function(x, which, horiz = FALSE, ...){

  fit <- x

  if(which %in% 1:5 == FALSE){warning("Choose one of the five plots available.")}

  if(which == 1){
    x <- fit$mxc$x
    mxcd <- fit$mxc[,-c(1:2,ncol(fit$mxc))]
    mxwo <- fit$mxc$mxwo
    mx <- fit$mxc$mx
    col <- fit$par$col

    matplot(x, log(mxcd), type = "l", lty = 1, ylim = c(min(log(mxwo)),max(log(mx))), xlab = "age", ylab = "mx (log scale)", las = 1, col = col, axes = FALSE,...)
    lines(x, log(mxwo), col = "grey", lwd = 2)
    lines(x, log(mx), lwd = 2)
    legend(title = "Cause-deleted", x = "topleft", legend = names(fit$par$typ), col = fit$par$col, lty = 1, cex = 0.7)
    legend(title = "Remaining", x = "bottomright", legend = c("All-cause","Other causes"), lwd = 2, col = c(1,8), cex = 0.7)
    axis(1)
    axis(2,las = 1, at = log(10^(-7:0)), labels = 10^(-7:0))
    box()
    title("Age-specific death rates")
  }


  if(which == 2){

    typ <- fit$par$typ
    fits <- fit$start
    x <- fit$mxc$x
    mxc <- fit$mxc$mx - fit$mxc[,-c(1:2,ncol(fit$mxc))]
    ncod <- length(typ)
    
    par(mfrow = disp(ncod, horiz = horiz))

    for(i in 1:length(fits)){

    plot(x, log(fit$mxc[,i+2]), col = 8, pch = 16, xlab = "age", ylab = "mx (log scale)", las = 1, axes = FALSE,...)
    lines(x, fits[[i]]$eta, col=2, lwd=2)
    lines(x, fits[[i]]$eta1, col=3)
    lines(x, fits[[i]]$eta2, col=4)
    #legend("top", legend = c("obs","total","hump","sen."), col=c(8,2:4),lty=c(-1,1,1,1), pch=c(16,-1,-1,-1),ncol=2)
    title(names(typ)[i])
    axis(1)
    axis(2,las = 1, at = log(10^(-7:0)), labels = 10^(-7:0))
    box()

      }
    par(mfrow = c(1,1))
  }

  if(which == 3){

    dif <- fit$dif
    neg <- fit$neg

    par(mfrow = c(1,2))
    # speed of convergence
    plot(1:sum(dif!=0), dif[1:sum(dif!=0)], ylim = c(0,0.1), las = 1, xlab = "iteration", ylab = "max relative change of eta",...)
    abline(h = 1e-3, col = 2)
    title("Speed of convergence")


    # elimination of negative contributions
    neg <- neg[neg != -Inf & neg != Inf]
    plot(1:sum(neg!=0), - neg[1:sum(neg!=0)] * 100, ylim = c(0,-min(neg*100)), las = 1, xlab = "iteration", ylab = "% of negative values in the contributions",...)
    abline(h = 0, col = 2)
    title("Elimination of negative values")

   par(mfrow = c(1,1))
  }

 if(which == 4){

    # comparison of all-cause and cause-deleted fits
    typ <- fit$par$typ 
    x <- fit$mxc$x
    col <- fit$par$col
    ncod <- length(typ)
    fits <- fit$fits
    mxcd <- fit$mxc[,3:ncol(fit$mxc)]

    par(mfrow = disp(ncod, horiz = horiz))
    for(i in 1:ncod){
      fit.i <- fits[[i]]

      plot(x, log(fit$mxc[,2]), pch = 16, col = 1, xlab = "age", ylab = "mx (log scale)", ylim = c(-10,-1), axes = FALSE,...)
      points(x, log(mxcd[,i]), pch = 16, col = col[i])
      lines(x, fit$all$eta)
      lines(x, fit$all$eta1, lty = 2)
      lines(x ,fit$all$eta2, lty = 3)
      lines(x, fit.i$eta, col = col[i])
      lines(x, fit.i$eta1, col = col[i], lty = 2)
      lines(x, fit.i$eta2, col = col[i], lty = 3)
      polygon(x = c(x,rev(x)), y = c(fit$all$eta1,rev(fit.i$eta1)), col = paste(col[i],60,sep=""), border = NA)
      axis(1)
      axis(2,las = 1, at = log(10^(-7:0)), labels = 10^(-7:0))
      box()
      title(names(typ)[i])
      }

    par(mfrow=c(1,1))
  }

  if(which == 5){
    
    typ <- fit$par$typ
    x <- fit$mxc$x
    col <- fit$par$col
    decomp <- fit$decomp
    decomp[decomp < 0] <- 0

    # cause- and age-specific contributions to the hump
    plot(x,fit$all$gamma1 / sum(fit$all$gamma1), xlim = c(10,60), type = "n", xlab = "age", ylab = "pdf", las = 1)
    polygon(c(x,rev(x)),c(rep(0,length(x)),rev(decomp[,1])) / sum(fit$all$gamma1),col = col[1], border = NA)
    if(ncol(decomp) > 1){
    for(i in 2:length(typ)){
      polygon(c(x,rev(x)),c(rowSums(as.matrix(decomp[,1:(i-1)])),rev(rowSums(as.matrix(decomp[,1:i])))) / sum(fit$all$gamma1),col = col[i], border = NA)
    }}
    lines(x,fit$all$gamma1 / sum(fit$all$gamma1),lwd = 2)

    #     axis(1,at = modeT, label = expression(M^T), lwd = 2, font = 2, line = 0.7)
    #     axis(1,at = modeA, label = expression(M^A), lwd = 3, col = col[1], font = 2, line = 0.7)
    #     axis(1,at=modeC,label=expression(M^C),lwd=3,col=col[3],font=2,line = 0.7)
    #     segments(x0=modeT,y0=0,x1=modeT,y1=max(fit$all$gamma1 / sum(fit$all$gamma1)),lty=2,lwd=2)
    #     segments(x0=modeA,y0=0,x1=modeA,y1=(fit$all$gamma1 / sum(fit$all$gamma1))[modeA-9],lty=2,col=col[1],lwd=2)
    #     segments(x0=modeC,y0=0,x1=modeC,y1=max(fit$decomp$C / sum(fit$all$gamma1)),lty=2,col=col[3],lwd=2)
    #     legend("topright", box.lty = 0, x.intersp = 0.3,
    #           legend = c("","area",expression(Delta^{kappa}),expression(bar(x)^kappa),expression(M^kappa),expression(sigma^kappa),
    #                      "A",paste(round(SA*100),"%"),paste(round(e10HA,2),"yr"),round(meanA,1),modeA,round(sdA,1),"C",paste(round(SC*100),"%"),paste(round(e10HC,2),"yr"),round(meanC,1),modeC,round(sdC,1)))

    legend("topright", fill = col, legend = names(typ))
    title("Contributions to the hump")

    #   negs <- pos <- decomp
    #   negs[negs > 0] <- 0
    #   pos[pos < 0] <- 0
    #   bp <- barplot(t(pos),space=0,col=col,border="grey",las=1,ylim=c(min(rowSums(negs)),max(rowSums(pos))))
    #   par(new = T)
    #   bp <- barplot(t(negs),space=0,col=col,border="grey",las=1,ylim=c(min(rowSums(negs)),max(rowSums(pos))),axes=F)
    #   #bp <- barplot(t(matrix(rep(mu1T,ncod),nrow=n) - do.call(cbind,muhat1)),space=0,col=col,border="grey",las=1)
    #   lines(bp-0.5, mu1T, lwd=2, type="s")
    #   axis(1,labels = x, at = x - min(x) + 0.5)
    #   title("Contributions to the hump: constrained")
    #   legend("topright",fill=col,legend=lab,border="grey")

  }

}
