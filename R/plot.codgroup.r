#' Plot diagnostic for a codgroup object
#'
#' @description Four plots (selectable by \code{which}) are currently available. They provide details on the clustering process that led to the suggested grouping.
#'
#' @param x codgroup object resulting from a call to the codgroup function
#' @param which which type of plot to be produced (see details)
#' @param ... other parameters to be passed through to plotting functions
#' @details
#' Four plots are curently available.
#'
#' \describe{
#'   \item{which = 1}{Dendrogram from the hierarchical clustering algorithm}
#'
#'   \item{which = 2}{Plot a list of different criteria for each possible values of \code{k}.}
#'
#'   \item{which = 3}{Plot the clusters on the plan between the first two dimensions of a Principal Component Analysis.}
#'
#'   \item{which = 4}{Plot the first difference of the force of mortality pooled by cluster.}
#'
#'   }
#'
#'
#' @examples
#'
#' data(USA2000m)
#'
#' grouping <- codgroup(USA2000m, k = "HGSD", x.range = 10:34)
#'
#' plot(grouping, which = 3)
#'
#' @seealso
#'
#' \link{codgroup}
#'
#' @export
#'
#' @import WeightedCluster
#' @import RColorBrewer
#' @method plot codgroup

plot.codgroup <- function(x, which,...){

  fit <- x

  if(which %in% 1:4 == FALSE){warning("Choose one of the four plots available.")}

  cluster <- fit$cluster
  k <- fit$k
  gr <- fit$groups
  x.range <- fit$x.range
  rxc <- as.data.frame(apply(fit$data$mxc,2,diff))
  x <- fit$data$x
  lab <- as.data.frame(fit$data$lab)
  d <- dist(x = t(rxc[which(x > min(x.range) & x < max(x.range)),]),method = "euclidian")
  opt <- as.clustrange(object = cluster, diss = d, ncluster = ncol(rxc))

  pca <- prcomp(t(rxc[which(x > min(x.range) & x < max(x.range)),]))

  if(which == 1){

    plot(cluster, main = "dendrogram", xlab = "causes", ylab = "", cex = 0.7, sub = NA, las = 1, axes = FALSE,...)

    }

  if(which == 2){

    par(xpd = TRUE)
    plot(opt, main = "cluster quality", stat = c("ASW","PBC","HGSD"), las = 1, withlegend = FALSE,...)
    cols <- brewer.pal(n = 12, "Set3")
    text(x = nrow(opt$stats),y = opt$stats[nrow(opt$stats),1],col = cols[1], "PBC", font = 2, pos = 1)
    text(x = nrow(opt$stats),y = opt$stats[nrow(opt$stats),3],col = cols[4], "HGSD", font = 2, pos = 1)
    text(x = nrow(opt$stats),y = opt$stats[nrow(opt$stats),4],col = cols[5], "ASW", font = 2, pos = 1)
    par(xpd = FALSE)
    abline(v = k, col = 8)

  }

  if(which == 3){

    par(xpd = TRUE)
    suppressWarnings(cols <- brewer.pal(n = k, name = "Set1"))
    pvar <- summary(pca)$importance
    xlim <- range(pca$x[,1]) * 1.2
    ylim <- range(pca$x[,2]) * 1.2
    plot(pca$x, pch = 16, col = cols[gr], xlab = paste(round(pvar[2,1]*100,2),"%"), ylab = paste(round(pvar[2,2]*100,2),"%"), axes = F, xlim = xlim, ylim = ylim,...)
    title("PCA")
    text(x = pca$x[,1], y = pca$x[,2], labels = lab$short, col = cols[gr], cex = 0.5, pos = 4)
    par(xpd = FALSE)
    abline(h = 0, col = 8)
    abline(v = 0, col = 8)
    par(xpd = TRUE)
    legend(x = min(pca$x[,1]), y = min(ylim), title = "clusters", legend = c(1:k), col = cols, pch = 16, ncol = k, cex = 0.8, bg = NA)
    par(xpd = FALSE)
  }

  if(which == 4){

    suppressWarnings(cols <- brewer.pal(n = k, name = "Set1"))
    rxc.gr <- matrix(NA, nrow = nrow(rxc), ncol = k)
    for(i in 1:k){if(sum(gr == i) > 1){rxc.gr[,i] <- rowMeans(rxc[,gr == i])}else{rxc.gr[,i] <- rxc[,gr == i]}}
    ylim <- range(rxc.gr[which(x > min(x.range) & x < max(x.range)),], na.rm = TRUE)
    matplot(x = rowMeans(cbind(x[-1],x[-length(x)])), y = rxc.gr, type = "l", lty = 1, xlab = "age", ylab = "",
            main = paste("First difference of mx by cluster (k = ",k,")",sep=""), las = 1, ylim = ylim, col = cols, xlim = c(0,max(x.range)+10),...)
    abline(v = range(x.range), col = 8, lwd = 3)
    abline(h = 0, col = 8)
    arrows(x0 = min(x.range), y0 = ylim[2], x1 = max(x.range), y1 = ylim[2], length = 0.1, code = 3,col="grey")
    text(x = mean(x.range), y = ylim[2], labels = "age range for clustering", pos = 1, cex = 0.8, col = "grey")
    legend("topright", legend = 1:k, lty = 1, col = cols, title = "Clusters",cex = 0.7, ncol = 1, x.intersp = 0.5)

  }

}
