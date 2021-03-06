#' Semi-automatized grouping of causes of deaths
#'
#' @description Using the age- and cause-specific mortality rates, this function applies hierarchical clustering to suggest possible cause-of-death groupings that isolate the causes that are suceptible to contribute to the young adult mortality hump.
#'
#' @param data list produced with \code{HCD2MH} or similarly structured
#' @param x.range age range to consider for the analysis
#' @param k either a fixed number of cluster or the name of the criterion to use for selection
#'
#' @details
#' This function is designed to help selecting the causes of death that contribute to the young adult mortality hump. It procedes in steps.
#'
#' \enumerate{
#'   \item Compute the first derivative of the force of mortality (\code{rxc}) for each cause, in order to focuse on the amount of deviation instead of the absolute death rate.
#'
#'   \item Using the provided age-range, compute the euclidian distance between each couple of causes.
#'
#'   \item Based on this distance, run a hierchical clustering method (\code{"complete"} algorithm of the \link{hclust} function).
#'
#'   \item If \code{k} is numerical, it is taken as the chosen number of clusters and each cause of death is assigned to one of the k groups.
#'         Alternatively, \code{k} can indicate one of the selection criteria available in the \code{WeightedCluster} package.
#'         Among the most interesting options, is the Average Silhouette Width (\code{ASW}) that compares the average distance of an observation from the other members of its group and its average distance from the closest group. The \code{ASW} is computed for each number of groups \code{k}, and the one maximising the \code{ASW} is selected for its ability to maximise the homogeneity within the groups and the heterogeneity between the groups.
#'
#'   }
#'
#'
#' @return codgroup returns a list of six elements containing
#'
#' \describe{
#'  \item{cluster}{An object of class \link{hclust} on which additional analysis can be performed.}
#'  \item{groups}{A membership vector indicating a group number for each of the causes of death.}
#'  \item{k}{The number of groups chosen.}
#'  \item{typ}{A list of groups of causes as needed by the function \link{codhump}.}
#'  \item{data}{The original data stored in a list produced with \code{HCD2MH} or similarly structured}
#'  \item{x.range}{A vector indicating the age range to consider for the analysis}
#' }
#'
#' @examples
#'
#' data(USA2000m)
#'
#' grouping <- codgroup(USA2000m, k = "ASW", x.range = 10:35)
#'
#' @seealso
#'
#' \code{\link{codhump}}, \code{\link{HCD2MH}}, \code{\link[WeightedCluster]{as.clustrange}}
#'
#' @export
#'
#' @import WeightedCluster

codgroup <- function(data, x.range = 10:35, k = "ASW"){

  x <- data$x
  mxc <- as.data.frame(data$mxc)
  lab <- as.data.frame(data$lab)
  names(mxc) <- data$lab$short
  names(mxc)[is.na(names(mxc))] <- LETTERS[1:sum(is.na(names(mxc)))]

  rxc <- as.data.frame(apply(mxc,2,diff))
  names(rxc) <- names(mxc)
  d <- dist(x = t(rxc[which(x > min(x.range) & x < max(x.range)),]),method = "euclidian")
  cluster <- hclust(d)

  opt <- as.clustrange(object = cluster, diss = d, ncluster = ncol(rxc))
  if(!is.numeric(k)){
    k <- summary(opt)[which(rownames(summary(opt)) == k),1]
  }

  pca <- prcomp(t(rxc[which(x > min(x.range) & x < max(x.range)),]))

  gr <- cutree(cluster, k = k)

  typ <- list()
  for(i in 1:k){typ[[i]] <- which(gr == i)}
  names(typ) <- LETTERS[1:k]
  typ <- typ[-which.max(lapply(typ,length))]

  out <- list(cluster = cluster, groups = gr, k = k, typ = typ, data = data, x.range = x.range)

  class(out) <- "codgroup"

  return(out)

}
