
#' @title Decompose the contributions to the young adult mortality hump by cause
#'
#' @description
#' \code{codhump} is used to fit a constrained non-parametric estimation of the cause- and age-specific contributions to the young adult mortality hump.
#'
#' @param data list produced with \code{HCD2MH} or similarly structured
#' @param typ list containing the cause-of-death typology to use (see details)
#' @param x.range age range to consider for the analysis
#' @param x.hump assumed end of the hump used to estimate the starting values (see details)
#' @param maxit maximum number of iterations of the optimization algorithm
#' @param lambda smoothing parameter for the senesence component
#' @param start.correct logical value to automatically correct incoherent starting values
#'
#' @details
#'
#' The estimation uses simultaneous constrained splines to estimate a Sum of Smooth Exponentials (SSE) model on cause-deleted forces of mortality.
#' Briefly, the SSE model describes the observed force of mortality over age \emph{mu} as the sum of three vectors \emph{mu1, mu2, mu3} over age.
#' In other words, it assumes that deaths are realizations from Poisson distributions with mean composed of three parts: infant, early-adulthood and old-age mortality.
#' For more information on the SSE model, see Camarda et al. (2016) and \link{sse.fit}. For the purpose of the study of the young adult mortality hump, the SSE model is here reduced to two components capturing the hump and the senescence parts of the force of mortality.
#'
#' In order to decompose the hump by cause of death, this model uses a constrained approach on cause-deleted forces of mortality that can be summarized in four steps.
#'
#' \enumerate{
#'    \item{Identify manually those causes of death that contribute to the young adult hump component}
#'    \item{Estimate an SSE model on the overall mortality in order to separate the senescent and young adult hump components}
#'    \item{Construct cause-deleted datasets by removing separately deaths from each cause that was identified in step 1.}
#'    \item{Simultaneously estimate SSE models for each of these cause-deleted datasets, interpreting the diminution of each component as the contribution of this cause to that component,
#'          and constraining the sum of all these contributions to be equal to the components estimated in step 2.}
#' }
#'
#' The SSE model on which this algorithm is based is more adaptive to specific mortality schedules than parametric models such as the Heligman-Pollard. It is thus designed to converge to meaningful results in the majority of cases.
#' It however sometimes needs some fine tuning in order to reach coherent results. Please pay special attention to the following parameters in order to maximize the chance to get a meaningful result.
#'
#' The \code{typ} argument defines the typology of causes of death that are assumed to contribute to the young adult mortality hump.
#' Each element of the list is a vector containing one or several numerical values corresponding to the columns of the \code{mxc} data frame from the \code{data} object.
#' If an element of \code{typ} contains only one cause, this cause will be considered on its own. If several causes are mentioned in the same element of the \code{typ} object, these causes will be grouped and considered together in the analysis.
#' The names of the elements of the \code{typ} object are recycled as the names of the causes of death in the typology.
#' The choice of the causes of death included in the \code{typ} argument has a profound influence on the results, and should therefore be made with caution.
#' As each case calls for a specific selection of causes of death, there is no general rule as for which causes should be assumed to contribute to the hump. However, the model assumes that the list accounts for all of the hump.
#' A good way to test this assumption is to plot the force of mortality after removing all of the deaths from the causes included in the list and check that this \emph{leftover} category does not display a hump.
#' Failure to include sufficient causes will result in a probable underestimation of the hump, as well as a difficulty for the algorithm to converge since the all-cause hump will not be totally accounted for by the cause-specific contributions to the hump.
#' On the other hand, specifying too many causes as contributing to the hump may result in difficulties to estimate cause-specific contributions that are not based on sufficient empirical evidence.
#' Typically, the \code{typ} object will include somewhere between 2 and 6 causes (or groups of causes) of death depending on the context and the number of causes available in the dataset.
#'
#' The \code{x.range} argument defines the age range on which the analysis takes place. The bottom boundary should correspond approximately to the lowest observed age-specific death rate, which is usually situated around 10 years of age.
#' The top of the range should include most of the adult years, but should also avoid the so-called mortality plateau often observed after age 100. Therefore, the \code{x.range} argument is defined by default as \code{10:90}, but the user may want to change these values if necessary.
#'
#' The \code{maxit} argument defines the maximum number of iterations used for the step 4 of the algorithm. This step uses a Penalized Composite Link Model (PCLM) along with an iterative re-weighted least squares (IRWLS) optimization procedure.
#' The \code{maxit} argument will therefore determine the maximum number of iterations used in the IRWLS loop, but also the speed of convergence since the re-weighting defines updated solution as \deqn{new = old * (1 - it/maxit) + new * (it/maxit)}
#' The \code{maxit} argument is defined at 200 by default, but it can be increased in case of an absence of convergence, even if the algorithm stopped before reaching \code{maxit} number of iterations.
#'
#' The \code{x.hump} argument is used to compute the starting values of the independant SSE models. More specifically, it is used to determined the apriori age range to be considered for the young adult mortality hump.
#' In case of a very narrow hump, it is advised to reduce the value of \code{x.hump} to 25 or even 20 (but not lower), wherease in case of an especially wide hump it may be useful to consider larger values for \code{x.hump}.
#' Inadequate values of \code{x.hump} may result in incoherent starting values for the cause-deleted SSE models and a lack of convergence in the step 4 of the algorithm (see \code{start.correct}).
#'
#' The \code{lambda} parameter controls the amount of smoothing imposed on the senescence component of the SSE model. Typically, a large value is advisable since this is the part of the force of mortality on which the model is not focused.
#' However, in some cases, especially when using abridged datasets, it may be useful to consider smaller values of \code{lambda} such as 10 or 100. A bad choice of \code{lambda} may result in poor starting values for the SSE and a lack of convergence in the step 4 of the algorithm.
#'
#' The \code{start.correct} argument is conceived as a safeguard against misspecified starting values of the SSE model. Specifically, while estimating the SSE model on the cause-deleted forces of mortality, if the hump component peaks after the middle of the \code{x.range} interval, the starting values are replaced with the all-cause components.
#' This parameter is designed as an attempt to salvage a bad choice in other arguments, especially \code{typ}, \code{maxit}, \code{x.hump} and \code{lambda}, but remains an emergency safeguard. In case of a lack of convergence, it is advised to change the values of the other parameters instead on relying on the \code{start.correct} argument to guarantee convergence.
#'
#' @return Returns an object of class codhump that includes
#'
#'
#' \describe{
#'  \item{mxc}{Data frame containing for each age the overall mortality rates, the cause-deleted rates and the rates of causes that do not contribute to the hump.}
#'  \item{all}{Fit of the \code{sse} model on the all-cause mortality.}
#'  \item{start}{Fit of the \code{sse} model for each cause-deleted mortality rates before constraining.}
#'  \item{fits}{Fit of the \code{sse} model for each cause-deleted mortality rates after constraining.}
#'  \item{decomp}{Age- and cause-specific contributions to the hump.}
#'  \item{neg}{Percents of negative contributions after each iteration.}
#'  \item{dif}{Maximum relative change among all coefficients of the components of the \code{sse} model after each iteration.}
#'  \item{par}{List of parameters provided to fit the model.}
#'  }
#'
#' @examples
#'
#' data("USA2000m")
#' typ <- list()
#' typ$tac <- 89
#' typ$sui <- 96
#' typ$hom <- 97
#' typ$poi <- c(93,94)
#' typ$oac <- c(95,98,100)
#' fit <- codhump(data = USA2000m, typ = typ, start.correct = TRUE)
#'
#' @references
#'
#' Camarda, C. G., Eilers, P. H. C., & Gampe, J. (2016). Sums of smooth exponentials to decompose complex series of counts. Statistical Modelling.
#'
#' Remund, A., Riffe, T., & Camarda, C. (2016). A Non-Parametric Approach to Decompose the Young Adult Mortality Hump by Causes of Death. Poster presented at the 2016 annual meeting of the Population Association of America. March 2016, Washington DC.
#'
#' @seealso
#'
#' \link{sse.fit}, \link{summary.codhump}, \link{plot.codhump}, \link{HCD2MH}
#'
#' @export
#'

codhump <- function(data, typ, x.range = 10:90, maxit = 200, x.hump = 35, start.correct = FALSE, lambda = 10^4){

  if(sum(duplicated(unlist(typ))) > 0){warning("At least one cause appears twice in the typology.")}

  Mxc <- data$mxc
  exp <- data$nx
  x <- data$x

  ncod <- length(typ)
  suppressWarnings(col <- brewer.pal(n = ncod, "Set2"))
  cod <- lapply(lapply(typ,function(t){-t}),cod.data,Mxc,exp)
  ys <- lapply(cod,function(mat){mat$d})
  ms <- lapply(cod,function(mat){mat$m})
  m <- rowSums(Mxc)
  y <- m * exp
  mB <- cod.data(which(!(1:ncol(Mxc) %in% unlist(typ))),Mxc,exp)$m
  yB <- cod.data(which(!(1:ncol(Mxc) %in% unlist(typ))),Mxc,exp)$d

  x.keep <- which(x > min(x.range) & x < max(x.range))
  x <- x[x.keep]
  e <- exp[x.keep]
  m <- m[x.keep]
  y <- y[x.keep]
  mB <- mB[x.keep]
  yB <- yB[x.keep]
  ys <- lapply(ys,function(vec){vec[x.keep]})
  ms <- lapply(ms,function(vec){vec[x.keep]})
  n <- length(x.keep)
  mxcd <- do.call(cbind,lapply(cod,function(x){x$m}))
  mxcd <- mxcd[x.keep,]

  ## fit SSE on Total deaths
  SSEfitT <- sse2.fit(x = x, d = y , e = e,
                        lambda1 = 1, lambda2 = lambda,
                        mon = FALSE, ridge = 0, x1 = x.hump, x2 = 50,
                        plotIT = FALSE, plotFIT = FALSE)
  eta1T <- SSEfitT$eta1
  eta2T <- SSEfitT$eta2
  mu1T <- exp(eta1T)
  mu2T <- exp(eta2T)
  yhat1T <- exp(eta1T)*e
  yhat2T <- exp(eta2T)*e

  ## fit single component to deaths from causes assumed w/o hump
  fitB <- Mort1Dsmooth(x = x, y = yB, offset = log(e),
               method = 3, lambda = lambda,
               ndx = ncol(SSEfitT$B), deg = 3, pord = 2)
  etaB <- fitB$logmortality
  muB <- exp(etaB)
  yhatB <- exp(etaB)*e

  ## fitting SSE on each "Total-cause" deaths, constraining the sum of the
  ## differences with all-cause to be equal to all-cause

  ## common B-splines
  B <- SSEfitT$B
  nb <- ncol(B)
  ## number of components
  nx <- 2
  ## number of CODs
  ncod <- ncod
  ## labels for CODs
  lab <- names(typ)

  ## starting values fitting independently
  SSEfits <- lapply(ys, sse2.fit, x = x, e = e,
                    lambda1 = 1, lambda2 = lambda,
                    mon = FALSE, ridge = 1e-2, max.it = 100, x2 = 50, x1 = x.hump,
                    plotIT = FALSE, plotFIT = FALSE)

  ## starting coeff
  coef.st <- lapply(SSEfits,function(fit){fit$coef})
  #  in case the hump disppears after removing this cause and the starting values are unreasonible, use overall as starting
  if(start.correct == TRUE){for(i in 1:ncod){if(which.max(coef.st[[i]][1:nb]) > nb/2){coef.st[[i]] <- SSEfitT$coef}}}
  coef.st <- unlist(coef.st)

  ## difference matrices for accident-hump
  D1 <- diff(diag(nb), diff=3)
  tDD1 <- t(D1) %*% D1
  ## difference matrices for accident-hump
  D2 <- diff(diag(nb), diff=2)
  tDD2 <- t(D2) %*% D2

  ## final design matrix
  U <- kronecker(diag(nx*ncod), B)

  ## final number of coefficients
  np <- ncol(U)

  ## smoothing parameters
  lambda1 <- 1 ; if(length(typ) == 1){lambda1 <- 50}
  lambda2 <- lambda

  ## component specific penalty terms
  P1 <- P2 <- list()
  for(i in 1:ncod){
    P1[[i]] <- lambda1 * tDD1
    P2[[i]] <- lambda2 * tDD2
  }
  ## overall penalty term
  P <- matrix(0,np,np)
  for(i in 1:ncod){
    P[1:nb+2*(i-1)*nb,1:nb+2*(i-1)*nb] <- P1[[i]]
    P[1:nb+(2*(i-1)+1)*nb,1:nb+(2*(i-1)+1)*nb] <- P2[[i]]
  }

  ## final response
  ## deaths from Mx - Mxc, deaths from total 1st comp, deaths from total 2nd comp
  yT <- c(unlist(ys), yhat1T * (ncod - 1), yhat2T * (ncod - 1) + yhatB)

  ## shape penalties
  kappa <- 10^6
  D1 <- W1 <- D2 <- W2 <- list()
  for(i in 1:ncod){
    D1[[i]] <- diff(diag(nb), diff=2)
    W1[[i]] <- diag(rep(0, nrow(D1[[i]])))
    D2[[i]] <- diff(diag(nb), diff=1)
    W2[[i]] <- diag(rep(0, nrow(D2[[i]])))
  }

  ## composite matrixlent
  C0 <- kronecker(matrix(1,1,nx), diag(e))
  C1 <- kronecker(diag(ncod), C0)
  C2 <- kronecker(matrix(1,1,ncod), diag(rep(e,nx)))
  C <- rbind(C1,C2)

  wei <- rep(1, length(yT))
  wei[n*ncod+ 1:(nx*n)] <- 1e5

  ## starting linear predictor
  coef <- coef.st
  eta <- U %*% coef
  Pr <- 10^-2*diag(np)

  dif <- neg <- vector(length = maxit)
  for(it in 1:maxit){

    ds <- 1 - it/maxit

    ## penalty for the shape constraints
    P1 <- P2 <- list()
    for(i in 1:ncod){
      P1[[i]] <- kappa * t(D1[[i]]) %*%W1[[i]]%*% D1[[i]]
      P2[[i]] <- kappa * t(D2[[i]]) %*%W2[[i]]%*% D2[[i]]
    }


    Psha <- matrix(0,np,np)
    for(i in 1:ncod){
      Psha[1:(nb*(2*i-1)),1:(nb*(2*i-1))] <- P1[[i]]
      Psha[1:(nb*2*i),1:(nb*2*i)]         <- P2[[i]]
    }

    gamma   <- exp(eta)
    mu      <- C %*% gamma
    X       <- (C * ((1 / mu) %*% t(gamma)) ) %*% U
    w       <- as.vector(wei*mu)
    r       <- wei* (yT - mu + C %*% (gamma * eta))
    G       <- t(X) %*% (w * X)
    GpP     <- G + P + Psha + Pr
    tXr     <- t(X) %*% r
    beta    <- solve(GpP, tXr)
    eta.old <- eta
    eta     <- U%*%beta
    eta     <- ds*eta.old + (1-ds)*eta
    ## cod and comp specific coeff
    beta1 <- beta2 <- list()
    for(i in 1:ncod){
      beta1[[i]] <- beta[1:nb+2*(i-1)*nb]
      beta2[[i]] <- beta[1:nb+(2*(i-1)+1)*nb]
    }

    ## update weights for shape constraints for each component and each cod
    W1.old <- W1
    W2.old <- W2
    for(i in 1:ncod){
      W1[[i]] <- diag(diff(beta1[[i]], diff=2) >= 0)
      W2[[i]] <- diag(diff(beta2[[i]]) <= 0)
    }

    ## stop when the coefficients stop changing
    dif.eta <- max(abs((eta - eta.old)/eta.old))
    if(dif.eta < 1e-3 & it > 4) break
    dif[it] <- dif.eta

    ## stop in case negative values start increasing
    etahat1 <- list()
    for(i in 1:ncod){etahat1[[i]] <- eta[1:n+2*(i-1)*n]}
    muhat1 <- lapply(etahat1,exp)
    decomp <- matrix(rep(mu1T,ncod),nrow=n) - do.call(cbind,muhat1)
    neg[it] <- sum(decomp[decomp < 0]) / sum(decomp[decomp > 0])
    if(it > 5){if(neg[it] < neg[it-1]) break}
    }

  etahat1 <- etahat2 <- list()
  for(i in 1:ncod){
    etahat1[[i]] <- eta[1:n+2*(i-1)*n]
    etahat2[[i]] <- eta[1:n+(2*(i-1)+1)*n]
  }
  muhat1 <- lapply(etahat1,exp)
  muhat2 <- lapply(etahat2,exp)
  decomp <- matrix(rep(mu1T,ncod),nrow=n) - do.call(cbind,muhat1)
  decomp <- as.data.frame(decomp)
  names(decomp) <- lab

  finalfits <- list()
  for(i in 1:ncod){
    fitcod <- list()
    fitcod$beta <- beta[1:n+2*(i-1)*n]
    fitcod$eta1 <- etahat1[[i]]
    fitcod$eta2 <- etahat2[[i]]
    fitcod$eta <- log(exp(etahat1[[i]]) + exp(etahat2[[i]]))
    fitcod$gamma1 <- muhat1[[i]]
    fitcod$gamma2 <- muhat2[[i]]
    fitcod$mu <- muhat1[[i]] + muhat2[[i]]
    fitcod$it <- it
    fitcod$B <- B
    finalfits[[i]] <- fitcod
  }
  names(finalfits) <- lab

  fit <- list(mxc = data.frame(x = x, mx = m, mxcd = mxcd, mxwo = mB), all = SSEfitT, start = SSEfits, fits = finalfits, decomp = decomp, neg = neg[1:it], dif = dif[1:it], par = list(data = data, Mxc = Mxc, exp = exp, typ = typ, x.range = x.range, x.hump = x.hump, maxit = maxit, start.correct = start.correct, coef.st = coef.st, col = col))

  class(fit) <- "codhump"

  return(fit)

  }

