#' @title Fit a Sum of Smooth Exponentials
#'
#' @description This function estimates a three-component non-parametric model of the age-specific mortality rates.
#' Its design resembles the parametric Heligman-Pollard model but it uses \emph{P}-splines instead of parametric functions to describe the force of mortality.
#'
#'
#' @param data data frame produced with \code{HMD2MH} or similarly structured
#' @param lambda.hump smoothing parameter for the hump component
#' @param lambda.sen smoothing parameter for the senescence part
#' @param x.hump upper age to compute the starting values for the hump component
#' @param x.sen lower age to compute the starting values for the senescence component
#' @param kappa ridge penalty
#' @param maxit maximum number of iterations of the optimization algorithm
#' @param verbose should the function print information during the optimization
#'
#'
#' @details Three non-parametric components are estimated for each phase of the force of mortality.
#' A first component estimates the decreasing trend observed during the first years of life.
#' A second component estimates the concave pattern observed during adolescence and early adulthood.
#' A third component estimates the approximately exponential increase observed during late adulthood.
#'
#' The function first estimates separate \emph{P}-splines for each components. Then, it uses a re-weighted iterative least square algorithm to adjust each component so that they sum up to the observed death rates.
#' The algorithm is relatively sensitive to starting values (i.e. the independent fit of each component). In case of a failed convergence, it is thus advised to reconsider the following parameters.
#'
#' The \code{x.hump} and \code{x.sen} arguments are used to compute the starting values of each component. More specifically, they are used to determined the age range on which the young adult mortality hump and the senescence components are fitted.
#' The values of \code{x.hump} and \code{x.sen} define in a way the interval where the hump and the senescence components overlap.
#' In case of a very narrow hump, it is advised to reduce the value of \code{x.hump} and/or of \code{x.sen}, wherease in case of an especially wide hump it may be useful to consider larger values for \code{x.hump} and \code{x.sen}.
#' Inadequate values of \code{x.hump} may result in incoherent starting values for the cause-deleted SSE models and a lack of convergence.
#'
#' The \code{lambda.hump} and \code{lambda.sen} parameters control the amount of smoothing imposed on the hump and senescence component respectively. Typically, the former is smaller than the latter, since this is the more instable part of the force of mortality.
#' However, in some cases, especially when using abridged datasets, it may be useful to consider smaller values of \code{lambda.sen}.
#'
#' The \code{maxit} argument defines the maximum number of iterations used for the adjustment step. This step uses a Penalized Composite Link Model (PCLM) along with an iterative re-weighted least squares (IRWLS) optimization procedure.
#' The \code{maxit} argument will therefore determine the maximum number of iterations used in the IRWLS loop, but also the speed of convergence since the re-weighting defines updated solution as \deqn{new = old * (1 - it/maxit) + new * (it/maxit)}
#' The \code{maxit} argument is defined at 200 by default, but it can be increased in case of an absence of convergence, even if the algorithm stopped before reaching \code{maxit} number of iterations.
#'
#' The \code{kappa} parameter regulates the force of the ridge penalty, which enforces the shape constraints, i.e. that the hump component is concave and the senescence component is monotonically increasing.
#'
#' @return
#' \describe{
#'    \item{data}{The original data frame produced with \code{HMD2MH} or similarly structured. Used later for the plot and summary methods.}
#'    \item{type}{The type of model used: non-parametric. To use in the summary method.}
#'    \item{method}{The algorithm used to fit the model: IRWLS. To use in the summary method.}
#'    \item{it}{The number of iterations of the IRWLS algorithm used to fit the model.}
#'    \item{maxit}{The maximum number of iterations of the optimization algorithm set before the estimation.}
#'    \item{deg2}{The degrees of freedom used for the \emph{P}-spline of the senescence component.}
#'    \item{deg3}{The degrees of freedom used for the \emph{P}-spline of the hump component.}
#'    \item{coef}{The concanated fitted coefficients of the \emph{P}-splines.}
#'    \item{XX}{List containing the base spline functions for each of the three components.}
#'    \item{mhat}{List containing three vectors containing each the age-specific mortality contributions of the estimated components. The sum of \code{mhat1}, \code{mhat2} and \code{mhat3} gives the overall fitted force of mortality.}
#' }
#'
#'
#' @examples
#'
#' data("CHE2010m")
#'
#' fit <- sse.fit(data = CHE2010m)
#'
#' @references
#'
#' Camarda, C. G., Eilers, P. H. C., & Gampe, J. (2016). Sums of smooth exponentials to decompose complex series of counts. Statistical Modelling.
#'
#' @seealso
#' \link{morthump}
#'
#' @export
#'
#' @import MortalitySmooth


sse.fit <- function(data, lambda.hump = 1, lambda.sen = 10, x.hump = 35, x.sen = 50, kappa = 10^5, maxit = 200, verbose = FALSE){

    x <- data$x
    n <- length(x)
    x <- 1:n
    y <- data$d
    e <- data$n
    lmx <- log(y/e)

    ## FITTING SECTION

    staINF <- 1
    endINF <- n
    staACC <- 1
    endACC <- n
    staAGI <- 1
    endAGI <- n

    ## let each component work in a particular age-range
    ## infant mortality
    p1 <- c(staINF, endINF)
    x1 <- x[x>=p1[1] & x<=p1[2]]
    w1 <- as.numeric(x%in%x1)
    y1 <- y[which(w1==1)]
    e1 <- e[which(w1==1)]
    ## aging mortality
    p2 <- c(staAGI, endAGI)
    x2 <- x[x>=p2[1] & x<=p2[2]]
    w2 <- as.numeric(x%in%x2)
    y2 <- y[which(w2==1)]
    e2 <- e[which(w2==1)]
    ## accident mortality
    p3 <- c(staACC, endACC)
    x3 <- x[x>=p3[1] & x<=p3[2]]
    w3 <- as.numeric(x%in%x3)
    y3 <- y[which(w3==1)]
    e3 <- e[which(w3==1)]

    ## Design matrix for the infant part
    B1 <- cbind(1,1/x1)
    nb1 <- ncol(B1)
    ## difference matrices
    D1 <- diff(diag(nb1), diff=2)
    tDD1 <- t(D1) %*% D1

    ## B-splines for the aging part
    deg2 <- floor(length(x)/5)
    xl2 <- min(x)
    xr2 <- max(x)
    xmax2 <- xr2 + 0.01 * (xr2 - xl2)
    xmin2 <- xl2 - 0.01 * (xr2 - xl2)
    B2 <- MortSmooth_bbase(x, xmin2, xmax2, deg2, 3)
    nb2 <- ncol(B2)
    ## difference matrices
    D2 <- diff(diag(nb2), diff=2)
    tDD2 <- t(D2) %*% D2

    ## B-splines for the accident-hump
    deg3 <- floor(length(x3)/5)
    xl3 <- min(x3)
    xr3 <- max(x3)
    xmax3 <- xr3 + 0.01 * (xr3 - xl3)
    xmin3 <- xl3 - 0.01 * (xr3 - xl3)
    B3 <- MortSmooth_bbase(x3, xmin3, xmax3, deg3, 3)
    nb3 <- ncol(B3)
    ## difference matrices
    D3 <- diff(diag(nb3), diff=3)
    tDD3 <- t(D3) %*% D3

    ## complete model matrix as a list
    ## in their order and with dimensions
    XX <- list(X1=B1, X2=B2, X3=B3)
    nx <- length(XX)
    nc <- unlist(lapply(XX, ncol))
    ind2 <- cumsum(nc)
    ind1 <- c(1, ind2[1:(nx-1)]+1)
    ind <- NULL # loop to avoid !
    for(i in 1:nx){
      ind[[i]] <- ind1[i]:ind2[i]
    }
    ## indicators for the fitted values
    indF <- cbind(w1, w2, w3)

    ## fitting a Poisson-GLM for infant mortality
    ww1 <- rep(0,length(y1))
    ww1[x1<=10] <- 1
    y1a <- c(y)
    y1a[which(y1a<=0)] <- 0
    y1 <- y1a[which(w1==1)]
    fit1 <- glm(y1~I(1/x1), offset=log(e1),
                weights=ww1,
                family=poisson)
    y1.st0 <- exp(XX[[1]] %*% fit1$coef)*e1
    lhaz1.st0 <- log(y1.st0/e1)

    y1.st <- rep(0,n)
    y1.st[which(w1==1)] <- y1.st0
    lhaz1.st <- log(y1.st/e)


    ## fitting P-splines over all x for the aging component
    y2a <- c(y-y1.st)
    y2a[which(y2a<=0)] <- 0
    y2 <- y2a[which(w2==1)]
    ww2 <- rep(0,length(y2))
    ww2[x2>=x.sen] <- 1
    fit2 <- suppressWarnings(Mort1Dsmooth(x=x2, y=y2, offset=log(e2),
                         w=ww2,
                         method=3, lambda=1000,
                         ndx = deg2, deg = 3, pord = 2))
    y2.st0 <- exp(XX[[2]] %*% fit2$coef)*e2
    lhaz2.st0 <- log(y2.st0/e2)
    y2.st <- rep(0,n)
    y2.st[which(w2==1)] <- y2.st0
    lhaz2.st <- log(y2.st/e)

    ## fitting P-splines over all x for the accident component
    # if y1+y2 > y => y3 < 0 => no fit / solution = only if possible
    if(all(y1.st + y2.st > y)){
     y3a <- c(y-y1.st-y2.st)
     y3a[which(y3a<=0)] <- 0
     y3 <- y3a[which(w3==1)]}
    ww3 <- rep(0, length(x3))
    ww3[x3>15 & x3<x.hump] <- 1
    fit3 <- suppressWarnings(Mort1Dsmooth(x=x3, y=y3, offset=log(e3),
                         w=ww3,
                         method=3, lambda=100000,
                         ndx = deg3, deg = 3, pord = 3))
    y3.st0 <- exp(XX[[3]] %*% fit3$coef)*e3
    lhaz3.st0 <- log(y3.st0/e3)
    y3.st <- rep(0,n)
    y3.st[which(w3==1)] <- y3.st0
    lhaz3.st <- log(y3.st/e)

    ## starting values for the whole fuction
    y.st <- y1.st + y2.st + y3.st
    lhaz.st <- log(y.st/e)

    ## concatenating the starting values
    coef.st <- as.vector(c(fit1$coef,
                           fit2$coef,
                           fit3$coef))
    coef <- coef.st

    ## penalty term just for the second component
    P2 <- lambda.sen * tDD2
    P3 <- lambda.hump * tDD3
    P <- bdiag(list(diag(0,nc[1]), P2, P3))

    ## including monotonicy for infant part
    D1mon <- diff(diag(nb1), diff=1)
    w1mon <- rep(0, nrow(D1mon))
    W1mon <- diag(w1mon)
    ## including monotonicy for aging part
    D2mon <- diff(diag(nb2), diff=1)
    w2mon <- rep(0, nrow(D2mon))
    W2mon <- diag(w2mon)
    ## including concaveness for accident-hump
    D3con <- diff(diag(nb3), diff=2)
    w3con <- rep(0, nrow(D3con))
    W3con <- diag(w3con)

    #plot(coef)

    ## iteration
    for(it in 1:maxit){

      d <- 1 - it/maxit

      ## penalty for the shape constraints
      P2mon <- kappa * t(D2mon) %*% W2mon %*% D2mon
      P3con <- kappa * t(D3con) %*% W3con %*% D3con
      Psha <- bdiag(list(diag(rep(0,nc[1])),
                         P2mon,
                         P3con))

      eta <- numeric(nx*n)
      for(i in 1:nx){
        eta0 <- rep(0, n)
        eta0[which(indF[,i]==1)] <- XX[[i]] %*% coef[ind[[i]]]
        eta[1:n+(i-1)*n] <- eta0
      }
      gamma <- exp(eta)*c(indF)

      mu <- numeric(n)
      for(i in 1:nx){
        mu <- (e * gamma[1:n+(i-1)*n]) + mu
      }

      w <- mu

      U <- matrix(NA, n, sum(nc))
      for(i in 1:nx){
        u <- gamma[1:n+(i-1)*n]/mu * e
        XXi <- matrix(0, nrow=n, ncol=nc[i])
        XXi[which(indF[,i]==1), ] <- XX[[i]]
        U0 <- u * XXi
        U[,ind[[i]]] <- U0
      }

      tUWU <- t(U) %*% (w * U)
      tUWUpP <- tUWU + P + Psha
      r <- y - mu
      tUr <- t(U) %*% r

      coef.old <- coef
      coef <- solve(tUWUpP, tUr + tUWU %*% coef)
      coef <- d * coef.old + (1 - d) * coef

      ## update weights for shape constraints
      ## aging
      W2mon.old <- W2mon
      coef2 <- coef[ind[[2]]]
      W2mon <- diag(diff(coef2) <= 0)
      ## accident-hump
      W3con.old <- W3con
      coef3 <- coef[ind[[3]]]
      W3con <- diag(diff(coef3, diff=2) >= 0)
      ## convergence criterion for the concaveness
      conv.con <- all(W2mon.old==W2mon,
                      W3con.old==W3con)

      dif.coef <- max(abs(coef.old - coef))/max(abs(coef))

      if(dif.coef < 1e-04 & it > 4 & conv.con) break
      if(verbose == T){cat(it, dif.coef, conv.con, "\n")}
    }


    coef.hat <- coef
    ## fitted values
    etas <- NULL
    for(i in 1:nx){
      etas[[i]] <- XX[[i]] %*% coef.hat[ind[[i]]]
    }
    eta1.hat <- etas[[1]]
    eta2.hat <- etas[[2]]
    eta3.hat <- etas[[3]]

    eta.hat <- numeric(nx*n)
    for(i in 1:nx){
      eta0 <- rep(0, n)
      eta0[which(indF[,i]==1)] <- XX[[i]]%*%coef.hat[ind[[i]]]
      eta.hat[1:n+(i-1)*n] <- eta0
    }

    gamma1.hat <- exp(eta1.hat)
    gamma2.hat <- exp(eta2.hat)
    gamma3.hat <- exp(eta3.hat)
    gamma.hat <- exp(eta.hat)*c(indF)

    mu1.hat <- exp(eta1.hat)*e1
    mu2.hat <- exp(eta2.hat)*e
    mu3.hat <- exp(eta3.hat)*e3
    mu.hat <- numeric(n)
    for(i in 1:nx){
      mu.hat <- (e * gamma.hat[1:n+(i-1)*n]) + mu.hat
    }

    fit <- list(data = data, type = "non-parametric", method = "IRWLS", it = it, maxit = maxit, deg2 = deg2, deg3 = deg3, coef = coef.hat, XX = XX, mhat = list(mhat1 = gamma1.hat, mhat2 = gamma2.hat, mhat3 = gamma3.hat))

    return(fit)

    }
