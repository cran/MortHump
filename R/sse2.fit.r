
# a function that fits a Sum of Smooth exponentials with only two components (senescence and hump)

# NB: 'data' should be truncated at the lowest point of the age-specific death rates to avoid the infant mortality component

sse2.fit <- function(x, d, e, lambda1, lambda2, kappa = 10^6, deg = NULL,
                     plotIT = TRUE, plotFIT = TRUE, lab = " ", mon = TRUE,
                     max.it=200, ridge=10^-4,
                     x1 = 35, x2 = 50){

  n <- length(x)
  y <- d
  lmx <- log(y/e)

  wei <- w1 <- w2 <- rep(1,n)
  wei[e==0] <- 0
  lhaz.act <- log(y/e)

  ## B-splines stuff common to both component
  if(is.null(deg)){
    deg <- floor(length(x)/5)
  }else{
    deg=deg
  }
  xl <- min(x)
  xr <- max(x)
  xmax <- xr + 0.01 * (xr - xl)
  xmin <- xl - 0.01 * (xr - xl)
  B <- MortSmooth_bbase(x, xmin, xmax, deg, 3)
  nb <- ncol(B)
  ## difference matrices for accident-hump
  D1 <- diff(diag(nb), diff=3)
  tDD1 <- t(D1) %*% D1
  ## difference matrices for accident-hump
  D2 <- diff(diag(nb), diff=2)
  tDD2 <- t(D2) %*% D2

  ## complete model matrix as a list
  ## in their order and with dimensions
  XX <- list(X1=B, X2=B)
  nx <- length(XX)
  nc <- unlist(lapply(XX, ncol))
  ## indicators for the coefficients of each basis
  ind2 <- cumsum(nc)
  ind1 <- c(1, ind2[1:(nx-1)]+1)
  ind <- NULL
  for(i in 1:nx){
    ind[[i]] <- ind1[i]:ind2[i]
  }
  ## indicators for the fitted values
  indF <- cbind(w1, w2)
  ## number of total coefficients
  np <- sum(nc)

  ## starting values

  ## fitting P-splines for the aging
  ## taking only ages "x.sen+" and extrapolating
  ## for the previous ages
  y2 <- y[which(w2==1)]
  ww2 <- rep(0,length(y2))
  ww2[x>=x2] <- 1
  fit2 <- suppressWarnings(Mort1Dsmooth(x=x, y=y2, offset=log(e),
                                        w=ww2*wei[w2==1],
                                        method=3, lambda=10^2,
                                        ndx = deg, deg = 3, pord = 2))
  y2.st0 <- exp(XX[[2]] %*% fit2$coef)*e
  y2.st <- rep(0,n)
  y2.st[which(w2==1)] <- y2.st0

  ## fitting P-splines for the accident component
  # y1a <- c(y-y2.st) # deleted this line because sometimes y2.st > y => y1a < 0 => error
  y1a <- y[which(w1 == 1)]
  y1a[which(y1a<=0)] <- 0
  y1 <- y1a[which(w1==1)]
  ww1 <- rep(0, length(x))
  ww1[x>3 & x<x1] <- 1
  wei1 <- wei[which(w1==1)]
  fit1 <- suppressWarnings(Mort1Dsmooth(x=x, y=y1, offset=log(e),
                                        w=ww1*wei1,
                                        method=3, lambda=10^4,
                                        ndx = deg, deg = 3, pord = 3))
  y1.st0 <- exp(XX[[1]] %*% fit1$coef)*e
  y1.st <- rep(0,n)
  y1.st[which(w1==1)] <- y1.st0

  lhaz1.st <- log(y1.st/e)
  lhaz2.st <- log(y2.st/e)
  y.st <- y1.st + y2.st
  lhaz.st <- log(y.st/e)

  ## ## plotting starting values
  ## plot(x, lhaz.act, ylim=c(-12,2))
  ## lines(x, lhaz.st, col=2, lwd=3)
  ## lines(x, lhaz1.st, col=3, lwd=2, lty=3)
  ## lines(x, lhaz2.st, col=4, lwd=2, lty=3)

  ## concatenating starting coefficients
  coef.st <- as.vector(c(fit1$coef,
                         fit2$coef))

  ## shape penalties
  kappa <- 10^6
  ## including log-concaveness for accident-hump
  D1con <- diff(diag(nb), diff=2)
  w1con <- rep(0, nrow(D1con))
  W1con <- diag(w1con)
  ## including monotonicy for aging part
  D2mon <- diff(diag(nb), diff=1)
  w2mon <- rep(0, nrow(D2mon))
  W2mon <- diag(w2mon)

  ## component specific penalty terms
  P1 <- lambda1 * tDD1
  P2 <- lambda2 * tDD2
  ## overall penalty term
  P <- matrix(0,np,np)
  P[1:nb,1:nb] <- P1
  P[1:nb+nb,1:nb+nb] <- P2

  coef <- coef.st

  ## plotting actual log-mortality
  if(plotIT) plot(x, lhaz.act, ylim=c(floor(min(log(y/e)[y != 0])),0))

  Pr <- ridge*diag(np)

  ## iterations for a given lambdas-combination
  for(it in 1:max.it){

    ds <- 1 - it/max.it

    ## penalty for the shape constraints
    P1con <- kappa * t(D1con) %*%W1con%*% D1con
    P2mon <- kappa * t(D2mon) %*%W2mon%*% D2mon
    Psha <- matrix(0,np,np)
    Psha[1:nb,1:nb] <- P1con
    Psha[1:nb+nb,1:nb+nb] <- P2mon


    ## linear predictor
    eta <- numeric(nx*n)
    for(i in 1:nx){
      eta0 <- rep(0, n)
      eta0[which(indF[,i]==1)] <- XX[[i]] %*% coef[ind[[i]]]
      eta[1:n+(i-1)*n] <- eta0
      if(plotIT){ ## plotting each component
        lines(x[which(indF[,i]==1)],
              eta0[which(indF[,i]==1)], col=i+2)
      }
    }
    ## components
    gamma <- exp(eta)*c(indF)
    ## expected values
    mu <- numeric(n)
    for(i in 1:nx){
      mu <- (e * gamma[1:n+(i-1)*n]) + mu
    }
    ## plotting overall log-mortality
    if(plotIT) lines(x, log(mu/e), col=2, lwd=4)
    ## weights for the IWLS
    w <- mu
    ## modified model matrix for a CLM
    U <- matrix(NA, n, sum(nc))
    for(i in 1:nx){
      u <- gamma[1:n+(i-1)*n]/mu * e
      XXi <- matrix(0, nrow=n, ncol=nc[i])
      XXi[which(indF[,i]==1), ] <- XX[[i]]
      U0 <- u * XXi
      U[,ind[[i]]] <- U0
    }
    U[is.nan(U)] <- 0
    ## regression parts for the P-CLM
    tUWU <- t(U) %*% (w*wei * U)
    tUWUpP <- tUWU + P + Psha + Pr
    r <- y - mu
    tUr <- t(U) %*% r
    ## updating coefficients with a d-step
    coef.old <- coef
    coef <- solve(tUWUpP, tUr + tUWU %*% coef)
    coef <- ds*coef.old + (1-ds)*coef
    ## update weights for shape constraints
    ## accident-hump, log-concaveness
    W1con.old <- W1con
    coef1 <- coef[ind[[1]]]
    W1con <- diag(diff(coef1, diff=2) >= 0)
    ## aging, monotonicity
    W2mon.old <- W2mon
    coef2 <- coef[ind[[2]]]
    W2mon <- diag(diff(coef2) <= 0)
    ## convergence criterion for coefficients
    dif.coef <- max(abs(coef.old-coef))/max(abs(coef))
    ## stopping loop at convergence
    if(dif.coef < 1e-04 & it > 4) break
    ## check convergence
    if(mon) cat(it, dif.coef, "\n")
    ## locator(1)
  }
  ## compute deviance
  yy <- y
  yy[y==0] <- 10^-8
  mumu <- mu
  mumu[mu==0] <- 10^-8
  dev <- 2*sum(y * log(yy/mumu))
  ## effective dimensions
  H <- solve(tUWUpP, tUWU)
  diagH <- diag(H)
  ed1 <- sum(diagH[ind[[1]]])
  ed2 <- sum(diagH[ind[[2]]])
  ## BIC
  bic <- dev + log(n)*sum(diagH)

  ## fitted coefficients
  coef.hat <- coef
  ## linear predictor for each component
  etas <- NULL
  for(i in 1:nx){
    etas[[i]] <- XX[[i]] %*% coef.hat[ind[[i]]]
  }
  eta1.hat <- etas[[1]]
  eta2.hat <- etas[[2]]

  ## linear predictor in a matrix over the whole x
  ETA.hat <- matrix(NA, n, nx)
  for(i in 1:nx){
    ETA.hat[which(indF[,i]==1),i] <- etas[[i]]
  }
  ## linear predictor for overall mortality
  eta.hat <- log(apply(exp(ETA.hat), 1, sum, na.rm=TRUE))

  ## fitted values for each component
  gamma1.hat <- exp(eta1.hat)
  gamma2.hat <- exp(eta2.hat)

  ## expected values
  mu.hat <- exp(eta.hat)*e

  ## deviance residuals
  res1 <- sign(y - mu.hat)
  res2 <- sqrt(2 * (y * log(y/mu.hat) - y + mu.hat))
  res <- res1 * res2

  if(plotFIT){
    ## plotting actual and fitted log-mortality
    testofit <- c(expression(paste("ln(",y, "/e)")),
                  expression(paste("ln(", mu, "/e)")),
                  expression(paste("ln(",gamma[1],")")),
                  expression(paste("ln(",gamma[2],")")))
    plot(x, lhaz.act, col=8, pch=16, ylim=c(floor(min(log(y/e)[y!=0])), 0),
         xlab="age",
         ylab="log-mortality",las=1)
    lines(x, eta.hat, col=2, lwd=2)
    lines(x, eta1.hat, col=3)
    lines(x, eta2.hat, col=4)
    legend("top", legend=testofit, col=c(8,2:4),
           lty=c(-1,1,1,1), pch=c(16,-1,-1,-1),ncol=2)
    title(lab)
  }
  ## return object
  out <- list(data=data, coef=coef.hat, XX=XX,
              eta1=eta1.hat, eta2=eta2.hat, eta=eta.hat,
              gamma1=gamma1.hat, gamma2=gamma2.hat,
              mu=mu.hat, resid=res,
              dev=dev, ed1=ed1, ed2=ed2, bic=bic, it=it, maxit=max.it,
              B=B,deg=deg)

  return(out)
}
