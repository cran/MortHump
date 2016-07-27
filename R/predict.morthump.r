# Predicted mortality rates based on the model used in the fitting of the morthump object
#' @method predict morthump

predict.morthump <- function(object, x = NULL, which = "fitted", tol = NULL,...){

  fit <- object

  model <- attr(fit, "model")

  if(is.null(x)){x <- fit$data$x}

  if(model == "hp"){

    if(which == "fitted"){y <- hp(x = x, par = coef(fit))}
    if(which == "hump" ){
      par <- coef(fit)
      par[-(4:6)] <- 0
      y <- hp(x = x, par = par)}

  }

  if(model == "hps"){

    if(which == "fitted"){y <- hps(x = x, par = coef(fit))}
    if(which == "hump" ){y <- rep(0, length(x))}

  }

  if(model == "hpk"){

    if(which == "fitted"){y <- hpk(x = x, par = coef(fit))}
    if(which == "hump" ){
      par <- coef(fit)
      par[-c(4:6,9)] <- 0
      y <- hpk(x = x, par = par)}

  }

  if(model == "sse"){

      decim <- Vectorize(function(x){min(which(x*10^(0:20) == floor(x*10^(0:20)))) - 1})

      if(all(x == round(x))){k <- 0}else{k <- max(decim(x))} # number of decimals required for the given input ages

      if(is.null(tol) == FALSE){k <- min(c(k,tol))} # limited if required

      n <- length(x)

      if(n == 1){x <- x + c(- 10^(-k), 0 , 10^(-k))} # for some reasons, the XX %*% coef step won't take a single value

      outer1 <- min(fit$data$x) # limits need to be identical to the original model to preserve the base
      outer2 <- max(fit$data$x)
      inner1 <- min(x)
      inner2 <- max(x)

      x1 <- unique(c(outer1, seq(inner1, inner2, 10^-k), outer2))

      xl <- min(x1)
      xr <- max(x1)
      xmax <- xr + 0.01 * (xr - xl)
      xmin <- xl - 0.01 * (xr - xl)

      deg2 <- fit$deg2
      deg3 <- fit$deg3

      B1 <- cbind(1,1/x1)
      B2 <- MortSmooth_bbase(x1, xmin, xmax, deg2, 3)
      B3 <- MortSmooth_bbase(x1, xmin, xmax, deg3, 3)
      XX <- list(X1=B1, X2=B2, X3=B3)

      XX <- lapply(XX, function(mat){mat[which(round(x1,k) %in% round(x,k)),]}) # only keep the ages we are interested in

      coef <- fit$coef
      nx <- length(XX)
      nc <- unlist(lapply(XX, ncol))
      ind2 <- cumsum(nc)
      ind1 <- c(1, ind2[1:(nx-1)]+1)
      ind <- NULL
      for(i in 1:nx){
        ind[[i]] <- ind1[i]:ind2[i]
      }

      etas <- NULL
      for(i in 1:nx){
        etas[[i]] <- XX[[i]] %*% coef[ind[[i]]]
      }
      eta1.hat <- etas[[1]]
      eta2.hat <- etas[[2]]
      eta3.hat <- etas[[3]]

      mhat <- lapply(etas, exp)

      if(n == 1){mhat <- lapply(mhat, function(mat){mat[2,]})}

    if(which == "fitted"){
      y <- rowSums(do.call(cbind,mhat))}
    if(which == "hump"){
      y <- mhat[[3]]
    }

  }

  return(as.vector(y))

}
