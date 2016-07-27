
# hpk.fit : a function that fits Kostaki's extension of the Heligman-Pollard model of mortality

hpk.fit <- function(data, method = "port", w = 1/data$m, start = NULL){

    # according to Heligman and Pollard (1980), weights should be 1/q
    # according to Brillinger (1986), weights should be 1/d
    # according to Heligman and Pollard (1980), the response variable should be q or q/(1-q), but it is defined here as m


    if(data$x[1] == 0){data$x[1] <- 1e-5}

    form <- as.formula(m ~ A ^ (x + B) ^ C + D * exp(-E * (x <= F) * (log(x) - log(F)) ^ 2) * exp(-k * E * (x > F) * (log(x) - log(F)) ^ 2) + G * H ^ x / (1 + G * (H ^ x)))

    if(is.null(start)){

    start <- list(A = 0.001, B = 0.005, C = 0.11, D = 0.0015, E = 8, F = 20, G = 0.00003, H = 1.105, k = 0.5)

    lower <- c(0.0001, 0.000001, 0.0001, 0, 1, 16, 0.0000001, 0.5, 0)

    upper <- c(0.1, 0.5, 1, 0.01, 50, 22, 0.01, 1.5, 1)

    }else{
      lower <- start$lower
      upper <- start$upper
      start <- start$start
    }

    if(method %in% c("port","lm","gnm","bayes") == F){warning("The method must be one of 'port, 'lm' or 'gnm'.")}

    if(method == "port"){

      fit <- nls(formula = form, data = data, start = start, lower = lower, upper = upper,
                 algorithm = "port", weights = w, control = list(maxiter=1000))
      fit$coef <- coef(fit)
    }

    if(method == "lm"){

      fit <- nlsLM(formula = form, data = data, start = start, lower = lower, upper = upper,
                   weights = w, control = list(maxiter=1000))
      fit$coef <- coef(fit)
    }

    if(method == "gnm"){

      # see Currie 2014 and Debon et al. 2005

      warning("Sorry, the estimation of the Kostaki model as a generalized linear model is not implemented yet.")

    }

    if(method == "bayes"){

      # fun <- function(m){rnorm(n = 8e3, mean = m, sd = m/10)}

      # prior <- do.call(cbind,lapply(start,fun))

      # fit <- hp.bm.imis(prior = prior, nrisk = data$exp, ndeath = data$d, K = 10)

      warning("Sorry, the estimation of the Kostaki model using Bayesian statistics is not implemented yet.")

    }

    if(data$x[1] < 1){data$x[1] <- 0}

    fit$data <- data
    fit$method <- method
    fit$w <- w
    fit$type <- "parametric"


    return(fit)

}
