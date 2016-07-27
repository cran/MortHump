
#' @title Isolate the young adult mortality hump from a set of age-specific mortality rates
#'
#' @description This function estimates a model of mortality on the provided set of age-specific death rates. Both parametric and non-parametric models are available.
#' They are all designed to estimate the size of the young adult mortality hump, i.e. the deviation in the force of mortality often observed during adolescence and early adulthood.
#'
#' @param data data frame produced with \code{HMD2MH} or similarly structured
#' @param model the name of the mortality model to be fitted
#' @param method the optimization algorithm to be used
#' @param w weights for the weighted least square estimation (parametric models only)
#' @param start optional list of three elements giving the starting, lower and upper values of the parameters (parametric models only)
#' @param maxit maximum number of iterations (\link{sse.fit})
#' @param x1 upper age limit for the hump component (\link{sse.fit})
#' @param x2 lower age limit for the senescence component (\link{sse.fit})
#' @param lambda.hump smoothing parameter for the hump component (\link{sse.fit})
#' @param lambda.sen smoothing parameter for the senescence part (\link{sse.fit})
#'
#' @details
#' Mortality models can be broadly divided into parametric and non-parametric types.
#'
#' \bold{Parametric models}
#'
#' Dozens of parametric models have been published (see Wunsch et al. 2002 for a review), but only three models are implemented here.
#' They are all based on the structure of the Heligman-Pollard model and are conceived as nested models that offer a range of ways to address the specificity of young adult mortality.
#' The \code{"hps"} model is the simplest and does not take the hump into account. The \code{"hp"} model assumes a symetrical hump with varying height (\code{D}), spread (\code{E}) and location (\code{F}).
#' The \code{"hpk"} model relaxes the assumption of symetry by introducing a different spread before and after the peak of the hump.
#' Algebraically, \code{"hp"} is equal to \code{"hps"} iff \eqn{E = 0}, and \code{"hpk"} is equal to \code{"hp"} iff \eqn{E1 = E2}.
#'
#' \code{"hp"} : Published by Heligman and Pollard (1980), this model contains eight parameters in three additive terms.
#' \deqn{A ^ ((x + B) ^ C) + D * exp(-E * (log(x) - log(F)) ^ 2) + (G * (H ^ x))/(1 + G * (H ^ x))}
#' where (1) A, B and C describe the infant mortality component, (2) D, E and F describe the hump component, and (3) G and H describe the senescence component.
#'
#' \code{"hpk"} : Published by Kostaki (1992), this extension of \code{"hp"} contains nine parameters in three additive terms.
#' \deqn{A ^ ((x + B) ^ C) + D * exp(-E1 * (log(x) - log(F)) ^ 2) + (G * (H ^ x))/(1 + G * (H ^ x)) for x \le F}
#' \deqn{A ^ ((x + B) ^ C) + D * exp(-E2 * (log(x) - log(F)) ^ 2) + (G * (H ^ x))/(1 + G * (H ^ x)) for x > F}
#' In contrast to the \code{"hp"} model, it allows the hump to be asymetrical by differenciating the spread of the hump before (\code{E1}) and after (\code{E2}) its peak.
#'
#' \code{"hps"} : This model is an adaptation of the \code{"hp"} model that does not account for the young adult hump.
#' It resembles in this respect the model published by Siler (1979), but is nested in the \code{"hp"} model.
#' \deqn{A ^ ((x + B) ^ C) + D  + (G * (H ^ x))/(1 + G * (H ^ x))}
#' In contrast to the \code{"hp"} model, it only has one parameter for the hump component (\code{D}) that has the same function as the constant term of the Makeham model (1860).
#'
#' \bold{Non-parametric model}
#'
#' The age-specific mortality rates can alternatively be fitted with non-parametric models, but very few of them adopt the multiple-component approach that is often found in parametric models.
#' In order to isolate the young adult mortality hump from the rest of the force of mortality, a non-parametric model should also include the notion of additive components. This is what makes the so-called Sum of Smooth Exponentials model an attractive alternative to parametric models.
#'
#' Briefly, the SSE model describes the observed force of mortality over age \eqn{\bm{\mu}}{\mu} as the sum of three vectors \eqn{[\bm{\gamma}_{1}:\bm{\gamma}_{2}:\bm{\gamma}_{3}]}{<\gamma_1:\gamma_2:\gamma_3>} over age. The subscripts (1 to 3) denote which mortality component each \eqn{\bm{\gamma}_{j}}{\gamma_j} refers to: infant, midlife and old-ages mortality, respectively.
#' In other words it assumes that deaths are realizations from Poisson distributions with mean composed of three parts:
#'
#' \deqn{\bm{y} \sim \mathcal{P}(\bm{e} \,\bm{\mu} = \bm{C} \,\bm{\gamma})\,}{y ~ P(e \mu = C \gamma) }
#'
#' Given \eqn{m} age groups, the matrix \eqn{\bm{C}}{C} is a \eqn{m} by \eqn{3} matrix repeating the exposures in each column, which, when multiplied with the \eqn{\bm{\gamma}_{j}}{\gamma_j}, sums the components and transforms the result into death counts. In matrix notation, \eqn{\bm{C}}{C} is given by
#'
#' \deqn{\bm{C} = \bm{1}_{1,3} \otimes diag (\bm{e})\,}{C = 1_{1,3} x diag(e)}
#'
#' where \eqn{\bm{1}_{1,3}}{1_{1,3}} is a \eqn{1\times 3}{1 by 3} matrix of ones, \eqn{diag(\bm{e})}{diag(e)} is the diagonal matrix of the exposure population and \eqn{\otimes}{x} denotes the Kronecker product.
#'
#' Unlike parametric models, in the \code{SSE} model there is no need to make strong assumptions about the functional form of each component.
#' For each component a discrete sequence is assumed and the exponential function is used to ensure non-negative elements:
#'
#' \deqn{\bm{\gamma}_{j} = \exp(\bm{X}_{j} \bm{\beta}_{j}) \, , \quad j=1,2,3\,}{\gamma_j = exp(X_j \beta_j),  j = 1,2,3}
#'
#' In other words, each component has to be described by a linear combination of a model matrix \eqn{\bm{X}_{j}}{X_j} and associated coefficients \eqn{\bm{\beta}_{j}}{\beta_j}.
#' The design matrices \eqn{\bm{X}_{j}}{X_j} can represent parametric or non-parametric structures. In this way the composite force of mortality \eqn{\bm{\mu}}{\mu} can be viewed as sum of \eqn{3} exponential components, which potentially can be smooth.
#' Furthermore, the \code{SSE} model allows incorporating shape constraints to enforce senescent and young-adult components to be monotonically increasing and log-concave, respectively.
#' A Penalized Composite Link Model (PCLM) has been proposed to estimate the SSE model. For more details, see Camarda et al. (2016).
#'
#'
#' \bold{Further Notes}
#'
#' All models use the age-specific mortality rates (\code{mx}) as the response variable in order to preserve comparability.
#' This is despite the original definitions that often use the mortality quotients (\code{qx}) or one of its transforms.
#'
#' By default, the parametric models are estimated with the "port" algorithm from the \link{nls} function. Is some rare cases, it may become stuck into local minimums. If this happens, try switching to the Levenberg-Marcquart algorithm.
#' The non-parametric models use a tailored type of optimization algorithm (see \link{sse.fit}).
#'
#' The optimization algorithm uses starting values, as well as lower and upper boundaries for each parameter. Most of the time, these values will be good enough to get a convergence.
#' If you are not satisfied with the fitted values (for instance, if one of the fitted parameter is a round number there is a good chance that the algorithm encountered a boundary value), you can specify your own starting, lower and upper values.
#' To change these values, use the \code{start} argument and construct it as a list of three objects: \code{start}, \code{lower} and \code{upper}. The first argument, start, defines the starting values for each of the parameters and is itself a named list (see example).
#' The other two elements defined the lower and upper boundaries of the parameters and are defined as vectors.
#'
#' @return Returns an object of class \code{morthump} containing the arguments used to fit the model as well as the estimated coefficients.
#'
#' @examples
#'
#' data("CHE2010m")
#'
#' # fits the Heligman-Pollard model (parametrical)
#' fit <- morthump(data = CHE2010m, model = "hp")
#'
#' # use custom starting, lower and upper values for the parameters
#' sv <- list(start = list(A=0.001, B=0.005, C=0.11, D=0.0015, E=8, F=20, G=3e-5, H=1.105),
#'            lower = c(0.0001, 0.000001, 0.0001, 0, 1, 16, 0.0000001, 0.5),
#'            upper = c(0.1, 0.5, 1, 0.01, 50, 30, 0.01, 1.5))
#'
#' fit <- morthump(data = CHE2010m, model = "hp", start = sv)
#'
#' # fits a Sum of Smooth Exponentials (non-parametrical)
#' fit <- morthump(data = CHE2010m, model = "sse")
#'
#' @references
#'
#' Camarda, C. G., Eilers, P. H. C., & Gampe, J. (2016). Sums of smooth exponentials to decompose complex series of counts. Statistical Modelling.
#'
#' Wunsch, G., Mouchart, M., & Duchene, J. (2002). The life table : modelling survival and death. Dordrecht ; London: Kluwer Academic.
#'
#' Heligman, L., & Pollard, J. H. (1980). The age pattern of mortality. Journal of the Institute of Actuaries, 107.
#'
#' Kostaki, A. (1992). A nine-parameter version of the Heligman-Pollard formula. Mathematical Population Studies, 3(4), 277-288.
#'
#' Siler, W. (1979). A Competing-Risk Model for Animal Mortality. Ecology, 60(4), 750-757.
#'
#' @seealso
#' \link{sse.fit}, \link{summary.morthump}, \link{plot.morthump}
#'
#' @export
#'
#' @import MortalitySmooth
#' @import Matrix
#'
#' @importFrom graphics abline arrows axis box legend lines matplot par plot points polygon segments text title
#' @importFrom stats as.formula coef cutree dist fitted glm hclust integrate loess loess.control nls optimize pf poisson prcomp predict smooth.spline spline
#' @importFrom utils data read.csv

morthump <- function(data, model, method = "port", w = 1/data$m, start = NULL, maxit = 500, x1 = 30, x2 = 50, lambda.sen = 100, lambda.hump = 1){

  if(model == "hp"){

    fit <- hp.fit(data, method = method, w = w, start = start)
    attr(fit,"model") <- "hp"

    }

  if(model == "hpk"){

    fit <- hpk.fit(data, method = method, w = w, start = start)
    attr(fit,"model") <- "hpk"

  }

  if(model == "hps"){

    fit <- hps.fit(data, method = method, w = w, start = start)
    attr(fit,"model") <- "hps"

  }

  if(model == "sse"){

    fit <- sse.fit(data, maxit = maxit, x.hump = x1, x.sen = x2, lambda.sen = lambda.sen, lambda.hump = lambda.hump)
    attr(fit,"model") <- "sse"

  }

  class(fit) <- "morthump"

  return(fit)
}
