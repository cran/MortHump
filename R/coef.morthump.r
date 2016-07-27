#' Extract the coefficients from a morthump object
#'
#' @param object an object of class morthump, resulting from the fit of a given mortality model
#' @param ... other arguments
#'
#' @details The nature of the coefficients depends on the type of mortality model used to fit the mortality schedule.
#' \itemize{
#'  \item
#'  If the model was specified as parametrical,
#'  the coefficients are the parameters as defined by the models' authors. For instance, in the case of the model defined by Heligman and Pollard (1980), the parameters are named A to H and each has a particular interpretation.
#'  \item
#'  If the model was specified as non-parametrical, the coefficients correspond to the B-spline coefficients.
#'  }
#'
#' @return fitted values for the mortality model defined in \link{morthump}
#'
#' @examples
#'
#' data("CHE2010m")
#'
#' # fits the Heligman-Pollard model (parametrical)
#' fit <- morthump(data = CHE2010m, model = "hp")

#' # extract the estimates for the eight coefficients
#' coef(fit)
#'
#' @seealso \code{\link{morthump}}
#'
#' @export
#' @method coef morthump


coef.morthump <- function(object,...){

  fit <- object

  par <- unlist(fit$coef)

  return(par)

}
