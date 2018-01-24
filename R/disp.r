#' Generate a grid layout for panel plots
#'
#' @description For a given number of plots to be displayed in a panel, this function suggests inputs for the \code{par(mfrow)} graphic option.
#'
#' @param n number of plots to be displayed
#' @param horiz should horizontal growth be favored instead of vertical (default)
#'
#' @details
#' This simple function is designed to generate the "squarest" layout grid given the number of plots to be displayed in the same panel.
#' If \code{n} has an integer square root, the dimensions are simply \code{sqrt(n) x sqrt(n)}. Otherwise, the function starts from the next smaller square grid and adds progressively more rows.
#'
#' @return disp returns a vector of length two that indicates the number of rows and columns to be passed to \code{par(mfrow)}.
#'
#' @examples
#' # Example with a random number of plots
#' n <- sample(2:16, 1)
#' par(mfrow = disp(n), mar = c(2,2,2,2))
#' for(i in 1:n){
#'   plot(rnorm(100),rnorm(100))
#' }
#' par(mfrow = c(1,1), mar = c(5,4,4,2) + .1)
#'
#' # Display grid for up to 20 plots
#' par(mfrow = disp(20), mar = c(2,2,2,2))
#' for(i in 1:20){
#'   mat <- matrix(NA,nrow = disp(i)[1], ncol = disp(i)[2])
#'   image(mat, axes = FALSE, main = i)
#'   grid(nx = disp(i)[1], ny = disp(i)[2], col = 1)
#'   box()
#' }
#' par(mfrow = c(1,1), mar = c(5,4,4,2) + .1)
#'
#'
#'
#' @export

disp <- function(n, horiz = FALSE){

  r <- c <- floor(sqrt(n))
  
  if(horiz == FALSE){
    while(r*c < n){r <- r + 1}}else{
    while(r*c < n){c <- c + 1}}

  return(c(r,c))
}
