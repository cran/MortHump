
# The Heligman-Pollard model of mortality

hp <- function(x, par){
  
  A <- par[1]
  B <- par[2]
  C <- par[3]
  D <- par[4]
  E <- par[5]
  F <- par[6]
  G <- par[7]
  H <- par[8]
  
  y <- A ^ ((x + B) ^ C) + D * exp(-E * (log(x) - log(F)) ^ 2) + (G * (H ^ x))/(1 + G * (H ^ x))
  
  names(y) <- NULL

  return(y)
}