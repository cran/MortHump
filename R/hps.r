
# An adaptation of the Siler model designed as a specific case of the Heligman-Pollard model of mortality with no hump

hps <- function(x, par){
  
  A <- par[1]
  B <- par[2]
  C <- par[3]
  D <- par[4]
  G <- par[5]
  H <- par[6]
  
  y <- A ^((x + B) ^ C) + D  + (G * (H ^ x))/(1 + G * (H ^ x))

  return(y)
}