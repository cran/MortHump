# Kostaki's extension of the Heliman-Pollard model of mortality

hpk <- function(x, par){

  A <- par[1]
  B <- par[2]
  C <- par[3]
  D <- par[4]
  E <- par[5]
  F <- par[6]
  G <- par[7]
  H <- par[8]
  k <- par[9]

  y <- A ^ (x + B) ^ C + D * exp(-E * (x <= F) * (log(x) - log(F)) ^ 2) * exp(-k * E * (x > F) * (log(x) - log(F)) ^ 2) + G * H ^ x / (1 + G * (H ^ x))

  return(y)
}
