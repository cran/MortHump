# function to isolate and/or group cod from the Mxc matrix
cod.data <- function(k, Mxc, exp){

  if(!(all(k < 0) | all(k > 0))){break ; warning("Cod must be all included or all excluded")}

  if(length(k) == 1){if(k > 0){m <- Mxc[,k]}else{m <- rowSums(Mxc[,k])}}else{
    m <- rowSums(Mxc[,k])}

  if(!is.null(dim(m))){warning("Error")}

  d <- m * exp

  return(data.frame(m,d,exp))
}
