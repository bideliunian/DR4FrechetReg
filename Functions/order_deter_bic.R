#############################################################################
#################  BIC-type criteria for order determination #############

bic <- function(x, y, criterion, method){
  # Para:
  #   x: A n by p design matrix
  #   y: n by n Gram matrix of response objects.
  #   criterion: 'lal' or 'zmp'
  #   method: "fols", "fopg", "fiht", "fsir", "fsave", "fdr", "fopg"
  # Return:
  #   A list containing estimated dimension, and values of loss function 
  x <- as.matrix(x)
  n <- dim(x)[1]
  p <- dim(x)[2]
  d0 <- 2
  
  f <- get(method)
  mat <- f(x=x, y=y, d=1)$candidate.matrix # candidate matrix
  eigen_mat <- eigen(mat)
  lam <- eigen_mat$values
  
  if(criterion=="lal"){
    gn <- c(0)
    for(k in 1:p){
      gn <- c(gn, sum(lam[1:k])-(2*lam[1])*n^(-1/2)*(log(n))^(1/2)*k)
    }
  }
  if(criterion=="zmp"){
    gn <- numeric()
    for(k in 0:(p-1)){
      c1 = (lam[1]/3)*(0.5* log(n)+0.1* n^(1/3))/(2*n)
      c2 = k*(2*p-k+1)
      gn = c(gn, sum(log(lam[(k+1):p]+1)-lam[(k+1):p])-c1*c2)
    }
    gn = c(gn,-c1*p*(2*p-p+1))
  }
  
  return(list(rhat = which.max(gn)-1, rcurve = gn))
}
