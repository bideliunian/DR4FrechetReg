######################################
## kernel pca with gram matrix input 
###################################

kpca <- function(x, s=5){
  # para: 
  #   x is n by n gram matrix
  #   s is the dimension of pc's
  # return:
  #   n by s matrix containing the first s pc's
  x <- as.matrix(x)
  n <- dim(x)[1]
  one <- rep(1,n)
  Q <- diag(n)-one%*%t(one)/n
  Gx <- Q%*%x%*%Q
  v <- eigen(Gx)$vectors[,1:s]
  u <- mppower(Gx, -1/2, 10^(-8))%*%v
  return(Gx%*%u)
}
