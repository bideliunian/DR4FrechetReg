############################# Weighted Inverse Regression Ensemble ###############################
# Chao Ying and Zhou Yu, Fr´echet Sufficient Dimension Reduction for Random Objects

wire <- function(x, y, d, isGram=TRUE, rho=1){
  # Para:
  #   x: A n by p matrix
  #   y: Response objects:
  #        A n by m matrix for distributional objects;
  #        A n by d by d array (a list of q by q matrices) for SPD matrices objects;
  #        A n by 1 matrix for spherical objects.
  #     OR n by n Gram matrix of response objects.
  #   d: dimension of SDR space
  #   isGram: if true, y is gram matrix
  #   rho: tuning parameter in the reproducing kernel
  # Return:
  #   WIRE estimator (p by d matrix)
  x <- as.matrix(x)
  z <- scale(x, center = T, scale = F)
  n <- dim(z)[1]
  p <- dim(z)[2]
  
  if(!(isGram)) {
    y <- gramw(y, complexity=rho)
  }
  
  g <- expand.grid(row = 1:n, col = 1:n)
  dist <- function(i, k){
      return(as.matrix(y[i, k]*(z[i,])%*%t(z[k,])))
  }
  k <- mapply(dist, i=g[,1], k=g[,2], SIMPLIFY = FALSE)
  kmat <- Reduce('+',k)/(n*(n-1))
  beta.final <- eigen(solve(var(x))%*%kmat)$vectors[,1:d]
  return(list(beta=beta.final, candidate.matrix=kmat))
}


