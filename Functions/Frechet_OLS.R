############################################################################
## Ordinary Least Squares and Frechet Ordinary Least Squares##################

# Ordinary Least Squares to estimate central mean subspace
ols <- function(x, y) {
  # Para:
  #   x: A n by p matrix
  #   y: A n by 1 matrix
  # Return:
  #   (Sigma_XX)^inv(Sigma_XY)
  y <- as.matrix(y)
  x <- as.matrix(x)
  return(solve(var(x))%*%cov(x,y))
}


# Frechet OLS
fols <- function(x, y, d, isGram=TRUE, rho=1){
  # Para:
  #   x: A n by p matrix
  #   y: Response objects:
  #        A n by m matrix for distributional objects; 
  #        A n by d by d array (a list of q by q matrices) for SPD matrices objects;
  #        A n by 1 matrix for spherical objects.
  #     or n by n Gram matrix of response objects.
  #   d: dimension of central space
  #   isGram: if true, y is gram matrix
  #   rho: tuning parameter in the reproducing kernel
  # Return:
  #   A list containing Frechet OLS estimator (p by d matrix), and Candidate Matrix (p by p matrix)
  x <- as.matrix(x)
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  if(!(isGram)) {
    y <- gramw(y, complexity=rho)
  }
  y <- as.matrix(y)
  
  # standardize x
  sig.nrt <- mpower(var(x), -1/2)
  x.std <- scale(x, center = T, scale = F)%*%sig.nrt
  
  beta <- apply(y, MARGIN = 2, cov, x=x.std)
  lambda <- beta%*%t(beta)/n
  B <- eigen(lambda)$vectors[,1:d]
  beta.final <- sig.nrt%*%B
  
  return(list(beta=beta.final, candidate.matrix=lambda))
}
