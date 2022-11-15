############################################################################
## Iterative Hessian Transformation and Frechet Iterative Hessian Transformation##################


# Function to implement IHT
iht <- function(z, y){
  # Para:
  #   z: A n by p matrix
  #   y: A n by 1 matrix
  # Return:
  #   IHT candidate matrix, p by p matrix
  z <- as.matrix(z)
  y <- as.vector(y)
  p <- dim(z)[2]
  n <- dim(z)[1]
  szy <- cov(z, y)
  szzy <- (t(z * (y-mean(y))) %*% z) /n
  imat <- szy
  for(i in 1:(p-1)) {
    imat=cbind(imat, mpower(szzy,i)%*%szy)
  }
  return(imat%*%t(imat))
}

# Frechet IHT
fiht <- function(x, y, d, isGram=TRUE, rho=1){
  # Para:
  #   x: A n by p matrix
  #   y: Response objects:
  #        A n by m matrix for distributional objects; 
  #        A n by d by d array (a list of q by q matrices) for SPD matrices objects;
  #        A n by 1 matrix for spherical objects.
  #     OR n by n Gram matrix of response objects.
  #   d: dimension of sdr space
  #   isGram: if true, y is gram matrix
  #   rho: tuning parameter in the reproducing kernel
  # Return:
  #   A list containing Frechet IHT estimator (p by d matrix), and Candidate Matrix (p by p matrix)
  signrt <- mpower(var(x),-1/2)
  z <- scale(x, center = T, scale = F)%*%signrt
  n <- dim(x)[1]
  p <- dim(x)[2]
  if(!(isGram)) {
    y <- gramw(y, complexity=rho) 
  }
  y <- split(y, col(y))
  beta <- lapply(y, iht, z=z) 
  lambda <- Reduce('+', beta)/n
  B <- eigen(lambda)$vectors[,1:d]
  beta.final <- signrt%*%B
  return(list(beta = beta.final, candidate.matrix=lambda))
}


