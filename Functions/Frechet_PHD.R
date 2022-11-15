############################################################################
## Principal Hessian Direction and Frechet PHD##################


phd <- function(x, y){
  # Para:
  #   x: A n by p matrix
  #   y: A n by 1 matrix
  # Return:
  #   IHT candidate matrix, p by p matrix
  x <- as.matrix(x)
  y <- as.vector(y)
  n <- length(y)
  return(t(x*(y-mean(y)))%*%x%*%t(t(x*(y-mean(y)))%*%x)/n)
}

fphd <- function(x, y, d, isGram=TRUE, rho=1){
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
  #   A list containing Frechet PHD estimator (p by d matrix), and Candidate Matrix (p by p matrix)
  signrt <- mpower(var(x),-1/2)
  x.std <- scale(x,center = T, scale = F)%*%signrt
  n <- dim(x.std)[1]
  p <- dim(x.std)[2]
  if(!(isGram)) {
    y <- gramw(y, complexity=rho)
  }
  beta <- apply(y, MARGIN = 2, phd, x=x.std)
  y.list <- split(y, col(y))
  beta <- lapply(y.list, phd, x=x.std)
  lambda <- Reduce('+', beta)/n
  B <- eigen(lambda)$vectors[,1:d]
  beta.final <- signrt%*%B
  return(list(beta=beta.final, candidate.matrix=lambda))
}



