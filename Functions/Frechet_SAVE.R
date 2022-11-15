############################################################################
## Sliced Average Variance Estimate and Frechet SAVE##################

save.sdr <- function(x, y, h=NULL){
  # Para:
  #   x: A n by p matrix, standardized predictor matrix
  #   y: A n by 1 matrix
  # Return:
  #   SAVE candidate matrix p by p matrix
  x <- as.matrix(x)
  n <- dim(x)[1]
  p <- dim(x)[2]
  if (is.null(h)){
    h <- max(round(n/(6*p)), 3) 
  } # take the number of slices as larger one between [n/6p] and 3
  
  y.dis <- discretize(y, h)
  y.label <- unique(y.dis)
  prob <- numeric()
  vxy = array(0, c(p, p, h))
  
  for(i in 1:h) {
    prob <- c(prob,length(y.dis[y.dis==y.label[i]])/n)
    vxy[,,i] = var(x[y.dis==y.label[i],])
  }
  savemat <- 0
  for(i in 1:h){
    savemat <- savemat + prob[i] * (diag(p) - vxy[,,i])%*%(diag(p) - vxy[,,i]) / h
  }
  return(savemat)
}


fsave <- function(x, y, d, h=NULL, isGram=TRUE, rho=1){
  # Para:
  #   x: A n by p matrix
  #   y: Response objects:
  #        A n by m matrix for distributional objects; 
  #        A n by d by d array (a list of q by q matrices) for SPD matrices objects;
  #        A n by 1 matrix for spherical objects.
  #     OR n by n Gram matrix of response objects.
  #   d: dimension of central space
  #   isGram: if true, y is gram matrix
  #   rho: tuning parameter in the reproducing kernel
  # Return:
  #   A list containing Frechet SAVE estimator (p by d matrix), and Candidate Matrix (p by p matrix)
  signrt <- mpower(var(x),-1/2)
  z <- scale(x, center = T, scale = F)%*%signrt
  n <- dim(z)[1]
  p <- dim(z)[2]
  if (is.null(h)){
    h <- max(round(n/(6*p)), 3) 
  } # take the number of slices as larger one between [n/6p] and 3
  
  if(!(isGram)) {
    y <- gramw(y, complexity=rho) 
  }
  y.list <- split(y, col(y)) 
  beta <- lapply(y.list, save.sdr, x=z) 
  lambda <- Reduce('+',beta) / n
  B <- eigen(lambda)$vectors[,1:d]
  beta.final <- signrt%*%B
  return(list(beta=beta.final, candidate.matrix=lambda))
}


###################################################################
## test for save #############################
# Example 5.1 in textbook
# n <- 200
# p <- 10
# x <- matrix(rnorm(n*p), nrow=n, ncol = p)
# y <- x[,1]^2 + 2*sin(x[,2]) + rnorm(n,mean=0, sd=0.2)
# signrt <- mpower(var(x),-1/2)
# z <- scale(x,center = T, scale = F)%*%signrt
# signrt%*%eigen(save.sdr(z, y, h=3))$vectors[,1:2]

