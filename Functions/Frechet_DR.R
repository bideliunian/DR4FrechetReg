############################################################################
## Directional Regression and Frechet DR##################

## Directional Regression for Central Space Estimation
dr <- function(x, y, h=NULL){
  # Para:
  #   x: A n by p matrix, standardized design matrix
  #   y: A n by 1 matrix
  # Return:
  #   DR candidate matrix p by p matrix
  x <- as.matrix(x)
  p <- ncol(x)
  n <- nrow(x)
  if (is.null(h)){
    h <- max(round(n/6/p), 3)
  }
  
  y.dis <- discretize(y, h)
  y.label <- unique(y.dis)
  prob <- numeric()
  
  for(i in 1:h) prob <- c(prob, length(y.dis[y.dis==y.label[i]])/n)
  vxy <- array(0, c(p, p, h))
  exy <- numeric()
  for(i in 1:h) {
    vxy[,,i] <- var(x[y.dis==y.label[i],])
    exy <- rbind(exy, apply(x[y.dis==y.label[i],], 2, mean))
  }
  mat1 <- matrix(0, p, p)
  mat2 <- matrix(0, p, p)
  for(i in 1:h){
    mat1 <- mat1 + prob[i] * (vxy[,,i] + exy[i,]%*%t(exy[i,]))%*%
      (vxy[,,i] + exy[i,]%*%t(exy[i,]))
    mat2 <- mat2 + prob[i] * exy[i,]%*%t(exy[i,])
  }
  mat <- 2*mat1 + 2*mat2%*%mat2 + 2*sum(diag(mat2))*mat2 - 2*diag(p)
  return(mat)
}

## Frechet DR
fdr <- function(x, y, d, h=NULL, isGram=TRUE, rho=1){
  # Para:
  #   x: A n by p matrix
  #   y: Response objects:
  #        A n by m matrix for distributional objects; 
  #        A n by d by d array (a list of q by q matrices) for SPD matrices objects;
  #        A n by 1 matrix for spherical objects.
  #     OR n by n Gram matrix of response objects.
  #   d: dimension of SDR space
  #   h: number of slices
  #   isGram: if true, y is gram matrix
  #   rho: tuning parameter in the reproducing kernel
  # Return:
  #   A list containing Frechet SIR estimator (p by d matrix), and Candidate Matrix (p by p matrix)
  signrt <- mpower(var(x),-1/2)
  z <- scale(x,center = T, scale = F)%*%signrt
  
  n <- dim(z)[1]
  p <- dim(z)[2]
  if(!(isGram)) {
    y = gramw(y, complexity=rho)
  }
  y <- split(y, col(y))
  beta <- lapply(y, dr, x=z) 
  lambda <- Reduce('+', beta)/n
  B <- eigen(lambda)$vectors[,1:d]
  beta.final <- signrt%*%B
  return(list(beta=beta.final, candidate.matrix=lambda))
}


# Example 5.1 in textbook
# n <- 200
# p <- 10
# x <- matrix(rnorm(n*p), nrow=n, ncol = p)
# y <- x[,1]^2 + 2*sin(x[,2]) + rnorm(n,mean=0, sd=0.2)
# signrt <- mpower(var(x),-1/2)
# z <- scale(x,center = T, scale = F)%*%signrt
# print('dr estimator')
# signrt%*%eigen(dr(z, y))$vectors[,1:2]
# 
# y_dist <- as.matrix(dist(y))
# fdr(x, y_dist, d=2)$beta
