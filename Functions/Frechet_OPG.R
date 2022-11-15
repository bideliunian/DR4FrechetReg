############################################################################
## Outer Product of Gradient and Frechet OPG##############################

########## simultaneous weighted least squares ###############
swls <- function(hmat, kmat, x, i){
  wi <- diag(kmat[,i])
  xdi <- t(t(x)-x[i,])
  xdi1 <- cbind(1,xdi)
  abmat <- solve(t(xdi1)%*%wi%*%xdi1)%*%t(xdi1)%*%wi%*%hmat
  return(list(a = abmat[1,], b = abmat[-1,], ab = abmat))
}

# Frechet OPG
fopg <- function(x, y, d, isGram=TRUE, rho=1, niter=5, h=NULL){
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
  #   niter: number of iteration for refined opg 
  #   h: tuning parameter for the smooth kernel of x
  # Return:
  #   A list containing Frechet OPG estimator (p by d matrix), and Candidate Matrix (p by p matrix)
  sig <- diag(var(x))
  z <- apply(x, 2, scale, center = TRUE, scale = TRUE)
  n <- dim(z)[1]
  p <- dim(z)[2]
  c0 <- 2.34
  p0 <- max(p, 3)
  rn <- n^(-1/(2*(p0+6)))
  if(is.null(h)){
    h <- c0*n^(-(1/(p0+6))) 
  }
  # if(is.null(b)){
  #   b <- c0*n^(-(1/(p0+5))) 
  # }
  if(!(isGram)) {
    y <- gramw(y, complexity=rho)
  }
  B <- diag(p)
  for(iit in 1:niter){
    kmat <- kern(z%*%B, h, type = 'Gaussian')
    bmat <- numeric()
    for(i in 1:n){
      bmat <- cbind(bmat, swls(y, kmat, z, i)$b)
    } 
    mat <- bmat%*%t(bmat)/(n^2)
    B <- eigen(mat)$vectors[,1:d]
    h <- max(rn*h, c0*n^((-1/(d+4))))
  }
  beta.final <- diag(sig^(-1/2))%*%B
  return(list(beta=beta.final, candidate.matrix=mat))
}

