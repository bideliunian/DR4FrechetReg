############################################################################
## Minimum Average Variance Estimate and Frechet MAVE##############################


fmave <- function(x, y, d, initial, isGram=TRUE, rho=1, niter=5, h=NULL){
  # Para:
  #   x: A n by p design matrix
  #   y: Response objects:
  #        A n by m matrix for distributional objects; 
  #        A n by d by d array (a list of q by q matrices) for SPD matrices objects;
  #        A n by 1 matrix for spherical objects.
  #     OR n by n Gram matrix of response objects.
  #   d: dimension of central space
  #   isGram: if true, y is gram matrix
  #   rho: tuning parameter in the reproducing kernel
  #   niter: number of iteration for refined mave
  #   initial: initial point of mave, take the solution of opg by default
  #   h: tuning parameter for the smooth kernel of x
  # Return:
  #   A list containing Frechet MAVE estimator (p by d matrix), and Candidate Matrix (p by p matrix)
  sig <- diag(var(x))
  z <- apply(x, 2, scale, center = TRUE, scale = TRUE)
  n <- dim(z)[1]
  p <- dim(z)[2]
  
  c0 <- 2.34
  p0 <- max(p,3)
  rn <- n^(-1/(2*(p0+6)))
  
  if(is.null(h)){
    h <- c0*n^(-(1/(p0+6)))
  }
  if(!(isGram)) {
    y <- gramw(y, complexity=rho)
  }
  # by default, use opg estimator as initial
  if(is.null(initial)) {
    bbt <- fopg(x, y, d, isGram, rho, 1)$beta 
  }
  else {
    bbt <- initial 
  }
  for(iit in 1:niter){
    kmat <- kern(z%*%bbt, h, type='Gaussian')
    mkmat <- apply(kmat, 2, mean)
    kmat <- t(kmat*(1/mkmat))
    acarray <- array(0, c(1+d,n,n))
    for(i in 1:n) {
      acarray[,,i] <- swls(y, kmat, z%*%bbt, i)$ab
    }
    vecb1 <- 0
    vecb2 <- 0
    for(i in 1:n) for(j in 1:n) {
      tmp <- t(t(z) - z[i,])
      tmp1 <- kmat[,i]
      tmp3 <- acarray[2:(1+d), j, i]
      vecb1 <- vecb1 + kronecker(t(tmp*tmp1)%*%tmp, tmp3%*%t(tmp3))
      tmp4 <- kmat[,i]*(y[,j] - acarray[1, j, i])
      vecb2 <- vecb2 + kronecker(t(tmp), tmp3)%*%tmp4
    }
    bbt <- t(matrix(solve(vecb1)%*%vecb2, d, p))
    h <- max(rn*h,c0*n^((-1/(d+4))))
  }
  beta.final <- diag(sig^(-1/2))%*%bbt
  return(beta.final)
}

