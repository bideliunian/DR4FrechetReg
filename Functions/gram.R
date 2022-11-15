##### compute the Gaussian kernel/ Epanechnikov kernel Gram matrix ###########
kern <- function(x, h, type='Gaussian'){
  # Para: 
  #   x: design matrix
  #   h: bindwidth for kernel smoother
  #   type: 'Gaussian' or 'Epan'
  x <- as.matrix(x)
  n <- dim(x)[1]
  k2 <- x%*%t(x)
  k1 <- t(matrix(diag(k2),n,n))
  k3 <- t(k1)
  k <- k1-2*k2+k3
  
  if(type=='Gaussian')
    return((1/h)*exp(-(1/(2*h^2))*k))
  
  else if(type=='Epan')
    return((1/h)*3/4*pmax(1-k/(h^2),0))
}

########################## compute distance matrix of spd matrices ########################
dist_spd <- function(x, metric){
  # Para: 
  #   x: list of spd matrices
  #   metric: "Frobenius", "log_Frobenius", or "affine_invariant"
  # Return:
  #   distance matrix
  n <- length(x)
  g <- expand.grid(row = 1:n, col = 1:n)
  ind <- g[upper.tri(diag(n), diag = FALSE), ]
  
  METRICS <- c("Frobenius", "log_Frobenius", "affine_invariant")
  metric <- pmatch(metric, METRICS)
  
  if (metric == 1) {
    distance <- function(i, k){return(frobenius(x[[i]]-x[[k]]))}
  }
  if (metric == 2) {
    logx <- lapply(x, mlog)
    distance <- function(i, k){return(frobenius(logx[[i]]-logx[[k]]))}
  }
  if (metric == 3) {
    distance <- function(i, k){
      s2 <- mpower(x[[k]], -1/2)
      s <- s2%*%x[[i]]%*%s2
      return(frobenius(mlog(s)))
    }
  }
  
  dmatrix <- matrix(0, nrow = n, ncol = n)
  dmatrix[upper.tri(dmatrix,diag = FALSE)] <- mapply(distance, i=ind[,1], k=ind[,2])
  dmatrix <- dmatrix+t(dmatrix)
  
  return(dmatrix)
}

###################### compute wasserstein distance matrix of distributional data ###########
dist_dist <- function(x){
  # Para: 
  #   x: n by m matrix, each row is iid observations from a distribution 
  # Return:
  #   wasserstein distance matrix
  x <- as.matrix(x)
  n <- dim(x)[1]
  g <- expand.grid(row = 1:n, col = 1:n)
  ind <- g[upper.tri(diag(n), diag = FALSE), ]
  
  distance <- function(i,k){
    return(wasserstein1d(x[i,], x[k,]))
  }
  
  dmatrix <- matrix(0,nrow = n, ncol = n)
  dmatrix[upper.tri(dmatrix,diag = FALSE)] <- mapply(distance, i=ind[,1], k=ind[,2])
  dmatrix <- dmatrix+t(dmatrix)
  
  return(dmatrix)
}


###################### compute spherical distance matrix of spherical data ###########
dist_sphere <- function(x){
  # Para: 
  #   x: n by d matrix, each row is iid observations on d-1 dimensional sphere 
  # Return:
  #   distance matrix
  x <- as.matrix(x)
  n <- dim(x)[1]
  g <- expand.grid(row = 1:n, col = 1:n)
  ind <- g[upper.tri(diag(n), diag = FALSE), ]
  
  distance <- function(i, k){
    return(acos(pmin(pmax(t(x[i,])%*%x[k,],-1.0), 1.0)))
  }
  
  dmatrix <- matrix(0, nrow = n, ncol = n)
  dmatrix[upper.tri(dmatrix, diag = FALSE)] <- mapply(distance, i=ind[,1], k=ind[,2])
  dmatrix <- dmatrix + t(dmatrix)
  
  return(dmatrix)
}


############# compute the Gaussian kernel Gram matrix for random element ###################### 
gram_matrix <- function(x, complexity, type=NULL, kernel='Gaussian'){
  # Para: 
  #   x: random objects:
  #        A n by m matrix for distributional objects; 
  #        A n by d by d array (a list of q by q matrices) for SPD matrices objects;
  #        A n by 1 matrix for spherical objects.
  #   complexity: tuning parameter in the reproducing kernel
  #   type: 'distribution', 'sphere', 'spd', 'euclidean'
  #   kernel: 'Gaussian' or 'Laplacian'
  # Return:
  #   n by n kernel Gram matrix
  x <- as.matrix(x)
  n <- dim(x)[1]
  # compute distance matrix
  if(type=='distribution'){
    k <- dist_dist(x)
  }
  else if(type=='sphere'){
    k <- dist_sphere(x)
  }
  else if(type=='spd'){
    k <- dist_spd(x, metric = 'Frobenius')
  }
  else if(type=='euclidean'){
    k <- as.matrix(dist(x))
  }
  
  if (kernel=='Gaussian'){
    sigma2 <- sum(k^2) / choose(n,2)
    gamma <- complexity / (sigma2)
    return(exp(-gamma*k^2))
  }
  else if (kernel=='Laplacian'){
    sigma <- sum(k) / choose(n,2)
    gamma <- complexity / (sigma)
    return(exp(-gamma*k)) 
  }
}

gram_matrix_new <- function(x, x.new, complexity, type=NULL, kernel='Gaussian'){
  # Para: 
  #   x: random objects:
  #        A n1 by m matrix for distributional objects; 
  #        A n1 by d by d array (a list of q by q matrices) for SPD matrices objects;
  #        A n1 by 1 matrix for spherical objects.
  #   x.new: random objects with same type of x and sample size n2
  #   complexity: tuning parameter in the reproducing kernel
  #   type: 'distribution', 'sphere', 'spd', 'euclidean'
  #   kernel: 'Gaussian' or 'Laplacian'
  # Return:
  #   n1 by n2 kernel Gram matrix
  x <- as.matrix(x)
  x.new <- as.matrix(x.new)
  n1 <- dim(x)[1]
  n2 <- dim(x.new)[1]
  # compute n1 by n2 distance matrix
  if(type=='distribution'){
    k <- outer( 
      1:n1, 1:n2,
      Vectorize(function(i,k) wasserstein1d(x[i,], x.new[k,]))
    )
  }
  else if(type=='sphere'){
    k <- outer( 
      1:n1, 1:n2,
      Vectorize(function(i,k) acos(pmin(pmax(t(x[i,])%*%x.new[k,],-1.0), 1.0)))
    )
  }
  else if(type=='spd'){
    k <- outer( 
      1:n1, 1:n2,
      Vectorize(function(i,k) frobenius(x[[i]]-x.new[[k]]))
    )
  }
  else if(type=='euclidean'){
    k <- outer( 
      1:n1, 1:n2,
      Vectorize(function(i,k) sqrt((x[i,] - x.new[k,])^2))
    )
  }
  
  if (kernel=='Gaussian'){
    sigma2 <- sum(k^2) / (n1*n2)
    gamma <- complexity / (sigma2)
    return(exp(-gamma*k^2))
  }
  else if (kernel=='Laplacian'){
    sigma <- sum(k) / (n1*n2)
    gamma <- complexity / (sigma)
    return(exp(-gamma*k)) 
  }
}

gram_wasserstein <- function(x, complexity){
  x <- as.matrix(x)
  n <- dim(x)[1]
  ##upper triangular index
  g <- expand.grid(row = 1:n, col = 1:n)
  ind <- g[upper.tri(diag(n), diag = FALSE), ]
  dist <- function(i,k){
    return(wasserstein1d(x[i,],x[k,]))
  }
  kupper <- mapply(dist, i=ind[,1], k=ind[,2])
  k <- matrix(0, nrow = n, ncol = n)
  k[upper.tri(k,diag = FALSE)] <- kupper^2
  k <- k+t(k)
  sigma2 <- sum(k)/choose(n,2)
  gamma <- complexity/(sigma2)
  return(exp(-gamma*k))
}
