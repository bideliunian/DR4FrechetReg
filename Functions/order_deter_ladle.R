#############################################################################
#################  Ladle estimator for order determination #############


ladle <- function(x, y, ymat, nboot, method){
  # Para:
  #   x: n by p matrix
  #   y: Response objects, 
  #     A n by m matrix for distributional objects, or 
  #     A n by d by d array (a list of q by q matrices) for SPD matrices objects, or
  #     A n by 1 matrix for spherical objects.
  #   ymat: n by n gram matrix.
  #   nboot: number of times for bootstrap
  #   method: "fols", "fopg", "fiht", "fsir", "fsave", "fdr", "fopg"
  x <- as.matrix(x)
  n <- dim(x)[1]
  p <- dim(x)[2]
  
  # kmax is the upper bound of candidate orders
  if(p <= 10) {
    kmax <- p-2
  }
  else {
    kmax <- round(p/log(p)) 
  }

  # compute phi_n
  phi <- function(kmax, eigen_val){
    eigen_all <- 1 + sum(eigen_val[1:(kmax+1)]) # adding 1 for robustness;
    # warning: adding 1 may affect the relative magnitude between phi and f parts, need more careful design
    return(eigen_val[1:(kmax+1)]/eigen_all)
  }
  
  f <- get(method)
  candidate_matrix <- f(x=x, y=ymat, d=1)$candidate.matrix # candidate matrix
  eigen_val_old <- eigen(candidate_matrix)$values  # eigenvalues
  eigen_vec_old <- eigen(candidate_matrix)$vectors  # eigenvectors
  
  pn <- phi(kmax, eigen_val_full)
  
  f <- function(kmax, evec1, evec2){
    out <- numeric()
    for(k in 0:kmax){
      if(k == 0) out <- c(out, 0)
      if(k == 1) out <- c(out, 1 - abs(t(evec1[,1])%*%evec2[,1]))
      if(k!=0 & k!=1) out <- c(out, 1 - abs(det(t(evec1[,1:k])%*%evec2[,1:k])))
    }
    return(out)
  }
  
  fn0 <- 0
  for(iboot in 1:nboot){
    bootindex <- round(runif(n, min=-0.5, max=n+0.5))
    xs <- x[bootindex,]
    ys <- y[bootindex]
    ysmat <- gramw(ys, complexity = 1)
    
    mat <- f(xs, ysmat, d=1)$candidate.matrix 
    eigen_val <- eigen(mat)$values
    eigen_vec <- eigen(mat)$vectors
    fn0 <- fn0 + f(kmax, eigen_vec_old, eigen_vec) / nboot
  }
  fn <- fn0 / (1 + sum(fn0)) # again, here adding 1 for robustness
  gn <- pn + fn
  rhat <- which.min(gn) - 1
  
  return(list(rset=(0:kmax), gn=gn, rhat=rhat))
}
