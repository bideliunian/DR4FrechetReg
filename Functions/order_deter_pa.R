#############################################################################
#################  Predictor augmentation for order determination #############

## compute candidate matrix
# candmat <- function(x, y, method, d0=2){
#   #   x: n by p matrix
#   #   y: n by n gram matrix
#   #   d0: order of reduced dimension, by defult 2 
#   #   method: "fols", "fopg", "fiht", "fsir", "fsave", "fdr", "fopg"
#   if(method == "fols") mat <- fols(x=x, y=y, d=d0)$candidate.matrix
#   if(method == "fphd") mat <- fphd(x=x, y=y, d=d0)$candidate.matrix
#   if(method == "fiht") mat <- fiht(x=x, y=y, d=d0)$candidate.matrix
#   if(method == "fsir") mat <- fsir(x=x, y=y, d=d0)$candidate.matrix
#   if(method == "fsave") mat <- fsave(x=x, y=y, d=d0)$candidate.matrix
#   if(method == "fdr") mat <- fdr(x=x, y=y, d=d0)$candidate.matrix
#   if(method == "fdr") mat <- fdr(x=x, y=y, d=d0)$candidate.matrix
#   if(method == "fopg") mat <- fopg(x=x, y=y, nit=1, d=d0)$candidate.matrix
#   
#   return(mat)
# }

# predictor augmentation 
pred_aug <- function(x, y, s, r, method){
  # Para:
  #   x: n by p matrix
  #   y: n by n gram matrix
  #   s: number of times of augmentation
  #   r: size of augmentation
  #   method: "fols", "fopg", "fiht", "fsir", "fsave", "fdr", "fopg", "wire"
  #   
  # Return:
  #   A list containing possible orders, objective function values, 
  #   and order selected by predictor augmentation
  x <- as.matrix(x)
  n <- dim(x)[1]
  p <- dim(x)[2]

  x_aug <- matrix(rnorm(n * r), n, r) # augmented predictors N(0, I)
  x_new <- cbind(x, x_aug)
  f <- get(method)
  candidate_matrix <- f(x=x_new, y=y, d=1)$candidate.matrix # candidate matrix
  eigen_val_full <- eigen(candidate_matrix)$values  # eigenvalues
  eigen_vec_full <- eigen(candidate_matrix)$vectors  # eigenvectors
  
  ## compute phi_n
  phi = function(eigen_val, p){
    eigen_all = 0.00001 + cumsum(eigen_val[1:(p+1)])
    return(eigen_val[1:(p+1)]/eigen_all)
  }
  phi_n = phi(eigen_val_full, p)
  
  f_n <- 0
  for(i in 1:s){
    # predictor augmentation
    x_add <- matrix(rnorm(n * r), n, r)
    x_all <- cbind(x, x_add)
    
    mat <- f(x=x_all, y=y, d=1)$candidate.matrix
    sub_evec <- eigen(mat)$vectors[(p+1):(p+r),1:p]
    sub_evec_norm <- apply(sub_evec, MARGIN=2, FUN = function(x) sum(t(x)%*%x))
    fn0 <- c(0,cumsum(sub_evec_norm))
    f_n <- f_n + fn0/s
  }
  
  g_n = phi_n + f_n
  rhat = which.min(g_n) - 1
  return(list(rset=(0:p), object_values=g_n, rhat=rhat))
}

# if (rho != 0) {
#   x_aug <- matrix(NA, n, r)
#   x_aug[, 1] <- rnorm(n)
#   for (i in 2:r) {
#     x_aug[, i] = rho * x_aug[, i - 1] + sqrt(1 - rho^2) * rnorm(n)
#   }
# }
# else {
#   x_aug <- matrix(rnorm(n * r), n, r)
# }
