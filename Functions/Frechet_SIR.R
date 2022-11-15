############################################################################
## Slice Inverse Regression and Frechet Slice Inverse Regression##################

# Function to implement SIR for SDR
sir <- function(x, y, h=NULL){
  # Para:
  #   x: A n by p matrix, standardized predictors
  #   y: A n by 1 matrix
  #   h: number of slices
  # Return:
  #   SIR candidate matrix p by p matrix
  x <- as.matrix(x)
  n <- dim(x)[1]
  p <- dim(x)[2]
  if (is.null(h)){
    h <- max(round(n/(2*p)), 8) 
  } # take the number of slices as larger one between [n/2p] and 8

  y_dis <- discretize(y, h)
  y_label <- unique(y_dis)
  prob <- numeric()
  exy <- numeric()
  for(i in 1:h) {
    prob <- c(prob,length(y_dis[y_dis==y_label[i]]) / n)
    exy <- rbind(exy, apply(x[y_dis==y_label[i],], 2, mean)) 
  }
  mat <- t(exy) %*% diag(prob) %*% exy
  return(mat)
}


# Frechet SIR
fsir <- function(x, y, d, isGram=TRUE, rho=1){
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
    y <- gramw(y, complexity=rho)
  }
  y <- split(y, col(y)) 
  beta <- lapply(y, sir, x=z) 
  lambda <- Reduce('+',beta)/n
  B <- eigen(lambda)$vectors[,1:d]
  beta.final <- signrt%*%B
  return(list(beta=beta.final, candidate.matrix=lambda))
}


# Function to implement SIR for SDR
# sir.old <- function(x, y){
#   # Para:
#   #   x: A n by p matrix
#   #   y: A n by 1 matrix
#   # Return:
#   #   SIR candidate matrix p by p matrix
#   x <- as.matrix(x)
#   n <- dim(x)[1]
#   p <- dim(x)[2]
#   h <- round(n / (2*p)) # take the number of slices as [n/2p]
#   
#   y_dis <- discretize(y, h)
#   y_less <- y_dis
#   y_label <- numeric()
#   for(i in 1:n) {
#     if(var(y_less)!=0) {
#       y_label <- c(y_label, y_less[1])
#       y_less <- y_less[y_less!=y_less[1]]
#     }
#   }
#   y_label <- c(y_label, y_less[1])
#   prob <- numeric()
#   exy <- numeric()
#   for(i in 1:h) {
#     prob <- c(prob,length(y_dis[y_dis==y_label[i]]) / n)
#     exy <- rbind(exy, apply(x[y_dis==y_label[i],], 2, mean)) 
#   }
#   mat <- t(exy) %*% diag(prob) %*% exy
#   return(mat)
# }
