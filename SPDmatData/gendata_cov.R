############# data generation for SPD matrix ##################

##### generate symmetric matrix variate normal distribution ##########
smvnormal <- function(d = 3, M = NULL, Sigma=NULL){
  # Para:
  #   d: dimension of spd matrix
  #   M: d by d mean matrix
  #   Sigma: Covariance matrix
  # Return:
  #   a d by d spd matrix, follows a normal distribution
  if(is.null(M)) M=matrix(0,d,d)
  if(is.null(Sigma)) Sigma=diag(d)
  if(any(abs(M-t(M)) > 1e-8)){
    stop("Mean covariate must be symmetric")
  }
  if(any(eigen(Sigma)$values<0)){
    stop("Cov must be symmetric positive definite")
  }
  x <- matrix(0, nrow = d, ncol = d)
  diag(x) <- rnorm(d)
  x[lower.tri(x)] <- rnorm(d*(d-1)/2, mean = 0, sd = 1/sqrt(2))
  x <- x + t(x) - diag(diag(x))
  G <- mpower(Sigma, 1/2)
  x <- t(G)%*%x%*%G+M
  return(x)
}

####################### generate data ###############################
gendata_cov <- function (n, p, rho=0, non_ellip=0, model, error) {
  # Para:  
  #   p the dim of covariates
  #   n sample size
  #   rho: generate predictor x from AR1(rho); if rho==-1, generate x from uniform distribution
  #   non_ellip: if 0, x from normal; if 1, x from a non-elliptical distribution by marginal transformation
  #   if 2, x from a non-elliptical distribution by joint transformation
  #   model: 'model1' or 'model2'
  #   error: 'Rlog_normal' or 'log_normal
  # Return: 
  #   A list containing 
  #   n: sample size; p: dimension of predictors;  
  #   d0: true order of dimension; b0: true sufficient directions;
  #   x: n*p design matrix; y: list with length n, with each element a d by d spd matrix 
  if(!(p>5)){
    stop("dimension p must greater or equal to 6")
  }
  if (!(error %in% c('Rlog_normal', 'log_normal'))) {
    stop("error must be one of 'Rlog_normal' 'log_normal'")
  }
  if (!(model %in% c('model1', 'model2','model3','model4'))) {
    stop("model must be one of 'model1' 'model2' ")
  }
  
  # generate AR1(rho) predictors
  if (rho > 0) {
    x <- matrix(NA, n, p)
    x[, 1] <- rnorm(n)
    for (i in 2:p) {
      x[, i] = rho * x[, i - 1] + sqrt(1 - rho^2) * rnorm(n)
    }
  }
  else {
    x <- matrix(rnorm(n * p), n, p)
  }
  # generate non-elliptical predictors if non_ellip != 0
  if (non_ellip == 1){
    x[,1] <- sin(x[,1])
    x[,2] <- abs(x[,2])
  }
  else if (non_ellip == 2){
    x[,1] <- abs(x[,3] + x[,4]) + rnorm(n=p, mean = 0,sd = abs(x[,3]))
    x[,2] <- sqrt(abs(x[,3] + x[,4])) + rnorm(n=p, mean=0, sd=abs(x[,4]))
  }
  else if (non_ellip == -1){
    x <- pnorm(x)
  }
  
  if (model=='model1') {
    d <- 2
    dFun <- function(t1) {
      dx <- diag(2)
      dx[1,2] = dx[2,1] <- (exp(t1)-1)/(exp(t1)+1)
      return(dx)
    }
    beta_1 <- c(1,1,rep(0,p-2))
    z1 <- x%*%beta_1
    dx <- sapply(X = z1, FUN = dFun, simplify = FALSE)
    beta <- beta_1
    d0 <- 1
    eps <- 0.3
  }
  else if (model=='model2') {
    d <- 3
    dFun <- function(t1,t2) {
      dx <- matrix(0,3,3)
      diag(dx) <- c(1, 1, 1)
      dx[2,3] = dx[3,2] = dx[1,2] = dx[2,1] <- 0.4*(exp(t1)-1)/(exp(t1)+1)
      dx[1,3] = dx[3,1] <- 0.4*sin(t2)
      #=0.4*(t3^2-1)/(t3^2+1)
      return(dx)
    }
    beta_1 <- c(1, 1, rep(0,p-2))
    beta_2 <- c(rep(0,p-2), 1, 1)
    beta <- cbind(beta_1, beta_2)
    z1 <- x %*% beta_1
    z2 <- x %*% beta_2
    dx <- mapply(FUN = dFun, z1, z2, SIMPLIFY = FALSE)
    d0 <- 2
    eps <- 0.1
  }
  if (error=='Rlog_normal'){
    geny <- function(dx){
      dxrt <- mpower(dx, 1/2)
      y <- t(dxrt) %*% mexp(smvnormal(d, Sigma = eps*diag(d))) %*% dxrt
      y <- round((y + t(y)) / 2, 7)
    }
    y <- lapply(X = dx, FUN = geny)
  }
  else if (error=='log_normal'){
    geny <- function(dx){
      logy <- smvnormal(d = d, M = mlog(dx), Sigma = eps*diag(d))
      y <- mexp(logy)
      y <- round((y + t(y))/2, 7)
    }
    y <- lapply(X = dx, FUN = geny)
  }
  data <- list()
  data$x <- x
  data$y <- y
  data$b0 <- beta
  data$d0 <- d0
  return(data)
}
#data = gendata_cov(n=100, p=10, model = 'model2', error = 'log_normal')
