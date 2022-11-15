############# data generation ##################
########### Preparation ###################

### truncated gamma distribution
rgammat <- function(n, range=c(0.1,10), shape, scale = 1) {
  # Para: 
  #   n sample size
  #   range: truncated range
  #   shape: shape parameter >0
  #   scale: scale parameter >0
  F.a <- pgamma(min(range), shape = shape, scale = scale)
  F.b <- pgamma(max(range), shape = shape, scale = scale)
  
  u <- runif(n, min = F.a, max = F.b)
  qgamma(u, shape = shape, scale = scale)
}

sigmoid <- function(x){
  return(1/(1+exp(-x)))
}


gendata_dist <- function(n, p, m, rho=0, non_ellip=0, model=NULL){
  # Para: 
  #   n: sample size;
  #   p: dimension of data
  #   m: number of discrete observations from a y
  #   rho: generate predictor x from AR1(rho); if rho==-1, generate x from uniform distribution
  #   non_ellip: if 0, x from normal; if 1, x from a non-elliptical distribution by marginal transformation
  #   if 2, x from a non-elliptical distribution by joint transformation, if -1, x from a uniform distribution
  # Return: 
  #   A list containing 
  #   n: sample size; p: dimension of predictors; m: number of discrete observations from a y; 
  #   d0: true order of dimension; b0: true sufficient directions
  #   x: n*p design matrix; y: n*m matrix for distributional data
  data <- list()
  
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
  data$x <- x
  
  # coefficients
  beta_1 <- c(1, 1, rep(0,p-2))
  beta_2 <- c(rep(0,p-2), 1, 1)
  beta_3 <- c(1, 2, rep(0,p-3), 2)
  beta_4 <- c(0, 0, 1, 2, 2, rep(0,p-5))
  
  # Model1: mu = N(exp(beta1*x), 1); sigma = 1
  if(model=='distex1'){
    d1 <- rowSums(x[,1:2])
    mu <- rnorm(n, mean = exp(d1), sd = 1)
    sigma <- 1
    data$y <- t(mapply(rnorm, n=m, mu, sigma))
    data$d0 <- 1
    data$b0 <- beta_1
  }
  # Model2: mu = N(exp(beta1*x), 1); sigma = exp(beta2*x)
  else if(model=='distex2'){
    d1 <- rowSums(x[,1:2])
    d2 <- (x[,p-1]+x[,p])
    mu <- rnorm(n, mean = exp(d1), sd = 1)
    sigma <- exp(d2)
    sigma[sigma > 10] <- 10
    sigma[sigma < 0.1] <- 0.1
    data$y <- t(mapply(rnorm, n=m, mu, sigma))
    data$d0 <- 2
    data$b0 <- cbind(beta_1, beta_2)
  }
  # Model 3: mu = N((beta3*x)^2, 0.5^2), 
  # sigma = tGamma(exp(2+beta4*x)/nu, nu/sqrt(exp(2+beta4*x))), with range (0.1, 10)
  else if(model=='distex3'){
    d1 <- (x[,1] + 2*x[,2] + 2*x[,p])
    d2 <- 2*(x[,3] + 2*x[,4] + 2*x[,5]) + 2
    mu <- rnorm(n, mean = 3*d1, sd = 0.5)
    sigma <- rgammat(n, shape = (d2)^2/0.5, scale = 0.5/abs(d2))
    data$y <- t(mapply(rnorm, n=m, mu, sigma))
    data$d0 <- 2
    data$b0 <- cbind(beta_3, beta_4)
  }
  # Model 4: mu = N(3*sin(beta3*x), 0.5^2), 
  # sigma = tGamma((2+beta4*x)^2/nu, nu/(2+beta4*x)), with range (0.1, 10)
  else if(model=='distex4'){
    d1 <- (x[,1] + 2*x[,2] + 2*x[,p])
    d2 <- 2*(x[,3] + 2*x[,4] + 2*x[,5]) + 2
    mu <- rnorm(n, mean = 3*sin(d1), sd = 0.5)
    sigma <- rgammat(n,shape = (d2)^2/0.5, scale = 0.5/abs(d2))
    data$y <- t(mapply(rnorm, n=m, mu, sigma))
    data$d0 <- 2
    data$b0 <- cbind(beta_3, beta_4)
  }
  else if(model=='distex5'){
    trans <- function(x) {
      k <- sample(c(-2,-1,1,2),1)
      x-sin(k*x)/abs(k)
    }
    d1 <- rowSums(x[,1:2])
    d2 <- abs(2 + x[,p-1] + x[,p])
    mu <- rnorm(n, mean = d1, sd = 1)
    sigma <- rgammat(n, shape = (d2)^2/0.5, scale = 0.5/abs(d2))
    data$y <- t(apply(as.matrix(mvrnorm(m, mu, diag(sigma))), MARGIN=c(1,2), trans))
    data$d0 <- 2
    data$b0 <- cbind(beta_1, beta_2)
  }
  
  data$n <- n
  data$p <- p
  data$m <- m
  return(data)
}

#################################################
### data visulization
##########################################
# data <- gendata_dist(100, 10, m=100, rho=0.5, non_ellip=2, mode = 'distex1')
# dens <- apply(data$y, 1, density)
# plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")))
# mapply(lines, dens, col=1:length(dens))

# view the shape of gamma distribution
# k <- c(1, 4, 9, 100)
# theta <- 1/sqrt(k)
# plot(0, 0, xlim = c(0, 10), ylim = c(0, 1), type = "n")
# for(i in seq_along(k)) {
#   curve(dgamma(x, shape = k[i], scale = theta[i]), from = 0, to = 10, col = i, add = TRUE) 
# }
