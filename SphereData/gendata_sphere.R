############# data generation for Spherical Data ##################

gendata_sphere <- function(n, p, rho=0, model){
  # Para:  
  #   p the dim of covariates
  #   n sample size
  #   rho: generate predictor x from AR1(rho) then transform to uniform distribution
  #   model: 'model1', 'model2' or 'model3'
  # Return: 
  #   A list containing 
  #   n: sample size; p: dimension of predictors;  
  #   d0: true order of dimension; b0: true sufficient directions;
  #   x: n by p design matrix; y: n by m matrix, m = 2 or 3
  data <- list()
  if(!(p>5)){
    stop("dimension p must greater or equal to 6")
  }
  if (!(model %in% c('model1', 'model2','model3'))) {
    stop("model must be one of 'model1' 'model2' 'model3 ")
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
  # generate predictors
  x <- pnorm(x)

  if(model=='model1'){
    delta <- rnorm(n, 0, 0.2)
    m_x <- list(cos(pi*x[,1]), sin(pi*x[,1]))
    eps <- list(-delta*sin(pi*x[,1]), delta*cos(pi*x[,1]))
    geny <- function(y_1, y_2){
      cos(abs(delta)) * y_1 + sin(abs(delta)) / abs(delta)*y_2
    }
    y <- mapply(geny, m_x, eps)
    data$x <- x
    data$y <- y
    data$d0 <- 1
    data$b0 <- c(1,rep(0,p-1))
  }
  
  else if(model=='model2'){
    # generate one sample
    genone <- function(x){
      delta <- rnorm(n = 2, 0, 0.2)
      m_x <- as.matrix(c(sqrt(1-x[2]^2) * cos(pi*x[1]), sqrt(1-x[2]^2) * sin(pi*x[1]),x[2]))
      ## identify two orthogonal basis of tangent space
      v_1 <- orth_mat(c(-sqrt(1-x[2]^2)*sin(pi*x[1]), sqrt(1-x[2]^2)*cos(pi*x[1]),0))## identify first basis 
      v_2 <- orth_mat(pracma::cross(v_1,m_x))## cross prod of v_1 and normal vector
      eps <- delta[1]*v_1+delta[2]*v_2## normal error on tangent plane
      y <- cos(sqrt(sum(eps^2)))*m_x+sin(sqrt(sum(eps^2)))/sqrt(sum(eps^2))*eps
      return(c(x,y))
    }
    data_full <- matrix(NA, nrow = n, ncol = p+3)
    for (i in 1:n) {
      data_full[i,] <- genone(x[i,])
    }
    data$x <- data_full[,1:p]
    data$y <- data_full[,(p+1):(p+3)]
    data$d0 <- 2
    data$b0 <- cbind(c(1,rep(0,p-1)),c(0,1,rep(0,p-2)))
  }
  
  else if(model=='model3'){
    eps_1 <- rnorm(n, 0, 0.2)
    eps_2 <- rnorm(n, 0, 0.2)
    m_x <- list(sin(pi*x[,1])*sin(pi*x[,2]),
                sin(pi*x[,1])*cos(pi*x[,2]),cos(pi*x[,1]))
    data$y <- cbind(sin(pi*x[,1]+eps_1)*sin(pi*x[,2]+eps_2),
                    sin(pi*x[,1]+eps_1)*cos(pi*x[,2]+eps_2),cos(pi*x[,1]+eps_1))
    data$x <- x
    data$d0 <- 2
    data$b0 <- cbind(c(1,rep(0,p-1)),c(0,1,rep(0,p-2)))
  }
  return(data)
}
