################### K-fold Cross Validation #####################
##################################################################

kfoldCV <- function(n, k){
  # Paras:
  #   x is random objects:
  #     if type == 'distribution',  a n by m matrix for distributional objects; 
  #     if type == 'spd', a n by d by d array (a list of q by q matrices) for SPD matrices objects;
  #     if type == 'sphere', a n by 1 matrix for spherical objects.
  #   k is the number of folds.
  # Returns:
  #   A list with CV splited random objects. 
  
  if ((n %% k) != 0){
    stop("The number of samples cannot be divied by the number of folds")
  } 
  else {
    n_divided <- n / k
  }
  
  index_list <- list()
  index_all <- c(1:n)
  
  for (i in 1:k){
    index_selected <- sample(1:((k - i + 1) * n_divided), n_divided, replace=FALSE)
    index_list[[i]] <- index_selected
  }
  
  return(index_list)
}

################## cv to choose bindwidth and kernel type ##################
cv4kernel <- function(x, y, index, type, complexity, kernel_type, method='fols', d=5) {
  # Para:
  #   x: A n by p design matrix
  #   y: 
  #     if type == 'distribution', a n by m matrix for distributional objects; 
  #     if type == 'spd', a n by d by d array (a list of q by q matrices) for SPD matrices objects;
  #     if type == 'sphere', a n by 1 matrix for spherical objects.
  #   index: index of the fold
  #   type: type of response y
  #   complexity: tuning parameter in the rk
  #   kernel_type: 'Gaussian' or 'Laplacian'
  #   
  # Return:
  #   
  n <- dim(x)[1]
  n_test <- length(index)
  n_train <- n - n_test
  x_train <- x[-index,]
  x_test <- x[index,]
  train_index <- setdiff(c(1:n), index)
  ygram <- gram_matrix(y, complexity=complexity, type=type, kernel=kernel_type)
  if (type == 'spd'){
    d <- dim(y[[1]])[1]
    y_train <- y[-index]
    y_test <- y[index]
  }
  else {
    y_train <- y[-index,]
    y_test <- y[index,]
  }
  ygram_train <- ygram[-index, -index]
  ygram_test <- ygram[index, index]
  
  # sdr on training set
  f <- get(method)
  bhat <- f(x=x_train, y=ygram_train, d=d)$beta
  x_suff_test <- x_test %*% bhat # sufficient predictors on testing set
  x_suff_train <- x_train %*% bhat # sufficient predictors on training set
  
  if (type == 'distribution'){
    qSup <- seq(0,1,0.02)
    qpred <- GloDenReg(xin=x_suff_train, yin=y_train, xout=x_suff_test, 
                       optns = list(qSup = qSup))$qout
    qtest <- t(apply(X=y_test, MARGIN=1, FUN=quantile, probs=qSup))
    pred_error <- sum((qpred - qtest)^2) / n_test  
  }
  
  else if (type == 'spd'){
    y_pred = GloCovReg(x=x_suff_train, M=y_train, xout=x_suff_test, 
                       optns=list(corrOut=FALSE,metric="frobenius"))$Mout
    pred_error <- sum((unlist(y_pred) - unlist(y_test))^2) / n_test
  }
  
  else if (type == 'sphere'){
    y_pred <- GloSpheReg(xin=x_suff_train, yin=y_train, xout=x_suff_test)$yout
    pred_error <- sum((y_pred - y_test)^2) / n_test
  }
  
  return(pred_error)
}