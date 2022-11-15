############### matrix power #####################
# matpower <- function(a, alpha){
#   # Para:
#   #   a: A symmetric matrix
#   #   alpha: index of power
#   # Return:
#   #   a^alpha
#   a <- (a + t(a)) / 2
#   tmp <- eigen(a)
#   return(tmp$vectors%*%diag((tmp$values)^alpha)%*%t(tmp$vectors))
# }

# ######## matrix power
# "%^%" <- function(x, n) 
#   with(eigen(x), vectors %*% (values^n * t(vectors)))

################ standadize matrix ##############
standmatrix <- function(x){
  mu <- apply(x,2,mean)
  sig <- var(x)
  signrt <- matpower(sig,-1/2)
  return(t(t(x)-mu)%*%signrt)
}

############### discretize response y ################
discretize <- function(y, h){
  # Para: 
  #     y: response variable
  #     h: number of slices
  # Return:
  #     discretized y
  if(is.vector(y)){
    n <- length(y) 
  }
  else{
    n <- nrow(y)
  }
  m <- floor(n/h)
  y <- y + .00001 * mean(y) * rnorm(n)
  y_sort = sort(y)
  divpt <- numeric()
  for(i in 1:(h-1)) {
    divpt  <-  c(divpt,y_sort[i*m+1]) 
  }
  y_discrete <- rep(0,n)
  y_discrete[y<divpt[1]] <- 1
  y_discrete[y>=divpt[h-1]] <- h
  for(i in 2:(h-1)) {
    y_discrete[(y>=divpt[i-1])&(y<divpt[i])]=i
  }
  return(y_discrete)
}

############## matrix exponential/logrithm/power #############
mexp <- function(x)  {
  x = round((x + t(x))/2,7)
  return(with(eigen(x), vectors %*% (exp(values) * t(vectors))))
}

mlog <- function(x)  {
  x = round((x + t(x))/2,7)
  return(with(eigen(x), vectors %*% (log(values) * t(vectors))))
}

mpower <- function(x, alpha)  {
  x = round((x + t(x))/2,7)
  return(with(eigen(x), vectors %*% ((values)^alpha * t(vectors))))
}


############# frobenius norm ###########################
frobenius <- function(X) {return(sqrt(sum(diag(t(X)%*%X))))}
tradeindex12 = function(k, n){
  j = ceiling(k/n)
  i = k - (j-1)*n
  return(c(i,j))
}


######## get mode ##################
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

####### projection matrix ###############
proj_matrix <- function(x) {
  x <- as.matrix(x)
  return(tcrossprod(x%*%mpower(crossprod(x),-1/2)))
}


###### orthogonal directions ###############
orth_mat <- function(x) {
  x <- as.matrix(x)
  return(x%*%mpower(crossprod(x),-1/2))
}

###### mp power of matrix ###############
mppower <- function(matrix, power, ignore){
  eig <- eigen(matrix)
  eval <- eig$values
  evec <- eig$vectors
  m <- length(eval[(eval)>ignore])
  tmp <- evec[,1:m]%*%diag(eval[1:m]^power)%*%t(evec[,1:m])
  return(tmp)
}
