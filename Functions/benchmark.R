#### benchmark ###############################
# The benchmark is the expectation of distance between the true subspace
# and a (Gaussian) randomly generated subspace of the same dimension
# Note: this benchmark may change with the dimension of true subspace (p and d)

bchmk <- function(x, N){
  # Para: 
  #   x: the true dimension reduction direction with dim p by d
  #   N: the number of MC simulations
  # Return: 
  x <- as.matrix(x)
  p <- dim(x)[1]
  d <- dim(x)[2]
  
  benchmark <- numeric()
  for (i in 1:N) {
    random_proj <- matrix(rnorm(p*d), nrow=p, ncol=d) 
    benchmark[i] <- eval(x, random_proj)
  }
  return(list(mean=mean(benchmark), sd=sd(benchmark)))
}
