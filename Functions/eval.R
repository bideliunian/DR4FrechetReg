########### evaluation the estimation of subspace ###############
# For two matrices B1 and B2, we compute
# || P_B1 - P_B2 || = ||B1(t(B1)B1)^-1B1 - B2(t(B2)B2)^-1B2 ||

eval <- function(x, y){
  # Para: x, y two matrix with same dimension
  # Return: Frobenius norm between the projection matrix
  x <- as.matrix(x)
  y <- as.matrix(y)
  error <- proj_matrix(x)-proj_matrix(y)
  return(norm(error, type = 'F'))
}


