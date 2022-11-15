############# trace correlation between two matrices ####################
tracecorr <- function(x, y, d){
  x <- as.matrix(x)
  y <- as.matrix(y)
  px <- x%*%solve(t(x)%*%x)%*%t(x)
  py <- y%*%solve(t(y)%*%y)%*%t(y)
  return(sum(as.vector(px)*as.vector(py))/d)
}