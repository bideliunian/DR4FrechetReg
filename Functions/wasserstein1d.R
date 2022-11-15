####### compute wasserstein distance between two 1 dimensional vector ##############

wasserstein1d <- function (a, b, p = 2) {
  # Para: a, b two one dimensional vector with same length; p is the index of Wasserstein p distance
  # Return: Wasserstein p distance
  m <- length(a)
  n <- length(b)
  stopifnot(m > 0 && n > 0 && m == n)
  return(mean(abs(sort(b) - sort(a))^p)^(1/p))
}