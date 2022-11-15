main_dist <- function(times, n, p, m, rho, non_ellip, model, data_type, kernel_type='Gaussian', 
                      complexity=1, methods=NULL){
  # Para:
  #   times: integer, repeat times of experiments
  #   n: sample size
  #   p: dimension of predictors
  #   m: number of observations for each distributions
  #   rho: AR1(rho) for generating predictors
  #   non_ellip: if 0, x from normal; if 1, x from a non-elliptical distribution by marginal transformation
  #     if 2, x from a non-elliptical distribution by joint transformation, if -1, x from a uniform distribution
  #   model: 'distex1', 'distex2', 'distex3', 'distex4'
  #   kernel_type: 'Gaussian' or 'Laplacian'
  #   complexity: tuning parameter in reproducing kernel, default 1
  #   methods: SDR methods, by default c('fols','fphd','fiht','fsir','fsave','fdr','fopg','wire')
  #
  # Return:
  #   list of length order_methods, with each element the percentage of correct order determination
  
  print(paste("-----------------------n=", n, "; p=", p, "; m=", m , "; model=", model,"----------------------"))
  if (is.null(methods)) {
    methods = c('fols','fphd','fiht','fsir','fsave','fdr','fopg','wire')
  }
  result <- matrix(NA, nrow = times, ncol=length(methods))
  result_order <- matrix(NA, nrow = times, ncol=length(methods))
  for (i in 1:times) {
    # generate data
    data <- gendata_dist(n=n, p=p, m=m, rho=rho, non_ellip=non_ellip, model=model)
    d0 <- data$d0
    b0 <- data$b0
    ygram <- gram_matrix(data$y, complexity = complexity, type=data_type, kernel=kernel_type)
    for (j in 1:length(methods)) {
      method <- methods[j]
      f <- get(method)
      print(paste("--------------------------", method, ":",  i, "th experiment----------------------------------"))
      tic()
      dhat <- pred_aug(x=data$x, y=ygram, s= 10, r =floor(p/5)+1, method=method)$rhat
      bhat <- f(x=data$x, y=ygram, d=dhat)$beta
      #bhat$fmave <- cal(x=data$x, ygram=ygram, d0=d0, method='fmave', initial=bhat$fopg)
      error <- eval(x=bhat, y=b0)
      result[i, j] <- error
      result_order[i, j] <- as.numeric(d0 == dhat)
      toc()
    }
  }
  return(list(est_error=result, order_acc=result_order))
}


# cal <- function(x, ygram, d0, method, initial=NULL){
#   print(method)
#   tic()
#   if(method=='fols'){
#     bhat <- fols(x=x, y=ygram, d=d0)$beta
#   }
#   else if(method=='fphd'){
#     bhat <- fphd(x=x, y=ygram, d=d0)$beta
#   }
#   else if(method=='fiht'){
#     bhat <- fiht(x=x, y=ygram, d=d0)$beta
#   }
#   else if(method=='fsir'){
#     bhat <- fsir(x=x, y=ygram, d=d0)$beta
#   }
#   else if(method=='fsave'){
#     bhat <- fsave(x=x, y=ygram, d=d0)$beta
#   }
#   else if(method=='fdr'){
#     bhat <- fdr(x=x, y=ygram, d=d0)$beta
#   }
#   else if(method=='wire'){
#     bhat <- wire(x=x, y=ygram, d=d0)
#   }
#   else if(method=='fopg'){
#     bhat <- fopg(x=x, y=ygram, d=d0)$beta
#   }
#   else if(method=='fmave'){
#     bhat <- fmave(x=x, y=ygram, d=d0, initial = initial)
#   }
#   toc()
#   return(bhat)
# }

