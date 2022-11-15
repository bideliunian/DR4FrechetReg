main_cov <- function(times, n, p, rho, non_ellip, model, data_type, error_type='log_normal', kernel_type='Gaussian', 
                      complexity=1, methods=NULL){
  # Para:
  #   times: integer, repeat times of experiments
  #   n: sample size
  #   p: dimension of predictors
  #   rho: AR1(rho) for generating predictors
  #   non_ellip: if 0, x from normal; if 1, x from a non-elliptical distribution by marginal transformation
  #     if 2, x from a non-elliptical distribution by joint transformation, if -1, x from a uniform distribution
  #   model: 'model1', 'model2'
  #   error_type: 'log_normal' or 'Rlog_normal'
  #   kernel_type: 'Gaussian' or 'Laplacian'
  #   complexity: tuning parameter in reproducing kernel, default 1
  #   methods: SDR methods, by default c('fols','fphd','fiht','fsir','fsave','fdr','fopg')
  #
  # Return:
  #   list of length order_methods, with each element the percentage of correct order determination
  
  print(paste("------------------n=", n, "; p=", p, "; error=", error_type , "; model=", model,"-----------------"))
  if (is.null(methods)) {
    methods = c('fols','fphd','fiht','fsir','fsave','fdr','fopg','wire')
  }
  result <- matrix(NA, nrow = times, ncol=length(methods))
  result_order <- matrix(NA, nrow = times, ncol=length(methods))
  for (i in 1:times) {
    # generate data
    data <- gendata_cov(n=n, p=p, rho=rho, non_ellip=non_ellip, model=model, error=error_type)
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

