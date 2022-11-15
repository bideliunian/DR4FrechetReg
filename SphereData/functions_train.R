main_sphere <- function(times, n, p, rho, model, data_type='sphere', kernel_type='Gaussian', 
                      complexity=1, methods=NULL){
  # Para:
  #   times: integer, repeat times of experiments
  #   n: sample size
  #   p: dimension of predictors
  #   rho: AR1(rho) for generating predictors, then transform to uniform (0,1)
  #   model: 'model1', 'model2' or 'model3'
  #   data_type: sphere
  #   kernel_type: 'Gaussian' or 'Laplacian'
  #   complexity: tuning parameter in reproducing kernel, default 1
  #   methods: SDR methods, by default c('fols','fphd','fiht','fsir','fsave','fdr','fopg', 'wire)
  #
  # Return:
  #   list of length methods, with each element the percentage of correct order determination
  #   and the estimation error  
  print(paste("------------------n=", n, "; p=", p, "; model=", model,"-----------------"))
  if (is.null(methods)) {
    methods = c('fols','fphd','fiht','fsir','fsave','fdr','fopg','wire')
  }
  result <- matrix(NA, nrow = times, ncol=length(methods))
  result_order <- matrix(NA, nrow = times, ncol=length(methods))
  for (i in 1:times) {
    # generate data
    data <- gendata_sphere(n=n, p=p, rho=rho, model=model)
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
