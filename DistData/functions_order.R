# main function
order_deter_dist <- function(times, n, p, m, rho, non_ellip, model, data_type, kernel_type='Gaussian', 
                        complexity=1, methods=NULL, order_methods=NULL){
  # Para:
  #   times: integer, repeat times of experiments
  #   n: sample size
  #   p: dimension of predictors
  #   m: number of observations for each distributions
  #   rho: AR1(rho) for generating predictors
  #   non_ellip: if 0, x from normal; if 1, x from a non-elliptical distribution by marginal transformation
  #     if 2, x from a non-elliptical distribution by joint transformation, if -1, x from a uniform distribution
  #   model: 'distex1', 'distex2', 'distex3', 'distex4'
  #   data_type: 'distribution
  #   kernel_type: 'Gaussian' or 'Laplacian'
  #   complexity: tuning parameter in reproducing kernel, default 1
  #   methods: SDR methods, by default c('fols','fphd','fiht','fsir','fsave','fdr','fopg','wire')
  #   order_methods: 'bic' or 'pa'(predictor augmentation)
  #
  # Return:
  #   list of length order_methods, with each element the percentage of correct order determination
  
  if (is.null(methods)) {
    methods = c('fols','fphd','fiht','fsir','fsave','fdr','fopg','wire')
  }
  if (is.null(order_methods)) {
    order_methods = c("bic", "pa")
  }
  result_pa = result_bic = matrix(data = NA, nrow = times, ncol = length(methods))
  for (i in 1:times) {
    cat(sprintf("-------------%i th experiment---------------", i))
    ### generate data
    data <- gendata_dist(n=n, p=p, m=m, rho=rho, non_ellip=non_ellip, model = model)
    d0 <- data$d0
    b0 <- data$b0
    ygram <- gram_matrix(data$y, complexity = complexity, type=data_type, kernel=kernel_type)
    for(j in 1:length(methods)){
      result_bic[i,j] = bic(x=data$x, y=ygram, criterion='lal', method=methods[j])$rhat
      result_pa[i,j] = pred_aug(x=data$x, y=ygram, s=10, r=floor(p/5)+1, method=methods[j])$rhat
    }
  }
  colnames(result_bic) = colnames(result_pa) = methods
  output <- list()
  result_list <- list(result_bic, result_pa)
  for (i in 1:length(result_list)){
    output[[i]] <- colSums(result_list[[i]] == d0) / times
  }
  names(output) <- order_methods
  
  print(paste("---n=", n, ", p=", p, ", rho=", rho, ", non_ellip=", non_ellip, ", model=", model,"---"))
  print(output)
  
  return(output)
}
