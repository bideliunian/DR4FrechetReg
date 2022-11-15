###########################################################
############ read the result and box plot #################
###########################################################

save_path <- "D:/Research/DR4FR/Codes/SphereData/OrderResults"
# if use aci
#save_path <- "~/work/DR4FR/SphereData/OrderResults"

############################## Model 1&2##############################
# set global parameters
times_repeat <- 100
n_list <- c(200, 400)
p_list <- c(10, 20)
rho <- 0.5
Models <- c('model1','model2','model3')
Methods <-c('fols','fphd','fiht','fsir','fsave','fdr','fopg','wire') #'fmave'
kernel_type <- "Laplacian"
kernel_para <- 0.01
data_type <- "sphere"

# record parameters in a grid
grid_1 <- expand.grid(n = n_list[1], p=p_list[1], rho=rho,  model = Models, 
                      stringsAsFactors=FALSE)
grid_2 <- expand.grid(n = n_list[2], p=p_list[2], rho=rho, model = Models, 
                      stringsAsFactors=FALSE)
grid <- rbind(grid_1, grid_2)

order_result <- matrix(0, nrow=nrow(grid), ncol=2*length(Methods))
for (i in 1:times_repeat) {
  load(file =  paste(save_path, "/order.rho=", rho, 
                     ".result",i,".RData",sep=""))
  order_result <- order_result + order_result_table[,5:ncol(order_result_table)] / times_repeat
}

result <- cbind(grid, order_result)

# write in latex
library(xtable)
result <- result[order(result$model),]
result_prop_bic <- round(as.matrix(result[,5:12]), 2)
colnames(result_prop_bic) <- Methods
result_prop_pa <- round(as.matrix(result[,13:20]), 2)
colnames(result_prop_pa) <- Methods
result_all <- matrix(0, nrow = 2*nrow(result), ncol = ncol(result_prop_bic)+2)
for (i in 1:nrow(result)) {
  for (j in 1:ncol(result_all)) {
    if(j==1){
      result_all[(2*(i-1)+2),j] <- 0
      if(i%%2 == 1){
        result_all[(2*(i-1)+1),j] <- paste("(", 10,",",200,")", sep = "")
      }
      else{
        result_all[(2*(i-1)+1),j] <- paste("(", 20,",",400,")", sep = "")
      }
    }
    else if(j==2){
      result_all[(2*(i-1)+1),j] <- "BIC"
      result_all[(2*(i-1)+2),j] <- "PA"
    }
    else{
      result_all[(2*(i-1)+1),j] <- paste(result_prop_bic[i,j-2] * 100, "%", sep = "")
      result_all[(2*(i-1)+2),j] <- paste(result_prop_pa[i,j-2] * 100, "%", sep = "")
    }
  }
}
names(result_all) <- Methods
xtable(result_all)