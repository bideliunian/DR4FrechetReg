###########################################################
############ read the result and box plot #################
###########################################################

save_path <- "D:/Research/DR4FR/Codes/SPDmatData/OrderResults"
# if use aci
#save_path <- "~/work/DR4FR/SPDmatData/OrderResults"

############################## Model 1&2##############################
# set global parameters
times_repeat <- 100
n_list <- c(200, 400)
p_list <- c(10, 20)
error_type <- "log_normal"
non_elliptical <- 1
rho <- 0.5
Models <- c('model1','model2')
Methods <-c('fols','fphd','fiht','fsir','fsave','fdr','fopg','wire') #'fmave'
kernel_type <- "Gaussian"
kernel_para <- 1
data_type <- "spd"

# record parameters in a grid
# record parameters in a grid
grid_1 <- expand.grid(n = n_list[1], p=p_list[1], error=error_type, rho=rho, 
                      non_ellip=non_elliptical, model = Models, 
                      stringsAsFactors=FALSE)
grid_2 <- expand.grid(n = n_list[2], p=p_list[2], error=error_type, rho=rho, 
                      non_ellip=non_elliptical, model = Models, 
                      stringsAsFactors=FALSE)
grid_12 <- rbind(grid_1,grid_2)



order_result_12 <- matrix(0, nrow=nrow(grid_12), ncol=2*length(Methods))
for (i in 1:times_repeat) {
  load(file =  paste(save_path, "/order.error=", error_type,".rho=", rho, "nonellip=", non_elliptical,
                     ".result",i,".RData",sep=""))
  order_result_12 <- order_result_12 + order_result_table[,7:ncol(order_result_table)] / times_repeat
}

result_12 <- cbind(grid_12, order_result_12)

# write in latex
library(xtable)
result_12 <- result_12[order(result_12$model),]
result_prop_bic <- round(as.matrix(result_12[,7:14]), 2)
colnames(result_prop_bic) <- Methods
result_prop_pa <- round(as.matrix(result_12[,15:22]), 2)
colnames(result_prop_pa) <- Methods
result_all <- matrix(0, nrow = 2*nrow(result_12), ncol = ncol(result_prop_bic)+2)
for (i in 1:nrow(result_12)) {
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