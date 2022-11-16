###########################################################
############ read the result and box plot #################
###########################################################

save_path <- "~/DR4FrechetReg/DistData/OrderResults"
# if use aci
#save_path <- "~/work/DR4FR/DistData/OrderResults"

############################## Model 1&2##############################
# set global parameters
times_repeat <- 100
n_list <- c(200, 400)
p_list <- c(10, 20)
m <- 100
non_elliptical <- 0
rho <- 0
Models <- c('distex1','distex2')
Methods <-c('fols','fphd','fiht','fsir','fsave','fdr','fopg','wire') #'fmave'
kernel_type <- "Gaussian"
kernel_para <- 1
data_type <- "distribution"

# record parameters in a grid
grid_1 <- expand.grid(n = n_list[1], p=p_list[1], m=m, rho=rho, 
                      non_ellip=non_elliptical, model = Models, 
                      stringsAsFactors=FALSE)
grid_2 <- expand.grid(n = n_list[2], p=p_list[2], m=m, rho=rho, 
                      non_ellip=non_elliptical, model = Models, 
                      stringsAsFactors=FALSE)
grid_12 <- rbind(grid_1,grid_2)


order_result_12 <- matrix(0, nrow=nrow(grid_12), ncol=2*length(Methods))
for (i in 1:times_repeat) {
  load(file = paste(save_path, "/order.model12.rho=", rho,"nonellip=", non_elliptical,
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

################################# Model 3&4 ###############################
# set global parameters
times_repeat <- 100
n_list <- c(200, 400)
p_list <- c(10, 20)
m <- 100
non_elliptical <- -1
rho <- 0.5
Models <- c('distex3','distex4')
Methods <-c('fols','fphd','fiht','fsir','fsave','fdr','fopg','wire') #'fmave'
kernel_type <- "Gaussian"
kernel_para <- 1
data_type <- "distribution"

# record parameters in a grid
grid_3 <- expand.grid(n = n_list[1], p=p_list[1], m=m, rho=rho, 
                      non_ellip=non_elliptical, model = Models, 
                      stringsAsFactors=FALSE)
grid_4 <- expand.grid(n = n_list[2], p=p_list[2], m=m, rho=rho, 
                      non_ellip=non_elliptical, model = Models, 
                      stringsAsFactors=FALSE)
grid_34 <- rbind(grid_3, grid_4)

order_result_34 <- matrix(0, nrow=nrow(grid_34), ncol=2*length(Methods))
for (i in 1:times_repeat) {
  load(file = paste(save_path, "/order.model34.rho=", rho,"nonellip=", non_elliptical,
                    ".result",i,".RData",sep=""))
  order_result_34 <- order_result_34 + order_result_table[,7:ncol(order_result_table)] / times_repeat
}

result_34 <- cbind(grid_34, order_result_34)

# write in latex
library(xtable)
result_34 <- result_34[order(result_34$model),]
result_prop_bic <- round(as.matrix(result_34[,7:14]), 2)
colnames(result_prop_bic) <- Methods
result_prop_pa <- round(as.matrix(result_34[,15:22]), 2)
colnames(result_prop_pa) <- Methods
result_all <- matrix(0, nrow = 2*nrow(result_34), ncol = ncol(result_prop_bic)+2)
for (i in 1:nrow(result_34)) {
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

############################ box plot ####################################

library(ggplot2)

# box plot of model 1
error_model1_p10 <- as.data.frame(t(error_array_12[1,,]))
error_model1_p20 <- as.data.frame(t(error_array_12[3,,]))
names(error_model1_p10) <- Methods
names(error_model1_p20) <- Methods
error_model1_p10_stack <- stack(error_model1_p10)
error_model1_p10_stack$p <- 10
error_model1_p20_stack <- stack(error_model1_p20)
error_model1_p20_stack$p <- 20
error_model1_stack <- rbind(error_model1_p10_stack, error_model1_p20_stack)
ggplot(error_model1_stack, aes(x = ind, y = values, fill=p)) + 
  geom_boxplot()+labs(title="",x="methods", y = "errors")

# box plot of model 2
error_model2_p10 <- as.data.frame(t(error_array_12[2,,]))
error_model2_p20 <- as.data.frame(t(error_array_12[4,,]))
names(error_model2_p10) <- Methods
names(error_model2_p20) <- Methods
error_model2_p10_stack <- stack(error_model2_p10)
error_model2_p10_stack$p <- 10
error_model2_p20_stack <- stack(error_model2_p20)
error_model2_p20_stack$p <- 20
error_model2_stack <- rbind(error_model2_p10_stack, error_model2_p20_stack)
ggplot(error_model2_stack, aes(x = ind, y = values, fill=p)) + 
  geom_boxplot()+labs(title="",x="methods", y = "errors")