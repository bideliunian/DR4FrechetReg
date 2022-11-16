###########################################################
############ read the result and box plot #################
###########################################################

save_path <- "~/DR4FrechetReg/SphereData/Results"

# if use aci
#save_path <- "~/work/DR4FR/SphereData/Results"

############################## Model 1&2##############################
############ rho = 0 #############################
# set global parameters
times_repeat <- 100
n_list <- c(200, 400)
p_list <- c(10, 20)
rho <- 0
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

error_array <- array(NA, dim=c(nrow(grid), length(Methods), times_repeat))
order_array <- array(NA, dim=c(nrow(grid), length(Methods), times_repeat))
for (i in 1:times_repeat) {
  load(file = paste(save_path, "/rho=", rho, ".result", i, ".RData", sep=""))
  # estimation error
  error_array[,,i] <-  matrix(unlist(result[1, ]), nrow=nrow(grid), byrow=TRUE)
  # percentage of correct identification
  order_array[,,i] <- matrix(unlist(result[2, ]), nrow=nrow(grid), byrow=TRUE)
}

error_mean <- apply(error_array, MARGIN=c(1,2), mean)
error_sd <- apply(error_array, MARGIN=c(1,2), sd)
order_perc <- apply(order_array, MARGIN=c(1,2), mean)

colnames(error_mean) <- paste(Methods, "mean")
colnames(error_sd) <- paste(Methods, "sd")
colnames(order_perc) <- Methods

result_1 <- cbind(grid, error_mean, error_sd, order_perc)

############ rho = 0.5 #############################
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

error_array <- array(NA, dim=c(nrow(grid), length(Methods), times_repeat))
order_array <- array(NA, dim=c(nrow(grid), length(Methods), times_repeat))
for (i in 1:times_repeat) {
  load(file = paste(save_path, "/rho=", rho, ".result", i, ".RData", sep=""))
  # estimation error
  error_array[,,i] <-  matrix(unlist(result[1, ]), nrow=nrow(grid), byrow=TRUE)
  # percentage of correct identification
  order_array[,,i] <- matrix(unlist(result[2, ]), nrow=nrow(grid), byrow=TRUE)
}

error_mean <- apply(error_array, MARGIN=c(1,2), mean)
error_sd <- apply(error_array, MARGIN=c(1,2), sd)
order_perc <- apply(order_array, MARGIN=c(1,2), mean)

colnames(error_mean) <- paste(Methods, "mean")
colnames(error_sd) <- paste(Methods, "sd")
colnames(order_perc) <- Methods

result_2 <- cbind(grid, error_mean, error_sd, order_perc)


############ function to write result in latex
result_12 <- rbind(result_1, result_2)
result_12 <- result_12[order(result_12$model),]
result_mean <- round(as.matrix(result_12[,5:12]),3)
result_sd <- round(as.matrix(result_12[,13:20]),3)
result_prop <- round(as.matrix(result_12[,21:ncol(result_12)]), 2)
result_all <- matrix(0, nrow = 3*nrow(result_12), ncol = ncol(result_mean)+1)
for (i in 1:nrow(result_12)) {
  for (j in 1:ncol(result_all)) {
    if(j==1){
      result_all[(3*(i-1)+1),j] <- 0
      if(i %%2 == 1){
        result_all[(3*(i-1)+2),j] <- paste("(", 10,",",200,")", sep = "")   
      }
      else{
        result_all[(3*(i-1)+2),j] <- paste("(", 20,"," ,400,")", sep = "")   
      }
      result_all[(3*(i-1)+3),j] <- 0
    }
    else{
      result_all[(3*(i-1)+1),j] <- paste(result_prop[i,j-1] * 100, "%", sep = "")
      result_all[(3*(i-1)+2),j] <- result_mean[i,j-1]
      result_all[(3*(i-1)+3),j] <- paste("(", result_sd[i,j-1],")", sep = "")  
    }
  }
}
names(result_all) <- names(result_sd)
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
