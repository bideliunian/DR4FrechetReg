###########################################################
############ read the result and box plot #################
###########################################################

save_path <- "~/DR4FrechetReg/SPDmatData/Results"

# if use aci
#save_path <- "~/work/DR4FR/SPDmatData/Results"

############################## Model 1&2##############################
############# predictor setting (a) ##########################
# set global parameters
times_repeat <- 100
n_list <- c(200, 400)
p_list <- c(10, 20)
error_type <- "log_normal"
non_elliptical <- 0
rho <- 0
Models <- c('model1','model2')
Methods <-c('fols','fphd','fiht','fsir','fsave','fdr','fopg','wire') #'fmave'
kernel_type <- "Gaussian"
kernel_para <- 1
data_type <- "spd"

# record parameters in a grid
grid_1 <- expand.grid(n = n_list[1], p=p_list[1], error=error_type, rho=rho, 
                      non_ellip=non_elliptical, model = Models, 
                      stringsAsFactors=FALSE)
grid_2 <- expand.grid(n = n_list[2], p=p_list[2], error=error_type, rho=rho, 
                      non_ellip=non_elliptical, model = Models, 
                      stringsAsFactors=FALSE)
grid_12 <- rbind(grid_1,grid_2)

error_array_12 <- array(NA, dim=c(nrow(grid_12), length(Methods), times_repeat))
order_array_12 <- array(NA, dim=c(nrow(grid_12), length(Methods), times_repeat))
for (i in 1:times_repeat) {
  load(file = paste(save_path, "/error=", error_type,".rho=", rho,"nonellip=", non_elliptical,
                    ".result",i,".RData",sep=""))
  # estimation error
  error_array_12[,,i] <-  matrix(unlist(result[1, ]), nrow=nrow(grid_12), byrow=TRUE)
  # percentage of correct identification
  order_array_12[,,i] <- matrix(unlist(result[2, ]), nrow=nrow(grid_12), byrow=TRUE)
}

error_mean_12 <- apply(error_array_12, MARGIN=c(1,2), mean)
error_sd_12 <- apply(error_array_12, MARGIN=c(1,2), sd)
order_perc_12 <- apply(order_array_12, MARGIN=c(1,2), mean)

colnames(error_mean_12) <- paste(Methods, "mean")
colnames(error_sd_12) <- paste(Methods, "sd")
colnames(order_perc_12) <- Methods

result_12_a <- cbind(grid_12, error_mean_12, error_sd_12, order_perc_12)

################ predictor setting (b) ##########################
non_elliptical <- 1
rho <- 0.5

error_array_12 <- array(NA, dim=c(nrow(grid_12), length(Methods), times_repeat))
order_array_12 <- array(NA, dim=c(nrow(grid_12), length(Methods), times_repeat))
for (i in 1:times_repeat) {
  load(file = paste(save_path, "/error=", error_type,".rho=", rho,"nonellip=", non_elliptical,
                    ".result",i,".RData",sep=""))
  # estimation error
  error_array_12[,,i] <-  matrix(unlist(result[1, ]), nrow=nrow(grid_12), byrow=TRUE)
  # percentage of correct identification
  order_array_12[,,i] <- matrix(unlist(result[2, ]), nrow=nrow(grid_12), byrow=TRUE)
}

error_mean_12 <- apply(error_array_12, MARGIN=c(1,2), mean)
error_sd_12 <- apply(error_array_12, MARGIN=c(1,2), sd)
order_perc_12 <- apply(order_array_12, MARGIN=c(1,2), mean)

colnames(error_mean_12) <- paste(Methods, "mean")
colnames(error_sd_12) <- paste(Methods, "sd")
colnames(order_perc_12) <- Methods

result_12_b <- cbind(grid_12, error_mean_12, error_sd_12, order_perc_12)


############ function to write result in latex
result_12 <- rbind(result_12_a, result_12_b)
result_12 <- result_12[order(result_12$model),]
result_mean <- round(as.matrix(result_12[,7:14]),3)
result_sd <- round(as.matrix(result_12[,15:22]),3)
result_prop <- round(as.matrix(result_12[,23:ncol(result_12)]), 2)
result_all <- matrix(0, nrow = 3*nrow(result_12), ncol = ncol(result_mean)+1)
for (i in 1:nrow(result_12)) {
  for (j in 1:ncol(result_all)) {
    if(j==1){
      result_all[(3*(i-1)+1),j] <- 0
      if(i%%2 == 1){
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
