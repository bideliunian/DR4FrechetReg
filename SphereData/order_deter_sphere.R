####################################################################################
##   Simulations for Scenario III: Order determination for sphere data response  ##
##    including Model III-1, III-2   ########################################

############################# PART 1: Preparation #################################

function_path <- "~/DR4FrechetReg/Functions"
working_path <- "~/DR4FrechetReg/SphereData"
save_path <- "~/DR4FrechetReg/SphereData/OrderResults"
# 
# function_path <- "~/work/DR4FR/Functions"
# working_path <- "~/work/DR4FR/SphereData"
# save_path <- "~/work/DR4FR/SphereData/OrderResults"

# source all function scipts from the function path
function_sources <- list.files(function_path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function_sources, source, .GlobalEnv)
source(paste(working_path,"gendata_sphere.R", sep="/"))
source(paste(working_path,"functions_order.R", sep="/"))

# packages
library(tictoc)

# reading the arguments passed from command line
args<- as.numeric(commandArgs(trailingOnly=TRUE))

# set global parameters
n_list <- c(200, 400)
p_list <- c(10, 20)
rho <- 0
Models <- c('model1','model2','model3')
Methods <-c('fols','fphd','fiht','fsir','fsave','fdr','fopg','wire') #'fmave'
repeat_time <- 1
kernel_type <- "Laplacian"
kernel_para <- 0.01
data_type <- "sphere"

# record parameters in a grid

grid_1 <- expand.grid(n = n_list[1], p=p_list[1], rho=rho, model = Models, 
                      stringsAsFactors=FALSE)
grid_2 <- expand.grid(n = n_list[2], p=p_list[2], rho=rho, model = Models, 
                      stringsAsFactors=FALSE)
grid <- rbind(grid_1,grid_2)

# set seed
set.seed(2021*args)

########################### PART 2: Order Determination ##############################

time_start <- proc.time()
order_result <- mapply(order_deter_sphere, times = repeat_time, n=grid$n, p=grid$p, 
                       rho = grid$rho, model=grid$model, data_type=data_type, 
                       kernel_type=kernel_type, complexity=kernel_para)
time_end <- proc.time()
time_run <- (time_end - time_start)[3]

order_result_mat <- matrix(NA,  nrow = nrow(grid), ncol = 2 * length(Methods))
for (i in 1:nrow(order_result_mat)) {
  order_result_mat[i, ] <- unlist(order_result[, i])
}
colnames(order_result_mat) <- names(unlist(order_result[, 1]))
order_result_table <- cbind(grid, as.data.frame(order_result_mat))
print(time_run)

# save result
save(order_result_table, file = paste(save_path, "/order.rho=", rho, 
                                ".result",args,".RData",sep=""))

