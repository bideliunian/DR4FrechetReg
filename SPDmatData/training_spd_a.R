####################################################################################
##   Simulations for Scenario II: SDR for spd matrices response  ##
##    including Model II-1, II-2   ########################################

############################# PART 1: Preparation #################################

function_path <- "~/DR4FrechetReg/Functions"
working_path <- "~/DR4FrechetReg/SPDmatData"
save_path <- "~/DR4FrechetReg/SPDmatData/Results"

# source all function scipts from the function path
function_sources <- list.files(function_path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function_sources, source, .GlobalEnv)
source(paste(working_path,"gendata_cov.R", sep="/"))
source(paste(working_path,"functions_train.R", sep="/"))

# packages
library(tictoc)

# reading the arguments passed from command line
args<- as.numeric(commandArgs(trailingOnly=TRUE))

# set global parameters
n_list <- c(200, 400)
p_list <- c(10, 20)
non_elliptical <- 0
rho <- 0
Models <- c('model1','model2')
Methods <-c('fols','fphd','fiht','fsir','fsave','fdr','fopg','wire') #'fmave'
repeat_time <- 1
kernel_type <- "Gaussian"
kernel_para <- 0.01
data_type <- "spd"
error_type <- 'log_normal'
# record parameters in a grid

grid_1 <- expand.grid(n = n_list[1], p=p_list[1], error=error_type, rho=rho, non_ellip=non_elliptical, model = Models, 
                      stringsAsFactors=FALSE)
grid_2 <- expand.grid(n = n_list[2], p=p_list[2], error=error_type, rho=rho, non_ellip=non_elliptical, model = Models, 
                      stringsAsFactors=FALSE)
grid <- rbind(grid_1,grid_2)

# set seed
set.seed(2021*args)

########################### PART 2: Order Determination and Training ##############################

time_start <- proc.time()
result <- mapply(main_cov, times = repeat_time, n=grid$n, p=grid$p, rho = grid$rho, non_ellip=grid$non_ellip, 
                 model=grid$model, data_type=data_type, error_type=error_type, kernel_type=kernel_type, 
                 complexity=kernel_para)
time_end <- proc.time()
time_run <- (time_end - time_start)[3]
print(time_run)

# save result
save(result, file = paste(save_path, "/error=", error_type,".rho=", rho,"nonellip=", non_elliptical,
                                      ".result",args,".RData",sep=""))

