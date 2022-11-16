####################################################################################
##   Simulations for Scenario I: SDR for distributional response  ##
##    including Model I-1, I-2   ########################################

############################# PART 1: Preparation #################################

function_path <- "~/DR4FrechetReg/Functions"
working_path <- "~/DR4FrechetReg/DistData"
save_path <- "~/DR4FrechetReg/DistData/Results"

# source all function scipts from the function path
function_sources <- list.files(function_path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function_sources, source, .GlobalEnv)
source(paste(working_path,"gendata_dist.R", sep="/"))
source(paste(working_path,"functions_train.R", sep="/"))

# packages
library(tictoc)

# reading the arguments passed from command line
args<- as.numeric(commandArgs(trailingOnly=TRUE))

# set global parameters
# set global parameters
n_list <- c(200, 400)
p_list <- c(10, 20)
m <- 100
non_elliptical <- -1
rho <- 0.5
Models <- c('distex3','distex4')
Methods <-c('fols','fphd','fiht','fsir','fsave','fdr','fopg','wire') #'fmave'
repeat_time <- 1
kernel_type <- "Gaussian"
kernel_para <- 1
data_type <- "distribution"
# record parameters in a grid

grid_1 <- expand.grid(n = n_list[1], p=p_list[1], m=m, rho=rho, non_ellip=non_elliptical, model = Models, 
                      stringsAsFactors=FALSE)
grid_2 <- expand.grid(n = n_list[2], p=p_list[2], m=m, rho=rho, non_ellip=non_elliptical, model = Models, 
                      stringsAsFactors=FALSE)
grid <- rbind(grid_1,grid_2)

# set seed
set.seed(2021*args)

########################### PART 2: Order Determination ##############################

time_start <- proc.time()
result <- mapply(main_dist, times = repeat_time, n=grid$n, p=grid$p, m=grid$m, rho = grid$rho, 
                       non_ellip=grid$non_ellip, model=grid$model, data_type=data_type, kernel_type=kernel_type, 
                       complexity=kernel_para)
time_end <- proc.time()
time_run <- (time_end - time_start)[3]
print(time_run)

# save result
save(result, file = paste(save_path, "/model34.rho=", rho,"nonellip=", non_elliptical,
                                      ".result",args,".RData",sep=""))

