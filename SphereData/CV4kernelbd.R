####################################################################
##  Scenario I: Frechet SDR for distributions  ####################
##   Using CV to determine the bd and kernel type  ###############
###################################################################

########### PART 1: Preparation #################################

function_path <- "D:/Research/DR4FR/Codes/Functions"
working_path <- "D:/Research/DR4FR/Codes/SphereData"
save_path <- "D:/Research/DR4FR/Codes/SphereData/Results"

# function_path <- "~/work/DR4FR/Functions"
# working_path <- "~/work/DR4FR/SphereData"
# save_path <- "~/work/DR4FR/SphereData/Results"

# packages
library(frechet)

# source all function scipts from the function path
function_sources <- list.files(function_path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function_sources, source, .GlobalEnv)
source(paste(working_path,"gendata_sphere.R", sep="/"))

# set parameters
times <- 5
n <- 100
p <- 10
rho <- 0
non_ellip <- 0
num_fold <- 5
order <- 5
model <- 'model1'
Methods <- 'fols'
repeat_time <- 10
kernel_types <- c("Gaussian", "Laplacian")
kernel_paras <- c(10^-2, 10^-1, 1, 10)
data_type <- "sphere"
grid <- expand.grid(n=n, p=p, rho=rho, model = Models, stringsAsFactors=FALSE)

########### PART 2: k-fold CV to choose kernel type and kernel bindwidth ##############
kernel_chosen <- matrix(NA, nrow = times, ncol = 2)
for (k in 1:times) {
  data <- gendata_sphere(n = n, p = p, rho = rho, model = model)
  d0 <- data$d0
  b0 <- data$b0
  
  index_list <- kfoldCV(n = n, k = num_fold)
  cv_pred_error_avg <- matrix(NA, nrow = length(kernel_types), ncol = length(kernel_paras))
  
  for (i in 1:length(kernel_types)){
    ktype <- kernel_types[i]
    for (j in 1:length(kernel_paras)){
      gamma <- kernel_paras[j]
      pred_error_list <- lapply(X = index_list, FUN = cv4kernel, x = data$x, y = data$y, type = data_type, complexity = gamma, 
             kernel_type = ktype, method = 'fols', d = order)
      cv_pred_error_avg[i, j] <- mean(unlist(pred_error_list))
    }
  }
  index_min <- which(cv_pred_error_avg == min(cv_pred_error_avg), arr.ind = TRUE)
  kernel_chosen[k, ] <- index_min
}

kernel_type_index <- getmode(kernel_chosen[,1])
gamma_index <- getmode(kernel_chosen[,2])
kernel_type_chosen <- kernel_types[kernel_type_index]
kernel_para_chosen <- kernel_paras[gamma_index]

print(paste("Selected kernel type is: ", kernel_type_chosen, "; complexity is: ", kernel_para_chosen))

