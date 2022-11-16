############## evaluate benchmark error via Monte Carlo ################

function_path <- "~/DR4FrechetReg/Functions"
working_path <- "~/DR4FrechetReg/SpehreData"
# 
# function_path <- "~/work/DR4FR/Functions"
# working_path <- "~/work/DR4FR/SpehreData"

# source all function scipts from the function path
function_sources <- list.files(function_path,
                               pattern="*.R$", full.names=TRUE,
                               ignore.case=TRUE)
sapply(function_sources, source, .GlobalEnv)
source(paste(working_path,"gendata_sphere.R", sep="/"))

# set global parameters
n_list <- c(200, 400)
p_list <- c(10, 20)
rho <- 0
Models <- c('model1','model2','model3')
# record parameters in a grid
grid_1 <- expand.grid(n = n_list[1], p=p_list[1], model = Models, 
                      stringsAsFactors=FALSE)
grid_2 <- expand.grid(n = n_list[2], p=p_list[2], model = Models, 
                      stringsAsFactors=FALSE)
grid <- rbind(grid_1,grid_2)
grid <- grid[order(grid$model),]
b0_list <- list()
for (i in 1:nrow(grid)){
  data <- gendata_sphere(n=grid[i,1], p=grid[i,2], rho=rho, model=grid[i,3])
  b0 <- data$b0
  b0_list[[i]] <- b0
}

set.seed(2021)
bk_list <- lapply(b0_list, bchmk, N=1000)
bk_mat <- round(matrix(unlist(bk_list), nrow = nrow(grid), ncol = 2, byrow = TRUE),3)
grid_bk <- cbind(grid, bk_mat)
