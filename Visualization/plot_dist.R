####################################################################################
##   Plots for Scenario I: distributional response  ##
############################################################

############################# PART 1: Preparation #################################

function_path <- "~/DR4FrechetReg/Functions"
working_path <- "~/DR4FrechetReg/DistData"
save_path <- "~/DR4FrechetReg/Visualization"

# source all function scipts from the function path
function_sources <- list.files(function_path,
                               pattern="*.R$", full.names=TRUE,
                               ignore.case=TRUE)
sapply(function_sources, source, .GlobalEnv)
source(paste(working_path,"gendata_dist.R", sep="/"))

# color palette
palette <- colorRampPalette(c('#4A52A7', '#59AFEA', '#7BCEB8', '#A7DA64',
                              '#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131'))


# source('len_depth.R')
# source('angle_depth.R')
# source('kerpca.R')
########################################
## plotting
require('plotly')

##############################
### model3
#######################################

# set global parameters
n <- 200
p <- 10
m <- 100
non_elliptical <- -1
rho <- 0
model <- 'distex3'
kernel_type <- "Gaussian"
kernel_para <- 10
data_type <- "distribution"

set.seed(2022)

data <- gendata_dist(n=n, p=p, m=m, rho=rho, non_ellip=non_elliptical, mode = model)
d0 <- data$d0
b0 <- data$b0
ygram <- gram_matrix(data$y, complexity = kernel_para, type=data_type, kernel=kernel_type)
bhat.opg <- fopg(x=data$x, y=ygram, d=data$d0)$beta
pred = data$x%*%bhat.opg

# density
dens <- apply(data$y, 1, density, n=201)
plotdata <- data.frame(
  x = unlist(lapply(dens, "[", "x")),
  z = unlist(lapply(dens, "[", "y")),
  y = rep(data$x[,3], each = length(dens[[1]]$x)),
  csd1 = rep(pred[,1], each = length(dens[[1]]$x)),
  csd2 = rep(pred[,2], each = length(dens[[1]]$x)))

fig_chaos <- plot_ly(plotdata, x = ~x, y = ~y, z = ~z, type = 'scatter3d',
                     mode = 'lines', color = ~z, split = ~y, alpha = 1, colors = palette(100))%>% 
  layout(scene = list(
    xaxis = list(title = "", range = c(-20,40)),
    yaxis = list(title = list(text="X3", font=list(size=18))),
    zaxis = list(title = "")),showlegend = FALSE)


## vs first sufficient predictor
fig_csd <- plot_ly(plotdata, x = ~x, y = ~csd1, z = ~z, type = 'scatter3d', 
                   mode = 'lines', color = ~z, split = ~csd1, alpha = 1, colors = palette(100))%>% 
  layout(
    scene = list(
      xaxis = list(title = "", range = c(-20,40)),
      yaxis = list(title = list(text="1st SP", font=list(size=18))),
      zaxis = list(title = "")
    ),showlegend = FALSE)


## vs second sufficient predictor
fig_csd_sp2 <- plot_ly(plotdata, x = ~x, y = ~csd2, z = ~z, type = 'scatter3d', 
                       mode = 'lines', color = ~z, split = ~csd2, alpha=1,colors = palette(100))%>% 
  layout(
    scene = list(
      xaxis = list(title = "", range = c(-20,40)),
      yaxis = list(title = list(text="2nd SP", font=list(size=18))),
      zaxis = list(title = "")
    ),showlegend = FALSE)


