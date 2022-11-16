####################################################################################
##   Plots for Scenario II: spd matrix response  ##
############################################################

############################# PART 1: Preparation #################################

function_path <- "~/DR4FrechetReg/Functions"
working_path <- "~/DR4FrechetReg/SPDmatData"
save_path <- "~/DR4FrechetReg/Visualization"

# source all function scipts from the function path
function_sources <- list.files(function_path,
                               pattern="*.R$", full.names=TRUE,
                               ignore.case=TRUE)
sapply(function_sources, source, .GlobalEnv)
source(paste(working_path,"gendata_cov.R", sep="/"))

# color palette
palette <- colorRampPalette(c('#4A52A7', '#59AFEA', '#7BCEB8', '#A7DA64',
                              '#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131'))


########################################
## graph
require('plotly')
require('mixtools')

##############################
### model1
#######################################

# set global parameters
n <- 200
p <- 10
non_elliptical <- 0
rho <- 0
model <- 'model1'
method <- 'fopg'
kernel_type <- "Gaussian"
kernel_para <- 1
data_type <- "spd"
error_type <- 'log_normal'

set.seed(2022)

## plotting 
data <- gendata_cov(n=n, p=p, rho=rho, non_ellip=non_elliptical, model=model, error=error_type)
d0 <- data$d0 
b0 <- data$beta
ygram <- gram_matrix(data$y, complexity = kernel_para, type=data_type, kernel=kernel_type)
bhat.fopg <- fopg(x=data$x, y=ygram, d=data$d0)$beta
d <- dim(data$y[[1]])[1]
y <- lapply(data$y, function(x) {colnames(x)=letters[1:d];rownames(x)=LETTERS[1:d]; return(x)})
ellip <- lapply(X = y, ellipse, mu=c(0,0))
pred <- data$x%*%bhat.fopg

## vs first sufficient predictor
p <- plot_ly()
for (j in 1:200) {
  p <- p %>% 
    add_trace(type = 'scatter3d', mode = 'lines', x = ellip[[j]][,1], y = ellip[[j]][,2],
              z = rep(pred[j], dim(ellip[[1]])[1]), opacity=rank(pred)[j]/600+0.1, 
              line = list(color = 'rgb(22, 96, 167)', width = 3), 
              showlegend=FALSE)
}
p_sp <- p %>% layout(scene = list(xaxis = list(title = ''),
                          yaxis = list(title = ''),
                          zaxis = list(title = list(text="1st SP", font=list(size=18)))))

## chaos
pc <- plot_ly()
for (j in 1:200) {
  pc <- pc %>% 
    add_trace(type = 'scatter3d', mode = 'lines', x = ellip[[j]][,1], y = ellip[[j]][,2],
              z = rep(data$x[j,10], dim(ellip[[1]])[1]), opacity=rank(data$x[,10])[j]/600+0.1, 
              line = list(color = 'rgb(22, 96, 167)', width = 3), 
              showlegend=FALSE)
}
p_random <- pc%>% layout(scene = list(xaxis = list(title = ''),
                           yaxis = list(title = ''),
                           zaxis = list(title = list(text="X10", font=list(size=18)))))
