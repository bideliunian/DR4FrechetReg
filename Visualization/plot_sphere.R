####################################################################################
##   Plots for Scenario III: sphere response  ##
############################################################

############################# PART 1: Preparation #################################

function_path <- "~/PARENT_DIRECTORY/Functions"
working_path <- "~/PARENT_DIRECTORY/SphereData"
save_path <- "~/PARENT_DIRECTORY/Visualization"

# source all function scipts from the function path
function_sources <- list.files(function_path,
                               pattern="*.R$", full.names=TRUE,
                               ignore.case=TRUE)
sapply(function_sources, source, .GlobalEnv)
source(paste(working_path,"gendata_sphere.R", sep="/"))

# color palette
palette <- colorRampPalette(c('#4A52A7', '#59AFEA', '#7BCEB8', '#A7DA64',
                              '#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131'))

require('plotly')

##############################
### model 1
#######################################

# set global parameters
n <- 200
p <- 10
rho <- 0
model <- 'model1'
method <- 'fopg'
kernel_type <- "Laplacian"
kernel_para <- 1
data_type <- "sphere"

set.seed(2022)

cyl = sqrt(1-seq(-1,1,length.out = 201)^2)
cyl.mat = (replicate(201, cyl))

data <- gendata_sphere(n=n, p=p, rho=rho, model=model)
d0 <- data$d0 
b0 <- data$beta
ygram <- gram_matrix(data$y, complexity = kernel_para, type=data_type, kernel=kernel_type)
bhat.fopg <- fopg(x=data$x, y=ygram, d=data$d0)$beta
pred <- data$x%*%bhat.fopg

par.data <- data.frame(pred = pred, y1 = data$y[,1], y2 = data$y[,2], 
                       rand = runif(min(pred),max(pred), n=n))

## vs 1st sufficient predictor
sp_model1 <- plot_ly(x = seq(min(par.data$pred),max(par.data$pred), length.out = 201), y = seq(-1,1,0.01), 
                    z = cyl.mat)%>%add_surface(color = rep(1,dim(cyl.mat)[1]^2), opacity = 0.5, showscale = FALSE)%>%
  add_surface(z = -cyl.mat, color = rep(1,dim(cyl.mat)[1]^2), opacity = 0.5, showscale = FALSE)%>%
  add_trace(data = par.data, x = ~pred, y = ~y1, z = ~y2, mode = "markers", type = "scatter3d", 
            marker = list(size = 4, color = '#BF382A', symbol = 104), opacity = 0.8, showlegend = FALSE)%>%
  layout(scene = list(xaxis = list(title = list(text="1st SP", font=list(size=18))),
                      yaxis = list(title = list(text="y1", font=list(size=18))),
                      zaxis = list(title = list(text="y2", font=list(size=18)))))
## chaos plot
cp_model1 <- plot_ly(x = seq(min(par.data$pred),max(par.data$pred), length.out = 201), y = seq(-1,1,0.01), 
                    z = cyl.mat)%>%add_surface(color = rep(1,dim(cyl.mat)[1]^2), opacity = 0.5, showscale = FALSE)%>%
  add_surface(z = -cyl.mat, color = rep(1,dim(cyl.mat)[1]^2), opacity = 0.5, showscale = FALSE)%>%
  add_trace(data = par.data, x = ~rand, y = ~y1, z = ~y2, mode = "markers", type = "scatter3d", 
            marker = list(size = 4, color = '#BF382A', symbol = 104), opacity = 0.8, showlegend = FALSE)%>%
  layout(scene = list(xaxis = list(title = list(text="X10", font=list(size=18))),
                      yaxis = list(title = list(text="y1", font=list(size=18))),
                      zaxis = list(title = list(text="y2", font=list(size=18)))))


##############################
### model 1
#######################################

# set global parameters
n <- 200
p <- 10
rho <- 0
model <- 'model2'
method <- 'fopg'
kernel_type <- "Laplacian"
kernel_para <- 1
data_type <- "sphere"

dd <- transform(expand.grid(theta=seq(0,pi,length=500),
                            phi=seq(0,2*pi,length=1000)),
                x = sin(theta)*cos(phi),
                y = sin(theta)*sin(phi),
                z = cos(theta))

set.seed(2022)

data <- gendata_sphere(n=n, p=p, rho=rho, model=model)
d0 <- data$d0 
b0 <- data$beta
ygram <- gram_matrix(data$y, complexity = kernel_para, type=data_type, kernel=kernel_type)
bhat.fopg <- fopg(x=data$x, y=ygram, d=data$d0)$beta
pred <- data$x%*%bhat.fopg

par.data2 <- data.frame(pred1 = pred[,1], pred2 = pred[,2], x = data$y[,1],
                        y = data$y[,2], z = data$y[,3], rand = runif(min(pred), max(pred), n=n))
## vs 1st sufficient predictor
sp_model2 = plot_ly(data=dd, type='mesh3d',
                    x = ~x,
                    y = ~y,
                    z = ~z, intensity = rep(0, dim(dd)[2]), colors = "#FFFFFF",opacity = 0.2, showscale = FALSE,
                    showlegend = FALSE)%>%
  add_trace(data = par.data2, x = ~x, y = ~y, z = ~z, mode = "markers", type = "scatter3d", 
            marker = list(size = 4, color = ~pred1, colorscale = c('#BF382A', '#0C4B8E'), showscale = TRUE),
            opacity = 1, showlegend = FALSE)%>%
  layout(scene = list(xaxis = list(title = list(text="y1", font=list(size=18))),
                      yaxis = list(title = list(text="y2", font=list(size=18))),
                      zaxis = list(title = list(text="y3", font=list(size=18)))))

## chaos plot
cp_model2 = plot_ly(data=dd, type='mesh3d',
                    x = ~x,
                    y = ~y,
                    z = ~z, intensity = rep(0, dim(dd)[2]), colors = "#FFFFFF",opacity = 0.2, showscale = FALSE, 
                    showlegend = FALSE)%>%
  add_trace(data = par.data2, x = ~x, y = ~y, z = ~z, mode = "markers", type = "scatter3d", 
            marker = list(size = 4, color = ~rand, colorscale = c('#BF382A', '#0C4B8E'), showscale = TRUE),
            opacity = 1, showlegend = FALSE)%>%
  layout(scene = list(xaxis = list(title = list(text="y1", font=list(size=18))),
                      yaxis = list(title = list(text="y2", font=list(size=18))),
                      zaxis = list(title = list(text="y3", font=list(size=18)))))

## vs second sufficient predictor
sp_model2_2 = plot_ly(data=dd, type='mesh3d',
                    x = ~x,
                    y = ~y,
                    z = ~z, intensity = rep(0, dim(dd)[2]), colors = "#FFFFFF",opacity = 0.2, showscale = FALSE, 
                    showlegend = FALSE)%>%
  add_trace(data = par.data2, x = ~x, y = ~y, z = ~z, mode = "markers", type = "scatter3d", 
            marker = list(size = 4, color = ~pred2, colorscale = c('#BF382A', '#0C4B8E'), showscale = TRUE),
            opacity = 1, showlegend = FALSE)%>%
  layout(scene = list(xaxis = list(title = list(text="y1", font=list(size=18))),
                      yaxis = list(title = list(text="y2", font=list(size=18))),
                      zaxis = list(title = list(text="y3", font=list(size=18)))))

