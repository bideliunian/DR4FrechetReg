####################################################################################
##   Mortality Data Analysis: SDR for distributional response  ##
####################################################################

############################# PART 1: Preparation #################################

function_path <- "D:/Research/DR4FR/Codes/Functions"
working_path <- "D:/Research/DR4FR/Codes/MortalityData"
save_path <- "D:/Research/DR4FR/Codes/MortalityData"

# function_path <- "~/work/DR4FR/Functions"
# working_path <- "~/work/DR4FR/MortalityData"
# save_path <- "~/work/DR4FR/MortalityData"

library(frechet)
library(pracma)
library(dplyr)
library(ggplot2)
library(plotly)
library(RColorBrewer)

# source all function scipts from the function path
function_sources <- list.files(function_path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function_sources, source, .GlobalEnv)

# load data
load(paste(working_path, 'mortality_data.RData', sep="/"))
load(paste(working_path, 'density.RData', sep="/"))

# set parameters
gamma <- 20
method <- 'fopg'
palette <- colorRampPalette(c('#4A52A7', '#59AFEA', '#7BCEB8', '#A7DA64',
                              '#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131')) # color palette

#############################
# plot the densities
#############################

# p <- ggplot()+labs(title="",x="Age(0-100)", y = "Density")
# plot_density <- function(l){
#  df <- data.frame(x=l$x,y=l$y)
#  p <<- p + geom_line(data=df, aes(x,y), alpha=0.1)
#  print(p)
# }
# lapply(density, plot_density)

####################### PART 2: Estimating Sufficient Predictors ###########################

############################
# transform response
############################
var_names <- c('deaths_age','pop_density','sex_ratio','mean_age_child','gdp_capita','gva_arg',
               'cpi','unimp','health_exp','land_arable')#'migr_ratio','gva_ind','gdp_growth','life_exp'
Xlist <- lapply(mortality_data, `[[`, "Value")
Xlist$deaths_age <- NULL
X <- matrix(unlist(Xlist), ncol = length(var_names)-1, byrow = FALSE)
colnames(X) <- var_names[-1]
X <- apply(X,2,function(x) return((x-mean(x))/sd(x)))## standardize individually
pairs(X, pch = 19, cex = 0.5)


########################
# SDR using fopg
# gram matrix 
gram_wass <- function(y, complexity){
  n <- length(y)
  ##upper triangular index
  g <- expand.grid(row = 1:n, col = 1:n)
  ind <- g[upper.tri(diag(n), diag = FALSE), ]
  dist <- function(i,k){
    return(dist4den(d1=y[[i]],d2=y[[k]]))
  }
  kupper <- mapply(dist, i=ind[,1], k=ind[,2])
  k <- matrix(0, nrow = n, ncol = n)
  k[upper.tri(k,diag = FALSE)]=kupper^2
  k <- k+t(k)
  sigma2 <- sum(k)/choose(n,2)
  gamma <- complexity/(sigma2)
  print(gamma)
  return(exp(-gamma*k))
}

y_gram <- gram_wass(y=density, complexity=gamma)
p <- dim(X)[2]
set.seed(1)
dhat <- pred_aug(x = X, y = y_gram, s = 10, r = floor(p/2)+1, method = method)$rhat
print(paste("-------------- estimated order of dimension:", dhat, "-----------------"))
f <- get(method)
result <- f(x=X, y = y_gram, d=dhat)
eigen_value <- eigen(result[[2]])$values%>%round(4)
eigen_ratio <- eigen_value[-length(eigen_value)]/eigen_value[-1]
# eigen values 
# 0.0405 0.0210 0.0122 0.0055 0.0040 0.0031 0.0010 0.0006 0.0002
beta <- result$beta
#beta1  0.416 -0.424  0.114  0.770  0.027 -0.146 -0.020 0.130 -0.053
#beta2  0.186  0.155 -0.159 -0.135 -0.576 -0.714 -0.083 0.198  0.096
#beta3 -0.108  0.498  0.507  0.135  0.487 -0.403  0.139 0.211  0.045


########################## PART 3: Plotting ###############################

######################
# scree plot 
########################
var_explained_df <- data.frame(PC = factor(c(paste0("PC",1:dim(X)[2])),levels = c(paste0("PC",1:dim(X)[2]))),
                               Eigenvalues=eigen_value, ER=c(eigen_ratio,1))#c(paste0("PC",1:9), "PC 10",)
var_explained_df %>%
  ggplot(aes(x=PC,y=Eigenvalues, group=1))+
  geom_point(size=2)+
  geom_line()+
  labs(title="")

x0 <-seq(0,100,length.out=101) #outputGrid
csd <- as.matrix(X)%*%beta
chaos_index <- runif(length(density), 0, 1)
plot_list <- density

for (i in 1:length(plot_list)) {
  plot_list[[i]]['z']=plot_list[[i]]['y']
  plot_list[[i]]['y']=list(rep(chaos_index[i],length(x0)))
  plot_list[[i]]['csd1']=list(rep(csd[i,1],length(x0)))
}

plotdata <- data.frame(
  x = unlist(lapply(plot_list, "[[", "x")),
  y = unlist(lapply(plot_list, "[[", "y")),
  z = unlist(lapply(plot_list, "[[", "z")),
  csd1 = rep(csd[,1], each = length(plot_list[[1]]$x)),
  gdp = rep(X[,4], each = length(plot_list[[1]]$x)),
  csd2 = rep(csd[,2], each = length(plot_list[[1]]$x))
  )

## summary stat
plotdata_sum <- data.frame(
  csd1 = as.vector(csd[,1]),
  csd2 = as.vector(csd[,2]),
  mean = sapply(X = plot_list,FUN = function(l) sum(l[["z"]]*l[["x"]]))%>%as.vector(),
  mode = sapply(X = lapply(plot_list,"[[","z"), FUN = function(l) (which.max(l)-1))%>%as.vector(),
  mom2 = sapply(X = plot_list, FUN = function(l) sum(l[["z"]]*(l[["x"]])^2)),
  mom3 = sapply(X = plot_list, FUN = function(l) sum(l[["z"]]*(l[["x"]])^3))
)
plotdata_sum$var <- plotdata_sum$mom2 - plotdata_sum$mean^2
plotdata_sum$sd <- sqrt(plotdata_sum$var)
plotdata_sum$skew <- (plotdata_sum$mom3 - 3*plotdata_sum$mean*plotdata_sum$var - plotdata_sum$mean^3) / sqrt(plotdata_sum$var^{3})
plotdata_sum$wass <- sapply(X = density, FUN = function(l) dist4den(l, density[[1]])^2)
plotdata_sum$gdp <- X[,4]

# 3D plot, densities in random order 
fig_chaos <- plot_ly(plotdata, x = ~x, y = ~y, z = ~z, type = 'scatter3d',
                     mode = 'lines', color = ~z, split = ~y, alpha=1, colors = palette(100))%>% 
  layout(scene = list(
    xaxis = list(title = list(text="Age-at-death", font=list(size=18))),
    yaxis = list(title = ""),
    zaxis = list(title = "")),showlegend = FALSE)


# 3D plot, densities vs first sufficient predictor
fig_csd <- plot_ly(plotdata, x = ~x, y = ~csd1, z = ~z, type = 'scatter3d', 
                   mode = 'lines', color = ~z, split = ~csd1, alpha=1,colors = palette(100))%>% 
  layout(
    scene = list(
      xaxis = list(title = list(text="Age-at-death", font=list(size=18))),
      yaxis = list(title = list(text="1st SP", font=list(size=18)), range=c(-2.5,4)),
      zaxis = list(title = "")
    ),showlegend = FALSE)


# 3D plot, densities vs second sufficient predictors
fig_csd_sp2 <- plot_ly(plotdata, x = ~x, y = ~csd2, z = ~z, type = 'scatter3d', 
                   mode = 'lines', color = ~z, split = ~csd2, alpha=1,colors = palette(100))%>% 
  layout(
    scene = list(
      xaxis = list(title = list(text="Age-at-death", font=list(size=18))),
      yaxis = list(title = list(text="2nd SP", font=list(size=18)), range=c(-4,2.5)),
      zaxis = list(title = "")
    ),showlegend = FALSE)



# 3D plot, densities vs GDP 
fig_gdp <- plot_ly(plotdata, x = ~x, y = ~gdp, z = ~z, type = 'scatter3d',
                   mode = 'lines', color = ~z, split = ~gdp, alpha=1,colors = palette(100))%>%
  layout(
    scene = list(
      xaxis = list(title = list(text="Age-at-death", font=list(size=18))),
      yaxis = list(title = list(text="GDP", font=list(size=18))),
      zaxis = list(title = "")
    ),showlegend = FALSE)


# summary statistics vs first sufficient predictor
fig_mean <- plot_ly(data = plotdata_sum, x= plotdata_sum$csd1, y =plotdata_sum$mean, name = 'Mean', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = list(text="1st Sufficient Predictor", font=list(size=24))),
         legend = list(font = list(size = 24)))
fig_mode <- plot_ly(data = plotdata_sum, x= plotdata_sum$csd1, y =plotdata_sum$mode, name = 'Mode', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = list(text="1st Sufficient Predictor", font=list(size=24))),
         legend = list(font = list(size = 24)))
fig_var <- plot_ly(data = plotdata_sum, x= plotdata_sum$csd1, y =plotdata_sum$var,  name = 'Variance', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = list(text="1st Sufficient Predictor", font=list(size=24))),
         legend = list(font = list(size = 24)))
fig_skew <- plot_ly(data = plotdata_sum, x= plotdata_sum$csd1, y =plotdata_sum$skew,  name = 'Skewness', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = list(text="1st Sufficient Predictor", font=list(size=24))),
         legend = list(font = list(size = 24)))
fig_sd <- plot_ly(data = plotdata_sum, x= plotdata_sum$csd1, y =plotdata_sum$sd,  name = 'Sd', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = list(text="1st Sufficient Predictor", font=list(size=24))),
         legend = list(font = list(size = 24)))
fig_sum <- subplot(fig_mean,fig_mode,fig_var,fig_skew,nrows = 2, shareX = TRUE, titleX = TRUE)


# summary statistics vs second sufficient predictor
fig_mean_sp2 <- plot_ly(data = plotdata_sum, x= plotdata_sum$csd2, y =plotdata_sum$mean, name = 'Mean', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = list(text="2nd Sufficient Predictor", font=list(size=24))),
         legend = list(font = list(size = 24)))
fig_mode_sp2 <- plot_ly(data = plotdata_sum, x= plotdata_sum$csd2, y =plotdata_sum$mode, name = 'Mode', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = list(text="2nd Sufficient Predictor", font=list(size=24))),
         legend = list(font = list(size = 24)))
fig_var_sp2 <- plot_ly(data = plotdata_sum, x= plotdata_sum$csd2, y =plotdata_sum$var,  name = 'Variance', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = list(text="2nd Sufficient Predictor", font=list(size=24))),
         legend = list(font = list(size = 24)))
fig_skew_sp2 <- plot_ly(data = plotdata_sum, x= plotdata_sum$csd2, y =plotdata_sum$skew,  name = 'Skewness', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = list(text="2nd Sufficient Predictor", font=list(size=24))),
         legend = list(font = list(size = 24)))
fig_sd_sp2 <- plot_ly(data = plotdata_sum, x= plotdata_sum$csd2, y =plotdata_sum$sd,  name = 'Sd', type = 'scatter', mode = 'markers')%>%
  layout(xaxis = list(title = list(text="2nd Sufficient Predictor", font=list(size=24))),
         legend = list(font = list(size = 24)))
fig_sum_sp2 <- subplot(fig_mean_sp2,fig_mode_sp2,fig_var_sp2,fig_skew_sp2,nrows = 2, shareX = TRUE, titleX = TRUE)


# summary statistics vs the first two sufficient predictors
vnames <- c('mean','sd')
mypalette <- brewer.pal(3, "Dark2")
fig_sum_meansd_1st <- subplot(fig_mean, fig_sd, nrows = 2, shareX = TRUE, titleX = TRUE) %>% layout(colorway=mypalette[1:2], annotations=list(color=vnames))
fig_sum_meansd_2nd <- subplot(fig_mean_sp2, fig_sd_sp2, nrows = 2, shareX = TRUE, titleX = TRUE)%>% layout(colorway=mypalette[1:2], annotations=list(color=vnames))
fig_sum_meansd <- subplot(fig_sum_meansd_1st, style(fig_sum_meansd_2nd, showlegend = FALSE), nrows = 1, titleX = TRUE)


# summary statistics vs the first two sufficient predictors
library(latticeExtra)
palette_2 <- rev(colorRampPalette(brewer.pal(6, "RdYlGn"))(20))
fig_mean_2 <- levelplot(mean ~ csd1 * csd2, plotdata_sum, panel = panel.levelplot.points, cex = 1,
          xlab = '1st sufficient predictor', ylab = '2nd sufficient predictor')
fig_mode_2 <- levelplot(mode ~ csd1 * csd2, plotdata_sum, panel = panel.levelplot.points, cex = 1,
                      xlab = '1st sufficient predictor', ylab = '2nd sufficient predictor')
fig_var_2 <- levelplot(var ~ csd1 * csd2, plotdata_sum, panel = panel.levelplot.points, cex = 1,
                      xlab = '1st sufficient predictor', ylab = '2nd sufficient predictor')
fig_skew_2 <- levelplot(skew ~ csd1 * csd2, plotdata_sum, panel = panel.levelplot.points, cex = 1,
                     xlab = '1st sufficient predictor', ylab = '2nd sufficient predictor')
