####################################################################################
##   Mortality Data Analysis: SDR for distributional response  ##
####################################################################

############################# PART 1: Preparation #################################


function_path <- "~/DR4FrechetReg/Functions"
working_path <- "~/DR4FrechetReg/StrokeData"
save_path <- "~/DR4FrechetReg/StrokeData"

library(frechet)
library(dplyr)
library(plotly)
library(RColorBrewer)
library(WRI)



# source all function scipts from the function path
function_sources <- list.files(function_path, 
                               pattern="*.R$", full.names=TRUE, 
                               ignore.case=TRUE)
sapply(function_sources, source, .GlobalEnv)


# set parameters
gamma <- 1
method <- 'fopg'
palette <- colorRampPalette(c('#4A52A7', '#59AFEA', '#7BCEB8', '#A7DA64',
                              '#EFF121', '#F5952D', '#E93131', '#D70131', '#D70131')) # color palette

# read in data
data(strokeCTdensity)
strokeCTdensity$predictors <- strokeCTdensity$predictors[c("log_b_vol","b_shapInd","B_TimeCT","age","weight")]
X <- as.matrix(strokeCTdensity$predict)
Y <- split(strokeCTdensity$densityCurve,row(strokeCTdensity$densityCurve))%>%
 lapply(function(x) x <- list(y=x, x=strokeCTdensity$densitySupport))

## visualization of the response density curves
# p <- ggplot()+labs(title="",x="Density Support (0,1)", y = "Density")
# plot <- function(l){
#   df <- data.frame(x=l$x,y=l$y)
#   p <<- p + geom_line(data=df, aes(x,y),alpha=0.2)
# }
# lapply(Y[c(1:30)], plot)

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

y_gram <- gram_wass(y=Y, complexity=gamma)
p <- dim(X)[2]

set.seed(2022)
dhat <- pred_aug(x = X, y = y_gram, s = 10, r = floor(p/2)+1, method = method)$rhat
print(paste("-------------- estimated order of dimension:", dhat, "-----------------"))
f <- get(method)
result <- f(x=X, y = y_gram, d=dhat)
eigen_value <- eigen(result[[2]])$values%>%round(4)
eigen_ratio <- eigen_value[-length(eigen_value)]/eigen_value[-1]
beta <- result$beta
print(paste("-------------- estimated direction:"))
t(round(beta, 4))

#scree plot 
var_explained_df <- data.frame(PC= paste0("PC",1:dim(X)[2]),
                               Eigenvalues=eigen(result[[2]])$values)
var_explained_df %>%
  ggplot(aes(x=PC,y=Eigenvalues, group=1))+
  geom_point(size=2)+
  geom_line()+
  labs(title="")

######## y vs 1st sufficient predictor 
csd <- X%*%beta
plot_list <- Y
for (i in 1:length(plot_list)) {
  plot_list[[i]]['z']=plot_list[[i]]['y']
  plot_list[[i]]['y']=list(rep(0.01*i,101))
  plot_list[[i]]['csd']=list(rep(csd[i],101))
}

plotdata <- data.frame(
  x = unlist(lapply(X = plot_list, FUN = "[[", "x")),
  y = unlist(lapply(X = plot_list, FUN = "[[", "y")),
  z = unlist(lapply(X = plot_list, FUN = "[[", "z")),
  dens = unlist(lapply(X = plot_list, FUN = "[[", "z")),
  csd1 = rep(csd, each = length(plot_list[[1]]$x))
  )

plotdata_sum <- data.frame(
  csd1 = as.vector(csd),
  mean = sapply(X = plot_list,FUN = function(l) sum(l[["z"]]*l[["x"]])/(length(l[["z"]])-1))%>%as.vector(),
  mode = sapply(X = lapply(plot_list,"[[","z"), FUN = function(l) (which.max(l)-1)/100)%>%as.vector(),
  mom2 = sapply(X = plot_list, FUN = function(l) sum(l[["z"]]*(l[["x"]])^2))/100,
  mom3 = sapply(X = plot_list, FUN = function(l) sum(l[["z"]]*(l[["x"]])^3))/100
)

plotdata_sum$var <- plotdata_sum$mom2-plotdata_sum$mean^2
plotdata_sum$sd <- sqrt(plotdata_sum$var)
plotdata_sum$skew <- (plotdata_sum$mom3-3*plotdata_sum$mean*plotdata_sum$var-plotdata_sum$mean^3)/sqrt(plotdata_sum$var^{3})

# chaos plot
fig_chaos <- plot_ly(plotdata, x = ~x, y = ~y, z = ~dens, type = 'scatter3d',
                     mode = 'lines', color = ~dens, split = ~y, alpha=1, colors = palette(100))%>% #'RePu'
  layout(scene = list(
      xaxis = list(title = list(text="hema density", font=list(size=24))),
      yaxis = list(title = ""),
      zaxis = list(title = list(text="density", font=list(size=24)))), showlegend = FALSE)

htmlwidgets::saveWidget(
  widget = fig_chaos, 
  file = paste(working_path, "fig_chaos.html", sep="/"),
  selfcontained = TRUE
)

## versus first cs direction
fig_csd <- plot_ly(plotdata, x = ~x, y = ~csd1, z = ~dens, type = 'scatter3d', 
                   mode = 'lines',color = ~dens, split = ~csd1, alpha=1, colors = palette(100))%>% 
  layout(scene = list(
      xaxis = list(title = list(text="hema density", font=list(size=24))),
      yaxis = list(title = list(text="1st SP", font=list(size=24))),
      zaxis = list(title = list(text="density", font=list(size=24)))),
      showlegend = FALSE)

htmlwidgets::saveWidget(
  widget = fig_chaos, 
  file = paste(working_path, "fig_csd.html", sep="/"),
  selfcontained = TRUE
)

## summary stats plots
fig_mean <- plot_ly(data = plotdata_sum, x= plotdata_sum$csd1, y = plotdata_sum$mean, name = 'Mean', type = 'scatter', mode = 'markers')%>%
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
fig_sum <- subplot(fig_mean,fig_mode,fig_sd,fig_skew,nrows = 2, shareX = TRUE, titleX = TRUE) 

htmlwidgets::saveWidget(
  widget = fig_sum_meansd, 
  file = paste(working_path, "fig_sum.html", sep="/"),
  selfcontained = TRUE
)