####################################################################################
############################   Processing Mortality Data ######################
##################################################################################

function_path <- "D:/Research/DR4FR/Codes/Functions"
working_path <- "D:/Research/DR4FR/Codes/MortalityData"
save_path <- "D:/Research/DR4FR/Codes/MortalityData"

# function_path <- "~/work/DR4FR/Functions"
# working_path <- "~/work/DR4FR/MortalityData"
# save_path <- "~/work/DR4FR/MortalityData"

setwd(paste(save_path, "raw_data", sep = '/'))

# load packages
library(frechet)
library("readxl")
library(dplyr)


# read and clean data, predictors
my_data <- list()
# deaths age table
deaths_age <- read_excel("WPP2019_MORT_F17_1_ABRIDGED_LIFE_TABLE_BOTH_SEXES.xlsx",sheet = "ESTIMATES", na = '...')%>%
  filter(Type=='Country/Area'&Period=='2015-2020'&`Age (x)`!=0)%>%
  select(`Region, subregion, country or area *`,`Country code`,`Age (x)`,`Number of deaths d(x,n)`)%>%
  rename(country.or.area=`Region, subregion, country or area *`, Country.code=`Country code`)
# popolation density
pop_density <- read_excel("WPP2019_POP_F06_POPULATION_DENSITY.xlsx",sheet = "ESTIMATES", na = '...')%>%
  filter(Type=='Country/Area')%>%
  select(`Region, subregion, country or area *`,`Country code`,'2015')%>%
  rename(country.or.area=`Region, subregion, country or area *`, Country.code=`Country code`, Value='2015')
# sex ratio
sex_ratio <- read_excel("WPP2019_POP_F04_SEX_RATIO_OF_TOTAL_POPULATION.xlsx",sheet = "ESTIMATES", na = '...')%>%
  filter(Type=='Country/Area')%>%
  select(`Region, subregion, country or area *`,`Country code`,'2015')%>%
  rename(country.or.area=`Region, subregion, country or area *`, Country.code=`Country code`,Value='2015')
# migration rate
migr_ratio <- read_excel("WPP2019_MIGR_F01_NET_MIGRATION_RATE.xlsx",sheet = "ESTIMATES", na = '...')%>%
  filter(Type=='Country/Area')%>%
  select(`Region, subregion, country or area *`,`Country code`,'2015-2020')%>%
  rename(country.or.area=`Region, subregion, country or area *`, Country.code=`Country code`,Value='2015-2020')
# MEAN_AGE_CHILDBEARING
mean_age_child <- read_excel("WPP2019_SA1_FERT_F08_MEAN_AGE_CHILDBEARING.xlsx",sheet = "ESTIMATES AND MEDIUM VARIANT", 
                             na = '...')%>%
  filter(Type=='Country/Area')%>%
  select(`Type of aggregate, group, and constituents *`,`Country code`,'2015-2020')%>%
  rename(country.or.area=`Type of aggregate, group, and constituents *`, Country.code=`Country code`,Value='2015-2020')%>%
  distinct()
# GDP per Capita
gdp_capita <- read.csv("SYB63_230_202009_GDP and GDP Per Capita.csv",fileEncoding="UTF-8-BOM")%>%
  filter(Series=='GDP per capita (US dollars)'&Year==2015)%>%
  select(country.or.area,Country.code,Value)
# GDP growth
gdp_growth <- read.csv("SYB63_230_202009_GDP and GDP Per Capita.csv",fileEncoding="UTF-8-BOM")%>%
  filter(Series=='GDP real rates of growth (percent)'&Year==2015)%>%
  select(country.or.area,Country.code,Value)
#  Gross Value Added by Agriculture
gva_arg <- read.csv("SYB63_153_202009_Gross Value Added by Economic Activity.csv",fileEncoding="UTF-8-BOM")%>%
  filter(Series=='Agriculture, hunting, forestry and fishing (% of gross value added)'&Year==2015)%>%
  select(country.or.area,Country.code,Value)
# Gross Value Added by Industry
gva_ind <- read.csv("SYB63_153_202009_Gross Value Added by Economic Activity.csv",fileEncoding="UTF-8-BOM")%>%
  filter(Series=='Industry (% of gross value added)'&Year==2015)%>%
  select(country.or.area,Country.code,Value)
# Public Expenditure on Education
edu_exp <- read.csv("SYB63_245_202009_Public Expenditure on Education.csv",fileEncoding="UTF-8-BOM")%>%
  filter(Series=='Public expenditure on education (% of GDP)'&Year==2015)%>%
  select(country.or.area,Country.code,Value)
# Gross enrollment ratio - Secondary (male)
edu_erll_m <- read.csv("SYB63_309_202009_Education.csv",fileEncoding="UTF-8-BOM")%>%
  filter(Series=='Gross enrollment ratio - Secondary (male)'&Year==2015)%>%
  select(country.or.area,Country.code,Value)
# Gross enrollment ratio - Secondary (female)
edu_erll_f <- read.csv("SYB63_309_202009_Education.csv",fileEncoding="UTF-8-BOM")%>%
  filter(Series=='Gross enrollment ratio - Secondary (female)'&Year==2015)%>%
  select(country.or.area,Country.code,Value)
# CPI
cpi <- read.csv("SYB63_128_202009_Consumer Price Index.csv",fileEncoding="UTF-8-BOM")%>%
  filter(Series=='Consumer price index: General'&Year==2015)%>%
  select(country.or.area,Country.code,Value)
# Labour Force and Unemployment
unimp <- read.csv("SYB63_329_202009_Labour Force and Unemployment.csv",fileEncoding="UTF-8-BOM")%>%
  filter(Series=='Unemployment rate - Total'&Year==2015)%>%
  select(country.or.area,Country.code,Value)
#crime <- read.csv("SYB63_328_202009_Intentional Homicides and Other Crimes.csv",fileEncoding="UTF-8-BOM")%>%
#  filter(Series=='Intentional homicide rates per 100,000'&Year==2017)%>%
#  select(country.or.area,Country.code,Value)
# Current health expenditure (% of GDP)
health_exp <- read.csv("SYB63_325_202009_Expenditure on Health.csv",fileEncoding="UTF-8-BOM")%>%
  filter(Series=='Current health expenditure (% of GDP)'&Year==2015)%>%
  select(country.or.area,Country.code,Value)
# Gross domestic expenditure on R & D: as a percentage of GDP (%)
sci_exp <- read.csv("SYB63_286_202009_GDP on R&D.csv",fileEncoding="UTF-8-BOM")%>%
  filter(Series=='Gross domestic expenditure on R & D: as a percentage of GDP (%)'&Year==2015)%>%
  select(country.or.area,Country.code,Value)
# 'Safely managed drinking water sources, total (Proportion of population with access)'
water <- read.csv("SYB63_315_202009_Water and Sanitation Services.csv",fileEncoding="UTF-8-BOM")%>%
  filter(Series=='Safely managed drinking water sources, total (Proportion of population with access)'&Year==2017)%>%
  select(country.or.area,Country.code,Value)
# 'Life expectancy at birth for both sexes (years)'
life_exp <- read.csv("SYB62_246_201907_Population Growth Fertility and Mortality Indicators.csv",fileEncoding="UTF-8-BOM")%>%
  filter(Series=='Life expectancy at birth for both sexes (years)'&Year==2015)%>%
  select(country.or.area,Country.code,Value)
#sanitation <- read.csv("SYB63_315_202009_Water and Sanitation Services.csv",fileEncoding="UTF-8-BOM")%>%
#  filter(Series=='Safely managed sanitation facilities, total (Proportion of population with access)'&Year==2017)%>%
#  select(country.or.area,Country.code,Value)
# 'Arable land (% of total land area)'
land_arable <- read.csv("SYB63_145_202009_Land.csv",fileEncoding="UTF-8-BOM")%>%
  filter(Series=='Arable land (% of total land area)'&Year==2017)%>%
  select(country.or.area,Country.code,Value)

my_data <- list(deaths_age, pop_density,sex_ratio,mean_age_child,gdp_capita, gva_arg,
                cpi,unimp,health_exp,land_arable)
var_names <- c('deaths_age','pop_density','sex_ratio','mean_age_child','gdp_capita','gva_arg',
               'cpi','unimp','health_exp','land_arable')#'migr_ratio','gva_ind','gdp_growth','life_exp'
names(my_data) <- var_names
country_code <- lapply(my_data, `[[`, "Country.code")
common_country <- Reduce(intersect, country_code)
common_country <- sort(common_country)
mortality_data <- lapply(my_data, filter, Country.code%in%common_country)%>% lapply(arrange,Country.code)

# save predictors
save(mortality_data, file = paste(save_path, "/mortality_data.RData", sep=""))

#######################
### create mortality density
n <- length(common_country)
m <- nrow(my_data_common$deaths_age)/n
y <- as.matrix(my_data_common$deaths_age[["Age (x)"]][1:m])
for (i in common_country) {
  y <- cbind(y,filter(my_data_common$deaths_age, Country.code==i)[["Number of deaths d(x,n)"]])
  #y <- cbind(y,my_data_common$deaths_age[["Number of deaths d(x,n)"]][(1:m)+i*m])
}
rownames(y) <- y[,1]
y <- y[,-1]
ylist <- split(y,col(y)) 
colnames(y) <- common_country
## using frechet packege
x0 <-seq(0,100,length.out=101) #outputGrid
biny <- c(0,1,seq(5,100, by=5))
density <- lapply(ylist, FUN = CreateDensity, bin=biny, optns=list(outputGrid=x0),y = NULL, histogram = NULL)
density <- lapply(density, function(x) x <- x[-1])
names(density) <- unique(mortality_data$deaths_age$country.or.area)

# save response densities
save(density, file = paste(save_path, "/density.RData", sep=""))