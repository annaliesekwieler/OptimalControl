library(chirps)
library(tidyverse)
library(raster)

#### load proportion of cases treated data and transform to data frame #####

# load 2016 data
kenya_2016_prop_treated = raster("../data/Kenya_PropCases/202106_Global_Antimalarial_Effective_Treatment_KEN_2016.tiff")

# convert to a matrix
kenya_2016_prop_treated_matrix = rasterToPoints(kenya_2016_prop_treated)

# convert to data frame and give informative column names
kenya_2016_prop_treated = data.frame(kenya_2016_prop_treated_matrix)
names(kenya_2016_prop_treated) = c("x", "y", "prop_cases")

#### load the malaria incidence ####

# load 2016 data
kenya_2016_incidence = raster("../data/Kenya_Inc_1000/202206_Global_Pf_Incidence_Rate_KEN_2016.tiff")

# convert to a matrix
kenya_2016_incidence_matrix = rasterToPoints(kenya_2016_incidence)

# convert to data frame and give informative column names
kenya_2016_incidence = data.frame(kenya_2016_incidence_matrix)
names(kenya_2016_incidence) = c("x", "y", "incidence_1000")



#### prevalence data ####

# load 2016 data
kenya_2016_prevalence = raster("../data/Kenya_Prevalence/202206_Global_Pf_Parasite_Rate_KEN_2016.tiff")

# convert to a matrix
kenya_2016_prevalence_matrix = rasterToPoints(kenya_2016_prevalence)

# convert to data frame and give informative column names
kenya_2016_prevalence = data.frame(kenya_2016_prevalence_matrix)
names(kenya_2016_prevalence) = c("x", "y", "prevalence")

head(kenya_2016_prevalence)



#### combine data frames  ####

kenya_2016_prop_treated = kenya_2016_prop_treated %>%
  mutate(x = round(x, 2),
         y = round(y, 2))


kenya_2016 = kenya_2016_incidence %>%
  full_join(kenya_2016_prevalence, by=c("x", "y")) %>%
  mutate(x = round(x, 2),
         y = round(y, 2)) %>%
  full_join(kenya_2016_prop_treated, by=c("x", "y"))

kenya_2016 = kenya_2016[complete.cases(kenya_2016),]

#### get temperature data from 2016 ####


# first 5 rasters
load("../results/Kenya_2016_Temperatures_First365.RData")
# 
# # store the daily temperature patterns here
# # kenya_2016_temperatures = list()
# 
# # for(i in 1:nrow(kenya_2016)){
# for(i in 301:365){
#   print(i)
#   lonlat_kenya = kenya_2016[i, c("x", "y")]
#   dates_2016 = seq(from=as.Date("2016-01-01"),  to=as.Date("2016-12-31"), by=1)
# 
#   min_temp = rep(NA, length(dates_2016))
# 
#   for(j in 1:(length(dates_2016)-1)){
#     min_temp[j] = get_chirts(lonlat_kenya, dates = seq(from = as.Date(dates_2016[j]), to = as.Date(dates_2016[j]) + 1,
#                                                                                                   by=1),
#                              server="CHC", var="Tmin")$chirts[1]
#   }
# 
#   min_temp = min_temp[1:365]
# 
#   max_temp = rep(NA, length(dates_2016))
# 
#   for(j in 1:(length(dates_2016)-1)){
#     max_temp[j] = get_chirts(lonlat_kenya, dates = seq(from = as.Date(dates_2016[j]), to = as.Date(dates_2016[j]) + 1,
#                                                        by=1),
#                              server="CHC", var="Tmax")$chirts[1]
#   }
# 
#   max_temp = max_temp[1:365]
# 
#   average_temp = (min_temp + max_temp) / 2
# 
#   kenya_2016_temperatures[[i]] = data.frame(x = lonlat_kenya$x, y = lonlat_kenya$y, daily_temp = average_temp)
# 
# 
# }

#### load functionsand set default parameter values ####

source("../code/functions.R")


# eggs laid per gonotrophic cycle
epsilon = 50

# egg-adult sensitivity to surface moisture
a.R = .5
# standardizd anomaly of cumulative weekly rainfall- NO DEFAULT, IS DATA
c.RD = 1
# mean anomaly level at which 0.5 survival is expected
b.R = 0 

PHM = 0.125

#### population data ####

# load Un population data
population_data = read.csv("../data/UN_Pop_Data.csv")

# convert to numeric data type and take birth and death rates per capita

as.numeric(population_data$Birth_Rate)
population_data_kenya_2016 = population_data %>%
  filter(Year == 2016, Region == "Kenya") %>%
  mutate(Birth_Rate = as.numeric(Birth_Rate),
         Death_Rate = as.numeric(Death_Rate)) %>%
  mutate(birth_rate = Birth_Rate / 1000,
         death_rate = Death_Rate / 1000) %>%
  dplyr::select(Year, Region, birth_rate, death_rate)

#### load 2015 prevalences and prop treated to begin the model ####

# load previous year's malaria prevalence to initiate
kenya_2015_prevalence = raster("../data/Kenya_Prevalence/202206_Global_Pf_Parasite_Rate_KEN_2015.tiff")
kenya_2015_prevalence_matrix = rasterToPoints(kenya_2015_prevalence)
kenya_2015_prevalence = data.frame(kenya_2015_prevalence_matrix)
names(kenya_2015_prevalence) = c("x", "y", "prevalence")

# load previous year's proportion of cases treated
kenya_2015_prop_treated = raster("../data/Kenya_PropCases/202106_Global_Antimalarial_Effective_Treatment_KEN_2015.tiff")
kenya_2015_prop_treated_matrix = rasterToPoints(kenya_2015_prop_treated)
kenya_2015_prop_treated = data.frame(kenya_2015_prop_treated_matrix)
names(kenya_2015_prop_treated) = c("x", "y", "prop_cases")

#### make options for location-varying parameters ####

PMH_minimum = 1e-20
PMH_maximum = 1

rs_minimum = 1/3000
rs_maximum = 1/5

mu_m_minimum = 1/100
mu_m_maximum = 1/5

#### run the model on one raster ####

SEI_model = function(time, y, params){
  # get rid of width , unpack parameters from vector instead of 
  # reference by name
  
  # drop DeSolve- use odein (much faster)
  
  with(as.list(c(y, params)), {
    
    X = Ic + Is
    
    temp = temps[floor(time)]
    
    if(temp > 31.7){
      temp=31.7
    }
    
    FOI = mu_SE(temp, X, PHM, epsilon, mu_m, PMH)
    
    # rate of change
    dS = -FOI*S - deltaT*S + mut + rc*Ic + rs*Is
    dE = FOI*S - tau*E - deltaT*E
    dIc = -rc*Ic + c * tau * E - deltaT*Ic
    dIs = -rs * Is + (1-c) * tau * E - deltaT*Is
    
    return(list(c(dS, dE, dIc, dIs)))
    
  })
  
}

# params is a vector c(rs, PHM, PMH, mu_m)
# index is the raster we are referring to
predicted_prevalence = function(params, index){
  
  state = c(S = (1 - (.001 + kenya_2015_prevalence$prevalence[index]*
                        kenya_2015_prop_treated$prop_cases[index] + 
                        kenya_2015_prevalence$prevalence[index]*
                        (1-kenya_2015_prop_treated$prop_cases[index]))),
            E = .001, Ic = kenya_2015_prevalence$prevalence[index]*
              kenya_2015_prop_treated$prop_cases[index], 
            Is = kenya_2015_prevalence$prevalence[index]*
              (1-kenya_2015_prop_treated$prop_cases[index]))
  
  times = seq(1, 365, by =1/8)
  
  parms = c(
    rc = 1/11,
    rs = params[1],
    tau = 1/10,
    c = kenya_2016$prop_cases[index],
    deltaT = population_data_kenya_2016$death_rate / 365,
    mut = population_data_kenya_2016$birth_rate / 365,
    PHM = PHM,
    PMH = params[2],
    epsilon = epsilon,
    mu_m = params[3]
  )
  
  temps = kenya_2016_temperatures[[index]]$daily_temp
  out = data.frame(ode(y = state, times = times, func = SEI_model,
                       parms = parms))
  
  return(mean(out$Ic + out$Is))
  
}
predicted_prevalence(params = c(.1, .5, .1), index = 2)


##### calibrate each raster to the true incidence and FOI (incidence) ####

outputs_given_params = list()
for(i in 1:100){
  print(i)
  outputs_given_params[[i]] = predicted_prevalence(c(.1, .5, .5, .1), index = i)
}


outputs_given_params = unlist(outputs_given_params)

# calibrate for a single raster

# mean squared error function
mse_prev_inc = function(params, index){
  model_prev = predicted_prevalence(params, index)
  
  temp = mean(kenya_2016_temperatures[[index]]$daily_temp)
  
  model_inc = mu_SE(temp, X = model_prev, PHM , epsilon, mu_m, PMH)
}

optim()

