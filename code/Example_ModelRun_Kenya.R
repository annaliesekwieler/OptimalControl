# This script shows how the SEIT model can be run on each raster in Kenya

#### load packages ####
library(raster)
library(tidyverse)
library(chirps)
library(zoo)
library(deSolve)
library(terra)
library(lubridate)

#### load proportion of cases treated data and transform to data frame #####

# load 2016 data
kenya_2020_prop_treated = raster("../data/Kenya_PropCases/202106_Global_Antimalarial_Effective_Treatment_KEN_2020.tiff")

# convert to a matrix
kenya_2020_prop_treated_matrix = rasterToPoints(kenya_2020_prop_treated)

# convert to data frame and give informative column names
kenya_2020_prop_treated = data.frame(kenya_2020_prop_treated_matrix)
names(kenya_2020_prop_treated) = c("x", "y", "prop_cases")

#### load the malaria incidence ####

# load 2016 data
kenya_2020_incidence = raster("../data/Kenya_Inc_1000/202206_Global_Pf_Incidence_Rate_KEN_2020.tiff")

# convert to a matrix
kenya_2020_incidence_matrix = rasterToPoints(kenya_2020_incidence)

# convert to data frame and give informative column names
kenya_2020_incidence = data.frame(kenya_2020_incidence_matrix)
names(kenya_2020_incidence) = c("x", "y", "incidence_1000")



#### prevalence data ####

# load 2016 data
kenya_2020_prevalence = raster("../data/Kenya_Prevalence/202206_Global_Pf_Parasite_Rate_KEN_2020.tiff")

# convert to a matrix
kenya_2020_prevalence_matrix = rasterToPoints(kenya_2020_prevalence)

# convert to data frame and give informative column names
kenya_2020_prevalence = data.frame(kenya_2020_prevalence_matrix)
names(kenya_2020_prevalence) = c("x", "y", "prevalence")

head(kenya_2020_prevalence)



#### combine data frames  ####

kenya_2020_prop_treated = kenya_2020_prop_treated %>%
  mutate(x = round(x, 2),
         y = round(y, 2))
  

kenya_2020 = kenya_2020_incidence %>%
  full_join(kenya_2020_prevalence, by=c("x", "y")) %>%
  mutate(x = round(x, 2),
         y = round(y, 2)) %>%
  full_join(kenya_2020_prop_treated, by=c("x", "y"))

kenya_2020 = kenya_2020[complete.cases(kenya_2020),]

# take only the first row
single_raster = kenya_2020[1,]

save(kenya_2020, file="../data/kenya_2020.RData")

#### monthly temperature averages ####


# CAN'T GET THIS TO WORK
# monthly_temperature = worldclim_tile("tavg", lat = single_raster$y, lon = single_raster$x, path=".")
# monthly_temperature_averages = as.data.frame(monthly_temperature)
# monthly_temperature[["names"]]

lonlat_kenya = kenya_2020[1, c("y", "x")]
dates_2016 = seq(from=as.Date("2016-01-01"),  to=as.Date("2016-12-31"), by=1)

load("../results/Kenya_2016_Temperatures_First365.RData")
temperature_list = kenya_2016_temperatures[[1]]


#### download rainfall data ####



load("../results/100_Rasters_Rainfall.RData")
rainfall_list = kenya_2020_rainfall[[1]]$chirps

# also get the rainfall for the past 6 days
last_days_2019 = get_chirps(lonlat_kenya, dates = c("2019-12-26", "2019-12-31"), server="CHC")$chirps
rainfall_list = c(last_days_2019, rainfall_list)


# get the historic rainfall data 
# dat = get_chirps(lonlat_kenya, dates = c("2005-01-01", "2019-12-31"),
#            server="CHC")
load(file="../results/historical_rainfall_first_raster.RData")
# convert to weekly estimates
dat$date= as.Date(dat$date)

# round the day down to the week
dat$date = round_date(dat$date, "week")
  
weekly_total_rainfall = dat %>%
  group_by(date) %>%
  summarize(total_rainfall = sum(chirps))

# mean and sd of rainfall
weekly_mean_rainfall = mean(weekly_total_rainfall$total_rainfall)
weekly_anomaly_rainfall = sd(weekly_total_rainfall$total_rainfall)


# find the anomalies for rainfall data by taking the past 7 days and finding the anomaly
rainfall_anomalies = rep(NA, 365)
for(i in 7:371){
  # sequence of past 7 days
  rainfall_anomalies[i-6] = (sum(rainfall_list[(i-6):i]) -weekly_mean_rainfall) / 
    weekly_anomaly_rainfall

}

#### support functions for model ####

# length of the gonotrophic cycle
GP = function(temp){
  1 / (.017 * temp - .165)
}

# eggs produced over the lifespan
B = function(temp, epsilon, mu_m){
  epsilon / (exp(GP(temp) * mu_m) - 1)
}


# egg-to-adult survivorship due to rainfall
PEAr = function(R, a.R, c.RD, b.R){
  1 / (1 + exp(-a.R*(c.RD-b.R)))
}

# water temperature relative to air temperature
Tw = function(temp){
  return(temp + 2)
}

# rate of immature survival
# temp and cRD are data
PEAt_PEAr = function(temp, cRD, aR, bR){
  TW = Tw(temp)
  PEAt = -.00924*(TW^2) + .453 * TW - 4.77

  PEAr = 1 / (1 + exp(-aR * (cRD - bR)))
  
  return(PEAt * PEAr)

}

# duration of egg-adult stage
tau_EA = function(temp){
  TW = Tw(temp)
  TW[TW > 34] = 34
  return(((.000111*TW*(TW - 14.7)) * sqrt(34 - TW))^(-1))
}

# mosquito biting rate
a = function(temp){
  0.017*temp - .165
}

# Sporogony
n = function(temp){
  (.0001122*temp*(temp - 15.384) * sqrt(35-temp))^(-1)
}

L = function(temp, epsilon, mu_m, cRD, aR, bR){
  B(temp, epsilon, mu_m) * PEAt_PEAr(temp, cRD, aR, bR) / tau_EA(temp)
}

# mosquito density
m = function(temp, epsilon, mu_m, cRD, aR, bR){
  L(temp, epsilon, mu_m, cRD, aR, bR) / mu_m
}


# primary function
EIR = function(temp, X, PHM, epsilon, mu_m, cRD, aR, bR){
  numerator = m(temp, epsilon, mu_m, cRD, aR, bR) * (a(temp)) ^ 2 * PHM * X * exp(-mu_m * n(temp))
  denominator = mu_m + a(temp) * PHM *X
  return(numerator / denominator)
}

# convert to FOI
mu_SE = function(temp, X, PHM, epsilon, mu_m, cRD, aR, bR, PMH){
  1 - exp(-PMH * EIR(temp, X, PHM, epsilon, mu_m, cRD, aR, bR))
}


#### code for model ####


SEI_model = odin::odin({
  
  # Derivatives 
  deriv(S) <- -FOI*S - deltaT*S + mut + rc*Ic + rs*Is
  deriv(E) <- FOI*S - tau*E - deltaT*E
  deriv(Ic) <- -rc*Ic + c * tau * E - deltaT*Ic
  deriv(Is) <- -rs * Is + (1-c) * tau * E - deltaT*Is
  
  # Initial conditions
  initial(S) <- S_init
  initial(E) <- E_init
  initial(Ic) <- Ic_init
  initial(Is) <- Is_init
  
  # make a spline of force of infection
  FOI = interpolate(tvec_spline, foi_in, 'spline')
  
  tvec_spline[] = user()
  foi_in[] = user()
  
  dim(tvec_spline) = 365
  dim(foi_in) = 365
  
  ## parameters
  deltaT = user()
  mut = user()
  rc = user()
  rs = user()
  tau = user()
  c = user()
  
  # Initial conditions
  S_init = user()
  E_init = user()
  Ic_init = user()
  Is_init = user()
  
})


#### run the model on a single location ####



# load previous year's malaria prevalence to initiate
kenya_2019_prevalence = raster("../data/Kenya_Prevalence/202206_Global_Pf_Parasite_Rate_KEN_2019.tiff")
kenya_2019_prevalence_matrix = rasterToPoints(kenya_2019_prevalence)
kenya_2019_prevalence = data.frame(kenya_2019_prevalence_matrix)
names(kenya_2019_prevalence) = c("x", "y", "prevalence")

# load previous year's proportion of cases treated
kenya_2019_prop_treated = raster("../data/Kenya_PropCases/202106_Global_Antimalarial_Effective_Treatment_KEN_2019.tiff")
kenya_2019_prop_treated_matrix = rasterToPoints(kenya_2019_prop_treated)
kenya_2019_prop_treated = data.frame(kenya_2019_prop_treated_matrix)
names(kenya_2019_prop_treated) = c("x", "y", "prop_cases")

state = c(S = (1 - (.001 + kenya_2019_prevalence$prevalence[1]*
                      kenya_2019_prop_treated$prop_cases[1] + 
                      kenya_2019_prevalence$prevalence[1]*
                      (1-kenya_2019_prop_treated$prop_cases[1]))), E = .001,
          Ic = kenya_2019_prevalence$prevalence[1]*
            kenya_2019_prop_treated$prop_cases[1], Is = kenya_2019_prevalence$prevalence[1]*
            (1-kenya_2019_prop_treated$prop_cases[1]))

times = seq(1, 365, by =1/8)

temps = kenya_2016_temperatures[[1]]$daily_temp

rainfall = rainfall_anomalies

#### set default parameter values #####

# eggs laid per gonotrophic cycle
epsilon = 50
# adult mosquito mortality (no intervention)
mu_m = -log(.98)

# epidemiological parameters

# transmission probabilities
PHM = 0.125
PMH = 0.5

# incubation rate in humans
tau = 1/10


# cumulative rainfall parameters
aR = 0.5

# recovery rates
rc = 1/5
rs = 1 / 100

##### make the force of infection ####


# force of infection for the full year
foi_in_vec = mu_SE(temps, X = kenya_2019_prevalence$prevalence[[1]], PHM = PHM, epsilon=epsilon, mu_m=mu_m,
                   rainfall, aR=aR, bR = -1, PMH=PMH)

t = seq(1, 365, by =1/8)
model = SEI_model$new(S_init = .95,
                      E_init=.01,
                      Ic_init = .02,
                      Is_init = .02,
                      deltaT = population_data_kenya_2020$death_rate / 365,
                      mut = population_data_kenya_2020$birth_rate / 365,
                      rc = rc,
                      rs = rs,
                      tau = tau,
                      c = kenya_2020$prop_cases[1],
                      foi_in = foi_in_vec,
                      tvec_spline = seq(1, 365))
out = model$run(t)
plot(out)

out=as.data.frame(out)
#### time to calibrate ####

true_prevalence = kenya_2020$prevalence[1]
model_prevalence = mean(out$Ic + out$Is)

predicted_prevalence = function(bR){
  
  foi_in_vec = mu_SE(temps, X = kenya_2019_prevalence$prevalence[[1]], PHM, epsilon, mu_m, rainfall, aR, bR=bR, PMH)
  
  model = SEI_model$new(S_init = .95,
                        E_init=.01,
                        Ic_init = .02,
                        Is_init = .02,
                        deltaT = population_data_kenya_2020$death_rate / 365,
                        mut = population_data_kenya_2020$birth_rate / 365,
                        rc = rc,
                        rs = rs,
                        tau = tau,
                        c = kenya_2020$prop_cases[1],
                        foi_in = foi_in_vec,
                        tvec_spline = seq(1, 365))
  
  out = model$run(t)
  
  out=as.data.frame(out)
  
  return(mean(out$Ic + out$Is))
  
}

# optimize the bR parameter

optimizing_function = function(parameter){
  bR = parameter[1]
  abs(predicted_prevalence(bR) - true_prevalence)
}

optimized_result = optim(par = c(0), optimizing_function, method="L-BFGS-B",
      lower=-100, upper = 100)





