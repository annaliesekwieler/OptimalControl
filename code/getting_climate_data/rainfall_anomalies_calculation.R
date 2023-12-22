##### Clean rainfall data by removing NA and calculating the rainfall anomalies

library(lubridate)
library(tidyverse)

#### Load data ####

# load the rainfall for 2020 and the last 6 days of 2019

# 2020
load("kenya_rainfall.RData")
rainfall = result_all

# 2019
load("last_6_days_rainfall_kenya.RData")
last_6_days = result

#### Exclude NA rasters ####

# rasters to exclude due to lack of rainfall data
load("rasters_to_exclude.RData")
rainfall = rainfall[-rasters_to_exclude]
last_6_days = last_6_days[-rasters_to_exclude]

#### Change to rainfall deviations ####

# combine these into a single list
kenya_rainfall = list()

for(i in 1:length(rainfall)){
  print(i)
  
  rainfall_vector = rainfall[[i]]
  last_6_days_vector = last_6_days[[i]]
  
  kenya_rainfall[[i]] = rbind(last_6_days_vector, rainfall_vector)
}

# round date to the week
rounded_dates = list()

for(i in 1:length(kenya_rainfall)){
  print(i)
  
  kenya_rainfall[[i]]$date = round_date(kenya_rainfall[[i]]$date, "week")
  
  rounded_dates[[i]] = kenya_rainfall[[i]]
  
}

# calculate the weekly total rainfall with a moving sum

weekly_sum = list()

for(i in 13377:length(kenya_rainfall)){
  print(i)
  
  weekly_sum[[i]] = kenya_rainfall[[i]] %>%
    group_by(date) %>%
    summarize(total_rainfall = sum(chirps))
  
  
}

# get the weekly mean and st of rainfall

weekly_mean_rainfall = weekly_anomaly_rainfall =  list()

for(i in 1:length(weekly_sum)){
  
  weekly_mean_rainfall[[i]]  = mean(weekly_sum[[i]]$total_rainfall)
  
  weekly_anomaly_rainfall[[i]]  = sd(weekly_sum[[i]]$total_rainfall)
  
  
}

# convert the daily rainfall into a weekly anomaly by taking the past 6 days
rainfall_anomalies = list()

for(j in 1:length(kenya_rainfall)){
  print(j)
  
  rainfall_anomaly = rep(NA, 365)
  
  for(i in 7:371){
    # sequence of past 7 days
    rainfall_anomaly[i-6] = (sum(kenya_rainfall[[j]]$chirps[(i-6):i]) -weekly_mean_rainfall[[j]]) / 
      weekly_anomaly_rainfall[[j]]
    
  }
  
  rainfall_anomalies[[j]] = rainfall_anomaly
}


save(rainfall_anomalies, file="../../data/kenya_climate_data/rainfall_anomalies.RData")

