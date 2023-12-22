library(rlist)
library(lubridate)

# load the parallelized RData files and consolidate into a single list of temperatures


rain_files = list.files("../../data/kenya_climate_data/", pattern = "rainfall_data_kenya*")

load(paste0("../../data/kenya_climate_data/", rain_files[[1]]))


index = 1
for(i in 1:length(result)){
  
  if(i == 1){
    result_new = list()
    result_new[[index]] = data.frame(result[[i]][[1]])
    index = index + 1
  } else if (length(result[[i]][[1]]) == 5){
    result_new[[index]] = data.frame(result[[i]][[1]])
    index = index + 1
  }
}

result_all = result_new

for(i in 2:length(rain_files)){
  
  print(i)
  
  load(paste0("../../data/kenya_climate_data/", rain_files[[i]]))
  index = 1
  for(j in 1:length(result)){
    
    if(j == 1){
      result_new = list()
      result_new[[index]] = data.frame(result[[j]][[1]])
      index = index + 1
    } else if (length(result[[j]][[1]]) == 5){
      result_new[[index]] = data.frame(result[[j]][[1]])
      index = index + 1
    }
  }
  
  result_all = c(result_all, result_new)
}

rainfall = result_all

length(rainfall)
# Check to make sure all values are numeric

num = rep(NA, length(rainfall))
for(i in 1:length(rainfall)){
  
  if(length(rainfall[[i]]) == 1){
    num[i] = 0
  } else{
    num[i] = 1
  }
}

# 99% of rasters have a value
mean(num)

#### temperature data ####

# load the parallelized RData files and consolidate into a single list of temperatures


temp_files = list.files("../../data/kenya_climate_data/", pattern = "temp_data_kenya*")

load(paste0("../../data/kenya_climate_data/", temp_files[[1]]))


index = 1
for(i in 1:length(result)){
  
  if(i == 1){
    result_new = list()
    result_new[[index]] = data.frame(result[[i]])
    index = index + 1
  } else if (length(result[[i]]) == 3){
    result_new[[index]] = data.frame(result[[i]])
    index = index + 1
  }
}

result_all = result_new

for(i in 2:length(temp_files)){
  
  print(i)
  
  load(paste0("../../data/kenya_climate_data/", temp_files[[i]]))
  index = 1
  for(j in 1:length(result)){
    
    if(j == 1){
      result_new = list()
    } 
    if (length(result[[j]]) == 3){
      result_new[[index]] = data.frame(result[[j]])
      index = index + 1
    }
  }
  
  result_all = c(result_all, result_new)
}

temp = result_all

length(temp)
# Check to make sure all values are numeric

num = rep(NA, length(temp))
for(i in 1:length(temp)){
  
  if(length(rainfall[[i]]) == 1){
    num[i] = 0
  } else{
    num[i] = 1
  }
}

# 99% of rasters have a value
mean(num)


#### consolidate and find rasters with both rainfall and temperature

lats_longs_temp = data.frame(x = NA, y = NA)

for(i in 1:length(temp)){
  if(length(temp[[i]]) == 3){
    lats_longs_temp = rbind(lats_longs_temp, head(temp[[i]][,c("x", "y")], 1))
  }
}

nrow(lats_longs_temp)
nrow(unique(lats_longs_temp))

lats_longs_rain = data.frame(lon = NA, lat = NA)

for(i in 1:length(rainfall)){
  if(length(rainfall[[i]]) == 5){
    lats_longs_rain = rbind(lats_longs_rain, head(rainfall[[i]][,c("lon", "lat")], 1))
  }
}

nrow(lats_longs_rain)
nrow(unique(lats_longs_rain))
lats_longs_rain = unique(lats_longs_rain)

lats_longs_temp = lats_longs_temp[-1,]

lats_longs_rain = lats_longs_rain[-1,]

head(lats_longs_temp)
head(lats_longs_rain)

load("../../data/kenya_data.RData")
still_needed_rain = kenya_data %>%
  anti_join(lats_longs_rain, by = c("x"="lon", "y"="lat")) %>%
  dplyr::select(x, y)

still_needed_temp = kenya_data %>%
  anti_join(lats_longs_temp, by = c("x", "y")) %>%
  dplyr::select(x, y)


# save(still_needed_rain, still_needed_temp, file = "../../data/kenya_climate_data/still_needed.RData")

# get latitude and longitude of rasters that have both rainfall and temperature
have_both = lats_longs_rain %>%
  rename(x = lon, y = lat) %>%
  semi_join(lats_longs_temp) %>%
  unique()

# get the rainfall and temperature lists in order according to this

retained_rainfall = list()
retained_temp = list()
for(i in 21903:nrow(have_both)){
  print(i)
  x_target = have_both[i,]$x
  y_target = have_both[i,]$y
  found_rain = F
  found_temp = F
  rain_ind = temp_ind = 1
  while(found_rain == F){
    if(length(rainfall[[rain_ind]]) != 5){
      rain_ind = rain_ind + 1
    } else if(unique(rainfall[[rain_ind]]$lon) == x_target &
              unique(rainfall[[rain_ind]]$lat) == y_target ){
      retained_rainfall = list.append(retained_rainfall, rainfall[[rain_ind]])
      found_rain = T
    } else{
      rain_ind = rain_ind + 1
    }
  }
  
  while(found_temp == F){
    if(length(temp[[temp_ind]]) != 3){
      temp_ind = temp_ind + 1
    } else if(unique(temp[[temp_ind]]$x) == x_target &
              unique(temp[[temp_ind]]$y) == y_target ){
      retained_temp = list.append(retained_temp, temp[[temp_ind]])
      found_temp = T
    } else{
      temp_ind = temp_ind + 1
    }
  }
}


# send out for the last 6 days of rainfall 
save(have_both, file = "../../data/kenya_climate_data/last_6_days_needed.RData")


# load the last 6 days of climate data
load("../../data/kenya_climate_data/last_6_days_rainfall_kenya.RData")

length(result)

for(i in 1:length(retained_rainfall)){
  retained_rainfall[[i]] = rbind(result[[i]], retained_rainfall[[i]])
}

kenya_rainfall = retained_rainfall
kenya_temps = retained_temp

# find the anomalies

rounded_dates = list()

for(i in 1:length(kenya_rainfall)){
  print(i)
  
  kenya_rainfall[[i]]$date = round_date(kenya_rainfall[[i]]$date, "week")
  
  rounded_dates[[i]] = kenya_rainfall[[i]]
  
}

weekly_sum = list()

for(i in 1:length(kenya_rainfall)){
  print(i)
  
  weekly_sum[[i]] = kenya_rainfall[[i]] %>%
    group_by(date) %>%
    summarize(total_rainfall = sum(chirps))
}

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

kenya_rainfall = rainfall_anomalies
kenya_temp = retained_temp

save(kenya_rainfall, kenya_temp, file = "../../data/kenya_climate_data/kenya_climate_data.RData")


