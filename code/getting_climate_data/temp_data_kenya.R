args = commandArgs(trailingOnly=TRUE)
j = as.numeric(args[1]) 

library(chirps)
library(doParallel)
library(tidyverse)




load("kenya_data.RData")
load("still_needed.RData")


kenya_data = kenya_data %>%
  semi_join(still_needed_temp)
  




kenya_2016_temperatures = list()
# 

dates_2016 = seq(from=as.Date("2016-01-01"),  to=as.Date("2016-12-31"), by=1)


get_temperature = function(row){
  
  min_temp = max_temp = rep(NA, length(dates_2016))
  
  lonlat_kenya = kenya_data[row, c("x", "y")]
  
  
  for(j in 1:(length(dates_2016)-1)){
    min_temp[j] = get_chirts(lonlat_kenya, dates = seq(from = as.Date(dates_2016[j]), to = as.Date(dates_2016[j]) + 1,
                                                       by=1),
                             server="CHC", var="Tmin")$chirts[1]
  }
  
  min_temp = min_temp[1:365]
  
  for(j in 1:(length(dates_2016)-1)){
    max_temp[j] = get_chirts(lonlat_kenya, dates = seq(from = as.Date(dates_2016[j]), to = as.Date(dates_2016[j]) + 1,
                                                       by=1),
                             server="CHC", var="Tmax")$chirts[1]
  }
  
  
  max_temp = max_temp[1:365]
  
  average_temp = (min_temp + max_temp) / 2
  
  return(data.frame(x = lonlat_kenya$x, y = lonlat_kenya$y, daily_temp = average_temp))
  
}



number_cores = 5
registerDoParallel(cores=number_cores)
cl=makeCluster(5)

beginning = 1

end = 5

result = foreach(i=beginning:end,.errorhandling = "pass") %dopar% get_temperature(i)

file_name = paste("temp_data_kenya_", j, ".RData", sep="")

save(result, file=file_name)