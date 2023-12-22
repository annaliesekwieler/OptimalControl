load("kenya_2020.RData")

library(chirps)


kenya_2016_temperatures = list()
# 

dates_2016 = seq(from=as.Date("2016-01-01"),  to=as.Date("2016-12-31"), by=1)

for(i in 1:nrow(kenya_2020)){
  
  lonlat_kenya = kenya_2020[i, c("x", "y")]
  
  
  min_temp = max_temp = rep(NA, length(dates_2016))

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

  kenya_2016_temperatures[[i]] = data.frame(x = lonlat_kenya$x, y = lonlat_kenya$y, daily_temp = average_temp)

  save(kenya_2016_temperatures, file="kenya_2016_temperatures.RData")

}