args = commandArgs(trailingOnly=TRUE)
j = as.numeric(args[1]) 


library(chirps)
library(doParallel)
library(tidyverse)



load("kenya_data.RData")
load("still_needed.RData")


kenya_data = kenya_data %>%
  semi_join(still_needed_rain)


# dates_2016 = seq(from=as.Date("2016-01-01"),  to=as.Date("2016-12-31"), by=1)


get_rainfall = function(row){
  
  lonlat_kenya = kenya_data[row, c("x", "y")]
  
  return(list(get_chirps(lonlat_kenya, dates = c("2020-01-01",
                                                 "2020-12-31"),
                         server="CHC"), lonlat_kenya))
  
}


number_cores = 20
registerDoParallel(cores=number_cores)
cl=makeCluster(20)

beginning = 1 + 100*(j-1)

end = ifelse(j == 13, 1217, (100 * j))


result = foreach(i=beginning:end,.errorhandling = "pass") %dopar% get_rainfall(i)

file_name = paste("rainfall_data_kenya_", j, ".RData", sep="")

save(result, file=file_name)