library(chirps)
library(doParallel)



load("last_6_days_needed.RData")


get_rainfall = function(row){
  
  lonlat_kenya = have_both[row, c("x", "y")]
  
  return(get_chirps(lonlat_kenya, dates = c("2019-12-26", "2019-12-31"),
                    server="CHC"))
  
}



number_cores = 20
registerDoParallel(cores=number_cores)
cl=makeCluster(20)
result = foreach(i=1:nrow(have_both),.errorhandling = "pass") %dopar% get_rainfall(i)

save(result, file="last_6_days_rainfall_kenya.RData")