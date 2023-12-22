# calibrates the parameter for every raster in Kenya
library(doParallel)
library(dplyr)


#### load data ####

# load the prevalence, incidence, and proportion of cases that are treated
load("kenya_data.RData")

# load kenya population data
load("population_data_kenya_2020.RData")

# load temperatures
load("kenya_climate_data.RData")

load("try_again_cal_kenya.RData")

kenya_data = kenya_data[try_again,]
kenya_temp= kenya_temp[try_again]
kenya_rainfall= kenya_rainfall[try_again]



#### run model ####

# load functions
source("functions.R")

# parallelizing, perform the bR optimizer for each rasters

number_cores = 20
registerDoParallel(cores=number_cores)
cl=makeCluster(20)
result = foreach(i=1:nrow(kenya_data),.errorhandling = "pass") %dopar% calibration_optimizer(temps = kenya_temp[[i]]$daily_temp,
                                                                                             rainfall = kenya_rainfall[[i]],
                                                                                             prop_cases_treated = kenya_data[i, "prop_cases"],
                                                                                             prevalence_2019 = kenya_data[i,"prevalence_2019"],
                                                                                             incidence_2020 = kenya_data[i,"incidence_2020"]/365,
                                                                                             prevalence_2020 = kenya_data[i, "prevalence_2020"])


save(result, file="calibrating_kenya_2.RData")


# result = list()
# 
# for(i in 1:nrow(kenya_data)) {
#   print(i)
#   result[[i]] = calibration_optimizer(temps = kenya_temp[[i]]$daily_temp,
#                         rainfall = kenya_rainfall[[i]],
#                         prop_cases_treated = kenya_data[i, "prop_cases"],
#                         prevalence_2019 = kenya_data[i,"prevalence_2019"],
#                         incidence_2020 = (kenya_data[i,"incidence_2020"] / 365),
#                         prevalence_2020 = kenya_data[i, "prevalence_2020"])
# }
