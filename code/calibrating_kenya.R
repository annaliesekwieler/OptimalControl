# calibrates the parameter for every raster in Kenya
library(doParallel)
library(dplyr)


#### load data ####

# load the prevalence, incidence, and proportion of cases that are treated
load("kenya_data.RData")


# load temperatures
load("kenya_climate_data.RData")

# kenya_temp = kenya_temp[indices]
# kenya_rainfall = kenya_rainfall[indices]

#### run model ####

# load functions
source("functions_SEISEI.R")

# parallelizing, perform the bR optimizer for each rasters

number_cores = 20
registerDoParallel(cores=number_cores)
cl=makeCluster(20)
result = foreach(i=1:nrow(kenya_data),.errorhandling = "pass") %dopar%  calibration_optimizer(temps = kenya_temp[[i]]$daily_temp,
                                                                                      rainfall = kenya_rainfall[[i]],
                                                                                      x = kenya_data[i, "x"],
                                                                                      y = kenya_data[i, "y"], 
                                                                                      prop_cases_treated = kenya_data[i, "prop_cases"],
                                                                                      prevalence_2019 = kenya_data[i,"prevalence_2019"],
                                                                                    incidence_2020 = kenya_data[i,"incidence_2020"]/365,
                                                                                      prevalence_2020 = kenya_data[i, "prevalence_2020"],
                                                                                    default_bR = 12)
stopCluster(cl)


save(result, file="calibrating_kenya.RData")
