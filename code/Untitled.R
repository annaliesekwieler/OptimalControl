# load data
load("../data/kenya_data.RData")

# load climate data
load("../data/kenya_climate_data/kenya_climate_data.RData")

# load population data
load("../data/population_data_kenya_2020.RData")

# take out rasters that weren't able to be calibrated
to_keep = which(!is.na(kenya_data$bR))
kenya_data = kenya_data[to_keep,]
kenya_temp = kenya_temp[to_keep]
kenya_rainfall = kenya_rainfall[to_keep]

# source functions with model code
source("functions.R")
