# load calibration results
load("../results/calibrating_kenya.RData")

# load prevalence and incidence
load("../data/kenya_2020.RData")
load("../data/rasters_to_exclude.RData")
kenya_2020= kenya_2020[-rasters_to_exclude,]

# get a vector of calibrated prevalence and incidence
pred_prev = pred_inc =  rep(NA, length(result))

for(i in 1:length(result)){
  pred_prev[i] = as.numeric(result[[i]][2])
  pred_inc[i] = as.numeric(result[[i]][3])
}
