library(tidyverse)


# load results from calibration and true prevalence/incidence

# calibration results
load("../results/calibrating_kenya.RData")

# make a vector of constants
predicted_prevalence = rep(NA, length(result))

bR_vals = convergence_vals = x_vals = y_vals = rep(NA, length(result))

for(i in 1:length(result)){

  if(length(result[[i]]) == 5){
    if(is.numeric(result[[i]][1][[1]])){
      predicted_prevalence[i] = result[[i]][1][[1]]
    }
    
 
    
    if(is.numeric(result[[i]][[2]][[1]])){
      bR_vals[i] = result[[i]][[2]][[1]]
    }
    
    
    
    if(is.numeric(result[[i]][[3]][[1]])){
      convergence_vals[i] = result[[i]][[3]][[1]]
    }
    
    if(is.numeric(result[[i]][[4]][[1]])){
      x_vals[i] = result[[i]][[4]][[1]]
    }
    
    if(is.numeric(result[[i]][[5]][[1]])){
      y_vals[i] = result[[i]][[5]][[1]]
    }
  }
  
  
  
  
}

# load the prevalence, incidence, and proportion of cases that are treated
load("../data/kenya_data.RData")


# load temperatures
load("../data/kenya_climate_data/kenya_climate_data.RData")

par(mfrow=c(1,2), mar=c(4,4,4,4))
# compare true prevalence and calibrated
plot(predicted_prevalence ~ 
       kenya_data$prevalence_2020[is.na(kenya_data$bR)],
     xlab="True Prevalence",
     ylab="Predicted Prevalence", xlim=c(0, 0.7), ylim=c(0, 0.7))
abline(a=0, b=1, lty="dashed")
abline(a=0, b=1, lty="dashed")



# find all of those that are within .01 of the true prevalence
ind_br = which(abs(predicted_prevalence - 
                     kenya_data$prevalence_2020[is.na(kenya_data$bR)]) < .01)

# mark those that don't currently have anything in the bR column
kenya_data$new_bR = NULL
kenya_data = kenya_data %>% 
  dplyr::left_join(data.frame(x = x_vals, y = y_vals, new_bR = bR_vals)[ind_br,], by = c("x", "y")) %>%
  mutate(bR = ifelse(is.na(bR), new_bR, bR))




save(kenya_data, file = "../data/kenya_data.RData")


# ensure that this worked
predicted_prevalence = rep(NA, nrow(kenya_data))
for(i in 1:nrow(kenya_data)){
  print(i)
  if(!is.na(kenya_data$bR[i])){
    predicted_prevalence[i] = pred_prevalence(temps = kenya_temp[[i]]$daily_temp,
                                              rainfall = kenya_rainfall[[i]],
                                              prop_cases_treated = kenya_data[i, "prop_cases"],
                                              prevalence_2019 = kenya_data[i,"prevalence_2019"],
                                              incidence_2020 = kenya_data[i,"incidence_2020"]/365,
                                              prevalence_2020 = kenya_data[i, "prevalence_2020"],
                                              bR = kenya_data$bR[i])
  }
 
}

plot(predicted_prevalence ~ kenya_data$prevalence_2020)
save(kenya_data, file = "../data/kenya_data.RData")


which(abs(predicted_prevalence - kenya_data$prevalence_2020) < .1)

plot(kenya_data$incidence_2020 / 365, pred_incidence, xlab="True incidence",
     ylab="Predicted incidence", xlim=c(0, .02), ylim=c(0,.02))
abline(a=0, b=1, lty="dashed")

# load data and code necessary to run
source("functions.R")
# load the prevalence, incidence, and proportion of cases that are treated


# plot the values of the parameters
par(mfrow=c(2,2))
hist(bR_vals, main="bR")
abline(v=0, col="red")

hist(epsilon_vals, main="epsilon")
abline(v=50, col="red")

hist(factor_vals, main="scaling factor")
abline(v=.5, col="red")

par(mfrow=c(2,2))
plot(bR_vals)
plot(epsilon_vals)
plot(factor_vals)


# find the rasters that have predicted values close to the true value

try_again = which(abs(pred_prevalence - kenya_data$prevalence_2020) > .02)

save(try_again, file = "../data/try_again_cal_kenya.RData")
