
library(raster)
library(chirps)
library(dplyr)
#### download raster info ####


# load 2020 data
kenya_2020_prop_treated = raster("../data/Kenya_PropCases/202106_Global_Antimalarial_Effective_Treatment_KEN_2020.tiff")

# convert to a matrix
kenya_2020_prop_treated_matrix = rasterToPoints(kenya_2020_prop_treated)

# convert to data frame and give informative column names
kenya_2020_prop_treated = data.frame(kenya_2020_prop_treated_matrix)
names(kenya_2020_prop_treated) = c("x", "y", "prop_cases")

#### load the malaria incidence ####

# load 2020 data
kenya_2020_incidence = raster("../data/Kenya_Inc_1000/202206_Global_Pf_Incidence_Rate_KEN_2020.tiff")

# convert to a matrix
kenya_2020_incidence_matrix = rasterToPoints(kenya_2020_incidence)

# convert to data frame and give informative column names
kenya_2020_incidence = data.frame(kenya_2020_incidence_matrix)
names(kenya_2020_incidence) = c("x", "y", "incidence_1000")



#### prevalence data ####

# load 2020 data
kenya_2020_prevalence = raster("../data/Kenya_Prevalence/202206_Global_Pf_Parasite_Rate_KEN_2020.tiff")

# convert to a matrix
kenya_2020_prevalence_matrix = rasterToPoints(kenya_2020_prevalence)

# convert to data frame and give informative column names
kenya_2020_prevalence = data.frame(kenya_2020_prevalence_matrix)
names(kenya_2020_prevalence) = c("x", "y", "prevalence")

head(kenya_2020_prevalence)



#### combine data frames  ####

kenya_2020_prop_treated = kenya_2020_prop_treated %>%
  mutate(x = round(x, 2),
         y = round(y, 2))


kenya_2020 = kenya_2020_incidence %>%
  full_join(kenya_2020_prevalence, by=c("x", "y")) %>%
  mutate(x = round(x, 2),
         y = round(y, 2)) %>%
  full_join(kenya_2020_prop_treated, by=c("x", "y"))

kenya_2020 = kenya_2020[complete.cases(kenya_2020),]

# #### get rainfall ####
# 
# store the daily temperature patterns here
# kenya_2020_rainfall = list()
load("../results/Rainfall_Kenya_First190.RData")
# for(i in 1:nrow(kenya_2016)){
for(i in 339:365){
  print(i)
  lonlat_kenya = kenya_2020[i, c("x", "y")]


  kenya_2020_rainfall[[i]] = get_chirps(lonlat_kenya, dates = c("2020-01-01", "2020-12-31"),
                                        server="CHC")

}

save(kenya_2020_rainfall, file="../results/100_Rasters_Rainfall.RData")
# 
# 


#### worldclim data ####




