library(raster)


# load Un population data
population_data = read.csv("../data/UN_Pop_Data.csv")

# convert to numeric data type and take birth and death rates per capita
population_data_kenya_2020 = population_data %>%
  filter(Year == 2020, Region == "Kenya") %>%
  mutate(Birth_Rate = as.numeric(Birth_Rate),
         Death_Rate = as.numeric(Death_Rate)) %>%
  mutate(birth_rate = Birth_Rate / 1000,
         death_rate = Death_Rate / 1000) %>%
  dplyr::select(Year, Region, birth_rate, death_rate)
save(population_data_kenya_2020, file = "../data/population_data_kenya_2020.RData")


##### incidence 2020 ####

kenya_2020_incidence = raster("../data/Kenya_Incidence/202206_Global_Pf_Incidence_Rate_KEN_2020.tiff")

kenya_2020_incidence_matrix = rasterToPoints(kenya_2020_incidence)

kenya_2020_incidence = data.frame(kenya_2020_incidence_matrix)
names(kenya_2020_incidence) = c("x", "y", "incidence")



#### prevalence 2020 ####


kenya_2020_prevalence = raster("../data/Kenya_Prevalence/202206_Global_Pf_Parasite_Rate_KEN_2020.tiff")

kenya_2020_prevalence_matrix = rasterToPoints(kenya_2020_prevalence)

kenya_2020_prevalence = data.frame(kenya_2020_prevalence_matrix)
names(kenya_2020_prevalence) = c("x", "y", "prevalence")


#### prevalence 2019 ####


kenya_2019_prevalence = raster("../data/Kenya_Prevalence/202206_Global_Pf_Parasite_Rate_KEN_2019.tiff")

kenya_2019_prevalence_matrix = rasterToPoints(kenya_2019_prevalence)

kenya_2019_prevalence = data.frame(kenya_2019_prevalence_matrix)
names(kenya_2019_prevalence) = c("x", "y", "prevalence")



#### proportion of cases treated ####

kenya_2020_prop_cases = raster("../data/Kenya_PropCases/202106_Global_Antimalarial_Effective_Treatment_KEN_2020.tiff")

kenya_2020_prop_cases_matrix = rasterToPoints(kenya_2020_prop_cases)

kenya_2020_prop_cases = data.frame(kenya_2020_prop_cases_matrix)
names(kenya_2020_prop_cases) = c("x", "y", "prop_cases")


#### temperature: only use rasters where temperature data exists ####

load("../data/kenya_temperatures.RData")
temp_data = result_all


#### check that all have same length ####

nrow(kenya_2020_prop_cases)
nrow(kenya_2020_incidence)
nrow(kenya_2020_prevalence)
nrow(kenya_2019_prevalence)

# round all coordinates
kenya_2020_prop_cases$x = round(kenya_2020_prop_cases$x, 2)
kenya_2020_prop_cases$y = round(kenya_2020_prop_cases$y, 2)

kenya_2020_incidence$x = round(kenya_2020_incidence$x, 2)
kenya_2020_incidence$y = round(kenya_2020_incidence$y, 2)

kenya_2020_prevalence$x = round(kenya_2020_prevalence$x, 2)
kenya_2020_prevalence$y = round(kenya_2020_prevalence$y, 2)

kenya_2019_prevalence$x = round(kenya_2019_prevalence$x, 2)
kenya_2019_prevalence$y = round(kenya_2019_prevalence$y, 2)

# combine into a single data frame

kenya_data = kenya_2020_prop_cases %>%
  inner_join(kenya_2020_incidence, by=c("x", "y")) %>%
  rename(incidence_2020 = incidence) %>%
  inner_join(kenya_2020_prevalence, by=c("x", "y")) %>%
  rename(prevalence_2020 = prevalence) %>%
  inner_join(kenya_2019_prevalence, by=c("x", "y")) %>%
  rename(prevalence_2019 = prevalence)

save(kenya_data, file="../data/kenya_data.RData")
  

# save for temperature and rainfall use
save


#### add temperature and rainfall data ####


#### save all of these ####

save(kenya_2019_prevalence, file="../data/kenya_2019_prevalence.RData")
save(kenya_2020_prop_cases, file="../data/kenya_2020_prop_cases.RData")
save(kenya_2020_prevalence, file="../data/kenya_2020_prevalence.RData")
save(kenya_2020_incidence, file="../data/kenya_2020_incidence.RData")


