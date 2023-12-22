# load the parallelized RData files and consolidate into a single list of temperatures

load("../../data/kenya_climate_data/kenya_2016_temperatures.RData")

length(kenya_2016_temperatures)

load("../../data/kenya_climate_data/temp_data_kenya_2.RData")
result_all = c(kenya_2016_temperatures, result)

load("../../data/kenya_climate_data/temp_data_kenya_3.RData")
result_all = c(result_all, result)

load("../../data/kenya_climate_data/temp_data_kenya_4.RData")
result_all = c(result_all, result)

load("../../data/kenya_climate_data/temp_data_kenya_5.RData")
result_all = c(result_all, result)

load("../../data/kenya_climate_data/temp_data_kenya_6.RData")
result_all = c(result_all, result)

load("../../data/kenya_climate_data/temp_data_kenya_7.RData")
result_all = c(result_all, result[-1])


print(length(result_all))

save(result_all, file="../../data/kenya_temperatures.RData")
