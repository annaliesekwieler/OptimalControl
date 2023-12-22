load("kenya_2020.RData")

library(chirps)

kenya_2020_rainfall = list()

for(i in 1:nrow(kenya_2020)){
  
  lonlat_kenya = kenya_2020[i, c("x", "y")]
  
   kenya_2020_rainfall[[i]] = get_chirps(lonlat_kenya, dates = c("2020-01-01",
                                                                "2020-12-31"),
                                        server="CHC")
  
}