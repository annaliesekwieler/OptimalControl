# This code does exploratory data analysis from VectorBase

#### packages ####

library(tidyverse)

#### load data ####

abundance_data = read.csv("../data/abundance_data.csv")

head(abundance_data)

#### analyze data ####

hist(abundance_data$Specimens.collected)

table(abundance_data$Collection.protocols)

abundance_data %>%
  filter(Collection.protocols == "man biting catch")

table(abundance_data$Collection.date.range)

table(abundance_data$)Collection.duration.days
