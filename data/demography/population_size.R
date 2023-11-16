data_popsize = fread("./data/demography/gambia_population_size.csv")
data_popsize = cbind(setAgeBreaks(data_popsize[, age]), data_popsize[, -"age"])
data_popsize