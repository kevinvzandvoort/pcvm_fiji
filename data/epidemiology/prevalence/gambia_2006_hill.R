data = fread("./data/epidemiology/prevalence/gambia_2006_hill.csv")
data = setAgeBreaks(c(1, 5, 15, 40)) %>%
  cbind(data[, -"age"] %>%
          .[, c("VT", "NVT", "S") := .(N*VT, N*NVT, N*S)])
data
