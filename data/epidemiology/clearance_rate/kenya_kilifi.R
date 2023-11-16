setAgeBreaks(as_units(c(0, 22, 41, 60), "month")) %>%
  .[, c("st", "value") := .("VT", 1/c(125, 67, 28, 22))] %>%
  rbind(setAgeBreaks(as_units(c(0, 22, 41, 60), "month")) %>%
          .[, c("st", "value") := .("NVT", 1/c(61, 50, 34, 24))])