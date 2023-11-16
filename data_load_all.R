#' All data is stored in the ./data subfolder

#' The first time getData is ran, it sources an R-script that returns a single object.
#' All subsequent times getData is called, it returns the object initially generated.
#' @param name string: path to the R-script, without the extension. Object will be stored as an RDS file with the same name, appended with "_processed".
#' @return returns the object created by the R-script
getData = function(name, force=FALSE){
  if(!force && file.exists(sprintf("%s_processed.RDS", name))){
    return(readRDS(sprintf("%s_processed.RDS", name)))
  } else {
    value = source(sprintf("%s.R", name))$value
    saveRDS(value, sprintf("%s_processed.RDS", name))
    return(value)
  }
}

#' Create one long list with all the data needed for the model
#' Nb, this is not required for the model to work, but I prefer to organize my data this way
data = list(
  epidemiology = list(
    prevalence = list(
      gambia_2006_hill = list(
        source = "Hill et al, 2006; https://doi.org/10.1086/506941", #I like to add the source where the data came from
        data = getData("./data/epidemiology/prevalence/gambia_2006_hill"))), #the file that pre-processes and reads the data (without the R extension)
    clearance_rate = list(
      source = "Lipsitch et al 2012; https://doi.org/10.1097/ede.0b013e31824f2f32",
      data = getData("./data/epidemiology/clearance_rate/kenya_kilifi"))),
  vaccine = list(
    efficacy_transmission = list(
      source = "https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(23)00061-0/fulltext#seccestitle170",
      data = setAgeBreaks(0) %>% .[, value := 0.5] %>% .[]),
    efficacy_total = list(
      source = "https://publications.aap.org/pediatrics/article/145/4/e20190377/36954/Efficacy-and-Effectiveness-of-the-PCV-10-and-PCV",
      data = setAgeBreaks(0) %>% .[, value := 0.80] %>% .[]),
    efficacy_duration = list(
      source = "https://pubmed.ncbi.nlm.nih.gov/26075814/",
      data = setAgeBreaks(0) %>% .[, value := 6] %>% .[])),
  demography = list(
    population_data = list(
        source = "UN World Population Prospects; https://population.un.org/wpp/Download/Standard/CSV",
        data = getData("./data/demography/population_size")),
      contact_data = list(
        source = "van Zandvoort, 2022; https://doi.org/10.1016/j.epidem.2022.100625",
        data = getData("./data/demography/contact_matrix"))))
