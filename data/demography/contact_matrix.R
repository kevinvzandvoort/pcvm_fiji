pacman::p_load(magrittr, socialmixr, data.table)

if(!file.exists("./data/demography/digaale_contact_data.RDS")){
  #' Get data from Zenodo
  socialmixr::get_survey("https://zenodo.org/record/5226281#.YR-TzlvTVH6") %>%
    saveRDS("./data/demography/digaale_contact_data.RDS")
}

#' The estimated population size in Digaale (for provided age groups)
#'  can manually be downloaded
if(!file.exists("./data/demography/digaale_survey_population.RDS")){
  digaale_survey_population =
    data.table::fread("https://zenodo.org/record/5226281/files/espicc_somaliland_digaale_survey_population.csv")
  saveRDS(digaale_survey_population, "./data/demography/digaale_survey_population.RDS")
}

#' Note that weekends fall on Fridays and Saturdays in Somaliland.
#' - The dayofweek variable provided in the dataset has been kept
#'   consistent with R defaults (0: Sunday to 6: Saturday)
digaale_contact_data = readRDS("./data/demography/digaale_contact_data.RDS")
digaale_contact_data$participants[, c("dayofweek", "dayofweek_name", "weekend")] %>%
      unique %>%
      setorder(dayofweek) %>%
      #' socialmixr currently assumes the weekend to fall on dayofweek
      #'  6 (Saturday) and 0 (Sunday)
      #' - dayofweek can be manually edited so that Fridays and Saturdays
      #'   are taken as the weekend, if you wish to weight contacts by
      #'   weekday
      .[, dayofweek := ifelse(dayofweek == 6, 0, dayofweek + 1)]
    
digaale_survey_population = readRDS("./data/demography/digaale_survey_population.RDS") %>%
  .[lower.age.limit >= 50, lower.age.limit := 50] %>%
  .[, .(population=sum(population)), by="lower.age.limit"]

#' The contact matrix can then be constructed as follows
#' - The provided survey_population can be used to construct a
#'   population representative matrix for Digaale IDP camp
#' - As the sample is not self-weighing (oversampling of young
#'   age groups), it is recommended to apply the survey_weight
#'   as weights
#' Note socialmixr's contact matrices show contactors in rows
#'  and contactees in columns
digaale_contact_matrix = digaale_contact_data %>%
  socialmixr::contact_matrix(survey.pop = digaale_survey_population,
                             age.limits = digaale_survey_population$lower.age.limit,
                             symmetric = TRUE, weights = "survey_weight", weigh.dayofweek = TRUE)

digaale_survey_population = cbind(setAgeBreaks(digaale_survey_population$lower.age.limit), digaale_survey_population[, -"lower.age.limit"])

#' put in right format for PCVm
digaale_contact_matrix = as.data.table(digaale_contact_matrix$matrix)
colnames(digaale_contact_matrix) = as.character(digaale_survey_population$name)
digaale_contact_matrix[, contactor_age_group := digaale_survey_population$name]
digaale_contact_matrix = melt(digaale_contact_matrix, measure.vars = digaale_survey_population$name, variable.name = "contactee_age_group")
digaale_contact_matrix[, contactee_age_group := factor(contactee_age_group, digaale_survey_population$name)] %>% .[]
setorder(digaale_contact_matrix, contactor_age_group, contactee_age_group)

list(contact_matrix = digaale_contact_matrix,
     population = digaale_survey_population)
