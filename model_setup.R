#' A - Generic constant values
MODEL_TIMESTEP = 1 #days
MODEL_START_DATE = as.Date("2022-01-01") #when does the model start?

#' Age groups that will be used in the model
age_groups_model = c(set_units(c(0:12), "month"), set_units(c(1:5, 6, 10, 15, 20, 30, 40, 50), "years")) %>%
  setAgeBreaks() %>% .[, age := .I]

#' A - Generic constant values
#' - Risk-ratio that will do nothing (TODO: refactor)
NULL_RR = c(0) %>% setAgeBreaks() %>% .[, value := 1] %>%
  combineAgeBreaks(x = age_groups_model, y = .) %>% .[, value]
populations = 3 #number of populations

#' B - Infection specific parameters
#' clearance rates per day
clearance_rates = data$epidemiology$clearance_rate$data
  
#' The competition parameter reflects the extent to which carrying one or more pneumococcal serotypes of one group (carrying
#' VT or NVT) protects against acquisition of serotypes in the other group (NVT or VT). I.e. if competition is 20%, those
#' who carry VT serotypes experience 0.8 times the rate of infection with NVTs compared to someone who is susceptible
#' NB - is serotype specific in the real world
pneumo_competition = 0.7 #could fit this, if data is available
  
#' C - Vaccine specific parameters
vaccine_efficacy_transmission = data$vaccine$efficacy_transmission$data
vaccine_efficacy_disease = vaccine_efficacy_transmission %>% copy %>% .[, value := 1 - (1 - data$vaccine$efficacy_total$data[, value])/value] %>% .[]
vaccine_waning_rates = data$vaccine$efficacy_duration$data %>% copy %>% .[, value := 1/(value*365)]
  
#' D - Population specific parameters
population_size_model = age_groups_model %>% combineAgeBreaks(data$demography$population_data$data, method = "sum", value.var = "total")
population_size_model = population_size_model[, value := total] %>% .[, -"total"]
contact_matrix_model = adjustContactMatrixAgeGroups(age_groups_model[, -"age"], data$demography$contact_data$data$contact_matrix,
                                                    data$demography$contact_data$data$population, population_size_model)
  
#' determines contact between populations
travel_matrix = diag(1, populations, populations)
  
#' determines migration between populations
migration_matrix = matrix(0, nrow = populations, ncol = populations)
  
#' Values repeated in vaccination_strategies
no_coverage = age_groups_model %>% getVaccineCoverage(set_units(0, "year"), 0)
ve = age_groups_model %>% combineAgeBreaks(vaccine_efficacy_transmission) %>% .[, value]
vwaning = age_groups_model %>% combineAgeBreaks(vaccine_waning_rates) %>% .[, value]

#' Specify populations to model
#' - Make sure to transpose the contact matrices so that columns denote contactees and rows contactors
#'   as this is required for the matrix multiplication when calculating the FOIs in the model
#'   i.e row i will show average number of contacts made by those of age i in all age groups, which will
#'   be multiplied with a column vector with prevalence in each contacted age groups and summed to give
#'   the average number of effective contacts with infectious individuals made by someone aged i (or j in
#'   contact matrix before transposing)
#' - Nb need to refactor start and stop of acq adjustments. Set for negative and very large timesteps now
#'   to make sure they are always active
model_populations = list(
  #one list for each parameter
  "unvaccinated" = list(
    "parameters" = list(
      #adjust acquisition of VTs and NVTs between the defined times (TODO refactor as optional)
      adjust_acq_vt = NULL_RR, adjust_acq_nvt = NULL_RR,
      adjust_acq_start = -1, adjust_acq_stop = 1e6,
      #the contact matrix, note that rows should be contactors and columns contactees
      betaVT = contact_matrix_model %>% t,
      betaNVT = contact_matrix_model %>% t,
      #the population size
      N = population_size_model[, value]),
    #specific parameters for each vaccine stratum in this population
    "arms" = list()),
  "3p+0" = list(
    "parameters" = list(
      adjust_acq_vt = NULL_RR, adjust_acq_nvt = NULL_RR,
      adjust_acq_start = -1, adjust_acq_stop = 1e6,
      betaVT = contact_matrix_model %>% t,
      betaNVT = contact_matrix_model %>% t,
      N = population_size_model[, value]),
    "arms" = list(
      "1p+0" = list(
        #no longer used. TODO refactor out
        "coverage" = no_coverage, "catchup_coverage" = no_coverage,
        #coverage for each arm
        #coverage for a catch-up campaign - these are applied to people who are in the defined age groups once
        "coverage_c" = list(
          list(value = getVaccineCoverage(age_groups_model, c(set_units(2, "months"), set_units(2, "years")), coverage = 0.80), time = 365*(1/12))),
        #coverage for routine vaccination - these are applied to people ageing into the defined age group(s) continuously
        "coverage_r" = list(
          list(value = getVaccineCoverage(age_groups_model, set_units(2, "months"), coverage = 0.80), time = 365*(1/12)),
          list(value = getVaccineCoverage(age_groups_model, set_units(2, "months"), coverage = 0.75), time = 365*1),
          list(value = getVaccineCoverage(age_groups_model, set_units(2, "months"), coverage = 0.85), time = 365*2)),
        "waning" = vwaning,
        "efficacy" = ve),
      "2p+0" = list(
        #no longer used. TODO refactor out
        "coverage" = no_coverage, "catchup_coverage" = no_coverage,
        #coverage for each arm
        "coverage_c" = list(list(value = no_coverage, time = 0)),
        "coverage_r" = list(
          list(value = getVaccineCoverage(age_groups_model, set_units(3, "months"), coverage = 0.95), time = 365*(1/12))),
        "waning" = vwaning,
        "efficacy" = ve),
      "3p+0" = list(#no longer used. TODO refactor out
        "coverage" = no_coverage, "catchup_coverage" = no_coverage,
        #coverage for each arm
        "coverage_c" = list(list(value = no_coverage, time = 0)),
        "coverage_r" = list(
          list(value = getVaccineCoverage(age_groups_model, set_units(4, "months"), coverage = 0.95), time = 365*(1/12))),
        "waning" = vwaning,
        "efficacy" = ve))),
  "1p+1" = list(
    "parameters" = list(
      adjust_acq_vt = NULL_RR, adjust_acq_nvt = NULL_RR,
      adjust_acq_start = -1, adjust_acq_stop = 1e6,
      betaVT = contact_matrix_model %>% t,
      betaNVT = contact_matrix_model %>% t,
      N = population_size_model[, value]),
    "arms" = list(
      "1p+0" = list(
        #no longer used. TODO refactor out
        "coverage" = no_coverage, "catchup_coverage" = no_coverage,
        #coverage for each arm
        "coverage_c" = list(list(value = no_coverage, time = 0)),
        "coverage_r" = list(
          list(value = getVaccineCoverage(age_groups_model, set_units(2, "months"), coverage = 0.80), time = 365*(1/12)),
          list(value = getVaccineCoverage(age_groups_model, set_units(2, "months"), coverage = 0.75), time = 365*1),
          list(value = getVaccineCoverage(age_groups_model, set_units(2, "months"), coverage = 0.85), time = 365*2)),
        "waning" = vwaning,
        "efficacy" = ve),
      "1p+1" = list(
        #no longer used. TODO refactor out
        "coverage" = no_coverage, "catchup_coverage" = no_coverage,
        #coverage for each arm
        "coverage_c" = list(list(value = no_coverage, time = 0)),
        "coverage_r" = list(
          list(value = getVaccineCoverage(age_groups_model, set_units(12, "months"), coverage = 0.95), time = 365*(1/12))),
        "waning" = vwaning,
        "efficacy" = ve))))

#' E - Generic parameters for all populations used in the model
params_vac = list(
  n_agrp = age_groups_model[, .N], #total number of age groups
  comp = 1 - pneumo_competition, #competition parameter
  clearVT = age_groups_model %>% combineAgeBreaks(clearance_rates[st == "VT"]) %>% .[, value], #clearance rates VTs
  clearNVT = age_groups_model %>% combineAgeBreaks(clearance_rates[st == "NVT"]) %>% .[, value], #clearance rates NVTs
  ageout = age_groups_model %>% .[, .(duration = (to - from) %>% set_units("days"))] %>% .[, 1/as.numeric(duration)], #rate at which people age (depends on defined age groups)
  trial_arms = model_populations, #populations that will be modelled (defined above)
  travel = travel_matrix, #travel matrix (contact between populations)
  migration = migration_matrix) #migration matrix (movement between populations)

#' For the unvaccinated scenario, all model population arms are emptied
params_unvac = params_vac
params_unvac$trial_arms = model_populations %>% lapply(function(p){p$arms = list(); return(p)})

#' F - Setup model parameters passed to the model
model_params = list(
  params_unvac = params_unvac,
  params_vac = params_vac,
  times_postvacc_eval = c(0:(365*5)), #for how long should the model be ran?
  state_prop_cluster = c( #initial states (note that we model proportions)
    rep(0.8, age_groups_model[, .N]), #S
    rep(0.1, age_groups_model[, .N]), #VT
    rep(0.1, age_groups_model[, .N]), #NVT
    rep(0, age_groups_model[, .N]) #B
  ) %>% rep(length(model_populations)))

#' Adjust for timestep
model_params = model_params %>% adjustForTimeStep()
contact_matrix_model = contact_matrix_model %>% adjustForTimeStep()

#' Add the (unadjusted) contact matrix to the model parameters
#' this will give one list (model_params) with all parameters and specifications used by the model
model_params$cm_unadjusted = contact_matrix_model

#' Priors for parameters to fit
#' - we will look at this later
#priors = rbindlist(list(
#  data.table(variable = "beta_1", min = 0, max = 1, plotmin = 1e-4, plotmax = 2e-2,
#             density = function(x) dbeta(x, shape1 = 0.1, shape2 = 10, log=TRUE),
#             sampler = function(n) rbeta(n, shape1 = 0.1, shape2 = 10)),
#  data.table(variable = "beta_2", min = 0, max = 1, plotmin = 1e-4, plotmax = 2e-2,
#             density = function(x) dbeta(x, shape1 = 0.1, shape2 = 10, log=TRUE),
#             sampler = function(n) rbeta(n, shape1 = 0.1, shape2 = 10)),
#  data.table(variable = "beta_3", min = 0, max = 1, plotmin = 1e-4, plotmax = 2e-2,
#             density = function(x) dbeta(x, shape1 = 0.1, shape2 = 10, log=TRUE),
#             sampler = function(n) rbeta(n, shape1 = 0.1, shape2 = 10)),
#  data.table(variable = "beta_VT_NVT_u5", min = 0.2, max = 5, plotmin = 0.2, plotmax = 2.5,
#             density = function(x) dlnorm(x, meanlog = 1, sdlog = 1, log=TRUE),
#             sampler = function(n) rlnorm(n, meanlog = 1, sdlog = 1)),
#  data.table(variable = "beta_VT_NVT_o5", min = 0.2, max = 5, plotmin = 0.2, plotmax = 2.5,
#             density = function(x) dlnorm(x, meanlog = 1, sdlog = 1, log=TRUE),
#             sampler = function(n) rlnorm(n, meanlog = 1, sdlog = 1))
#))

#' Function that will be used to replace the model parameters with new values
#' this function takes two parameters:
#' - a named vector with the new parameters that you are using
#' - the model_params object that you defined earlier
updateParameters = function(params, model_params){
  
  #' Assign betaVT to the correct age groups
  betaVT = set_units(c(0, 5, 15), "year") %>% setAgeBreaks() %>%
    .[, value := c(params["beta_1"]*(params["beta_VT_NVT_u5"]),
                   params["beta_2"]*(params["beta_VT_NVT_o5"]),
                   params["beta_3"]*(params["beta_VT_NVT_o5"]))] %>%
    combineAgeBreaks(x = age_groups_model, y = .) %>% .[, value]
  
  betaNVT = set_units(c(0, 5, 15), "year") %>% setAgeBreaks() %>%
    .[, value := c(params["beta_1"],
                   params["beta_2"],
                   params["beta_3"])] %>%
    combineAgeBreaks(x = age_groups_model, y = .) %>% .[, value]
  
  new_matrix_VT = sweep(model_params$cm_unadjusted, 2, betaVT, "*") %>% t
  new_matrix_NVT = sweep(model_params$cm_unadjusted, 2, betaNVT, "*") %>% t
  
  #' do this for both unvac and vac (if applicable)
  if(length(model_params[["params_vac"]]) > 0){
    scenarios = c("params_unvac", "params_vac")
  } else {
    scenarios = "params_unvac"
  }
  
  for(s in scenarios){
    for(p in names(model_params[[s]][["trial_arms"]])){
      model_params[[s]][["trial_arms"]][[p]][["parameters"]][["betaVT"]] = new_matrix_VT
      model_params[[s]][["trial_arms"]][[p]][["parameters"]][["betaNVT"]] = new_matrix_NVT
    }
  }
    
  return(model_params)
}
