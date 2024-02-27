# Github test
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
populations = 2 #number of populations

#---------------------------------------------
#' TO AMEND TO 1 population                   |
#' populations = 1 #number of populations     |
#' BUT BREAKS MODEL SINGLE RUNCHANGED TO 1    |
#---------------------------------------------

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

#' --------------------------------------------------------------------------------------------------------------------
#' The travel_matrix defines potential contacts or travel between different populations within the model.             |
#' Originally designed to accommodate a scenario with 3  populations, it uses the 'populations' variable              |
#' to dynamically set its dimensions. Updating  model to single, non-traveling population, populations adjusted to=1  |
#' travel_matrix is now a 1x1  matrix becuase (diag,1, populations, popualtions) is the same as diag(1,1,1) yes?      | 
#' implicitly reflecting the absence of inter-population travel or contact. So no need to change travel_matrix, as it |
#' automatically conforms to the single population structure by the adjustment of 'populations' to 1.                 |
#' -------------------------------------------------------------------------------------------------------------------
  
#' determines migration between populations
migration_matrix = matrix(0, nrow = populations, ncol = populations)


#' --------------------------------------------------------------------------------------------------------------------
#' migration_matrix is for potential migration events between different populations within the model
#' Updating populations to=1, inter-population migration becomes irrelevant. no modifiation of migration_matrix is 
#' needed with populations=1, correct? And this keeps the model's integrity and allows to change later if wanted?
#' --------------------------------------------------------------------------------------------------------------------
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
#' model_populations = list(
#'  #one list for each parameter
#'  "unvaccinated" = list(
#'    "parameters" = list(
#'      #adjust acquisition of VTs and NVTs between the defined times (TODO refactor as optional)
#'      adjust_acq_vt = NULL_RR, adjust_acq_nvt = NULL_RR,
#'      adjust_acq_start = -1, adjust_acq_stop = 1e6,
#'      #the contact matrix, note that rows should be contactors and columns contactees
#'      betaVT = contact_matrix_model %>% t,
#'      betaNVT = contact_matrix_model %>% t,
#'      #the population size
#'      N = population_size_model[, value]),
#'    #specific parameters for each vaccine stratum in this population
#'    "arms" = list(
#'    "unvaccinated" = list(
#'        #no longer used. TODO refactor out
#'        "coverage" = no_coverage, "catchup_coverage" = no_coverage,
#'        #coverage for each arm
#'        #coverage for a catch-up campaign - these are applied to people who are in the defined age groups once
#'        "coverage_c" = list(
#'          list(value = no_coverage, time = 0)),
#'        #coverage for routine vaccination - these are applied to people ageing into the defined age group(s) continuously
#'        "coverage_r" = list(
#'        list(value = no_coverage, time = 0)),
#'        "waning" = vwaning,
#'        "efficacy" = ve))),
#'  "3p+0" = list(
#'    "parameters" = list(
#'      adjust_acq_vt = NULL_RR, adjust_acq_nvt = NULL_RR,
#'      adjust_acq_start = -1, adjust_acq_stop = 1e6,
#'      betaVT = contact_matrix_model %>% t,
#'      betaNVT = contact_matrix_model %>% t,
#'      N = population_size_model[, value]),
#'    "arms" = list(
#'      "unvaccinated" = list(
#'        #no longer used. TODO refactor out
#'        "coverage" = no_coverage, "catchup_coverage" = no_coverage,
#'        #coverage for each arm
#'        #coverage for a catch-up campaign - these are applied to people who are in the defined age groups once
#'        "coverage_c" = list(
#'          list(value = no_coverage, time = 0),
#'          list(value = getVaccineCoverage(age_groups_model, c(set_units(2, "months"), set_units(5, "years")), coverage = 0.80),
#'               time = 365*(1/12),
#'               coverage_to = rep("1p+0", age_groups_model[, .N]))),
#'        #coverage for routine vaccination - these are applied to people ageing into the defined age group(s) continuously
#'        "coverage_r" = list(
#'          list(value = no_coverage, time = 0),
#'          list(value = getVaccineCoverage(age_groups_model, c(set_units(2, "months")), coverage = 0.80), time = 365*(1/12),
#'               coverage_to = rep("1p+0", age_groups_model[, .N]))),
#'        "waning" = vwaning,
#'        "efficacy" = ve),
#'      "1p+0" = list(
#'        #no longer used. TODO refactor out
#'        "coverage" = no_coverage, "catchup_coverage" = no_coverage,
#'        "coverage_c" = list(
#'          list(value = no_coverage, time = 0)),
#'        "coverage_r" = list(
#'          list(value = no_coverage, time = 0),
#'          list(value = getVaccineCoverage(age_groups_model, c(set_units(3, "months")), coverage = 0.95), time = 365*(1/12),
#'               coverage_to = rep("2p+0", age_groups_model[, .N]))),
#'        "waning" = vwaning,
#'        "efficacy" = ve),
#'      "2p+0" = list(
#'        #no longer used. TODO refactor out
#'        "coverage" = no_coverage, "catchup_coverage" = no_coverage,
#'        #coverage for each arm
#'        "coverage_c" = list(
#'          list(value = no_coverage, time = 0)),
#'        "coverage_r" = list(
#'          list(value = no_coverage, time = 0),
#'          list(value = getVaccineCoverage(age_groups_model, c(set_units(4, "months")), coverage = 0.90), time = 365*(1/12),
#'               coverage_to = rep("3p+0", age_groups_model[, .N]))),
#'         "waning" = vwaning,
#'          "efficacy" = ve),
#'      "3p+0" = list(#no longer used. TODO refactor out
#'        "coverage" = no_coverage, "catchup_coverage" = no_coverage,
#'       "waning" = vwaning, "efficacy" = ve))))
#' ---------------------------------------------------------------------------------------------------------------------
#'
#' Specify populations to model for PCV schedules "3+0" and "1+1"
#'
#' This section outlines the configuration for different pneumococcal conjugate vaccine (PCV) schedules,
#' specifically focusing on the "3+0" (three doses in the primary series without a booster) and "1+1"
#' (one dose in the primary series plus a booster). It excludes the "2+0" schedule.
#'
#' - Make sure to transpose the contact matrices so that columns denote contactees and rows contactors.
#'   This is crucial for the matrix multiplication when calculating the force of infection (FOI) in the model.
#'   For example, row i in the transposed matrix represents the average number of contacts made by individuals
#'   of age i across all age groups. This row is then multiplied by a column vector representing the prevalence
#'   of infection in each age group, which is summed to calculate the average number of effective contacts
#'   with infectious individuals made by someone aged i. The transposition ensures the correct orientation
#'   for these calculations.
#'
#' - Note: The start and stop times for adjusting acquisition (acq) need refactoring. Currently, they are set
#'   to negative and very large time steps, respectively, to ensure they are always considered active in the
#'   model's calculations. This approach ensures that the adjustments for acquisition are applied universally
#'   across the modeled time frame, but may be refined for more specific temporal dynamics in future revisions.
#'
#' This configuration excludes the "2+0" schedule to concentrate on evaluating the impacts and outcomes
#' of the "3+0" and "1+1" schedules within the modeled populations.

#' Define model populations for evaluating PCV schedules
model_populations = list(
  # Unvaccinated population configuration
  "unvaccinated" = list(
    "parameters" = list(
      # Placeholder values for adjusting acquisition rates for vaccine types (VT) and non-vaccine types (NVT)
      adjust_acq_vt = NULL_RR, adjust_acq_nvt = NULL_RR,
      # Time range for acquisition adjustment: active from the start and essentially indefinitely
      adjust_acq_start = -1, adjust_acq_stop = 1e6,
      # Transpose contact matrices for VT and NVT to align rows with contactors and columns with contactees
      betaVT = contact_matrix_model %>% t(),
      betaNVT = contact_matrix_model %>% t(),
      # Population size derived from model-specific parameter
      N = population_size_model[, value]),
    # Vaccination arm configurations for the unvaccinated population
    "arms" = list(
      "unvaccinated" = list(
        # Coverage settings indicate no vaccination coverage or catch-up campaigns
        "coverage" = no_coverage, "catchup_coverage" = no_coverage,
        # Definitions for catch-up (c) and routine (r) vaccination coverage, both set to none
        "coverage_c" = list(list(value = no_coverage, time = 0)),
        "coverage_r" = list(list(value = no_coverage, time = 0)),
        # Parameters for vaccine waning and efficacy, specific to this arm
        "waning" = vwaning,
        "efficacy" = ve))),
  # Configuration for the "3+0" PCV schedule
  "3p+0" = list(
    "parameters" = list(
      # Similar configuration for adjusting acquisition rates and time range as the unvaccinated group
      adjust_acq_vt = NULL_RR, adjust_acq_nvt = NULL_RR,
      adjust_acq_start = -1, adjust_acq_stop = 1e6,
      # Transposed contact matrices for calculating FOIs with adjusted population and vaccination status
      betaVT = contact_matrix_model %>% t(),
      betaNVT = contact_matrix_model %>% t(),
      N = population_size_model[, value]),
    # Vaccination arm configurations under the "3+0" schedule
    "arms" = list(
      "unvaccinated" = list(
        "coverage" = no_coverage, "catchup_coverage" = no_coverage,
        # Catch-up and routine vaccination coverage, with specific values for initiating "1p+1" coverage
        "coverage_c" = list(
          list(value = no_coverage, time = 0),
          list(value = getVaccineCoverage(age_groups_model, c(set_units(2, "months"), set_units(5, "years")), coverage = 0.80),
               time = 365*(1/12),
               coverage_to = rep("1p+1", age_groups_model[, .N]))),
        "coverage_r" = list(
          list(value = no_coverage, time = 0),
          list(value = getVaccineCoverage(age_groups_model, c(set_units(2, "months")), coverage = 0.80), time = 365*(1/12),
               coverage_to = rep("1p+1", age_groups_model[, .N]))),
        "waning" = vwaning,
        "efficacy" = ve),
      # Configuration for the "1+1" PCV schedule within the "3+0" context
      "1p+1" = list(
        "coverage" = no_coverage, "catchup_coverage" = no_coverage,
        # Only routine vaccination coverage specified for transitioning to "3p+0" with high coverage
        "coverage_c" = list(list(value = no_coverage, time = 0)),
        "coverage_r" = list(
          list(value = no_coverage, time = 0),
          list(value = getVaccineCoverage(age_groups_model, c(set_units(3, "months")), coverage = 0.95), time = 365*(1/12),
               coverage_to = rep("3p+0", age_groups_model[, .N]))),
        "waning" = vwaning,
        "efficacy" = ve),
      # Default parameters for the "3+0" vaccination arm, reflecting no direct coverage but inheriting global parameters
      "3p+0" = list(
        "coverage" = no_coverage, "catchup_coverage" = no_coverage,
        "waning" = vwaning,
        "efficacy" = ve))))
#' ---------------------------------------------------------------------------------------------------------------------

#' update coverage_to for populations
model_populations = lapply(model_populations, function(population){
  arm_names = names(population$arms)
  population$arms = lapply(arm_names, function(name, arms){
    arm = arms[[name]]
    arm$coverage_r = arm$coverage_r %>% lapply(function(x, arms, name){
      if(!is.null(x$coverage_to)){
        x$coverage_to = x$coverage_to %>% sapply(function(z, arms){
          which(arms == z) - 1
        }, arms)  
      } else {
        x$coverage_to = rep(which(arms == name), age_groups_model[, .N])
      }
      
      return(x)
    }, names(arms), name)
    
    arm$coverage_c = arm$coverage_c %>% lapply(function(x, arms, name){
      if(!is.null(x$coverage_to)){
        x$coverage_to = x$coverage_to %>% sapply(function(z, arms){
          which(arms == z) - 1
        }, arms)  
      } else {
        x$coverage_to = rep(which(arms == name), age_groups_model[, .N])
      }
      
      return(x)
    }, names(arms), name)
    
    return(arm)
  }, population$arms)
  
  names(population$arms) = arm_names
  
  return(population)
})

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
params_unvac$trial_arms = model_populations %>% lapply(function(p){
  p$arms = list(
    "unvaccinated" = list("coverage" = no_coverage, "catchup_coverage" = no_coverage,
                          "coverage_c" = list(), "coverage_r" = list()))
  return(p)
})
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
