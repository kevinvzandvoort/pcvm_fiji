#' load/update all libraries
pacman::p_load(data.table, Rcpp, RcppArmadillo, inline, deSolve, rootSolve, magrittr, units, ggplot2)

#' set arguments if R is running interactively, otherwise arguments should be passed to the Rscript command
#'  1: working directory of pcvmr (make sure to setwd manually if running interactively)
#'  2: working directory of pcvm
#'  3: name of output folder to store objects (e.g. date, or specific scenario name)
#'  4: integer identifying the chain that is ran (will also be used to set the seed)
.args = if(interactive()) c(getwd(), "model/pcvm", "20231102") else commandArgs(trailingOnly = TRUE)
.args = setNames(.args, c("pcvmr_dir", "pcvm_dir", "output_subdir"))
setwd(.args["pcvmr_dir"])

#' Create folder to store output
OUTPUT_FOLDER = sprintf("%s/output", .args["pcvmr_dir"])
if(!dir.exists(OUTPUT_FOLDER)) dir.create(OUTPUT_FOLDER)
OUTPUT_FOLDER = sprintf("%s/%s", OUTPUT_FOLDER, .args["output_subdir"])
if(!dir.exists(OUTPUT_FOLDER)) dir.create(OUTPUT_FOLDER)

#' Load PCVm
PCVM_FOLDER = .args["pcvm_dir"] #directory of pcvm
PCVM_COMPILE = FALSE #should pcvm be recompiled?
PCVM_VERSION = 1 #which pcvm version should be used? only use 1 for now
source(sprintf("%s/index.R", PCVM_FOLDER))

#' Set up project specific functions and parameters
source("./functions.R")
source("./data_load_all.R")

#' load the defined model_params and updateParameters function
source("./model_setup.R")

#' create a vector that will adjust some of the parameters
#' - we will fit these (or similar paramaters), but let's use these hard-coded values for now
#' - beta_1, beta_2, beta_3: these are applied to children aged <5y, 5-14y, and 15+y and represent
#'   the susceptibility of people, or the probability that contact would result in transmission (effective contact)
#' - beta_VT_NVT_u5 and beta_VT_NVT_o5: these are risk ratios for the relative transmissability of vaccine-types compared to non-vaccine types
adjust_params = c("beta_1" = 1e-2, "beta_2" = 1e-2, "beta_3" = 1e-3, "beta_VT_NVT_u5" = 1, "beta_VT_NVT_o5" = 0.8)

#' use these values to adjust the model parameters
current_model_params = updateParameters(adjust_params, model_params)

#' we can now run the model
#' let's first run the model for two years, and plot some results
results = runModel(initial_states = current_model_params$state_prop_cluster,
         model_params = current_model_params$params_unvac,
         times = c(0:365))

#' we can plot the prevalence over time, for the different compartment, for some age groups
results[age_group %in% c("[0m, 1m)", "[1y, 2y)", "[10y, 15y)")] %>%
  ggplot(aes(x=time, y=value, colour=compartment, linetype=age_group))+
  facet_grid(dose~population)+
  geom_line(linewidth=1)+
  theme_bw()

#' note how it takes a while before the compartments reach an equilibrium
#' we shouldn't start any vaccination before an equilibrium is reached
#' to do so, we can set the steady_state option to TRUE
results_steady_state = runModel(initial_states = current_model_params$state_prop_cluster,
                                model_params = current_model_params$params_unvac,
                                steady_state = TRUE)

#' we can plot the prevalence of different compartments at steady state
results_steady_state[age_group %in% c("[0m, 1m)", "[1y, 2y)", "[5y, 6y)", "[10y, 15y)", "[40y, 50y)") & population == "unvaccinated"] %>%
  ggplot(aes(x=age_group, y=value, colour=compartment, fill=compartment))+
  facet_grid(dose~population)+
  geom_col(position = "dodge")+
  theme_bw()

#' if we use the steady state values of this model as initial states, the model
#' will remain at steady state indefinitely
results = runModel(initial_states = results_steady_state[, value],
                   model_params = current_model_params$params_unvac,
                   times = c(0:365))

results[age_group %in% c("[0m, 1m)", "[1y, 2y)", "[10y, 15y)")] %>%
  ggplot(aes(x=time, y=value, colour=compartment, linetype=age_group))+
  facet_grid(dose~population)+
  geom_line(linewidth=1)+
  theme_bw()

#' to add vaccination, we need to add the states of the other vaccinated strata to our steady states
initial_states = results_steady_state[, -"age_group"] %>%
  setNames(c("age", "variable", "cluster", "dose", "compartment", "time", "value")) %>%
  eqStatesVaccinate(model_params$params_vac)

#' we can now run the model with vaccination
results = runModel(initial_states = initial_states[, state],
                   model_params = current_model_params$params_vac, #note we changed this to the parameters with vaccination
                   times = c(0:c(365*5)))

#' we can plot the prevalence again, but now vaccination occurs in the population
results[age_group %in% c("[0m, 1m)", "[1y, 2y)", "[10y, 15y)")] %>%
  ggplot(aes(x=time, y=value, colour=compartment, linetype=age_group))+
  facet_grid(dose~population)+
  geom_line(linewidth=1)+
  theme_bw()

#' we can see the results a bit better if we aggregate some groups
#' Nb - can be programmed more efficiently
N = model_params$params_unvac$trial_arms %>% seq_along %>%
  lapply(function(x, clusters) data.table(N = clusters[[x]]$parameters$N, cluster = names(clusters)[x]) %>% .[, age := 1:.N],
         model_params$params_unvac$trial_arms) %>% rbindlist
results = results %>% merge(N, by.x = c("population", "age"), by.y = c("cluster", "age")) %>% .[, value := value * N]

results_aggregated = combineAgeBreaks(setAgeBreaks(c(1,5,10,15)),
                 results[, .(value = sum(value), N=mean(N)), by=c("age", "population", "compartment", "time")] %>%
                 merge(age_groups_model, by="age"), value.var = c("value", "N"), method = "mean", by = c("population", "compartment", "time"))

results_aggregated %>% 
  .[name %in% c("[0y, 1y)", "[1y, 5y)", "[10y, 15y)")] %>%
  ggplot(aes(x=time, y=value/N, colour=compartment, linetype=name))+
  facet_grid(name~population)+
  geom_line(linewidth=1)+
  theme_bw()