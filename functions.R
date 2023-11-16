#' Function to run lapply and return a named list
lapplyNamed = function(x, fun, ...){
  if(is.null(names(x))) stop("Element has no names")
  #if(length(formals(fun)) != 1) stop("Function does not only have a single argument")
  if(is.list(x)){
    lapply(names(x), function(x, vals = x){ x = setNames(list(vals[[x]]), x); fun(x, ...) }, x)
  } else {
    lapply(names(x), function(x, vals = x){ x = vals[x]; fun(x, ...) }, x) 
  }
}

#' runs the model with the defined initial states, model parameters, and returns output on specified times
runModel = function(initial_states, model_params, times = c(0:365), steady_state = FALSE){
  if(steady_state){
    result = runsteady(
      y=initial_states, func = "derivs", initpar = model_params,
      dllname = ifelse(PCVM_VERSION == 2, "pcvm2", "pcvm"), nout = 1,
      initfunc = "rt_initmod", outnames = "output", jactype = "fullint")
    
    #' Get states at steady state
    result = data.table(time = attr(result, "time")) %>%
      cbind(data.table(i = 1:length(result$y), val = result$y)) %>%
      dcast(time~i, value.var="val") %>%
      reshapeModelOutput(model_params)
  } else {
    events = sort(unique(unlist(sapply(model_params$trial_arms, function(cluster)
      sapply(cluster$arms, function(arm) sapply(arm$coverage_c, "[[", "time"))))))
    
    if(!is.null(events)){
      result = lsoda(
        y=initial_states, times=times, func = "derivs", parms = model_params,
        dllname = ifelse(PCVM_VERSION == 2, "pcvm2", "pcvm"), nout = 1,
        initfunc = "initmod", outnames = "output",
        events=list(func="vaccineCampaignEvent",
                    time=sort(unique(unlist(sapply(model_params$trial_arms, function(cluster)
                      sapply(cluster$arms, function(arm) sapply(arm$coverage_c, "[[", "time"))))))))  
    } else {
      result = lsoda(
        y=initial_states, times=times, func = "derivs", parms = model_params,
        dllname = ifelse(PCVM_VERSION == 2, "pcvm2", "pcvm"), nout = 1,
        initfunc = "initmod", outnames = "output")
    }
    
    
    result = result %>% as.data.table %>% .[, -"output"] %>%
      reshapeModelOutput(model_params)  
  }
  
  result = result %>% merge(age_groups_model[, c("age", "name")], by="age")
  colnames(result) = c("age", "variable", "population", "dose", "compartment", "time", "value", "age_group")
  
  setorder(result, population, dose, compartment, age)
  
  return(result)  
}