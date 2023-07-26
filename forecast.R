forecast_simulation <- function(events, T_max = Inf, N_max = Inf,
                                  kernel, parameters,
                                  mu_fn = mu_none, mu_fn_diff = mu_diff_none,
                                  mu_t_max = NULL, imported = F,
                                  print_level = 1) {
  
  if (!exists("delay", parameters)) {
    parameters$delay <- 0
  }
  
  # Initialize current_time to the maximum time in the historical data
  current_time <- max(events)
  # Store the maximum time in the original historical data
  max_original_event <- current_time
  
  # Initialize an empty vector to store the new events
  new_events <- numeric()
  
  while (length(new_events) < N_max && current_time < T_max) {
    # Generate a new event
    current_time <- simulate_next_event(time = current_time,
                                        events = events, T_max = T_max,
                                        num_children = 1, kernel = kernel,
                                        parameters = parameters, mu_fn = mu_fn,
                                        mu_fn_diff = mu_fn_diff, mu_t_max = mu_t_max,
                                        imported = imported, print_level = print_level)
    
    if (length(current_time) == 0 || current_time > T_max ||
        is.infinite(current_time)) {
      break
    }
    
    # Add the new event to the new_events vector
    new_events <- c(new_events, current_time)
    
    # Add the new event to the history
    events <- sort(c(events, current_time))
  }
  
  return(c(max_original_event,new_events)) #return the first time point from which we forecast
  # along with the forecasted points
}

