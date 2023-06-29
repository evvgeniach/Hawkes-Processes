
# Creating a data frame for new_times
df_new_times <- data.frame(t = new_times, N = seq(1, length(new_times)), type = "Real Data")

# Initialize a list to hold all simulated event data frames
list_df_events <- list()

for(i in 1:N_runs){
  
  events <- list_events[[i]]
  
  # Creating a data frame for events
  df_events <- data.frame(t = events, N = seq(1, length(events)), type = paste("Simulated Data ", i))
  
  # Store each data frame in the list
  list_df_events[[i]] <- df_events
}

## CUMULATIVE PLOT

# Combine all the data into one dataframe
all_simulations_zika <- do.call(rbind, list_df_events)

# Convert Simulation column into factor to help with plotting
all_simulations_zika$Simulation <- as.factor(all_simulations_zika$type)

df_new_times$Simulation <- as.factor(df_new_times$type)
# Add real data to the same dataframe
all_data_zika <- rbind(all_simulations_zika, df_new_times)

# Create the plot
plot_zika<- ggplot() +
  theme_bw() +
  labs(title = "Comparison of Simulated and Real Data for Zika",
       x = "Time (days)",
       y = expression(N(t))) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# Add the simulated data lines
plot_zika <- plot_zika +
  geom_path(data = subset(all_data_zika, type != "Real Data"), 
            aes(x = t, y = N, group = Simulation), 
            color = "black", alpha = 0.3, size = 0.3)

# Add the real data line
plot_zika <- plot_zika +
  geom_path(data = subset(all_data_zika, type == "Real Data"), 
            aes(x = t, y = N), color = "red", size = 1)

print(plot_zika)

## NON cumulative plot

# Initialize an empty list to hold the data frames for each simulation
list_df_non_cum_zika <- list()
list_df_non_cum_zika_days <- list()
# Loop over the list of simulations
for(i in 1:N_runs){
  # Convert simulation events from seconds to dates
  non_cum <- as.POSIXct(list_events[[i]] * 86400, origin = "2015-10-19")
  
  # Convert the filtered events to weeks from the origin
  non_cum1 <- floor(as.numeric(difftime(non_cum, "2015-10-19", units = "weeks")))
  tbl_zika <- factor(non_cum1, levels = 0:13) ## there are 14 weeks
  all_dates_zika <- as.Date(seq(as.Date("2015-10-19"), max(as.Date(zika_girardot_2015$date)), by = "day"))
  non_cum_days_zika <- as.data.frame(table(as.Date(non_cum)))
  names(non_cum_days_zika) <- c("date", "cases")
  non_cum_days_zika$date <- as.Date(non_cum_days_zika$date)
  ## Fill the missing dates with a value of 0
  data_all_zika <- merge(data.frame(date = all_dates_zika), non_cum_days_zika, by = "date", all.x = TRUE)
  # Replace NA values with 0
  data_all_zika$cases[is.na(data_all_zika$cases)] <- 0
  # Create a data frame for the non-cumulative counts for this simulation
  non_cum_df_zika <- data.frame(
    Week = as.numeric(names(table(tbl_zika))), 
    Simulated_Count = as.vector(table(tbl_zika)),
    True_count = week_dat_3$cases
  )
  non_cum_df_zika_days <- data.frame(
    Days = all_dates_zika,
    True_days = zika_girardot_filled$cases,
    Simulated_days = data_all_zika$cases
  )
  # Add the data frame to the list
  list_df_non_cum_zika[[i]] <- non_cum_df_zika
  list_df_non_cum_zika_days[[i]] <- non_cum_df_zika_days
}

# Prepare the simulated data
for(i in 1:N_runs){
  list_df_non_cum_zika[[i]]$type <- paste("Simulation ", i)
  list_df_non_cum_zika_days[[i]]$type <- paste("Simulation ", i)
}

# Combine all the data into one dataframe
all_simulations_non_cum_zika <- do.call(rbind, list_df_non_cum_zika)

# Convert Simulation column into factor to help with plotting
all_simulations_non_cum_zika$Simulation <- as.factor(all_simulations_non_cum_zika$type)

# Create the plot for WEEKS
plot_non_cum_zika <- ggplot() +
  theme_bw() +
  labs(title = "Counts of Simulated Events and True Observations of Zika",
       x = "Time (Weeks)",
       y = "Count") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# Add the simulated data lines
plot_non_cum_zika <- plot_non_cum_zika +
  geom_path(data = all_simulations_non_cum_zika, 
            aes(x = Week, y = Simulated_Count, group = Simulation), 
            color = "black", alpha = 0.1, size = 0.3)

# Add the real data line
plot_non_cum_zika <- plot_non_cum_zika +
  geom_path(data = list_df_non_cum_zika[[1]], 
            aes(x = Week, y = True_count), color = "red", size = 1)

print(plot_non_cum_zika)
#ggsave("Zika_non_cum_plot_weeks.pdf", plot = plot_non_cum_zika, dpi = 600)

### DAYS

# Create a combined data frame with all the simulation results
all_simulations_df_zika <- do.call(rbind, list_df_non_cum_zika_days)
all_simulations_df_zika$Simulation <- as.factor(all_simulations_df_zika$type)
all_simulations_df_zika$Num_Days <- as.numeric(difftime(all_simulations_df_zika$Days, min(all_simulations_df_zika$Days), units = "days"))

# Create the plot for DAYS
plot_days_zika <- ggplot() +
  theme_bw() +
  labs(title = "Counts of Simulated Events and True Observations of Zika",
       x = "Time (Days)",
       y = "Count") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# Add the simulated data lines
plot_days_zika <- plot_days_zika +
  geom_line(data = all_simulations_df_zika, aes(x = Num_Days, y = Simulated_days, group = Simulation), color = "black", alpha = 0.1)

# Add the real data line
plot_days_zika <- plot_days_zika +
  geom_line(data = all_simulations_df_zika, aes(x = Num_Days, y = True_days), color = "red", size = 0.3)

print(plot_days_zika)


# Reshape the data
library(tidyverse)
all_simulations_df_zika_long <- all_simulations_df_zika %>%
  pivot_longer(c(Simulated_days, True_days), names_to = "Type", values_to = "Count")
# Create the bar plot for DAYS
plot_days_zika <- ggplot(all_simulations_df_zika_long, aes(x = as.factor(Num_Days), y = Count, fill = Type)) +
  geom_col(position = position_dodge()) +
  labs(title = "Counts of Simulated Events and True Observations of Zika",
       x = "Time (Days)",
       y = "Count") +
  theme_bw() +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 12, angle = 90))

print(plot_days_zika)

## bar plot now

plot_days_zika <- ggplot() +
  theme_bw() +
  labs(title = "Counts of Simulated Events and True Observations of Zika",
       x = "Time (Days)",
       y = "Count") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# Add the simulated data bars
plot_days_zika <- plot_days_zika +
  geom_bar(data = all_simulations_df_zika, aes(x = Num_Days, y = Simulated_days, fill = Simulation), stat = "identity", alpha = 0.1, color = "black")

# Add the real data bar
plot_days_zika <- plot_days_zika +
  geom_bar(data = all_simulations_df_zika, aes(x = Num_Days, y = True_days), stat = "identity", fill = "red", color = "red")

print(plot_days_zika)

################## Intensities ################

# Calculate intensities for all simulations and real data
list_intensities <- list()
parameters = list(alpha = alpha_zika,
                  delta = delta_zika,
                  #A = A_zika,
                  #B = B_zika,
                  #C = C_zika,
                  M = M_zika,
                  N= N_zika,
                  #P = P_zika,
                  delay = 0)
for(i in 1:N_runs){
  events <- list_events[[i]]
  data <- compute_intensity_function(events = events, kernel = exp_kernel, 
                                     T_max = max(new_times), parameters = parameters, mu_fn = mu_fn, 
                                     N = 5000)
  mu_ts <- mu_fn(events, parameters =  parameters)
  event_intensities <- mu_ts + conditional_intensity_list(times = events +1e-10, 
                                                          events = events, 
                                                          kernel = exp_kernel, 
                                                          parameters = parameters)
  data_events <- data.frame(t = events, intensity = event_intensities, type = paste("Simulated Data ", i))
  list_intensities[[i]] <- data_events
}

# Calculate intensities for the real data
data <- compute_intensity_function(events = new_times, kernel =exp_kernel, 
                                   T_max = max(new_times), parameters = parameters, mu_fn = mu_fn, 
                                   N = 5000)
mu_ts <- mu_fn(new_times, parameters =  parameters)
event_intensities <- mu_ts + conditional_intensity_list(times = new_times +1e-10, 
                                                        events = new_times, 
                                                        kernel = exp_kernel, 
                                                        parameters = parameters)
df_new_times <- data.frame(t = new_times, intensity = event_intensities, type = "Real Data")


# plot the intensities
plot <- ggplot() +
  theme_bw() +
  labs(title = "Comparison of Simulated and Real Intensities",
       x = "Time (weeks)",
       y = expression(lambda(t))) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# Add the simulated data lines with alpha for transparency
for(i in 1:N_runs){
  plot <- plot + geom_line(data = list_intensities[[i]], aes(x = t, y = intensity), color = "black", alpha = 0.2)
}

# Add the real data line
plot <- plot + geom_line(data = df_new_times, aes(x = t, y = intensity), color = "red", alpha = 0.9)

# Display the plot
print(plot)




