## Ebola with first count of infection removed as outlier
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
library(outbreaks)
library(epihawkes)
library(ggplot2)
library(DEoptim)
ebola_kikwit_1995_mod <- ebola_kikwit_1995[-(1:59),]  # Remove the first row

total_days_mod <- as.integer(max(ebola_kikwit_1995_mod$date) - min(ebola_kikwit_1995_mod$date)) + 1
cumulative_cases_mod <- cumsum((ebola_kikwit_1995_mod$onset))

ggplot(data.frame(total_days = 1:total_days_mod, cumulative_cases = cumulative_cases_mod), aes(x = total_days, y = cumulative_cases)) +
  geom_step(color = "steelblue", linewidth = 1) +
  theme_bw()+
  labs(title = "Cumulative Plot of Cases for Ebola in Kikwit",
       x = "Total Days",
       y = "Cumulative Cases")

ggplot(ebola_kikwit_1995_mod, aes(x = 1:total_days_mod, y = onset)) +
  theme_bw()+
  geom_line(color = "red") +
  labs(title = "Line Plot of Onset",
       x = "Number of Days",
       y = "Count")

# Add a new column that is the number of weeks since the start date
ebola_kikwit_1995_mod$Week <- as.numeric(difftime(ebola_kikwit_1995_mod$date, min(ebola_kikwit_1995_mod$date), units = "weeks"))
ebola_kikwit_1995_mod$Week <- floor(ebola_kikwit_1995_mod$Week)

week_dat_mod <- ebola_kikwit_1995_mod %>%
  group_by(Week) %>%
  summarise(onset = sum(onset))

week_dat_mod$cumulative_cases <- cumsum(week_dat_mod$onset)

ggplot(week_dat_mod, aes(x = Week, y = onset)) +
  geom_line(color = "red") +
  labs(title = "Line Plot of Cases for Ebola in Kikwit",
       x = "Number of Weeks",
       y = "Count")

ggplot(week_dat_mod, aes(x = Week, y = cumulative_cases)) +
  geom_step(color = "steelblue", linewidth = 1) +
  labs(title = "Cumulative Plot of Cases for Ebola in Kikwit",
       x = "Total Weeks",
       y = "Cumulative Cases")

bar_plot_ebola_mod <- ggplot(week_dat_mod, aes(x = Week, y = onset)) +
  geom_col(fill = "steelblue", width = 0.5) +
  labs(title = "Bar Plot of Cases for Ebola in Kikwit",
       x = "Number of Weeks",
       y = "Count") +
  theme(
    plot.title = element_text(size = 20),  
    axis.title.x = element_text(size = 16),  
    axis.title.y = element_text(size = 16),  
    plot.subtitle = element_text(size = 14),  
    axis.text = element_text(size = 18)
  )

set.seed(52143)
generate_timestamps <- function(n) {
  sample(seq(from = 0, to = 24*60*60 - 1), size = n, replace = TRUE)
}

df3_mod <- ebola_kikwit_1995_mod %>%
  filter(onset > 0) %>%
  mutate(date2 = as.POSIXct(date))

df3_mod$case_hour <- purrr::map(df3_mod$onset, generate_timestamps)
df3_mod = tidyr::unnest(df3_mod, case_hour)
df3_mod$case_epoch = df3_mod$date2 + df3_mod$case_hour

df3_mod <- unnest(df3_mod, case_epoch)
df3_mod <- df3_mod[order(df3_mod$case_epoch),]

min_time_mod <- min(df3_mod$case_epoch)
df3_mod$case_epoch <- df3_mod$case_epoch - min_time_mod
new_times_ebola_mod <- c(as.integer(df3_mod$case_epoch))/86400 # turn seconds to days
mu_fn <- mu_quadratic
mu_diff_fn <- mu_diff_quadratic
mu_int_fn <- mu_int_quadratic

neg_log_likelihood_quadratic <- function(parameters, events, delay = 0, kernel, mu_fn = mu_none, 
                                         mu_diff_fn = mu_diff_none, mu_int_fn = mu_int_none, 
                                         print_level = 0) 
{
  names(parameters) <- c("alpha", "delta", "A", "B", "C")
  out <- neg_log_likelihood(parameters, events, delay, kernel, mu_fn, mu_diff_fn, mu_int_fn, print_level)
  return(out)
}

optim_ebola<- DEoptim(neg_log_likelihood_quadratic, lower = c(0,0,0,-10,-10), upper = c(10,10,10,10,10), events = new_times_ebola_mod,
                      kernel = ray_kernel, 
                      mu_fn = mu_fn, 
                      delay = 10,
                      mu_diff_fn = mu_diff_fn,
                      mu_int_fn = mu_int_fn, control = list(parallelType = "parallel"))
# Define the number of starting points
n_start_points <- 15

# Generate starting points
start_points <- as.list(replicate(n_start_points, list(
  alpha = log(sample(2:30, 1)), 
  delta = log(sample(2:30, 1)), 
  A = log(sample(2:30, 1)),
  B = log(sample(2:30, 1)),
  C = log(sample(2:30, 1))), simplify = FALSE))
outtt <- list()
for(i in 1:15){
  outtt[[i]]<-optimx(par = unlist(start_points[[i]]), fn = neg_log_likelihood, gr = ray_derivatives,
       method="BFGS",
       events = new_times_ebola_mod, 
       kernel = ray_kernel,
       delay = 10,
       mu_fn = mu_fn, 
       mu_diff_fn = mu_diff_fn,
       mu_int_fn = mu_int_fn)}
alphamod<- as.numeric(optim_ebola$optim$bestmem[1])
deltamod<-as.numeric(optim_ebola$optim$bestmem[2])
Amod<-as.numeric(optim_ebola$optim$bestmem[3])
Bmod<- as.numeric(optim_ebola$optim$bestmem[4])
Cmod <- as.numeric(optim_ebola$optim$bestmem[5])
library(doParallel)
library(foreach)
no_cores <- detectCores()
registerDoParallel(cores=no_cores)

list_events_ebola_mod <- foreach(i = 1:500, .packages = "epihawkes") %dopar% {
  events1 <- hawkes_simulation(events = c(0), kernel = ray_kernel, 
                               T_max = max(new_times_ebola_mod),
                               parameters = list(alpha = 0.03733332  ,delta = 0.08762098 , 
                                                 A = 0.1352101 ,
                                                 B = 0.09432507 ,
                                                 C = -0.0009531409 ,
                                                 delay = 10), 
                               mu_fn = mu_fn,
                               N_max = length(new_times_ebola_mod),
                               mu_fn_diff = mu_diff_fn,
                               print_level = print_level)
}
N_runs<-500
# Creating a data frame for new_times
df_new_times_ebola <- data.frame(t = new_times_ebola_mod, N = seq(1, length(new_times_ebola_mod)), type = "Real Data")

# Initialize a list to hold all simulated event data frames
list_df_events_ebola <- list()

for(i in 1:N_runs){
  
  events <- list_events_ebola_mod[[i]]
  
  # Creating a data frame for events
  df_events_ebola <- data.frame(t = events, N = seq(1, length(events)), type = paste("Simulated Data ", i))
  
  # Store each data frame in the list
  list_df_events_ebola[[i]] <- df_events_ebola
}

## CUMULATIVE PLOT

# Combine all the data into one dataframe
all_simulations <- do.call(rbind, list_df_events_ebola)

# Convert Simulation column into factor to help with plotting
all_simulations$Simulation <- as.factor(all_simulations$type)

df_new_times_ebola$Simulation <- as.factor(df_new_times_ebola$type)
# Add real data to the same dataframe
all_data <- rbind(all_simulations, df_new_times_ebola)

# Create the plot
plot_ebola <- ggplot() +
  theme_bw() +
  labs(x = "Time (days)",
       y = expression(N(t))) +
  theme(axis.title.x = element_text(size = 20, family="Calibri"),
        axis.title.y = element_text(size = 20, family="Calibri"),
        axis.text = element_text(size = 20, family="Calibri"))

# Add the simulated data lines
plot_ebola <- plot_ebola +
  geom_path(data = subset(all_data, type != "Real Data"), 
            aes(x = t, y = N, group = Simulation), 
            color = "black", alpha = 0.3, size = 0.3)

# Add the real data line
plot_ebola <- plot_ebola +
  geom_path(data = subset(all_data, type == "Real Data"), 
            aes(x = t, y = N), color = "red", size = 1)

print(plot_ebola)



## NON cumulative plot

#true_cum_ebola <- as.POSIXct(new_times_ebola_mod * 86400, origin = "1995-03-06")
#true_cum_ebola1<-floor(as.numeric(difftime(true_cum_ebola, "1995-03-06", units = "weeks")))
#true_tbl <- factor(true_cum_ebola1, levels = 0:18) ##there are 28 weeks
#all_dates <- as.Date(seq(as.Date("1995-03-06"), max(as.Date(ebola_kikwit_1995_mod$date)), by = "day"))
#true_non_cum_days <- as.data.frame(table(as.Date(true_cum_ebola)))
#names(true_non_cum_days) <- c("date", "onset")
#true_non_cum_days$date <- as.Date(true_non_cum_days$date)
## Fill the missing dates with a value of 0
#data_all_true <- merge(data.frame(date = all_dates), true_non_cum_days, by = "date", all.x = TRUE)
# Replace NA values with 0
#data_all_true$onset[is.na(data_all_true$onset)] <- 0

# Initialize an empty list to hold the data frames for each simulation
list_df_non_cum_ebola <- list()
list_df_non_cum_ebola_days <- list()
# Loop over the list of simulations
for(i in 1:N_runs){
  # Convert simulation events from seconds to dates
  non_cum_ebola <- as.POSIXct(list_events_ebola_mod[[i]] * 86400, origin = "1995-03-06")
  # Convert the filtered events to weeks from the origin
  non_cum_ebola1 <- floor(as.numeric(difftime(non_cum_ebola, "1995-03-06", units = "weeks")))
  tbl <- factor(non_cum_ebola1, levels = 0:18) ##there are 28 weeks
  all_dates <- as.Date(seq(as.Date("1995-03-06"), max(as.Date(ebola_kikwit_1995_mod$date)), by = "day"))
  non_cum_days <- as.data.frame(table(as.Date(non_cum_ebola)))
  names(non_cum_days) <- c("date", "onset")
  non_cum_days$date <- as.Date(non_cum_days$date)
  ## Fill the missing dates with a value of 0
  data_all <- merge(data.frame(date = all_dates), non_cum_days, by = "date", all.x = TRUE)
  # Replace NA values with 0
  data_all$onset[is.na(data_all$onset)] <- 0
  
  # Create a data frame for the non-cumulative counts for this simulation
  non_cum_df_ebola <- data.frame(
    Week = as.numeric(names(table(tbl))), 
    Simulated_Count = as.vector(table(tbl)),
    True_count = week_dat_mod$onset
  )
  non_cum_df_ebola_days <- data.frame(
    Days = all_dates,
    True_days = ebola_kikwit_1995_mod$onset,
    Simulated_days = data_all$onset
  )
  # Add the data frame to the list
  list_df_non_cum_ebola[[i]] <- non_cum_df_ebola
  list_df_non_cum_ebola_days[[i]] <- non_cum_df_ebola_days
}

# Prepare the simulated data
for(i in 1:N_runs){
  list_df_non_cum_ebola[[i]]$type <- paste("Simulation ", i)
  list_df_non_cum_ebola_days[[i]]$type <- paste("Simulation ", i)
}

# Combine all the data into one dataframe
all_simulations_non_cum <- do.call(rbind, list_df_non_cum_ebola)

# Convert Simulation column into factor to help with plotting
all_simulations_non_cum$Simulation <- as.factor(all_simulations_non_cum$type)

# Create the plot for WEEKS
plot_non_cum_ebola <- ggplot() +
  theme_bw() +
  labs(title = "Non-Cumulative Counts of Simulated Events and True Observations of Ebola",
       x = "Time (Weeks)",
       y = "Count") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# Add the simulated data lines
plot_non_cum_ebola <- plot_non_cum_ebola +
  geom_path(data = all_simulations_non_cum, 
            aes(x = Week, y = Simulated_Count, group = Simulation), 
            color = "black", alpha = 0.1, size = 0.3)

# Add the real data line
plot_non_cum_ebola <- plot_non_cum_ebola +
  geom_path(data = list_df_non_cum_ebola[[1]], 
            aes(x = Week, y = True_count), color = "red", size = 1)

print(plot_non_cum_ebola)
##ggsave("Ebola_non_cum_plot_weeks.pdf", plot = plot_non_cum_ebola, dpi = 600)

### DAYS

# Create a combined data frame with all the simulation results
all_simulations_df <- do.call(rbind, list_df_non_cum_ebola_days)
all_simulations_df$Simulation <- as.factor(all_simulations_df$type)
all_simulations_df$Num_Days <- as.numeric(difftime(all_simulations_df$Days, min(all_simulations_df$Days), units = "days"))
max_days <- max(all_simulations_df$Num_Days, na.rm = TRUE)
# Create the plot for DAYS
plot_days_ebola <- ggplot() +
  theme_bw() +
  labs(x = "Time (Days)",
       y = "Count") +
  theme(axis.title.x = element_text(size = 20, family="Calibri"),
        axis.title.y = element_text(size = 20, family="Calibri"),
        axis.text = element_text(size = 20, family="Calibri"))


# Add the simulated data lines
plot_days_ebola <- plot_days_ebola +
  geom_line(data = all_simulations_df, aes(x = Num_Days, y = Simulated_days, group = Simulation), color = "black", alpha = 0.1)

# Add the real data line
plot_days_ebola <- plot_days_ebola +
  geom_line(data = all_simulations_df, aes(x = Num_Days, y = True_days), color = "red", size = 0.8)

print(plot_days_ebola)
##ggsave("Ebola_non_cum_plot_days.pdf", plot = plot_days_ebola, dpi = 600)


### Intensities ####

# Calculate intensities for all simulations and real data
list_intensities_ebola_mod <- list()
for(i in 1:N_runs){
  events_ebola_mod <- list_events_ebola_mod[[i]]
  data1 <- compute_intensity_function(events = events_ebola_mod, kernel = ray_kernel, 
                                      T_max = max(new_times_ebola_mod), parameters = list(alpha = 0.03733332  ,delta = 0.08762098 , 
                                                                              A = 0.1352101 ,
                                                                              B = 0.09432507 ,
                                                                              C = -0.0009531409 ,
                                                                              delay = 10), mu_fn = mu_fn, 
                                      N = 5000)
  if (is.null(mu_fn(events_ebola_mod[1], parameters = list(alpha = 0.03733332  ,delta = 0.08762098 , 
                                                           A = 0.1352101 ,
                                                           B = 0.09432507 ,
                                                           C = -0.0009531409 ,
                                                           delay = 10)))) {
    mu_ts1 <- rep(0, length(events_ebola_mod))
  }
  else{
  mu_ts1 <- mu_fn(events_ebola_mod, parameters = list(alpha = 0.03733332  ,delta = 0.08762098 , 
                                                      A = 0.1352101 ,
                                                      B = 0.09432507 ,
                                                      C = -0.0009531409 ,
                                                      delay = 10))
  }
  event_intensities1 <- mu_ts1 + conditional_intensity_list(times = events_ebola_mod+1e-10, 
                                                            events = events_ebola_mod, 
                                                            kernel = ray_kernel, 
                                                            parameters = list(alpha = 0.03733332  ,delta = 0.08762098 , 
                                                                              A = 0.1352101 ,
                                                                              B = 0.09432507 ,
                                                                              C = -0.0009531409 ,
                                                                              delay = 10))
  data_events <- data.frame(t = events_ebola_mod, intensity = event_intensities1, type = paste("Simulated Data ", i))
  list_intensities_ebola_mod[[i]] <- data_events
}

# Calculate intensities for the real data
data_ebola_true_mod <- compute_intensity_function(events = new_times_ebola_mod, kernel = ray_kernel, 
                                               T_max = max(new_times_ebola_mod), parameters = list(alpha = 0.03733332  ,delta = 0.08762098 , 
                                                                                       A = 0.1352101 ,
                                                                                       B = 0.09432507 ,
                                                                                       C = -0.0009531409 ,
                                                                                       delay = 10), mu_fn = mu_fn, 
                                               N = 5000)

mu_ts_true1 <- mu_fn(new_times_ebola_mod, parameters = list(alpha = 0.03733332  ,delta = 0.08762098 , 
                                                            A = 0.1352101 ,
                                                            B = 0.09432507 ,
                                                            C = -0.0009531409 ,
                                                            delay = 10))

event_intensities_true1 <- mu_ts_true1 + conditional_intensity_list(times = new_times_ebola_mod +1e-10, 
                                                                    events = new_times_ebola_mod, 
                                                                    kernel = ray_kernel, 
                                                                    parameters = list(alpha = 0.03733332  ,delta = 0.08762098 , 
                                                                                      A = 0.1352101 ,
                                                                                      B = 0.09432507 ,
                                                                                      C = -0.0009531409 ,
                                                                                      delay = 10))
df_new_times_ebola_mod <- data.frame(t = new_times_ebola_mod, intensity = event_intensities_true1, type = "Real Data")

intensities_all_ebola_mod <- do.call(rbind, list_intensities_ebola_mod)
intensities_all_ebola_mod$Simulation <- as.factor(intensities_all_ebola_mod$type)

# plot the intensities
plot_ebola_mod <- ggplot() +
  theme_bw() +
  labs(x = "Time (days)",
       y = expression(lambda(t))) +
  theme(axis.title.x = element_text(size = 20, family="Calibri"),
        axis.title.y = element_text(size = 20, family="Calibri"),
        axis.text = element_text(size = 20, family="Calibri"))
# Add the simulated data lines
plot_ebola_mod <- plot_ebola_mod +
  geom_path(data = intensities_all_ebola_mod, 
            aes(x = t, y = intensity, group = Simulation), 
            color = "black", alpha = 0.1, size = 0.3)

# Add the real data line
plot_ebola_mod <- plot_ebola_mod +
  geom_path(data = df_new_times_ebola_mod, 
            aes(x = t, y = intensity), color = "red", size = 1)

print(plot_ebola_mod)



#### Goodness of fit ####


events_one_sim_new <- new_times_ebola_mod + 0.0001
new_times_ebola_new <- events_one_sim_new

cumulative_intensities_ebola_new <- sapply(new_times_ebola_new, function(t) {
  integral_intensity(events = new_times_ebola_new[new_times_ebola_new <= t], int_kernel = int_ray, 
                     parameters = list(alpha = 0.03733332,delta = 0.08762098 , 
                                       A = 0.1352101 ,
                                       B = 0.09432507 ,
                                       C = -0.0009531409 ,
                                       delay = 10),
                     mu_fn = mu_fn, mu_diff_fn = mu_diff_fn,
                     mu_int_fn = mu_int_fn)
})


ggplot(data.frame(x = 1:length(new_times_ebola_new) , y=cumulative_intensities_ebola_new), aes(x = x, y = y)) +
  geom_point(size=0.7) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed",size = 0.9) +
  xlab("i") + 
  ylab(expression(Lambda(t[i]))) + 
  ggtitle("Cumulative Intensity vs Event Index") + 
  theme_minimal()

#### intensities using Inter-arrival times for 1st day range
inter_arrival_int_ebola_new <-c()
for(i in 1:length(cumulative_intensities_ebola_new)-1){
  inter_arrival_int_ebola_new<-c(inter_arrival_int_ebola_new, cumulative_intensities_ebola_new[i+1] - cumulative_intensities_ebola_new[i])
}

uk<- 1-exp(-inter_arrival_int_ebola_new)

uk<-sort(uk)
bk<-c()
for(i in 1:length(uk)){
  bk<-c(bk,(i-(1/2))/length(uk))
}
## or could do qunif(ppoints(nrow(df)))

df <- data.frame(CalculatedIntensities = uk, UniformRandomData = bk)

# Calculate confidence intervals
confint_n <- length(df$CalculatedIntensities)
conf_int <- 1.36 / sqrt(confint_n)

# Add upper and lower confidence intervals
df$upperCI <- bk + conf_int
df$lowerCI <- bk - conf_int
# Plot the data using ggplot
library(extrafont)
library(ggplot2)
ggplot(df, aes(x = UniformRandomData, y = CalculatedIntensities)) +
  geom_point(size=0.8) +
  geom_abline(intercept = 0, slope = 1 , color = "red", size = 1) +
  geom_line(aes(y = upperCI), linetype = "dashed", color = "red") +
  geom_line(aes(y = lowerCI), linetype = "dashed", color = "red") +
  labs(x = "Quantiles", 
       y = "Cumulative Distribution Function") +
  theme_bw() +
  theme(text = element_text(size = 20, family = "Calibri"),
        axis.title = element_text(size = 20, family = "Calibri"),
        axis.text = element_text(size = 20, family = "Calibri")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), 
                     limits = c(0, 1))+
  scale_x_continuous(limits = c(0, max(df$UniformRandomData)))

#### QQ plot ####


df$lowerBetaCI <- qbeta(0.025, (1:confint_n), confint_n:1)
df$upperBetaCI <- qbeta(0.975, (1:confint_n), confint_n:1)

ggplot(df, aes(x = UniformRandomData, y = CalculatedIntensities)) +
  geom_point(size=0.8) +
  geom_abline(intercept = 0, slope = 1 , color = "red", size = 1) +
  geom_line(aes(y = upperBetaCI), linetype = "dashed", color = "red") +
  geom_line(aes(y = lowerBetaCI), linetype = "dashed", color = "red") +
  labs(x = "Theoretical Quantiles", 
       y = "Empirical Quantiles") +
  theme_bw() +
  theme(text = element_text(size = 20, family = "Calibri"),
        axis.title = element_text(size = 20, family = "Calibri"),
        axis.text = element_text(size = 20, family = "Calibri")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), 
                     limits = c(0, 1))+
  scale_x_continuous(limits = c(0, max(df$UniformRandomData)))



########### FORECASTING ###########

# Determine the split point for the data
# Forecast the last 5 events
split_point <- 286
# Create the training and test sets
train_times_ebola_mod <- new_times_ebola_mod[1:split_point]
test_times_ebola_mod <- new_times_ebola_mod[(split_point+1):length(new_times_ebola_mod)]


optim_ebola_train <- DEoptim(neg_log_likelihood_quadratic, lower = c(0,0,0,-10,-10), upper = c(10,10,10,10,10), events = train_times_ebola_mod,
                             kernel = ray_kernel, 
                             mu_fn = mu_fn, 
                             delay = 10,
                             mu_diff_fn = mu_diff_fn,
                             mu_int_fn = mu_int_fn, control = list(parallelType = "parallel"))

outtt_train <- list()
for(i in 1:15){
  outtt_train[[i]]<-optimx(par = unlist(start_points[[i]]), fn = neg_log_likelihood, gr = ray_derivatives,
                           method="BFGS",
                           events = train_times_ebola_mod, 
                           kernel = ray_kernel,
                           delay = 10,
                           mu_fn = mu_fn, 
                           mu_diff_fn = mu_diff_fn,
                           mu_int_fn = mu_int_fn)}

library(doParallel)
library(foreach)
no_cores <- detectCores()
registerDoParallel(cores=no_cores)

forecast_events_ebola_mod <- foreach(i = 1:1000, .packages = "epihawkes") %dopar% {
  events1 <- forecast_simulation(events = train_times_ebola_mod, kernel = ray_kernel, 
                               T_max = max(new_times_ebola_mod),
                               parameters = list(alpha = 0.128131,delta = 0.191355 , 
                                                 A =  0.263768 ,
                                                 B = 0.056011  ,
                                                 C = -0.000594 ,
                                                 delay = 10), 
                               mu_fn = mu_fn,
                               mu_fn_diff = mu_diff_fn,
                               print_level = print_level)
}

# Create data frame for real data
df_test <- data.frame(t = c(max(train_times_ebola_mod),test_times_ebola_mod), N = seq(286, 291), type = "Real Data")

# Initialize a list to hold all simulated event data frames
list_df_simulations <- list()

# Loop over each simulation
for(i in 1:length(forecast_events_ebola_mod)) {
  # Create a data frame for each simulation
  df_simulated <- data.frame(t = forecast_events_ebola_mod[[i]], N = seq(286, length(forecast_events_ebola_mod[[i]]) + length(train_times_ebola_mod) - 1), type = paste("Simulated Data", i))
  
  # Store each data frame in the list
  list_df_simulations[[i]] <- df_simulated
}

# Combine all the data into one dataframe
all_data <- do.call(rbind, list_df_simulations)

# Convert Simulation column into factor to help with plotting
all_data$Simulation <- as.factor(all_data$type)

df_test$Simulation <- as.factor(df_test$type)

all_data1 <- rbind(all_data, df_test)


## CUMULATIVE PLOT

# Create the plot
plot_ebola_forecast <- ggplot() +
  theme_bw() +
  labs(x = "Time (days)",
       y = expression(N(t))) +
  theme(text = element_text(size = 22, family = "Calibri"),
        axis.title = element_text(size = 22, family = "Calibri"),
        axis.text = element_text(size = 22, family = "Calibri"))


plot_ebola_forecast <- plot_ebola_forecast +
  geom_path(data = subset(all_data1, type != "Real Data"), 
            aes(x = t, y = N, group = Simulation), 
            color = "black", alpha = 0.1, size = 0.3) +
  geom_path(data = subset(all_data1, type == "Real Data"), 
            aes(x = t, y = N), color = "darkred", size = 1)+
  geom_point(data = subset(all_data1, type == "Real Data"), 
             aes(x = t, y = N), color = "red", shape=20, size=4)

print(plot_ebola_forecast)


# Initialize an empty vector to store RMSE values
RMSEs <- c()

for(i in 1:length(forecast_events_ebola_mod)){
  RMSEs<- c(RMSEs, sqrt(mean(((forecast_events_ebola_mod[[i]][-1])-df_test$t[-1])^2))) #remove the first observation as it is the last point from the training data
}

mean(RMSEs)

sd_RMSEs_ebola <- sd(RMSEs, na.rm = TRUE)


n_RMSEs_ebola <- sum(!is.na(RMSEs))

# Calculate the standard error of the mean RMSEs_covid
sem_RMSEs_ebola <- sd_RMSEs_ebola / sqrt(n_RMSEs_ebola)

# Calculate the mean RMSE
mean_RMSEs_ebola <- mean(RMSEs, na.rm = TRUE)

# Calculate the 95% confidence intervals
CI_lower_ebola <- mean_RMSEs_ebola - 1.96 * sem_RMSEs_ebola
CI_upper_ebola <- mean_RMSEs_ebola + 1.96 * sem_RMSEs_ebola

# Return the confidence interval
c(CI_lower_ebola, CI_upper_ebola)

# Generate R bootstrap replicates
set.seed(123) 
R <- 10000 
results_boot_ebola <- boot(data=na.omit(RMSEs), statistic=mean_fun, R=R)
results_boot_ebola
# Calculate the 95% confidence interval
boot.ci(results_boot_ebola, type="bca")

#### RATIOS ####

mu_divide_int <- mu_ts_true1 / event_intensities_true1
kernel_divide_int <- conditional_intensity_list(times = new_times_ebola_mod +1e-10,
                                                events = new_times_ebola_mod,
                                                kernel = ray_kernel,
                                                parameters = list(alpha = 0.03733332  ,
                                                                  delta = 0.08762098 ,
                                                                  A = 0.1352101 ,
                                                                  B = 0.09432507 ,
                                                                  C = -0.0009531409 ,
                                                                  delay = 10)) / event_intensities_true1



plot_df1 <- data.frame(time = new_times_ebola_mod,
                 mu_divide_int = mu_divide_int,
                 kernel_divide_int = kernel_divide_int)

df_long1 <- reshape2::melt(plot_df1, id.vars = "time")

# Create the plot
ggplot(df_long1, aes(x = time, y = value, color = variable)) +
  geom_line(size = 1) +  
  scale_color_manual(values = c("lightblue", "green")) +
  labs(x = "Time (days)", y = "Ratio") + 
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 22, family = "Calibri"),
        axis.title = element_text(size = 22, family = "Calibri"),
        axis.text = element_text(size = 22, family = "Calibri"))

############################## not included in report ####################################


# Initialize an empty list to hold the data frames for each simulation
#forecast_df_non_cum_ebola <- list()
forecast_df_non_cum_ebola_days <- list()
# Loop over the list of simulations
for(i in 1:1000){
  # Convert simulation events from seconds to dates
  non_cum_ebola <- as.POSIXct(forecast_events_ebola_mod[[i]] * 86400, origin = "1995-03-06")
  true_obs <-as.POSIXct(new_times_ebola_mod * 86400, origin = "1995-03-06")
  # Convert the filtered events to weeks from the origin
  #non_cum_ebola1 <- floor(as.numeric(difftime(non_cum_ebola, "1995-06-17", units = "weeks")))
  all_dates <- as.Date(seq(as.Date("1995-03-06"), max(as.Date(ebola_kikwit_1995_mod$date)), by = "day"))
  non_cum_days <- as.data.frame(table(as.Date(non_cum_ebola)))
  true_cum_days <- as.data.frame(table(as.Date(true_obs)))
  names(non_cum_days) <- c("date", "onset")
  names(true_cum_days) <- c("date", "onset")
  non_cum_days$date <- as.Date(non_cum_days$date)
  true_cum_days$date <- as.Date(true_cum_days$date)
  ## Fill the missing dates with a value of 0
  data_all <- merge(data.frame(date = all_dates), non_cum_days, by = "date", all.x = TRUE)
  # Replace NA values with 0
  data_all$onset[is.na(data_all$onset)] <- 0
  ## Fill the missing dates with a value of 0
  data_all_true <- merge(data.frame(date = all_dates), true_cum_days, by = "date", all.x = TRUE)
  # Replace NA values with 0
  data_all_true$onset[is.na(data_all_true$onset)] <- 0
  
  # Create a data frame for the non-cumulative counts for this simulation
  #non_cum_df_ebola <- data.frame(
   # Week = as.numeric(names(table(tbl))), 
  #  Simulated_Count = as.vector(table(tbl)),
   # True_count = week_dat_mod$onset
  #)
  non_cum_df_ebola_days <- data.frame(
    Days = all_dates,
    True_days = data_all_true$onset,
    Simulated_days = data_all$onset
  )
  # Add the data frame to the list
  #list_df_non_cum_ebola[[i]] <- non_cum_df_ebola
  forecast_df_non_cum_ebola_days[[i]] <- non_cum_df_ebola_days
}

# Prepare the simulated data
for(i in 1:1000){
  #list_df_non_cum_ebola[[i]]$type <- paste("Simulation ", i)
  forecast_df_non_cum_ebola_days[[i]]$type <- paste("Simulation ", i)
}



### DAYS

# Create a combined data frame with all the simulation results
all_forecasts_df <- do.call(rbind, forecast_df_non_cum_ebola_days)
all_forecasts_df$Simulation <- as.factor(all_forecasts_df$type)
all_forecasts_df$Num_Days <- as.numeric(difftime(all_forecasts_df$Days, min(all_forecasts_df$Days), units = "days"))
all_forecasts_df$Days <- as.Date(all_forecasts_df$Days)
# Create the plot for DAYS
plot_days_ebola_forecast <- ggplot() +
  theme_bw() +
  labs(title = "Counts of Simulated Events and True Observations of Ebola",
       x = "Time (Days)",
       y = "Count") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14, family = "Calibri"),
        axis.text = element_text(size = 12, family = "Calibri"))
library(lubridate)
all_forecasts_df<- all_forecasts_df %>% filter(Days >= as.Date("1995-06-18"))

# Add the simulated data lines
plot_days_ebola_forecast <- plot_days_ebola_forecast +
  geom_line(data = all_forecasts_df, aes(x = Num_Days, y = Simulated_days, group = Simulation), color = "black", alpha = 0.3)

# Add the real data line
plot_days_ebola_forecast <- plot_days_ebola_forecast +
  geom_line(data = all_forecasts_df, aes(x = Num_Days, y = True_days), color = "red", size = 0.5)

print(plot_days_ebola_forecast)

#####  Forecasts without N_max specified: #####

N_ebola <- c()
for(i in 1:1000){
  N_ebola <- c(N_ebola, length(forecast_events_ebola_mod[[i]])-1) #the first event is the last from the training data so exclude it
}

sd_N_ebola <- sd(N_ebola)



# Calculate the standard error of the mean
sem_N_ebola <- sd_N_ebola / sqrt(1000)

# Calculate the mean 
mean_N_ebola <- mean(N_ebola)

# Calculate the 95% confidence intervals
CI_lower_N_ebola <- mean_N_ebola - 1.96 * sem_N_ebola
CI_upper_N_ebola <- mean_N_ebola + 1.96 * sem_N_ebola

# Return the confidence interval
c(CI_lower_N_ebola, CI_upper_N_ebola)

# Generate R bootstrap replicates
set.seed(123) 
R <- 10000 
results_boot_ebola <- boot(data=na.omit(N_ebola), statistic=mean_fun, R=R)
results_boot_ebola
# Calculate the 95% confidence interval
boot.ci(results_boot_ebola, type="bca")

