#these libraries need to be loaded
library(utils)
library(dplyr)
library(tidyr)

#read the Dataset sheet into “R”. The dataset will be called "data".
covid <- read.csv("https://opendata.ecdc.europa.eu/covid19/nationalcasedeath_eueea_daily_ei/csv", na.strings = "", fileEncoding = "UTF-8-BOM")
countries<-unique(covid$countriesAndTerritories)
covid <- covid %>% dplyr::filter(countriesAndTerritories == "Malta")
covid$dateRep <- as.Date(covid$dateRep, format="%d/%m/%Y")
covid <- covid %>% dplyr::filter(dateRep >= as.Date("2022-09-01") & dateRep <= as.Date("2022-10-25"))
covid_days <- as.integer(max(covid$dateRep) - min(covid$dateRep)) + 1

all_dates_covid <- seq(min(covid$dateRep), max(covid$dateRep), by = "day")

covid_filled = tidyr::complete(covid, dateRep = all_dates_covid)
covid_filled$cases = abs(ifelse(is.na(covid_filled$cases), 0, covid_filled$cases)) # remove the NAs and turn negative value to positive
covid_filled = dplyr::group_by(covid_filled, dateRep)
covid_filled = dplyr::summarise(covid_filled, total_cases = sum(cases, na.rm = TRUE))


covid_filled <- as.data.frame(covid_filled)
library(extrafont)
library(ggplot2)
ggplot(covid_filled, aes(x = 1:covid_days, y = total_cases)) +
  geom_col(fill = "steelblue", width = 0.7) +
  labs(title = "Plot of cases for COVID-19",
       x = "Number of days",
       y = "Cases") +
  theme(text = element_text(size = 18, family = "Calibri"),  # Increase text size and set font
        axis.text = element_text(size = 18, family = "Calibri"))

ggplot(covid_filled, aes(x = 1:covid_days, y = total_cases)) +
  geom_line() +
  labs(title = "Plot of cases for COVID-19",
       x = "Number of days",
       y = "Cases") +
  theme(text = element_text(size = 18, family = "Calibri"),  # Increase text size and set font
        axis.text = element_text(size = 18, family = "Calibri"))


# Cumulative plot
cum_cases_covid <- cumsum((covid_filled$total_cases))
ggplot(data.frame(total_days = 1:covid_days, cum_cases_covid), aes(x = total_days, y = cum_cases_covid)) +
  geom_step(color = "steelblue", linewidth = 1) +
  labs(title = "Cumulative Plot of Cases for COVID-19",
       x = "Total Days",
       y = "Cumulative Cases") +
  theme(text = element_text(size = 18, family = "Calibri"),  # Increase text size and set font
        axis.text = element_text(size = 18, family = "Calibri"))

# Calculate the start date
start_date <- min(covid_filled$dateRep)

# Add a new column that is the number of weeks since the start date
covid_filled$Week <- floor(as.numeric(difftime(covid_filled$dateRep, start_date, units = "weeks")))

# Summarise the data by week
week_dat_covid <- covid_filled %>%
  dplyr::group_by(Week) %>%
  dplyr::summarise(total_cases  = sum(total_cases))

# Calculate the cumulative cases
week_dat_covid$cum_cases <- cumsum(week_dat_covid$total_cases)

# Plot of Weeks and Cases
ggplot(week_dat_covid, aes(x = Week, y = total_cases)) +
  geom_line(color = "red") +
  labs(title = "Plot of Weeks and Cases for COVID-19",
       x = "Weeks",
       y = "Cases") +
  theme(text = element_text(size = 18, family = "Calibri"),  # Increase text size and set font
        axis.text = element_text(size = 18, family = "Calibri"))

# Cumulative plot
ggplot(week_dat_covid, aes(x = Week, y = cum_cases)) +
  geom_step(color = "steelblue", linewidth = 1) +
  labs(title = "Cumulative Plot of Cases for COVID-19",
       x = "Total Weeks",
       y = "Cumulative Cases") +
  theme(text = element_text(size = 18, family = "Calibri"),  # Increase text size and set font
        axis.text = element_text(size = 18, family = "Calibri"))

ggplot(week_dat_covid, aes(x = Week, y = total_cases)) +
  geom_col(fill = "red", width = 0.3) +
  labs(title = "Plot of Weeks and Cases for COVID-19",
       x = "Weeks",
       y = "Cases") +
  theme(text = element_text(size = 18, family = "Calibri"),  # Increase text size and set font
        axis.text = element_text(size = 18, family = "Calibri"))


############### Optimizing #################


library(lubridate)
library(purrr)
set.seed(50)
# Create a function to generate a sequence of random hours
#generates a number of seconds from 0 (the start of the day) to 246060 - 1 (the end of the day). 
#We then add these seconds to the date in date2 to create a more detailed timestamp in case_epoch. 
#Now, the case_epoch column will have unique timestamps down to the second for each case on each date.
generate_timestamps <- function(n) {
  sample(seq(from = 0, to = 24*60*60 - 1), size = n, replace = TRUE)
}

df_covid = covid_filled
df_covid$date2 = as.POSIXct(df_covid$dateRep)
df_covid$case_time = purrr::map(df_covid$total_cases, generate_timestamps)
df_covid = tidyr::unnest(df_covid, case_time)
df_covid$case_epoch = df_covid$date2 + df_covid$case_time
head(df_covid)



library(tidyr)
df_covid<-unnest(df_covid, case_epoch)
df_covid<- df_covid[order(df_covid$case_epoch),]
library(epihawkes)
library(DEoptim)
min_time <- min(df_covid$case_epoch)
df_covid$case_epoch <- df_covid$case_epoch - min_time
new_times_covid <- c(as.integer(df_covid$case_epoch))/86400 # turn seconds to days

mu_fn <- mu_constant
mu_diff_fn <- mu_diff_constant
mu_int_fn <- mu_int_constant
#delay = 5, https://www.cambridge.org/core/journals/epidemiology-and-infection/article/estimation-of-the-incubation-period-and-generation-time-of-sarscov2-alpha-and-delta-variants-from-contact-tracing-data/29F10A5324C8E6EB1567496D3B154740
# with delay 0, parameters: alpha = 0.8404435, delta= 2.5316182, A = 14.5791273.
# test with delay 5 as well.
covid_optim1 <- DEoptim(neg_log_likelihood_constant, lower = c(0,0,0), 
                      upper = c(10,10,20), 
                      events = new_times_covid, 
                      kernel = exp_kernel, 
                      mu_fn = mu_fn, 
                      delay = 0,
                      mu_diff_fn = mu_diff_fn,
                      mu_int_fn = mu_int_fn, 
                      control = list(itermax = 200, parallelType = "parallel"))

# Define the number of starting points
n_start_points <- 20

# Generate starting points
start_points <- as.list(replicate(n_start_points, list(
  alpha = log(sample(2:30, 1)), 
  delta = log(sample(2:30, 1)), 
  A = log(sample(2:30, 1)),
  B = log(sample(2:30, 1))
  #C = log(sample(2:10, 1))
), simplify = FALSE))

# register the number of cores to use for parallel processing
cl <- makeCluster(detectCores())
registerDoParallel(cl)

npar <- 4 #we have 4 parameters

# rewrite the loop with foreach
outtt_covid <- foreach(i = 1:20, .packages=c('optimx','epihawkes')) %dopar% {
  result <- optimx(par = unlist(start_points[[i]]), fn = my_neg_log_likelihood, gr = transformed_gradients,
                   method="BFGS",
                   events = new_times_covid, 
                   kernel = ray_kernel,
                   delay = 5,
                   mu_fn = mu_fn, 
                   mu_diff_fn = mu_diff_fn,
                   mu_int_fn = mu_int_fn)
  result[1:npar] <-exp(result[1:npar])
  return(result)
}

# stop the cluster
stopCluster(cl)
## parameters 0.093282    0.000000    0.426823
mu_fn <- mu_constant
mu_diff_fn <- mu_diff_constant
mu_int_fn <- mu_int_constant
neg_log_likelihood_none <- function(parameters, events, delay = 0, kernel, mu_fn = mu_none, 
                                      mu_diff_fn = mu_diff_none, mu_int_fn = mu_int_none, 
                                      print_level = 0) 
{
  names(parameters) <- c("alpha", "delta")
  out <- neg_log_likelihood(parameters, events, delay, kernel, mu_fn, mu_diff_fn, mu_int_fn, print_level)
  return(out)
}
covid_optim2 <- DEoptim(neg_log_likelihood_constant, lower = c(0,0,0), 
                        upper = c(20,20,20), 
                        events = new_times_covid[356:3158], #day 50 to end
                        kernel = exp_kernel, 
                        mu_fn = mu_fn, 
                        delay = 0,
                        mu_diff_fn = mu_diff_fn,
                        mu_int_fn = mu_int_fn, 
                        control = list(itermax = 200, parallelType = "parallel"))



N_runs <-500
library(doParallel)
library(foreach)
T_max1_covid = max(new_times_covid)
T_max2_covid = max(new_times_covid)

set.seed(25) 
#alpha_covid <-0.075973
#delta_covid <- 0.045560
#A_covid <- 2.327589 
#B_covid <- -2.647020
alpha_covid1 <- as.numeric(covid_optim1$optim$bestmem[1])
delta_covid1 <- as.numeric(covid_optim1$optim$bestmem[2])
A_covid1<- as.numeric(covid_optim1$optim$bestmem[3])
#alpha_covid2 <- as.numeric(covid_optim2$optim$bestmem[1])
#delta_covid2 <- as.numeric(covid_optim2$optim$bestmem[2])
#A_covid2<- as.numeric(covid_optim2$optim$bestmem[3])
B_covid1<- as.numeric(covid_optim1$optim$bestmem[4])
#C_zika <-as.numeric(zika_optim$optim$bestmem[5])
#M_zika <-as.numeric(zika_optim$optim$bestmem[5])
#N_zika <- as.numeric(zika_optim$optim$bestmem[6])
#P_zika <- as.numeric(zika_optim$optim$bestmem[7])
# Register the parallel backend
no_cores <- detectCores()
registerDoParallel(cores=no_cores)

# Running the simulation in parallel
list_events_covid <- foreach(i = 1:N_runs, .packages = "epihawkes") %dopar% {
  events1 <- hawkes_simulation(events = c(0), kernel = exp_kernel, 
                              T_max = T_max1_covid,
                              parameters = list(alpha =  alpha_covid1,
                                                delta = delta_covid1,
                                                A = A_covid1,
                                                #B = B_covid1,
                                                #C = C_zika,
                                                #M = 18.404159,
                                                #N= -3.054688,
                                                #P = 98.929048,
                                                delay = 0), 
                              mu_fn = mu_fn,
                              mu_fn_diff = mu_diff_fn,
                              N_max = length(new_times_covid),
                              print_level = 1)
  #events2 <- hawkes_simulation(events = max(events1), kernel = exp_kernel, 
   #                            T_max = T_max2_covid,
    #                           parameters = list(alpha =  alpha_covid2,
     #                                            delta = delta_covid2,
      #                                           A = A_covid2,
                                                 #B = B_covid2,
                                                 #C = C_zika,
                                                 #M = 18.404159,
                                                 #N= -3.054688,
                                                 #P = 98.929048,
       #                                          delay = 0), 
        #                       mu_fn = mu_constant,
         #                      mu_fn_diff = mu_diff_constant,
          #                     N_max = length(new_times_covid[355:3158]),
           #                    print_level = 1)
  events1
}


# Creating a data frame for new_times_covid
df_new_times_covid <- data.frame(t = new_times_covid, N = seq(1, length(new_times_covid)), type = "Real Data")

# Initialize a list to hold all simulated event data frames
list_df_events_covid <- list()

for(i in 1:N_runs){
  
  events <- list_events_covid[[i]]
  
  # Creating a data frame for events
  df_events_covid <- data.frame(t = events, N = seq(1, length(events)), type = paste("Simulated Data ", i))
  
  # Store each data frame in the list
  list_df_events_covid[[i]] <- df_events_covid
}

## CUMULATIVE PLOT

# Combine all the data into one dataframe
all_simulations_covid <- do.call(rbind, list_df_events_covid)

# Convert Simulation column into factor to help with plotting
all_simulations_covid$Simulation <- as.factor(all_simulations_covid$type)

df_new_times_covid$Simulation <- as.factor(df_new_times_covid$type)
# Add real data to the same dataframe
all_data_covid <- rbind(all_simulations_covid, df_new_times_covid)

# Create the plot
plot_covid<- ggplot() +
  theme_bw() +
  labs(title = "Comparison of Simulated and Real Data for Zika",
       x = "Time (days)",
       y = expression(N(t))) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# Add the simulated data lines
plot_covid <- plot_covid +
  geom_path(data = subset(all_data_covid, type != "Real Data"), 
            aes(x = t, y = N, group = Simulation), 
            color = "black", alpha = 0.3, size = 0.3)

# Add the real data line
plot_covid <- plot_covid +
  geom_path(data = subset(all_data_covid, type == "Real Data"), 
            aes(x = t, y = N), color = "red", size = 1)

print(plot_covid)

## NON cumulative plot

# Initialize an empty list to hold the data frames for each simulation
list_df_non_cum_covid <- list()
list_df_non_cum_covid_days <- list()
# Loop over the list of simulations
for(i in 1:N_runs){
  # Convert simulation events from seconds to dates
  non_cum <- as.POSIXct(list_events_covid[[i]] * 86400, origin = "2022-09-01")
  
  # Convert the filtered events to weeks from the origin
  non_cum1 <- floor(as.numeric(difftime(non_cum, "2022-09-01", units = "weeks")))
  tbl_covid <- factor(non_cum1, levels = 0:7) ## there are 9 weeks
  all_dates_covid <- as.Date(seq(as.Date("2022-09-01"), max(as.Date(covid$dateRep)), by = "day"))
  non_cum_days_covid <- as.data.frame(table(as.Date(non_cum)))
  names(non_cum_days_covid) <- c("date", "cases")
  non_cum_days_covid$date<- as.Date(non_cum_days_covid$date)
  ## Fill the missing dates with a value of 0
  data_all_covid <- merge(data.frame(date = all_dates_covid), non_cum_days_covid, by = "date", all.x = TRUE)
  # Replace NA values with 0
  data_all_covid$cases[is.na(data_all_covid$cases)] <- 0
  # Create a data frame for the non-cumulative counts for this simulation
  non_cum_df_covid <- data.frame(
    Week = as.numeric(names(table(tbl_covid))), 
    Simulated_Count = as.vector(table(tbl_covid)),
    True_count = week_dat_covid$total_cases
  )
  non_cum_df_covid_days <- data.frame(
    Days = all_dates_covid,
    True_days = covid_filled$total_cases,
    Simulated_days = data_all_covid$cases
  )
  # Add the data frame to the list
  list_df_non_cum_covid[[i]] <- non_cum_df_covid
  list_df_non_cum_covid_days[[i]] <- non_cum_df_covid_days
}

# Prepare the simulated data
for(i in 1:N_runs){
  list_df_non_cum_covid[[i]]$type <- paste("Simulation ", i)
  list_df_non_cum_covid_days[[i]]$type <- paste("Simulation ", i)
}

# Combine all the data into one dataframe
all_simulations_non_cum_covid <- do.call(rbind, list_df_non_cum_covid)

# Convert Simulation column into factor to help with plotting
all_simulations_non_cum_covid$Simulation <- as.factor(all_simulations_non_cum_covid$type)

# Create the plot for WEEKS
plot_non_cum_covid <- ggplot() +
  theme_bw() +
  labs(title = "Counts of Simulated Events and True Observations of Zika",
       x = "Time (Weeks)",
       y = "Count") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# Add the simulated data lines
plot_non_cum_covid <- plot_non_cum_covid +
  geom_path(data = all_simulations_non_cum_covid, 
            aes(x = Week, y = Simulated_Count, group = Simulation), 
            color = "black", alpha = 0.1, size = 0.3)

# Add the real data line
plot_non_cum_covid <- plot_non_cum_covid +
  geom_path(data = list_df_non_cum_covid[[1]], 
            aes(x = Week, y = True_count), color = "red", size = 1)

print(plot_non_cum_covid)
#ggsave("Zika_non_cum_plot_weeks.pdf", plot = plot_non_cum_zika, dpi = 600)

### DAYS

# Create a combined data frame with all the simulation results
all_simulations_df_covid <- do.call(rbind, list_df_non_cum_covid_days)
all_simulations_df_covid$Simulation <- as.factor(all_simulations_df_covid$type)
all_simulations_df_covid$Num_Days <- as.numeric(difftime(all_simulations_df_covid$Days, min(all_simulations_df_covid$Days), units = "days"))
library(extrafont)
# Create the plot for DAYS
plot_days_covid <- ggplot() +
  theme_bw() +
  labs(x = "Time (Days)",
       y = "Count") +
  theme(axis.title = element_text(size = 18, family = "Calibri"),
        axis.text = element_text(size = 19, family = "Calibri"))

# Add the simulated data lines
plot_days_covid <- plot_days_covid +
  geom_line(data = all_simulations_df_covid, aes(x = Num_Days, y = Simulated_days, group = Simulation), color = "black", alpha = 0.1)

# Add the real data line
plot_days_covid <- plot_days_covid+
  geom_line(data = all_simulations_df_covid, aes(x = Num_Days, y = True_days), color = "red", size = 0.3)

print(plot_days_covid)



################## Intensities ################

# Calculate intensities for all simulations and real data
nass_covid<-c()
for(i in 1:N_runs){
  nass_covid<-c(nass_covid,max(length(list_events_covid[[i]])))}
min(nass_covid)
# Calculate intensities for all simulations and real data
list_intensities_covid <- list()

for(i in 1:N_runs){
  events_covid <- list_events_covid[[i]]
  length<- length(events_covid)
  data1 <- compute_intensity_function(events = events_covid[1:355], kernel = exp_kernel, 
                                      T_max = max(events_covid[1:355]), parameters = list(alpha = alpha_covid1,
                                                                                        delta = delta_covid1,
                                                                                        A = A_covid1,
                                                                                        delay = 0), mu_fn = mu_fn, 
                                      N = 5000)
  data2 <- compute_intensity_function(events = events_covid[356:length(events_covid)], kernel = exp_kernel, 
                                      T_max = max(events_covid[356:length(events_covid)]), parameters = list(alpha = alpha_covid2,
                                                                                         delta = delta_covid2,
                                                                                         A = A_covid2,
                                                                                         delay = 0), mu_fn = mu_fn, 
                                      N = 5000)
  mu_ts1 <- mu_fn(events_covid[1:355], parameters = list(alpha = alpha_covid1,
                                                         delta = delta_covid1,
                                                         A = A_covid1,
                                                         delay = 0))
  mu_ts2 <- mu_fn(events_covid[356:length(events_covid)], parameters = list(alpha = alpha_covid2,
                                                                            delta = delta_covid2,
                                                                            A = A_covid2,
                                                                            delay = 0))
  event_intensities1 <- mu_ts1 + conditional_intensity_list(times = events_covid[1:355]+1e-10, 
                                                            events = events_covid[1:355], 
                                                            kernel = exp_kernel, 
                                                            parameters = list(alpha = alpha_covid1,
                                                                              delta = delta_covid1,
                                                                              A = A_covid1,
                                                                              delay = 0))
  event_intensities2 <- mu_ts2 + conditional_intensity_list(times = events_covid[356:length(events_covid)] +1e-10, 
                                                            events = events_covid[356:length(events_covid)], 
                                                            kernel = exp_kernel, 
                                                            parameters = list(alpha = alpha_covid2,
                                                                              delta = delta_covid2,
                                                                              A = A_covid2,
                                                                              delay = 0))
  event_intensities <- c(event_intensities1, event_intensities2)
  data_events <- data.frame(t = events_covid, intensity = event_intensities, type = paste("Simulated Data ", i))
  list_intensities_covid[[i]] <- data_events
}

# Calculate intensities for the real data
data_covid_true1 <- compute_intensity_function(events = new_times_covid[1:355], kernel = exp_kernel, 
                                               T_max = T_max1_covid, parameters = list(alpha = alpha_covid1,
                                                                                 delta = delta_covid1,
                                                                                 A = A_covid1,
                                                                                 delay = 0), mu_fn = mu_fn, 
                                               N = 5000)
data_covid_true2 <- compute_intensity_function(events = new_times_covid[355:3158], kernel = exp_kernel, 
                                               T_max = T_max2_covid, parameters = list(alpha = alpha_covid2,
                                                                                 delta = delta_covid2,
                                                                                 A = A_covid2,
                                                                                 delay = 0), mu_fn = mu_fn, 
                                               N = 5000)

mu_ts_true1 <- mu_fn(new_times_covid[1:355], parameters = list(alpha = alpha_covid1,
                                                               delta = delta_covid1,
                                                               A = A_covid1,
                                                               delay = 0))
mu_ts_true2 <- mu_fn(new_times_covid[356:3158], parameters = list(alpha = alpha_covid2,
                                                                  delta = delta_covid2,
                                                                  A = A_covid2,
                                                                  delay = 0))

event_intensities_true1 <- mu_ts_true1 + conditional_intensity_list(times = new_times_covid[1:355] +1e-10, 
                                                                    events = new_times_covid[1:355], 
                                                                    kernel = exp_kernel, 
                                                                    parameters = list(alpha = alpha_covid1,
                                                                                      delta = delta_covid1,
                                                                                      A = A_covid1,
                                                                                      delay = 0))
event_intensities_true2 <- mu_ts_true2 + conditional_intensity_list(times = new_times_covid[356:3158] +1e-10, 
                                                                    events = new_times_covid[356:3158], 
                                                                    kernel = exp_kernel, 
                                                                    parameters = list(alpha = alpha_covid2,
                                                                                      delta = delta_covid2,
                                                                                      A = A_covid2,
                                                                                      delay = 0))
event_intensities_true <- c(event_intensities_true1, event_intensities_true2)
df_new_times_covid <- data.frame(t = new_times_covid, intensity = event_intensities_true, type = "Real Data")

intensities_all_covid <- do.call(rbind, list_intensities_covid)
intensities_all_covid$Simulation <- as.factor(intensities_all_covid$type)


# plot the intensities
plot_covid_int <- ggplot() +
  theme_bw() +
  labs(title = "Comparison of Simulated and Real Intensities",
       x = "Time (days)",
       y = expression(lambda(t))) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))

# Add the simulated data lines
plot_covid_int <- plot_covid_int +
  geom_path(data = intensities_all_covid, 
            aes(x = t, y = intensity, group = Simulation), 
            color = "black", alpha = 0.1, size = 0.3)

# Add the real data line
plot_covid_int <- plot_covid_int +
  geom_path(data = df_new_times_covid, 
            aes(x = t, y = intensity), color = "red", size = 0.71)

print(plot_covid_int)


############# Goodness of fit ##############

new_times_covid1 <- new_times_covid + 0.001
#new_times_covid2 <- new_times_covid1[1:355] 
cumulative_intensities_covid1 <- sapply(new_times_covid1, function(t) {
  integral_intensity(events = new_times_covid1[new_times_covid1 <= t], int_kernel = int_ray, 
                     parameters = list(alpha = alpha_covid1,
                                       delta = delta_covid1,
                                       A = A_covid1,
                                       #B = B_covid,
                                       #C = C_zika,
                                       #M = 18.404159,
                                       #N= -3.054688,
                                       #P = 98.929048,
                                       delay = 0), mu_fn = mu_fn, mu_diff_fn = mu_diff_fn,
                     mu_int_fn = mu_int_fn)
})

new_times_covid3 <- new_times_covid1[356:3158] 
cumulative_intensities_covid2 <- sapply(new_times_covid3, function(t) {
  integral_intensity(events = new_times_covid3[new_times_covid3 <= t], int_kernel = int_exp, 
                     parameters = list(alpha = alpha_covid2,
                                       delta = delta_covid2,
                                       A = A_covid2,
                                       #B = B_covid,
                                       #C = C_zika,
                                       #M = 18.404159,
                                       #N= -3.054688,
                                       #P = 98.929048,
                                       delay = 0), mu_fn = mu_fn, mu_diff_fn = mu_diff_fn,
                     mu_int_fn = mu_int_fn)
})



ggplot(data.frame(x = 1:length(new_times_covid1) , y = cumulative_intensities_covid1), aes(x = x, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  xlab("i") + 
  ylab(expression(Lambda(t[i]))) + 
  ggtitle("Cumulative Intensity vs Event Index") + 
  theme_minimal()

## good

#### intensities using Inter-arrival times 
inter_arrival_int <-c()
for(i in 1:length(cumulative_intensities_covid1)-1){
  inter_arrival_int<-c(inter_arrival_int, cumulative_intensities_covid1[i+1] - cumulative_intensities_covid1[i])
}


uk<- 1-exp(-inter_arrival_int)

uk<-sort(uk)
bk<-c()
for(i in 1:length(uk)){
  bk<-c(bk,(i-(1/2))/length(uk))
}
## or could do qunif(ppoints(nrow(df)))

# Create a data frame to hold your data
df <- data.frame(CalculatedIntensities = uk, UniformRandomData = bk)

# Calculate confidence intervals
confint_n <- length(df$CalculatedIntensities)
conf_int <- 1.36 / sqrt(confint_n)

# Add upper and lower confidence intervals to your data frame
df$upperCI <- bk + conf_int
df$lowerCI <- bk - conf_int
# Plot the data using ggplot
library(extrafont)
library(ggplot2)
ggplot(df, aes(x = UniformRandomData, y = CalculatedIntensities)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.6) +
  geom_line(aes(y = upperCI), linetype = "dashed", color = "red") +
  geom_line(aes(y = lowerCI), linetype = "dashed", color = "red") +
  labs(x = "Quantiles", 
       y = "Cumulative Distribution Function") +
  theme_minimal() +
  theme(text = element_text(size = 18, family = "Calibri"),
        axis.title = element_text(size = 18, family = "Calibri"),
        axis.text = element_text(size = 18, family = "Calibri"))

#### QQ plot ####


df$lowerBetaCI <- qbeta(0.025, (1:confint_n), confint_n:1 + 1)
df$upperBetaCI <- qbeta(0.975, (1:confint_n), confint_n:1 + 1)

ggplot(df, aes(x = UniformRandomData, y = CalculatedIntensities)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.6) +
  geom_line(aes(y = upperBetaCI), linetype = "dashed", color = "red") +
  geom_line(aes(y = lowerBetaCI), linetype = "dashed", color = "red") +
  labs(x = "Theoretical Quantiles", 
       y = "Empirical Quantiles") +
  theme_minimal() +
  theme(text = element_text(size = 18, family = "Calibri"),
        axis.title = element_text(size = 18, family = "Calibri"),
        axis.text = element_text(size = 18, family = "Calibri"))
