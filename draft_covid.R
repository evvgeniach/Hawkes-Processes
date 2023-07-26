#these libraries need to be loaded
library(utils)
library(dplyr)
library(tidyr)
library(readr)
library(lubridate)

# Parse the dates
date1 <- ymd("2022-09-01")
date2 <- ymd("2022-10-25")

# Calculate the week numbers
week1 <- isoweek(date1) #35
week2 <- isoweek(date2) # 43


variants <- read_csv("data.csv")
variants <- variants %>% dplyr::filter(country == "Malta")
variants <-variants %>% dplyr::filter(year_week >= "2022-35" & year_week <= "2022-43")

# Count the number of occurrences of each variant
variant_counts <- variants %>% count(variant)

most_frequent_variant <- variant_counts %>% arrange(desc(n))

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
ggplot(covid_filled, aes(x = as.Date(dateRep), y = total_cases)) +
  theme_bw() +
  geom_col(fill = "steelblue", color = "black", width = 0.6) +
  labs(x = "Date of case report",
       y = "Count") +
  #scale_x_date(date_labels = "%d-%b", date_breaks = "2 weeks", ) +
  theme(
    axis.title.x = element_text(size = 24, family="Calibri"),
    axis.title.y = element_text(size = 24, family="Calibri"),
    axis.text = element_text(size = 24, family="Calibri")
  )

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
  theme_bw() +
  geom_col(fill = "steelblue", color = "black", width = 0.5) +
  labs(x = "Weeks",
       y = "Count") +
  theme(
    axis.title.x = element_text(size = 24, family="Calibri"),
    axis.title.y = element_text(size = 24, family="Calibri"),
    axis.text = element_text(size = 24, family="Calibri")
  )

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
# constant mu with delay 0, parameters: alpha = 0.8404435, delta= 2.5316182, A = 14.5791273.
# test with delay 5 as well: parameters: alpha =  0.05750214, delta = 0.7811818 ,A = 20.31431, better goodness of fit.
covid_optim1 <- DEoptim(neg_log_likelihood_constant, lower = c(0,0,0), 
                      upper = c(10,100,25), 
                      events = new_times_covid, 
                      kernel = exp_kernel, 
                      mu_fn = mu_fn, 
                      delay = 5,
                      mu_diff_fn = mu_diff_fn,
                      mu_int_fn = mu_int_fn, 
                      control = list(itermax = 200, parallelType = "parallel"))

# Define the number of starting points
n_start_points <- 20

# Generate starting points
start_points <- as.list(replicate(n_start_points, list(
  alpha = log(sample(2:30, 1)), 
  delta = log(sample(2:30, 1)), 
  A = log(sample(2:30, 1))
  #B = log(sample(2:30, 1))
  #C = log(sample(2:10, 1))
), simplify = FALSE))

# register the number of cores to use for parallel processing
cl <- makeCluster(detectCores())
registerDoParallel(cl)

npar <- 3 #we have 3 parameters

# rewrite the loop with foreach
outtt_covid <- foreach(i = 1:20, .packages=c('optimx','epihawkes')) %dopar% {
  result <- optimx(par = unlist(start_points[[i]]), fn = neg_log_likelihood, gr = exp_derivatives,
                   method="BFGS",
                   events = new_times_covid, 
                   kernel = exp_kernel,
                   delay = 5,
                   mu_fn = mu_fn, 
                   mu_diff_fn = mu_diff_fn,
                   mu_int_fn = mu_int_fn)
  #result[1:npar] <-exp(result[1:npar])
  return(result)
}

# stop the cluster
stopCluster(cl)
## parameters 0.093282    0.000000    0.426823


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

# Register the parallel backend
no_cores <- detectCores()
registerDoParallel(cores=no_cores)

# Running the simulation in parallel
list_events_covid <- foreach(i = 1:N_runs, .packages = "epihawkes") %dopar% {
  events1 <- hawkes_simulation(events = c(0), kernel = exp_kernel, 
                              T_max = T_max1_covid,
                              parameters = list(alpha =  0.05750214 ,
                                                delta = 0.7811818 ,
                                                A = 20.31431 ,
                                                delay = 5), 
                              mu_fn = mu_fn,
                              mu_fn_diff = mu_diff_fn,
                              #N_max = length(new_times_covid),
                              print_level = 1)
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
  labs(x = "Time (days)",
       y = expression(N(t))) +
  theme(axis.title.x = element_text(size = 20, family="Calibri"),
        axis.title.y = element_text(size = 20, family="Calibri"),
        axis.text = element_text(size = 20, family="Calibri")) +
  scale_y_continuous(breaks = seq(0, 1500, by = 500), limits = c(0, 1500)) +
  scale_x_continuous(breaks = seq(0, 55, by = 10), limits = c(0, 55))

# Add the simulated data lines
plot_covid <- plot_covid +
  geom_path(data = subset(all_data_covid, type != "Real Data"), 
            aes(x = t, y = N, group = Simulation), 
            color = "black", alpha = 0.3, size = 0.3)

# Add the real data line
plot_covid <- plot_covid +
  geom_path(data = subset(all_data_covid, type == "Real Data"), 
            aes(x = t, y = N), color = "red", linewidth = 1)

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
  labs(x = "Time (days)",
       y = "Count") +
  theme(axis.title.x = element_text(size = 20, family="Calibri"),
        axis.title.y = element_text(size = 20, family="Calibri"),
        axis.text = element_text(size = 20, family="Calibri"))+
  scale_x_continuous(breaks = seq(0, 55, by = 10), limits = c(0, 55))
# Add the simulated data lines
plot_days_covid <- plot_days_covid +
  geom_line(data = all_simulations_df_covid, aes(x = Num_Days, y = Simulated_days, group = Simulation), color = "black", alpha = 0.1)

# Add the real data line
plot_days_covid <- plot_days_covid+
  geom_line(data = all_simulations_df_covid, aes(x = Num_Days, y = True_days), color = "red", size = 0.8)

print(plot_days_covid)



################## Intensities ################
# Calculate intensities for all simulations and real data
list_intensities_covid <- list()

for(i in 1:N_runs){
  events_covid <- list_events_covid[[i]]
  data1 <- compute_intensity_function(events = events_covid, kernel = exp_kernel, 
                                      T_max = T_max1_covid, parameters = list(alpha =  0.05750214 ,
                                                                              delta = 0.7811818 ,
                                                                              A = 20.31431 ,
                                                                              #B = B_covid1,
                                                                              #C = C_zika,
                                                                              #M = 18.404159,
                                                                              #N= -3.054688,
                                                                              #P = 98.929048,
                                                                              delay = 5), mu_fn = mu_fn, 
                                      N = 5000)
  mu_ts1 <- mu_fn(events_covid, parameters = list(alpha =  0.05750214 ,
                                                  delta = 0.7811818 ,
                                                  A = 20.31431 ,
                                                  #B = B_covid1,
                                                  #C = C_zika,
                                                  #M = 18.404159,
                                                  #N= -3.054688,
                                                  #P = 98.929048,
                                                  delay = 5))
  event_intensities1 <- mu_ts1 + conditional_intensity_list(times = events_covid+1e-10, 
                                                            events = events_covid, 
                                                            kernel = exp_kernel, 
                                                            parameters = list(alpha =  0.05750214 ,
                                                                              delta = 0.7811818 ,
                                                                              A = 20.31431 ,
                                                                              #B = B_covid1,
                                                                              #C = C_zika,
                                                                              #M = 18.404159,
                                                                              #N= -3.054688,
                                                                              #P = 98.929048,
                                                                              delay = 5))
  data_events <- data.frame(t = events_covid, intensity = event_intensities1, type = paste("Simulated Data ", i))
  list_intensities_covid[[i]] <- data_events
}

# Calculate intensities for the real data
data_covid_true1 <- compute_intensity_function(events = new_times_covid, kernel = exp_kernel, 
                                               T_max = T_max1_covid, parameters = list(alpha =  0.05750214 ,
                                                                                       delta = 0.7811818 ,
                                                                                       A = 20.31431 ,
                                                                                       #B = B_covid1,
                                                                                       #C = C_zika,
                                                                                       #M = 18.404159,
                                                                                       #N= -3.054688,
                                                                                       #P = 98.929048,
                                                                                       delay = 5), mu_fn = mu_fn, 
                                               N = 5000)

mu_ts_true1 <- mu_fn(new_times_covid, parameters = list(alpha =  0.05750214 ,
                                                        delta = 0.7811818 ,
                                                        A = 20.31431 ,
                                                        #B = B_covid1,
                                                        #C = C_zika,
                                                        #M = 18.404159,
                                                        #N= -3.054688,
                                                        #P = 98.929048,
                                                        delay = 5))

event_intensities_true1 <- mu_ts_true1 + conditional_intensity_list(times = new_times_covid +1e-10, 
                                                                    events = new_times_covid, 
                                                                    kernel = exp_kernel, 
                                                                    parameters = list(alpha =  0.05750214 ,
                                                                                      delta = 0.7811818 ,
                                                                                      A = 20.31431 ,
                                                                                      #B = B_covid1,
                                                                                      #C = C_zika,
                                                                                      #M = 18.404159,
                                                                                      #N= -3.054688,
                                                                                      #P = 98.929048,
                                                                                      delay = 5))
df_new_times_covid <- data.frame(t = new_times_covid, intensity = event_intensities_true1, type = "Real Data")

intensities_all_covid <- do.call(rbind, list_intensities_covid)
intensities_all_covid$Simulation <- as.factor(intensities_all_covid$type)


# plot the intensities
plot_covid_int <- ggplot() +
  theme_bw() +
  labs(x = "Time (days)",
       y = expression(lambda(t))) +
  theme(axis.title.x = element_text(size = 20, family="Calibri"),
        axis.title.y = element_text(size = 20, family="Calibri"),
        axis.text = element_text(size = 20, family="Calibri")) +
  scale_x_continuous(breaks = seq(0, 55, by = 10), limits = c(0, 55))
# Add the simulated data lines
plot_covid_int <- plot_covid_int +
  geom_path(data = intensities_all_covid, 
            aes(x = t, y = intensity, group = Simulation), 
            color = "black", alpha = 0.1, size = 0.3)

# Add the real data line
plot_covid_int <- plot_covid_int +
  geom_path(data = df_new_times_covid, 
            aes(x = t, y = intensity), color = "red", size = 1)

print(plot_covid_int)


############# Goodness of fit ##############

new_times_covid1 <- new_times_covid + 0.001
#new_times_covid2 <- new_times_covid1[1:355] 
cumulative_intensities_covid1 <- sapply(new_times_covid1, function(t) {
  integral_intensity(events = new_times_covid1[new_times_covid1 <= t], int_kernel = int_exp, 
                     parameters = list(alpha =  0.05750214 ,
                                       delta = 0.7811818 ,
                                       A = 20.31431 ,
                                       #B = B_covid1,
                                       #C = C_zika,
                                       #M = 18.404159,
                                       #N= -3.054688,
                                       #P = 98.929048,
                                       delay = 5), mu_fn = mu_fn, mu_diff_fn = mu_diff_fn,
                     mu_int_fn = mu_int_fn)
})





ggplot(data.frame(x = 1:length(new_times_covid1) , y = cumulative_intensities_covid1), aes(x = x, y = y)) +
  geom_point(size=0.7) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size=0.9) +
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
  geom_abline(intercept = 0, slope = 1 , color = "red", size = 1) +
  geom_line(aes(y = upperCI), linetype = "dashed", color = "red") +
  geom_line(aes(y = lowerCI), linetype = "dashed", color = "red") +
  geom_point(size=0.8) +
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
  geom_abline(intercept = 0, slope = 1 , color = "red", size = 1) +
  geom_line(aes(y = upperBetaCI), linetype = "dashed", color = "red") +
  geom_line(aes(y = lowerBetaCI), linetype = "dashed", color = "red") +
  geom_point(size=0.8) +
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
split_point <- 1188
# Create the training and test sets
train_times_covid <- new_times_covid[1:split_point]
test_times_covid <- new_times_covid[(split_point+1):length(new_times_covid)]


optim_covid_train <- DEoptim(neg_log_likelihood_constant, lower = c(0,0,0), upper = c(10,80,25), events = train_times_covid,
                             kernel = exp_kernel, 
                             mu_fn = mu_fn, 
                             delay = 5,
                             mu_diff_fn = mu_diff_fn,
                             mu_int_fn = mu_int_fn, control = list(parallelType = "parallel"))

outtt_train_covid <- list()
for(i in 1:20){
  outtt_train_covid[[i]]<-optimx(par = unlist(start_points[[i]]), fn = neg_log_likelihood, gr = exp_derivatives,
                           method="BFGS",
                           events = train_times_covid, 
                           kernel = exp_kernel,
                           delay = 5,
                           mu_fn = mu_fn, 
                           mu_diff_fn = mu_diff_fn,
                           mu_int_fn = mu_int_fn)}

library(doParallel)
library(foreach)
no_cores <- detectCores()
registerDoParallel(cores=no_cores)

forecast_events_covid <- foreach(i = 1:1000, .packages = "epihawkes") %dopar% {
  events1 <- forecast_simulation(events = train_times_covid, kernel = exp_kernel, 
                               T_max = max(new_times_covid),
                               parameters = list(alpha = 0.09239509 ,delta = 1.48763  , 
                                                 A =  20.46666  ,
                                                 delay = 5), 
                               mu_fn = mu_fn,
                               N_max = length(test_times_covid),
                               mu_fn_diff = mu_diff_fn,
                               print_level = print_level)
}

# Create data frame for real data
df_test <- data.frame(t = c(max(train_times_covid), test_times_covid), N = seq(1188, 1193), type = "Real Data")

# Initialize a list to hold all simulated event data frames
list_df_simulations <- list()

# Loop over each simulation
for(i in 1:length(forecast_events_covid)) {
  if(length(forecast_events_covid[[i]]) > 0){
    # Create a data frame for each simulation
    df_simulated <- data.frame(t = forecast_events_covid[[i]], N = seq(1188, length(forecast_events_covid[[i]]) + length(train_times_covid) - 1), type = paste("Simulated Data", i))
  } else {
    # Create an empty data frame with the appropriate columns if no new events occurred
    df_simulated <- data.frame(t = numeric(), N = numeric(), type = character())
  }
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
plot_covid_forecast <- ggplot() +
  theme_bw() +
  labs(x = "Time (days)",
       y = expression(N(t))) +
  theme(text = element_text(size = 22, family = "Calibri"),
        axis.title = element_text(size = 22, family = "Calibri"),
        axis.text = element_text(size = 22, family = "Calibri"))

plot_covid_forecast <- plot_covid_forecast +
  #geom_path(data = subset(all_data_covid, type == "Real Data"), 
  #          aes(x = t, y = N), color = "red", linewidth = 1) +
  geom_path(data = subset(all_data1, type != "Real Data"), 
            aes(x = t, y = N, group = Simulation), 
            color = "black", alpha = 0.1, size = 0.3) +
  geom_path(data = subset(all_data1, type == "Real Data"), 
            aes(x = t, y = N), color = "darkred", size = 1)+
  geom_point(data = subset(all_data1, type == "Real Data"), 
             aes(x = t, y = N), color = "red", shape=20, size=4)

print(plot_covid_forecast)



# Initialize an empty vector to store RMSE values
RMSEs_covid <- c()

for(i in 1:length(forecast_events_covid)){
  # Get the forecasted events, excluding the first one
  forecasted_events <- forecast_events_covid[[i]][-1]
  
  # Get the test events, excluding the first one as it is the last point from the training data
  test_events <- df_test$t[-1]
  
  # If there are fewer forecasted events than test events, 
  # append the last forecasted event to the forecasted events until they're the same length
  if (length(forecasted_events) < length(test_events)) {
    forecasted_events <- c(forecasted_events, rep(tail(forecasted_events, n = 1), length(test_events) - length(forecasted_events)))
  }
  # Now that forecasted_events and test_events are the same length, calculate the RMSE
  RMSEs_covid <- c(RMSEs_covid, sqrt(mean((forecasted_events - test_events)^2)))
  }

mean(RMSEs_covid, na.rm = TRUE)

# Calculate the standard deviation of RMSEs_covid
sd_RMSEs_covid <- sd(RMSEs_covid, na.rm = TRUE)

# Calculate the number of RMSEs_covid
n_RMSEs_covid <- sum(!is.na(RMSEs_covid))

# Calculate the standard error of the mean RMSEs_covid
sem_RMSEs_covid <- sd_RMSEs_covid / sqrt(n_RMSEs_covid)

# Calculate the mean RMSE
mean_RMSEs_covid <- mean(RMSEs_covid, na.rm = TRUE)

# Calculate the 95% confidence intervals
CI_lower <- mean_RMSEs_covid - 1.96 * sem_RMSEs_covid
CI_upper <- mean_RMSEs_covid + 1.96 * sem_RMSEs_covid

# Return the confidence interval
c(CI_lower, CI_upper)

# Load the boot package
library(boot)

# Define a function to calculate the mean
mean_fun <- function(data, indices) {
  return(mean(data[indices]))
}

# Generate R bootstrap replicates
set.seed(123)  # For reproducibility
R <- 10000  # Choose a number of bootstrap replicates
results_boot <- boot(data=na.omit(RMSEs_covid), statistic=mean_fun, R=R)

# Calculate the 95% confidence interval
boot.ci(results_boot, type="bca")


#####################################################
# Initialize an empty list to hold the data frames for each simulation
forecast_df_non_cum_covid_days <- list()
# Loop over the list of simulations
for(i in 1:1000){
  # Convert simulation events from seconds to dates
  non_cum_covid <- as.POSIXct(forecast_events_covid[[i]] * 86400, origin = "2022-09-01")
  true_obs <-as.POSIXct(test_times_covid * 86400, origin = "2022-09-01")
  # Convert the filtered events to weeks from the origin
  all_dates_covid <- as.Date(seq(as.Date("2022-09-01"), max(as.Date(covid$dateRep)), by = "day"))
  non_cum_days_covid <- as.data.frame(table(as.Date(non_cum_covid)))
  names(non_cum_days_covid) <- c("date", "cases")
  non_cum_days_covid$date<- as.Date(non_cum_days_covid$date)
  ## Fill the missing dates with a value of 0
  data_all_covid <- merge(data.frame(date = all_dates_covid), non_cum_days_covid, by = "date", all.x = TRUE)
  # Replace NA values with 0
  data_all_covid$cases[is.na(data_all_covid$cases)] <- 0
  # True observations now
  true_cum_days <- as.data.frame(table(as.Date(true_obs)))
  names(true_cum_days) <- c("date", "cases")
  non_cum_days_covid$date <- as.Date(non_cum_days_covid$date)
  true_cum_days$date <- as.Date(true_cum_days$date)
  ## Fill the missing dates with a value of 0
  data_all_true <- merge(data.frame(date = all_dates_covid), true_cum_days, by = "date", all.x = TRUE)
  # Replace NA values with 0
  data_all_true$cases[is.na(data_all_true$cases)] <- 0
  
  non_cum_df_covid_days <- data.frame(
    Days = all_dates_covid,
    True_days = data_all_true$cases,
    Simulated_days = data_all_covid$cases
  )
  forecast_df_non_cum_covid_days[[i]] <- non_cum_df_covid_days
}

# Prepare the simulated data
for(i in 1:1000){
  #list_df_non_cum_ebola[[i]]$type <- paste("Simulation ", i)
  forecast_df_non_cum_covid_days[[i]]$type <- paste("Simulation ", i)
}



### DAYS

# Create a combined data frame with all the simulation results
all_forecasts_df <- do.call(rbind, forecast_df_non_cum_covid_days)
all_forecasts_df$Simulation <- as.factor(all_forecasts_df$type)
all_forecasts_df$Num_Days <- as.numeric(difftime(all_forecasts_df$Days, min(all_forecasts_df$Days), units = "days"))
all_forecasts_df$Days <- as.Date(all_forecasts_df$Days)
# Create the plot for DAYS
plot_days_covid_forecast <- ggplot() +
  theme_bw() +
  labs(title = "Counts of Simulated Events and True Observations of Covid",
       x = "Time (Days)",
       y = "Count") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 14, family = "Calibri"),
        axis.text = element_text(size = 12, family = "Calibri"))
library(lubridate)
all_forecasts_df<- all_forecasts_df %>% filter(Days >= as.Date("2022-10-25"))

# Add the simulated data lines
plot_days_covid_forecast <- plot_days_covid_forecast +
  geom_line(data = all_forecasts_df, aes(x = Num_Days, y = Simulated_days, group = Simulation), color = "black", alpha = 0.3)

# Add the real data line
plot_days_covid_forecast <- plot_days_covid_forecast +
  geom_line(data = all_forecasts_df, aes(x = Num_Days, y = True_days), color = "red", size = 0.5)

print(plot_days_covid_forecast)


#### RATIOS #####
mu_divide_int_covid <- mu_ts_true1 / event_intensities_true1
kernel_divide_int_covid <- conditional_intensity_list(times = new_times_covid +1e-10, 
                                                     events = new_times_covid, 
                                                     kernel = exp_kernel, 
                                                     parameters = list(alpha =  0.05750214 ,
                                                                       delta = 0.7811818 ,
                                                                       A = 20.31431 ,
                                                                       delay = 5)) / event_intensities_true1



plot_df_covid <- data.frame(time = new_times_covid,
                           mu_divide_int = mu_divide_int_covid,
                           kernel_divide_int = kernel_divide_int_covid)

df_long_covid <- reshape2::melt(plot_df_covid, id.vars = "time")

# Create the plot
ggplot(df_long_covid, aes(x = time, y = value, color = variable)) +
  geom_line(size = 1) +  
  scale_color_manual(values = c("lightblue", "green")) +
  labs(x = "Time (days)", y = "Ratio") + 
  theme_bw() +
  theme(legend.position = "none", text = element_text(size = 22, family = "Calibri"),
        axis.title = element_text(size = 22, family = "Calibri"),
        axis.text = element_text(size = 22, family = "Calibri"))
