library(readr)
library(ggplot2) 
library(outbreaks)
library(tidyverse)
install.packages("remotes")
remotes::install_github("mrc-ide/epihawkes")
install.packages("DEoptim")
zika_days<-as.integer(max(zika_girardot_2015$date) - min(zika_girardot_2015$date))+1


all_dates <- seq(min(zika_girardot_2015$date), max(zika_girardot_2015$date), by = "day")

# Convert the dates in the original data frame to Date objects
zika_girardot_2015$date <- as.Date(zika_girardot_2015$date)

zika_girardot_filled <- zika_girardot_2015 %>%
  complete(date = all_dates) %>%
  mutate(cases = replace_na(cases, 0))
zika_girardot_filled<- as.data.frame(zika_girardot_filled)

ggplot(zika_girardot_filled, aes(x = as.Date(date), y = cases)) +
  theme_bw() +
  geom_col(fill = "steelblue", color = "black", width = 0.6) +
  labs(x = "Date of symptom onset",
       y = "Count") +
  scale_x_date(date_labels = "%d-%b", date_breaks = "2 weeks", 
  ) +
  theme(
    axis.title.x = element_text(size = 24, family="Calibri"),
    axis.title.y = element_text(size = 24, family="Calibri"),
    axis.text = element_text(size = 24, family="Calibri")
  )
## Cumulative plot now
cum_cases_2 <- cumsum((zika_girardot_filled$cases))
ggplot(data.frame(total_days = 1:zika_days, cum_cases_2), aes(x = total_days, y = cum_cases_2)) +
  geom_step(color = "steelblue", linewidth = 1) +
  labs(title = "Cumulative Plot of Cases for Zika virus",
       x = "Total Days",
       y = "Cumulative Cases")
# Convert date to Date class, if it's not already
zika_girardot_2015$date <- as.Date(zika_girardot_2015$date)

# Calculate the start date
start_date <- min(zika_girardot_2015$date)

# Add a new column that is the number of weeks since the start date
zika_girardot_2015$Week <- floor(as.numeric(difftime(zika_girardot_2015$date, start_date, units = "weeks")))

# Summarise the data by week
week_dat_3 <- zika_girardot_2015 %>%
  group_by(Week) %>%
  summarise(cases = sum(cases))

# Calculate the cumulative cases
week_dat_3$cum_cases <- cumsum(week_dat_3$cases)



# Cumulative plot
ggplot(week_dat_3, aes(x = Week, y = cum_cases)) +
  geom_step(color = "steelblue", linewidth = 1) +
  labs(title = "Cumulative Plot of Cases for Zika virus",
       x = "Total Weeks",
       y = "Cumulative Cases")

ggplot(week_dat_3, aes(x = Week, y = cases)) +
  theme_bw() +
  geom_col(fill = "steelblue", color = "black", width = 0.5) +
  labs(x = "Weeks",
       y = "Count") +
  theme(
    axis.title.x = element_text(size = 24, family="Calibri"),
    axis.title.y = element_text(size = 24, family="Calibri"),
    axis.text = element_text(size = 24, family="Calibri")
  )

library(dplyr)
library(lubridate)
library(purrr)
set.seed(20)
generate_timestamps <- function(n) {
  sample(seq(from = 0, to = 24*60*60 - 1), size = n, replace = TRUE)
}

df2 <- zika_girardot_2015 %>%
  mutate(date2 = as.POSIXct(date)) %>%
  mutate(case_hour = map(cases, generate_timestamps))

library(tidyr)
df2<-unnest(df2, case_hour)
df2$case_epoch <- df2$date2 + df2$case_hour
df2<- df2[order(df2$case_epoch),]
library(epihawkes)
library(DEoptim)
min_time <- min(df2$case_epoch)
df2$case_epoch <- df2$case_epoch - min_time
new_times <- c(as.integer(df2$case_epoch))/86400 # turn seconds to days
mu_fn <- mu_quadratic
mu_diff_fn <- mu_diff_quadratic
mu_int_fn <- mu_int_quadratic
## first attempt with sinusoidal with linearity with rayleigh: worked ok but not great results (delay 22)
## second attempt with sinusoidal with linearity and period with rayleigh: error: check mu function is not increasing (delay 22)
## third attempt with quadratic with rayleigh: same error (delay 22)
## fourth attempt with sinusoidal with linearity with exp: worked ok but not great results (delay 22) (-4066 val)
## fifth attempt with sinusoidal with linearity with exp: worked well results were fairly ok (delay 0) (-4469)
## parameters: alpha = 0.7184937,delta=0.7718912,A=0.8994892,B=0.0001325746,M=0.05607725,N=0.3833896
## sixth attempt with sinusoidal with linearity with exp: worked ok results were not ok (delay 5)
## seventh attempt with sinusoidal with linearity with exp:  not good  (delay 10)
## eigth attempt with sinusoidal with linearity with rayleigh (delay 0) (-4489), error:
##Error in { : task 1 failed - "lower < upper  is not fulfilled"
## ninth attempt: constant mu with rayleigh and delay 22: not good
## tenth attempt: linear mu with rayleigh and delay 22: error
## 11th attempt : quadratic mu with rayleigh and delay 22: same error
## 12th attempt: quadratic mu with rayleigh and delay 0: same error
## 13th attempt: sinusoidal with linearity with exp and delay 0: not very good but ok (-4473.404707 )
## parameters: alpha= 0.680943,delta= 0.714993 ,  A=0.007718, B=0.000185, M= 1.580220 , N=0.018217
## 14th attempt: sinusoidal with linearity with period with exp and delay 0: error: mu increasing

######### Negative constants now (i.e. A, B, M, N.) ############
## 15th: sinusoidal_linear_period with exp and delay 0: ok results
## parameters: alpha= 1.209857,delta= 1.280710 ,  A=2.491259, B=-0.587885, M= 18.653866 , N=19.816728, P = 0.018002
## 16th: sinusoidal_linear_period with rayleigh and delay 0: very bad results

## 16th attempt: sinusoidal with constant with exp and delay 0 : error mu increasing
## 17th attempt: sinusoidal with exp and delay 0: not that good but ok
## 18th attempt: quadratic with exp and delay 0 : error mu increasing
## 19th attempt: linear with exp and delay 0: not good
## 20th attempt: constant with exp and delay 0: not good
## 21st attempt: constant with exp and delay 10: not good, worse than before

## 22nd attempt: sinusoidal mu with rayleigh and delay 0: error mu
## 23rd attempt: sinusoidal mu with exponential and delay 10: not that good
## parameters: alpha = 0.2975260, delta = 0.6415258, M = 17.7091641, N = 0.6893677
## 24th attempt: sinusoidal mu with constant with exponential and delay 10: not that good
## parameters: alpha = 0.3080146, delta = 0.6612765, A=-4.981723, M = 19.99926, N = 5.595824
## 25th attempt: sinusoidal mu with linear with exponential and delay 10: not that good
## parameters: alpha = 0.3123655, delta = 0.6906676, A=-6.80306, B=-0.1451268, M = 19.77705, N = 18.62615
## 26th attempt: constant with exponential and delay 0: not good either
## 27th attempt: sinusoidal with exp and delay 0: not that good
## parameters: $alpha 0.6802088, $delta 0.713373, $M 1.577837,$N 0.01906901,$delay 0
## 28th attempt: sinusoidal with linearity with exp and delay 0: (-4593.755221) worst
## params: alpha=0.824368,    0.837784   18.783695   -2.268760   19.015402   17.876115
## 29th attempt: sinusoidal with linearity and period with exp and delay 0: error

## 30th attempt: sinusoidal with linearity and period with ray and delay 0: worst results
## 31st attempt: constant with ray and delay 0: good
## parameters: alpha = 2.464507, delta = 2.626202, A=1.276636, delay = 0.
## 32nd attempt: quadratic with rayleigh and delay 10: alpha=0.1412493, delta=0.2595605,A=0.1587731, B=1.86735,C=-0.03459858

neg_log_likelihood_quadratic <- function(parameters, events, delay = 0, kernel, mu_fn = mu_none, 
                                        mu_diff_fn = mu_diff_none, mu_int_fn = mu_int_none, 
                                        print_level = 0) 
{
  names(parameters) <- c("alpha", "delta", "A", "B", "C")
  out <- neg_log_likelihood(parameters, events, delay, kernel, mu_fn, mu_diff_fn, mu_int_fn, print_level)
  return(out)
}
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2819875/ , delay = 10
zika_optim <- DEoptim(neg_log_likelihood_quadratic, lower = c(0,0,0,-10,-10), 
                      upper = c(10,10,15,15,15), 
                    events = new_times,
                    kernel = ray_kernel, 
                    mu_fn = mu_fn, 
                    delay = 10,
                    mu_diff_fn = mu_diff_fn,
                    mu_int_fn = mu_int_fn, 
                    control = list(itermax = 200, parallelType = "parallel"))

n_start_points <- 15

start_points <- as.list(replicate(n_start_points, list(
  alpha = log(sample(2:30, 1)), 
  delta = log(sample(2:30, 1)), 
  A = log(sample(2:30, 1)),
  B = log(sample(2:30, 1)),
  C = log(sample(2:30, 1))
  ), simplify = FALSE))

my_neg_log_likelihood <- function(log_params, ...) {
  # Transform the parameters back to their original scale
  params <- exp(log_params)
  # Compute the objective function value
  val <- neg_log_likelihood(params, ...)
  return(val)
}

library(doParallel)
library(foreach)

# register the number of cores to use for parallel processing
cl <- makeCluster(detectCores())
registerDoParallel(cl)


# rewrite the loop with foreach
outtt_zika <- foreach(i = 1:20, .packages=c('optimx','epihawkes')) %dopar% {
  result <- optimx(par = unlist(start_points[[i]]), fn = my_neg_log_likelihood, gr = transformed_gradients,
                   method="BFGS",
                   events = new_times, 
                   kernel = ray_kernel,
                   delay = 10,
                   mu_fn = mu_fn, 
                   mu_diff_fn = mu_diff_fn,
                   mu_int_fn = mu_int_fn)
  result[1:npar] <-exp(result[1:npar])
  return(result)
}

# stop the cluster
stopCluster(cl)

N_runs <- 500
library(doParallel)
library(foreach)
T_max = max(new_times)
set.seed(25)
alpha_zika <- as.numeric(zika_optim$optim$bestmem[1])
delta_zika <- as.numeric(zika_optim$optim$bestmem[2])
A_zika <-as.numeric(zika_optim$optim$bestmem[3])
B_zika <-as.numeric(zika_optim$optim$bestmem[4])
C_zika <-as.numeric(zika_optim$optim$bestmem[5])
#M_zika <-as.numeric(zika_optim$optim$bestmem[5])
#N_zika <- as.numeric(zika_optim$optim$bestmem[6])
#P_zika <- as.numeric(zika_optim$optim$bestmem[7])
# Register the parallel backend
no_cores <- detectCores()
registerDoParallel(cores=no_cores)

# Running the simulation in parallel
list_events <- foreach(i = 1:N_runs, .packages = "epihawkes") %dopar% {
  events <- hawkes_simulation(events = c(0), kernel = ray_kernel, 
                              T_max = max(new_times),
                              parameters = list(alpha = alpha_zika  ,
                                                delta = delta_zika  ,
                                                A = A_zika ,
                                                B = B_zika,
                                                C = C_zika,
                                                #M = 18.404159,
                                                #N= -3.054688,
                                                #P = 98.929048,
                                                delay = 10), 
                              mu_fn = mu_fn,
                              mu_fn_diff = mu_diff_fn,
                              N_max = length(new_times),
                              print_level = 1)
  events
}


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
  labs(x = "Time (days)",
       y = expression(N(t))) +
  theme(axis.title.x = element_text(size = 20, family="Calibri"),
        axis.title.y = element_text(size = 20, family="Calibri"),
        axis.text = element_text(size = 20, family="Calibri"))

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

#save 8 by 6 dimensions landscape

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
  labs(x = "Time (Weeks)",
       y = "Count") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 18, family = "Calibri"),
        axis.text = element_text(size = 18, family = "Calibri"))

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
  labs(x = "Time (days)",
       y = "Count") +
  theme(axis.title.x = element_text(size = 20, family="Calibri"),
        axis.title.y = element_text(size = 20, family="Calibri"),
        axis.text = element_text(size = 20, family="Calibri"))

# Add the simulated data lines
plot_days_zika <- plot_days_zika +
  geom_line(data = all_simulations_df_zika, aes(x = Num_Days, y = Simulated_days, group = Simulation), color = "black", alpha = 0.1)

# Add the real data line
plot_days_zika <- plot_days_zika +
  geom_line(data = all_simulations_df_zika, aes(x = Num_Days, y = True_days), color = "red", size = 0.8)

print(plot_days_zika)


################## Intensities ################

# Calculate intensities for all simulations and real data
list_intensities <- list()
parameters = list(alpha = alpha_zika,
                  delta = delta_zika,
                  A = A_zika,
                  B = B_zika,
                  C = C_zika,
                  #M = 18.404159,
                  #N= -3.054688,
                  #P = 98.929048,
                  delay = 10)
for(i in 1:N_runs){
  events <- list_events[[i]]
  data <- compute_intensity_function(events = events, kernel = ray_kernel, 
                                     T_max = max(new_times), parameters = parameters, mu_fn = mu_fn, 
                                     N = 5000)
  mu_ts <- mu_fn(events, parameters =  parameters)
  event_intensities <- mu_ts + conditional_intensity_list(times = events +1e-10, 
                                                          events = events, 
                                                          kernel = ray_kernel, 
                                                          parameters = parameters)
  data_events <- data.frame(t = events, intensity = event_intensities, type = paste("Simulated Data ", i))
  list_intensities[[i]] <- data_events
}

# Calculate intensities for the real data
data <- compute_intensity_function(events = new_times, kernel =ray_kernel, 
                                   T_max = max(new_times), parameters = parameters, mu_fn = mu_fn, 
                                   N = 5000)
mu_ts <- mu_fn(new_times, parameters =  parameters)
event_intensities <- mu_ts + conditional_intensity_list(times = new_times +1e-10, 
                                                        events = new_times, 
                                                        kernel = ray_kernel, 
                                                        parameters = parameters)
df_new_times <- data.frame(t = new_times, intensity = event_intensities, type = "Real Data")


# plot the intensities
plot <- ggplot() +
  theme_bw() +
  labs(x = "Time (days)",
       y = expression(lambda(t))) +
  theme(axis.title.x = element_text(size = 20, family="Calibri"),
        axis.title.y = element_text(size = 20, family="Calibri"),
        axis.text = element_text(size = 20, family="Calibri"))

# Add the simulated data lines with alpha for transparency
for(i in 1:N_runs){
  plot <- plot + geom_line(data = list_intensities[[i]], aes(x = t, y = intensity), color = "black", alpha = 0.1, size = 0.3)
}

# Add the real data line
plot <- plot + geom_line(data = df_new_times, aes(x = t, y = intensity), color = "red", size = 1)

# Display the plot
print(plot)


############# Goodness of fit ##############

new_times1 <- new_times + 0.0001
cumulative_intensities1 <- sapply(new_times1, function(t) {
  integral_intensity(events = new_times1[new_times1 <= t], int_kernel = int_ray, 
                     parameters = parameters, mu_fn = mu_fn, mu_diff_fn = mu_diff_fn,
                     mu_int_fn = mu_int_fn)
})

compensator<-ggplot(data.frame(x = 1:length(new_times1) , y=cumulative_intensities1), aes(x = x, y = y)) +
  geom_point(size=0.7) +
  theme_bw()+
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed", size=0.9) +
  xlab("i") + 
  ylab(expression(Lambda(t[i])))+
  theme(text = element_text(size = 20, family = "Calibri"),
        axis.title = element_text(size = 20, family = "Calibri"),
        axis.text = element_text(size = 20, family = "Calibri"))

## Very good

#### intensities using Inter-arrival times 
inter_arrival_int <-c()
for(i in 1:length(cumulative_intensities1)-1){
  inter_arrival_int<-c(inter_arrival_int, cumulative_intensities1[i+1] - cumulative_intensities1[i])
}


uk<- 1-exp(-inter_arrival_int)

uk<-sort(uk)
bk<-c()
for(i in 1:length(uk)){
  bk<-c(bk,(i-(1/2))/length(uk))
}
## or could do qunif(ppoints(nrow(df)))

# Create a data frame to hold the data
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
split_point <- 1931
# Create the training and test sets
train_times_zika <- new_times[1:split_point]
test_times_zika <- new_times[(split_point+1):length(new_times)]

optim_zika_train <- DEoptim(neg_log_likelihood_quadratic, lower = c(0,0,0,-10,-10), 
                            upper = c(10,10,15,15,15), events = train_times_zika,
                            kernel = ray_kernel, 
                            mu_fn = mu_fn, 
                            delay = 10,
                            mu_diff_fn = mu_diff_fn,
                            mu_int_fn = mu_int_fn, control = list(parallelType = "parallel"))
## very similar with the original estimated parameters:
## 0.130588    0.235585    0.185065    1.904290   -0.035446
outtt_train_zika <- list()
for(i in 1:15){
  outtt_train_zika[[i]]<-optimx(par = unlist(start_points[[i]]), fn = neg_log_likelihood, gr = exp_derivatives,
                                method="BFGS",
                                events = train_times_zika, 
                                kernel = ray_kernel,
                                delay = 5,
                                mu_fn = mu_fn, 
                                mu_diff_fn = mu_diff_fn,
                                mu_int_fn = mu_int_fn)}

library(doParallel)
library(foreach)
no_cores <- detectCores()
registerDoParallel(cores=no_cores)

forecast_events_zika <- foreach(i = 1:1000, .packages = "epihawkes") %dopar% {
  events1 <- forecast_simulation(events = train_times_zika, kernel = ray_kernel, 
                                   T_max = max(new_times),
                                   parameters = list(alpha = 0.130588 ,delta = 0.235585  , 
                                                     A =  0.185065  , B = 1.904290, C=-0.035446,
                                                     delay = 10), 
                                   mu_fn = mu_fn,
                                   N_max = length(test_times_zika),
                                   mu_fn_diff = mu_diff_fn,
                                   print_level = print_level)
}



# Create data frame for real data
df_test <- data.frame(t = c(max(train_times_zika), test_times_zika), N = seq(1931, 1936), type = "Real Data")

# Initialize a list to hold all simulated event data frames
list_df_simulations <- list()

# Loop over each simulation
for(i in 1:length(forecast_events_zika)) {
  # Create a data frame for each simulation
  df_simulated <- data.frame(t = forecast_events_zika[[i]], N = seq(1931, length(forecast_events_zika[[i]]) + length(train_times_zika) - 1), type = paste("Simulated Data", i))
  
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
plot_zika_forecast <- ggplot() +
  theme_bw() +
  labs(x = "Time (days)",
       y = expression(N(t))) +
  theme(text = element_text(size = 22, family = "Calibri"),
        axis.title = element_text(size = 22, family = "Calibri"),
        axis.text = element_text(size = 22, family = "Calibri"))


plot_zika_forecast <- plot_zika_forecast +
  #geom_path(data = subset(all_data_zika, type == "Real Data"), 
  #          aes(x = t, y = N), color = "red", linewidth = 1) +
  geom_path(data = subset(all_data1, type != "Real Data"), 
            aes(x = t, y = N, group = Simulation), 
            color = "black", alpha = 0.1, size = 0.3) +
  geom_path(data = subset(all_data1, type == "Real Data"), 
            aes(x = t, y = N), color = "darkred", size = 1)+
  geom_point(data = subset(all_data1, type == "Real Data"), 
             aes(x = t, y = N), color = "red", shape=20, size=4)

print(plot_zika_forecast)

# Initialize an empty vector to store RMSE values
RMSEs_zika <- c()

for(i in 1:length(forecast_events_zika)){
  # Get the forecasted events, excluding the first one
  forecasted_events <- forecast_events_zika[[i]][-1]
  
  # Get the test events, excluding the first one as it is the last point from the training data
  test_events <- df_test$t[-1]
  
  # If there are fewer forecasted events than test events, 
  # append the last forecasted event to the forecasted events until they're the same length
  if (length(forecasted_events) < length(test_events)) {
    forecasted_events <- c(forecasted_events, rep(tail(forecasted_events, n = 1), length(test_events) - length(forecasted_events)))
  }
  # Now that forecasted_events and test_events are the same length, calculate the RMSE
  RMSEs_zika <- c(RMSEs_zika, sqrt(mean((forecasted_events - test_events)^2)))
}



sd_RMSEs_zika <- sd(RMSEs_zika, na.rm = TRUE)


n_RMSEs_zika <- sum(!is.na(RMSEs_zika))

# Calculate the standard error of the mean RMSEs_covid
sem_RMSEs_zika <- sd_RMSEs_zika / sqrt(n_RMSEs_zika)

# Calculate the mean RMSE
mean_RMSEs_zika <- mean(RMSEs_zika, na.rm = TRUE)

# Calculate the 95% confidence intervals
CI_lower_zika <- mean_RMSEs_zika - 1.96 * sem_RMSEs_zika
CI_upper_zika <- mean_RMSEs_zika + 1.96 * sem_RMSEs_zika

# Return the confidence interval
c(CI_lower_zika, CI_upper_zika)

# Generate R bootstrap replicates
set.seed(123) 
R <- 10000 
results_boot_zika <- boot(data=na.omit(RMSEs_zika), statistic=mean_fun, R=R)
results_boot_zika
# Calculate the 95% confidence interval
boot.ci(results_boot_zika, type="bca")

#### RATIOS #####
mu_divide_int_zika <- mu_ts / event_intensities
kernel_divide_int_zika <- conditional_intensity_list(times = new_times +1e-10, 
                                                     events = new_times, 
                                                     kernel = ray_kernel, 
                                                     parameters = parameters) / event_intensities



plot_df_zika <- data.frame(time = new_times,
                           mu_divide_int = mu_divide_int_zika,
                           kernel_divide_int = kernel_divide_int_zika)

df_long_zika <- reshape2::melt(plot_df_zika, id.vars = "time")

# Create the plot
ggplot(df_long_zika, aes(x = time, y = value, color = variable)) +
  geom_line(size = 1) +  
  scale_color_manual(values = c("lightblue", "green"), 
                     labels = c(TeX("$\\frac{\\hat{\\mu}(t)}{\\hat{\\lambda}(t)}$"), 
                                TeX("$\\frac{\\hat{\\varphi}(t - t_i)}{\\hat{\\lambda}(t)}$")),
                     name = "Variable") +
  labs(x = "Time (days)", y = "Ratio") + 
  theme_bw() +
  theme(legend.position = "right", text = element_text(size = 20, family = "Calibri"),
        axis.title = element_text(size = 22, family = "Calibri"),
        axis.text = element_text(size = 22, family = "Calibri"))

#### forecasts without N_max specified: #####

no_cores <- detectCores()
registerDoParallel(cores=no_cores)

forecast_events_zika <- foreach(i = 1:1000, .packages = "epihawkes") %dopar% {
  events1 <- forecast_simulation(events = train_times_zika, kernel = ray_kernel, 
                                 T_max = max(new_times),
                                 parameters = list(alpha = 0.130588 ,delta = 0.235585  , 
                                                   A =  0.185065  , B = 1.904290, C=-0.035446,
                                                   delay = 10), 
                                 mu_fn = mu_fn,
                                 #N_max = length(test_times_zika),
                                 mu_fn_diff = mu_diff_fn,
                                 print_level = print_level)
}

N_zika <- c()
for(i in 1:1000){
  N_zika <- c(N_zika, length(forecast_events_zika[[i]])-1) #the first event is the last from the training data so exclude it
}

sd_N_zika <- sd(N_zika)


# Calculate the standard error of the mean
sem_N_zika <- sd_N_zika / sqrt(1000)

# Calculate the mean 
mean_N_zika <- mean(N_zika)

# Calculate the 95% confidence intervals
CI_lower_N_zika <- mean_N_zika - 1.96 * sem_N_zika
CI_upper_N_zika <- mean_N_zika + 1.96 * sem_N_zika

# Return the confidence interval
c(CI_lower_N_zika, CI_upper_N_zika)

# Generate R bootstrap replicates
set.seed(123) 
R <- 10000 
results_boot_zika <- boot(data=na.omit(N_zika), statistic=mean_fun, R=R)
results_boot_zika
# Calculate the 95% confidence interval
boot.ci(results_boot_zika, type="bca")