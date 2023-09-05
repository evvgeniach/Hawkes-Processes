library(epihawkes)
library(DEoptim)
parameters<- list(alpha = 0.8,
                  delta= 1.2,
                  A = 4.5,
                  delay = 2)
mu_term<- "constant"
mu_fn <- mu_constant
mu_fn_diff <- mu_diff_constant
mu_int_fn <- mu_int_constant
N_runs <- 500
# Register the parallel backend
library(doParallel)
library(foreach)
no_cores <- detectCores()
registerDoParallel(cores=no_cores)

# Running the simulation in parallel
events_sim1 <- foreach(i = 1:N_runs, .packages = "epihawkes") %dopar% {
  events <- hawkes_simulation(events = c(0), kernel = ray_kernel, 
                              T_max = 50,
                              parameters = parameters, 
                              mu_fn = mu_fn,
                              mu_fn_diff = mu_diff_fn,
                              N_max = 400,
                              print_level = 1)
  events
}

neg_log_likelihood_constant <- function(parameters, events, delay = 0, kernel, mu_fn = mu_none, 
                                        mu_diff_fn = mu_diff_none, mu_int_fn = mu_int_none, 
                                        print_level = 0) 
{
  names(parameters) <- c("alpha", "delta", "A")
  out <- neg_log_likelihood(parameters, events, delay, kernel, mu_fn, mu_diff_fn, mu_int_fn, print_level)
  return(out)
}

experiment_optim<- list()
for(i in 1:N_runs){
  experiment_optim[[i]]<-DEoptim(neg_log_likelihood_constant, lower = c(0,0,0), upper = c(10,10,10), events = events_sim1[[i]], 
          kernel = ray_kernel, 
          mu_fn = mu_fn, 
          delay = 2,
          mu_diff_fn = mu_diff_fn,
          mu_int_fn = mu_int_fn, control = list(parallelType = "parallel"))
}

experiment_optim_params<- list()
for(i in 1:N_runs){
  experiment_optim_params[[i]]<- experiment_optim[[i]]$optim$bestmem
  names(experiment_optim_params[[i]]) <- c('alpha', 'beta', 'mu')
}


library(ggplot2)
library(tidyr)
library(dplyr)

# First create a dataframe from estimated parameters:
df_estimated <- as.data.frame(do.call(rbind, experiment_optim_params))
df_estimated$run <- 1:N_runs 
df_estimated$Type <- "Estimated"

df_true <- data.frame(run = 1:N_runs,
                      alpha = rep(0.8, N_runs),
                      beta = rep(1.2, N_runs),
                      mu = rep(4.5, N_runs))
df_true$Type <- "True"


df_experiment <- rbind(df_estimated, df_true)

df_long <- gather(df_experiment, Parameter, Value, -run, -Type)

library(extrafont)
# calculate the average estimated parameter values
df_avg <- df_estimated %>%
  summarise(across(c(alpha, beta, mu), mean)) %>%
  pivot_longer(everything(), names_to="Parameter", values_to="Avg_Value")

sim_params_plots <- ggplot(df_long, aes(x=run, y=Value, color=Type)) +
  geom_line(aes(size=Type)) +
  theme_bw() +
  facet_wrap(~Parameter, scales="free_y", labeller = label_parsed) +
  #geom_hline(data=df_avg, aes(yintercept=Avg_Value, color="Average"), linetype="dashed") +
  scale_color_manual(values=c("Estimated" = "darkgreen", "True" = "green"))+#, "Average" = "blue")) +
  scale_size_manual(values=c("Estimated" = 0.7, "True" = 1.5)) +
  labs(x="Simulation", y="Parameter Value", color="Type") +
  theme(text = element_text(size = 34, family = "Calibri"),
        axis.title = element_text(size = 34, family = "Calibri"),
        axis.text = element_text(size = 34, family = "Calibri"),
        legend.text = element_text(size = 34, family = "Calibri"),
        legend.title = element_text(size = 34, family = "Calibri"),
        strip.text = element_text(size = 38, family = "Calibri", margin = margin(14, 0, 14, 0)))

## boxplot
# Convert Parameter to a factor in df_true_hline
df_true_hline <- data.frame(Parameter = c("alpha", "beta", "mu"),
                            TrueValue = c(0.8,1.2,4.5))
df_true_hline$Parameter <- factor(df_true_hline$Parameter, levels = c("alpha", "beta", "mu"))

ggplot(df_long, aes(x=Parameter, y=Value)) +
  geom_boxplot(data=subset(df_long, Type == "Estimated"), alpha=0.6, fill = "lightgreen") +
  geom_hline(data=df_true_hline, aes(yintercept=TrueValue, color="True values"),size=1, linetype="dashed") +
  scale_color_manual(values=c("True values" = "darkred")) +
  facet_grid(. ~ Parameter, scales="free", labeller = label_parsed) +
  theme_bw() +
  labs(x="Parameter", y="Parameter Value", fill="Type", color="Line Type") +
  theme(text = element_text(size = 18, family = "Calibri"),
        axis.title.y = element_text(size = 18, family = "Calibri"),
        axis.text.y = element_text(size = 18, family = "Calibri"),
        strip.placement = "outside",
        axis.line.x = element_blank(),   # Remove x-axis line
        axis.text.x = element_blank(),   # Remove x-axis text
        axis.ticks.x = element_blank(),  # Remove x-axis ticks
        strip.text = element_text(size = 22, family = "Calibri", margin = margin(10, 0, 10, 0)))


# Calculate 95% confidence intervals for estimated parameters
conf_intervals_exp <- df_estimated %>%
  summarise(across(c(alpha, beta, mu), list(lower = ~quantile(., 0.025), upper = ~quantile(., 0.975))))

# Calculate standard deviations for each parameter
sd_alpha <- sd(df_estimated$alpha)
sd_beta <- sd(df_estimated$beta)
sd_mu <- sd(df_estimated$mu)

# Calculate the standard error of the mean for each parameter
sem_alpha <- sd_alpha / sqrt(N_runs)
sem_beta <- sd_beta / sqrt(N_runs)
sem_mu <- sd_mu / sqrt(N_runs)

# Calculate the mean for each parameter
mean_alpha <- mean(df_estimated$alpha)
mean_beta <- mean(df_estimated$beta)
mean_mu <- mean(df_estimated$mu)

# Calculate the 95% confidence intervals for alpha
CI_lower_alpha <- mean_alpha - 1.96 * sem_alpha
CI_upper_alpha <- mean_alpha + 1.96 * sem_alpha

# Calculate the 95% confidence intervals for beta
CI_lower_beta <- mean_beta - 1.96 * sem_beta
CI_upper_beta <- mean_beta + 1.96 * sem_beta

# Calculate the 95% confidence intervals for mu
CI_lower_mu <- mean_mu - 1.96 * sem_mu
CI_upper_mu <- mean_mu + 1.96 * sem_mu

# Return the confidence intervals
list(
  alpha = c(CI_lower_alpha, CI_upper_alpha),
  beta = c(CI_lower_beta, CI_upper_beta),
  mu = c(CI_lower_mu, CI_upper_mu)
)

# Generate R bootstrap replicates
set.seed(123) 
R <- 10000 
results_boot_exp_mu <- boot(data=na.omit(df_estimated$mu), statistic=mean_fun, R=R)
results_boot_exp_mu
# Calculate the 95% confidence interval
boot.ci(results_boot_exp_mu, type="bca")



MSEs_alpha <- c()
MSEs_beta <- c()
MSEs_mu <- c()
for(i in 1:N_runs) {
  estimated_params <- experiment_optim_params[[i]]
  true_params <- c(alpha = 0.8, beta = 1.2, mu = 4.5)
  
  MSE_alpha <- mean((estimated_params['alpha'] - true_params['alpha'])^2)
  MSE_beta <- mean((estimated_params['beta'] - true_params['beta'])^2)
  MSE_mu <- mean((estimated_params['mu'] - true_params['mu'])^2)
  
  MSEs_alpha <- c(MSEs_alpha, MSE_alpha)
  MSEs_beta <- c(MSEs_beta, MSE_beta)
  MSEs_mu <- c(MSEs_mu, MSE_mu)
}
mean_MSE_alpha <- mean(MSEs_alpha, na.rm = TRUE)
sd_MSE_alpha <- sd(MSEs_alpha, na.rm = TRUE)

mean_MSE_beta <- mean(MSEs_beta, na.rm = TRUE)
sd_MSE_beta <- sd(MSEs_beta, na.rm = TRUE)

mean_MSE_mu <- mean(MSEs_mu, na.rm = TRUE)
sd_MSE_mu <- sd(MSEs_mu, na.rm = TRUE)
sem_MSE_alpha <- sd_MSE_alpha / sqrt(N_runs)
sem_MSE_beta <- sd_MSE_beta / sqrt(N_runs)
sem_MSE_mu <- sd_MSE_mu / sqrt(N_runs)
CI_alpha <- c(mean_MSE_alpha - 1.96 * sem_MSE_alpha, mean_MSE_alpha + 1.96 * sem_MSE_alpha)
CI_beta <- c(mean_MSE_beta - 1.96 * sem_MSE_beta, mean_MSE_beta + 1.96 * sem_MSE_beta)
CI_mu <- c(mean_MSE_mu - 1.96 * sem_MSE_mu, mean_MSE_mu + 1.96 * sem_MSE_mu)
list(
  alpha = CI_alpha,
  beta = CI_beta,
  mu = CI_mu
)


### APPENDIX ####

## histogram
ggplot(df_long, aes(x=Value, fill=Type)) +
  geom_histogram(data=subset(df_long, Type == "Estimated"), bins=30, alpha=0.6) +
  geom_vline(data=subset(df_long, Type == "True"), aes(xintercept=Value), color="red", size=0.7, linetype = "dashed") +
  #geom_vline(data=df_avg, aes(xintercept=Avg_Value, color="Average"), linetype="dashed") +
  facet_wrap(~Parameter, scales="free_x") +
  scale_color_manual(values=c("Estimated" = "black", "True" = "red"))+#, "Average" = "blue")) +
  theme_minimal() +
  labs(x="Estimated Parameter Value", y="Count", fill="Type", color="Additional Lines") +
  theme(text = element_text(size = 18, family = "Calibri"),
        axis.title = element_text(size = 18, family = "Calibri"),
        axis.text = element_text(size = 18, family = "Calibri"),
        strip.text = element_text(size = 22, family = "Calibri"))



##########

# Initialize a list to hold all simulated event data frames
list_df_events_experiment <- list()

for(i in 1:N_runs){
  
  events <- events_sim1[[i]]
  
  # Creating a data frame for events
  df_events_experiment <- data.frame(t = events, N = seq(1, length(events)), type = paste("Simulated Data ", i))
  
  # Store each data frame in the list
  list_df_events_experiment[[i]] <- df_events_experiment
}

## CUMULATIVE PLOT

# Combine all the data into one dataframe
all_simulations_experiment <- do.call(rbind, list_df_events_experiment)

# Convert Simulation column into factor to help with plotting
all_simulations_experiment$Simulation <- as.factor(all_simulations_experiment$type)

# Create the plot
plot_experiment<- ggplot() +
  theme_bw() +
  labs(x = "Time (days)",
       y = expression(N(t))) +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.title = element_text(size = 28,family="Calibri"),
        axis.text = element_text(size = 28,family="Calibri"))

# Add the simulated data lines
plot_experiment <- plot_experiment +
  geom_path(data = all_simulations_experiment, 
            aes(x = t, y = N, group = Simulation), 
            color = "black", alpha = 0.3, size = 0.3)



print(plot_experiment)

###############################################################


########################################

new_times_experiment1 <- events_sim1[[1]] + 0.0001
cumulative_intensities_experiment1 <- sapply(new_times_experiment1, function(t) {
  integral_intensity(events = new_times_experiment1[new_times_experiment1 <= t], int_kernel = int_ray, 
                     parameters = parameters, mu_fn = mu_fn, mu_diff_fn = mu_diff_fn,
                     mu_int_fn = mu_int_fn)
})

ggplot(data.frame(x = 1:length(new_times_experiment1) , y=cumulative_intensities_experiment1), aes(x = x, y = y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  xlab("i") + 
  ylab(expression(Lambda(t[i]))) + 
  ggtitle("Cumulative Intensity vs Event Index") + 
  theme_minimal()

## Very good

#### intensities using Inter-arrival times 
inter_arrival_int_experiment <-c()
for(i in 1:length(cumulative_intensities_experiment1)-1){
  inter_arrival_int_experiment<-c(inter_arrival_int_experiment, cumulative_intensities_experiment1[i+1] - cumulative_intensities_experiment1[i])
}


uk_experiment<- 1-exp(-inter_arrival_int_experiment)

uk_experiment<-sort(uk_experiment)
bk_experiment<-c()
for(i in 1:length(uk_experiment)){
  bk_experiment<-c(bk_experiment,(i-(1/2))/length(uk_experiment))
}


# Create a data frame to hold data
df_experiment <- data.frame(CalculatedIntensities = uk_experiment, UniformRandomData = bk_experiment)

# Calculate confidence intervals
confint_n_experiment <- length(df_experiment$CalculatedIntensities)
conf_int_experiment <- 1.36 / sqrt(confint_n_experiment)

# Add upper and lower confidence intervals to data frame
df_experiment$upperCI <- bk_experiment + conf_int_experiment
df_experiment$lowerCI <- bk_experiment - conf_int_experiment
# Plot the data using ggplot
library(extrafont)
library(ggplot2)
ggplot(df_experiment, aes(x = UniformRandomData, y = CalculatedIntensities)) +
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


df_experiment$lowerBetaCI <- qbeta(0.025, (1:confint_n_experiment), confint_n_experiment:1 )
df_experiment$upperBetaCI <- qbeta(0.975, (1:confint_n_experiment), confint_n_experiment:1 )

ggplot(df_experiment, aes(x = UniformRandomData, y = CalculatedIntensities)) +
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


