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


# Define the number of starting points
n_start_points <- 30

# Generate starting points
start_points <- as.list(replicate(n_start_points, list(
  alpha = log(sample(2:30, 1)), 
  delta = log(sample(2:30, 1)), 
  A = log(sample(2:30, 1))), simplify = FALSE))


outtt <- list()
for(i in 1:30){
  npar <- 3 #we have 3 parameters
  outtt[[i]]<- optimx(par = unlist(start_points[[i]]), fn = my_neg_log_likelihood, gr = transformed_gradients_exp,
                      method="BFGS",
                      events = events_sim1, 
                      delay = 0,
                      kernel = exp_kernel, 
                      mu_fn = mu_fn, 
                      mu_diff_fn = mu_diff_fn,
                      mu_int_fn = mu_int_fn)
  outtt[[i]][1:npar] <-exp(outtt[[i]][1:npar])
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
df_estimated$run <- 1:N_runs  # assuming N_runs is defined somewhere
df_estimated$Type <- "Estimated"

# Then create a dataframe for true parameters

df_true <- data.frame(run = 1:N_runs,
                      alpha = rep(0.8, N_runs),
                      beta = rep(1.2, N_runs),
                      mu = rep(4.5, N_runs))
df_true$Type <- "True"


# Combine the two dataframes:
df_experiment <- rbind(df_estimated, df_true)

# Reshape to a long format:
df_long <- gather(df_experiment, Parameter, Value, -run, -Type)
library(extrafont)
# First, calculate the average estimated parameter values
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


#bias <- sapply(df_estimated[, c("alpha", "beta", "A")], mean) - c(alpha=0.8, beta=1.2, A=4.5)
#mse <- sapply(df_estimated[, c("alpha", "beta", "A")], function(est) mean((est - c(alpha=0.8, beta=1.2, A=4.5))^2))

#print(paste("Bias for alpha, beta, and A are: ", bias["alpha"], ", ", bias["beta"], ", and ", bias["A"], " respectively."))
#print(paste("MSE for alpha, beta, and A are: ", mse["alpha"], ", ", mse["beta"], ", and ", mse["A"], " respectively."))

#### Confidence intervals ####

# Set the confidence level
conf_level <- 0.95

# Calculate confidence intervals for each parameter
ci_alpha <- quantile(df_estimated$alpha, probs = c((1-conf_level)/2, 1-(1-conf_level)/2))
ci_beta <- quantile(df_estimated$beta, probs = c((1-conf_level)/2, 1-(1-conf_level)/2))
ci_A <- quantile(df_estimated$A, probs = c((1-conf_level)/2, 1-(1-conf_level)/2))

# Calculate coverage probability for each parameter
coverage_alpha <- ifelse(ci_alpha[1] <= 0.8 & ci_alpha[2] >= 0.8, 1, 0)
coverage_beta <- ifelse(ci_beta[1] <= 1.2 & ci_beta[2] >= 1.2, 1, 0)
coverage_A <- ifelse(ci_A[1] <= 4.5 & ci_A[2] >= 4.5, 1, 0)

# Combine into a data frame
coverage_df <- data.frame(
  Parameter = c("alpha", "beta", "A"),
  Lower_CI = c(ci_alpha[1], ci_beta[1], ci_A[1]),
  Upper_CI = c(ci_alpha[2], ci_beta[2], ci_A[2]),
  Coverage = c(coverage_alpha, coverage_beta, coverage_A)
)

print(coverage_df)

##### ADD confidence intervals to initial plot #####

coverage_df_long <- coverage_df %>%
  pivot_longer(cols = c("Lower_CI", "Upper_CI"), names_to = "Type", values_to = "Value") %>%
  mutate(run = if_else(Type == "Lower_CI", 1, N_runs))  # assuming N_runs is defined somewhere

ggplot(df_long, aes(x=run, y=Value, color=Type)) +
  geom_line() +
  facet_wrap(~Parameter, scales="free_y", labeller = label_parsed) +
  geom_hline(data=subset(coverage_df_long, Type == "Lower_CI"), aes(yintercept=Value), linetype="dashed", color="blue") +
  geom_hline(data=subset(coverage_df_long, Type == "Upper_CI"), aes(yintercept=Value), linetype="dashed", color="blue") +
  scale_color_manual(values=c("Estimated" = "lightblue", "True" = "red")) +
  theme_minimal() +
  labs(x="Simulation", y="Parameter Value", color="Type") +
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


run_optimization_ray <- function(start_params) {
  log_start_params <- log(start_params)
  print(start_params)
  print(log_start_params)
  out_new <- optimx(par = log_start_params, fn = my_neg_log_likelihood, gr = transformed_gradients, 
                    events = new_times_ebola, 
                    kernel = ray_kernel, 
                    mu_fn = mu_fn, mu_diff_fn = mu_diff_fn,
                    mu_int_fn = mu_int_fn)
  
  # Transform parameters back to their original scale
  out_new[1:5] <- exp(out_new[1:5])
  return(out_new)
}


results_ray <- lapply(start_points, run_optimization_ray)

DEoptim(fn = neg_log_likelihood_constant, events = events_sim1,
        lower=c(0,0,0,0), upper=c(10,10,10,10),
        kernel = ray_kernel, 
        mu_fn = mu_fn, mu_diff_fn = mu_diff_fn,
        mu_int_fn = mu_int_fn,
        print_level = 0)



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


# Create a data frame to hold your data
df_experiment <- data.frame(CalculatedIntensities = uk_experiment, UniformRandomData = bk_experiment)

# Calculate confidence intervals
confint_n_experiment <- length(df_experiment$CalculatedIntensities)
conf_int_experiment <- 1.36 / sqrt(confint_n_experiment)

# Add upper and lower confidence intervals to your data frame
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


