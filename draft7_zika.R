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

ggplot(zika_girardot_filled, aes(x = 1:zika_days, y = cases)) +
  geom_col(fill = "steelblue", width = 0.7) +
  labs(title = "Plot of cases for Zika virus in Girardot",
       x = "Number of days",
       y = "Cases")
## Cumulative plot now
cum_cases_2 <- cumsum((zika_girardot_filled$cases))
ggplot(data.frame(total_days = 1:zika_days, cum_cases_2), aes(x = total_days, y = cum_cases_2)) +
  geom_step(color = "steelblue", linewidth = 1) +
  labs(title = "Cumulative Plot of Cases for Zika virus",
       x = "Total Days",
       y = "Cumulative Cases")


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

neg_log_likelihood_linear <- function(parameters, events, delay = 0, kernel, mu_fn = mu_none, 
                                        mu_diff_fn = mu_diff_none, mu_int_fn = mu_int_none, 
                                        print_level = 0) 
{
  names(parameters) <- c("alpha", "delta", "A", "B")
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
# Define the number of starting points
n_start_points <- 20

# Generate starting points
start_points <- as.list(replicate(n_start_points, list(
  alpha = log(sample(2:30, 1)), 
  delta = log(sample(2:30, 1)), 
  A = log(sample(2:30, 1)),
  B = log(sample(2:10, 1)),
  C = log(sample(2:10, 1))
  ), simplify = FALSE))


library(doParallel)
library(foreach)

# register the number of cores to use for parallel processing
cl <- makeCluster(detectCores())
registerDoParallel(cl)

npar <- 5 #we have 4 parameters

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

plot_events(list_events[[2]], T_max = max(new_times))
plot_events(new_times, T_max = max(new_times))

