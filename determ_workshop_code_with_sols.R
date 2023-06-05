# WORKSHOP 1: CALIBRATING A DETERMINISTIC MODEL
# This tutorial offers an interactive introduction to the main functionality of the HMER package, 
# using a simple deterministic epidemiological model. 


###############################################  DEPENDENCIES  ################################################ 

library(hmer)
library(deSolve)
library(ggplot2)
library(reshape2)
library(purrr)
library(tidyverse)
library(lhs)
set.seed(123)


#############################################  HELPER FUNCTIONS  ############################################## 

# Define the helper function `ode_results`, to obtain the solution of the ODEs at any given parameter set.
ode_results <- function(parms, end_time = 365*2) {
  forcer = matrix(c(0, parms['beta1'], 100, parms['beta2'], 180, parms['beta2'], 270, parms['beta3']), 
                  ncol = 2, byrow = TRUE)
  force_func = approxfun(x = forcer[,1], y = forcer[,2], method = "linear", rule = 2)
  des = function(time, state, parms) {
    with(as.list(c(state, parms)), {
      dS <- b*(S+E+I+R)-force_func(time)*I*S/(S+E+I+R) + omega*R - mu*S
      dE <- force_func(time)*I*S/(S+E+I+R) - epsilon*E - mu*E
      dI <- epsilon*E - alpha*I - gamma*I - mu*I
      dR <- gamma*I - omega*R - mu*R
      return(list(c(dS, dE, dI, dR)))
    })
  }
  yini = c(S = 900, E = 100, I = 0, R = 0)
  times = seq(0, end_time, by = 1)
  out = deSolve::ode(yini, times, des, parms)
  return(out)
}

# Define the helper function `get_results` that acts as `ode_results`, but has the additional feature 
# of allowing us to decide which outputs and times should be returned
get_results <- function(params, times, outputs) {
  t_max <- max(times)
  all_res <- ode_results(params, t_max)
  actual_res <- all_res[all_res[,'time'] %in% times, c('time', outputs)]
  shaped <- reshape2::melt(actual_res[,outputs])
  return(setNames(shaped$value, paste0(shaped$Var2, actual_res[,'time'], sep = "")))
}

# Define the helper function `plot_runs` which takes a data frame containing 180 parameter sets and plots the trajectories for
# S, E, I, R for them.
plot_runs <- function(points){
  sol <- list()
  for (i in 1:nrow(points)) sol[[i]]<-ode_results(points[i,])
  par(mar = c(3, 3, 3, 3))
  do.call(plot, sol)
  }


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  
########################################  2. INTRODUCTION TO THE MODEL  #######################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Set some parameter values, run the model on them and plot the model output
example_params <- c(
  b = 1/(60*365),
  mu = 1/(76*365),
  beta1 = 0.2, beta2 = 0.1, beta3 = 0.3,
  epsilon = 0.13,
  alpha = 0.01,
  gamma = 0.08,
  omega = 0.003
)
solution <- ode_results(example_params)
par(mar = c(2, 2, 2, 2))
plot(solution)


### Solution to the task on familiarising yourself with the model ###

# Increase the value of the force of infection parameters (beta 1 and beta 3), run the model, and plot the new 
# solution on top of the previous one
higher_foi_params <- c(
b = 1/(60*365),
mu = 1/(76*365),
beta1 = 0.3, beta2 = 0.1, beta3 = 0.5,
epsilon = 0.13,
alpha = 0.01,
gamma = 0.08,
omega = 0.003
)
higher_foi_solution <- ode_results(higher_foi_params)
plot(solution, higher_foi_solution)

# Increase the rate of becoming infectious (sigma), run the model, and plot the new solution on top of the previous one
higher_epsilon_params <- c(
b = 1/(60*365),
mu = 1/(76*365),
beta1 = 0.3, beta2 = 0.1, beta3 = 0.5,
epsilon = 0.21,
alpha = 0.01,
gamma = 0.08,
omega = 0.003
)
higher_epsilon_solution <- ode_results(higher_epsilon_params)
plot(higher_foi_solution,higher_epsilon_solution)

# Decrease the recovery rate (gamma), run the model, and plot the new solution on top of the previous one
smaller_recovery_params <- c(
b = 1/(60*365),
mu = 1/(76*365),
beta1 = 0.3, beta2 = 0.1, beta3 = 0.5,
epsilon = 0.21,
alpha = 0.01,
gamma = 0.05,
omega = 0.003
)
smaller_recovery_solution <- ode_results(smaller_recovery_params)
plot(higher_epsilon_solution, smaller_recovery_solution)

### End of the solution ###



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############################  3. WAVE0 - PARAMETER RANGES, TARGETS AND DESIGN POINTS  ###########################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Define the parameter ranges 
ranges = list(
  b = c(1e-5, 1e-4), # birth rate
  mu = c(1e-5, 1e-4), # rate of death from other causes
  beta1 = c(0.2, 0.3), # infection rate from time t=0
  beta2 = c(0.1, 0.2), # infection rate from time t=100
  beta3 = c(0.3, 0.5), # infection rate from time t=270
  epsilon = c(0.07, 0.21), # rate of becoming infectious after infection
  alpha = c(0.01, 0.025), # rate of death from the disease
  gamma = c(0.05, 0.08), # recovery rate
  omega = c(0.002, 0.004) # rate at which immunity is lost following recovery
)

# Define the targets' mean values and standard deviations
targets <- list(
  I25 = list(val = 115.88, sigma = 5.79),
  I40 = list(val = 137.84, sigma = 6.89),
  I100 = list(val = 26.34, sigma = 1.317),
  I200 = list(val = 0.68, sigma = 0.034),
  I300 = list(val = 29.55, sigma = 1.48),
  I350 = list(val = 68.89, sigma = 3.44),
  R25 = list(val = 125.12, sigma = 6.26),
  R40 = list(val = 256.80, sigma = 12.84),
  R100 = list(val = 538.99, sigma = 26.95),
  R200 = list(val = 444.23, sigma = 22.21),
  R300 = list(val = 371.08, sigma = 15.85),
  R350 = list(val = 549.42, sigma = 27.47)
)

# The parameter set below was used to determine the targets. The model was run on it and the outputs taken to 
# be the mean value of the targets. The standard deviations were defined as 5% of the corresponding mean value.
chosen_params <- list(b = 1/(76*365), mu = 1/(76*365), beta1 = 0.214, 
                      beta2 = 0.107, beta3 = 0.428, epsilon = 1/7, alpha = 1/50, gamma = 1/14, omega = 1/365)

# Define two Latin hypercube designs through the function `maximinLHS`. This function assumes that each parameter 
# is distributed on [0,1]
initial_LHS_training <- maximinLHS(90, 9)
initial_LHS_validation <- maximinLHS(90, 9)
initial_LHS <- rbind(initial_LHS_training, initial_LHS_validation)

# Rescale the parameter ranges from [0,1] to the correct ranges, and add columns names to identify the parameters
initial_points <- setNames(data.frame(t(apply(initial_LHS, 1, 
                  function(x) x*unlist(lapply(ranges, function(x) x[2]-x[1])) + 
                  unlist(lapply(ranges, function(x) x[1]))))), names(ranges))

# Run the model on `initial_points` and add column names to identify the different targets 
initial_results <- setNames(data.frame(t(apply(initial_points, 1, get_results, 
                   c(25, 40, 100, 200, 300, 350), c('I', 'R')))), names(targets))

# Bind `initial_points` and the corresponding model outputs `initial_results` by column
wave0 <- cbind(initial_points, initial_results)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#################################################  4. EMULATORS  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Split the dataframe `wave0` into a training and a validation set
training <- wave0[1:90,]
validation <- wave0[91:180,]

# # # # # # # # # # # # # # # # # # # #  4.2 TRAINING EMULATORS  # # # # # # # # # # # # # # # # # # # # # # #

# Train the first set of emulators using the function `emulator_from_data`
ems_wave1 <- emulator_from_data(training, names(targets), ranges, 
                                c_lengths= rep(0.55,length(targets)))
# Show the emulator specification for the number of recovered people at t=200
ems_wave1$R200

# Plot the emulator expectation for the number of recovered people at t=200 in the (beta1,gamma)-plane
emulator_plot(ems_wave1$R200, params = c('beta1', 'gamma'))


###  Solution to the task on active variables ###

# Show what variables are active for the emulator of the number of recovered people at t=200
ems_wave1$R200$active_vars

# Shows whether beta_3, the fifth parameter, is active
ems_wave1$R200$active_vars[5]

# Loop through `ems_wave1` to look at the role of beta_3 in each of the emulators 
beta3_role <- logical()
for (i in 1:length(ems_wave1)) beta3_role[i] <- ems_wave1[[i]]$active_vars[5]
names(beta3_role) <- names(ems_wave1)
beta3_role

### End of the solution ###


# Produce a plot showing what variables are active for different emulators, using the function `plot_actives`
plot_actives(ems_wave1)

# Plot the emulation variance for the number of recovered people at t=200 in the (beta1,gamma)-plane
emulator_plot(ems_wave1$R200, plot_type = 'var', params = c('beta1', 'gamma'))

# Examine the adjusted R^2 squared of the regression hypersurface for the number of recovered people at t=200
summary(ems_wave1$R200$model)$adj.r.squared



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
##############################################  5. IMPLAUSIBILITY  ############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#Plot the emulator implausibility for the number of recovered people at t=200 in the (beta1,gamma)-plane
emulator_plot(ems_wave1$R200, plot_type = 'imp', 
              targets = targets, params = c('beta1', 'gamma'), cb=TRUE)

# Plot the emulator implausibility for all emulators in the (beta1,gamma)-plane
emulator_plot(ems_wave1, plot_type = 'imp', 
              targets = targets, params = c('beta1', 'gamma'), cb=TRUE)


### Solution to the task on the functionalities of `emulator_plot` ###

# Visualise the maximum implausibility passing all emulators to `emulator_plot` and setting `plot_type='nimp'`
emulator_plot(ems_wave1, plot_type = 'nimp', targets = targets, params = c('beta1', 'gamma'), cb=T)

# Visualise the second-maximum implausibility adding nth=2 to the previous function call
emulator_plot(ems_wave1, plot_type = 'nimp', targets = targets, params = c('beta1', 'gamma'), cb=T, nth=2)

# Focus on early times outputs (up to $t=200$), and produce implausibility plots for them
restricted_ems <- ems_wave1[c(1,2,3,4,7,8,9,10)]
emulator_plot(restricted_ems, plot_type = 'imp', targets = targets, params = c('beta1', 'gamma'), cb=T)

# Set the unshown parameters to be as in `chosen_params`:
emulator_plot(restricted_ems, plot_type = 'nimp', targets = targets[c(1,2,3,4,7,8,9,10)],
              params = c('beta1', 'gamma'),
              fixed_vals = chosen_params[!names(chosen_params) %in% c('beta1', 'gamma')], cb=T)+
              geom_point(aes(x=0.214, y=1/14), size=3)

### End of the solution ###



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###########################################  6. EMULATOR DIAGNOSTICS  #########################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Produce three diagnostics of the emulator for the number of recovered people at t=200 using 
#`validation_diagnostics`
vd <- validation_diagnostics(ems_wave1$R200, validation = validation, targets = targets, plt=TRUE)

# Double the value of sigma in the emulator for the number of recovered people at t=200
sigmadoubled_emulator <- ems_wave1$R200$mult_sigma(2)

# Produce the three diagnostics again
vd <- validation_diagnostics(sigmadoubled_emulator,
                             validation = validation, targets = targets, plt=TRUE)


### Solution to the task on varying the value of sigma  ###

# Set sigma to be 10 times smaller than its default value
tinysigma_emulator <- ems_wave1$R200$mult_sigma(0.1)

# Produce the three validation diagnostics 
vd <- validation_diagnostics(tinysigma_emulator, validation = validation, targets = targets, plt=TRUE)

# Set sigma to be ten times larger than its default value
hugesigma_emulator <- ems_wave1$R200$mult_sigma(10)

# Produce the three validation diagnostics
vd <- validation_diagnostics(hugesigma_emulator, validation = validation, targets = targets, plt=TRUE)

### End of the solution ###



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############################################  7. PROPOSING NEW POINTS  ########################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Generate 180 new parameter sets using `generate_new_design`
new_points_restricted <- generate_new_design(restricted_ems, 180, targets, verbose=TRUE)

# Plot the new parameter sets using `plot_wrap`
plot_wrap(new_points_restricted, ranges)


### Solution to the task on generating new points using all emulators ###

# Pass the list of wave 1 emulators to `generate_new_design` and plot them using `plot_wrap`
new_points <- generate_new_design(ems_wave1, 180, targets, verbose = TRUE)
plot_wrap(new_points, ranges)

### End of the solution ###


# Pass the list of wave 1 emulators to `space_removed` to quantify how much of the input space has been cut out.
# Here we set ppd to 2, to speed the process up. 
space_removed(ems_wave1, targets, ppd=2) + geom_vline(xintercept = 3, lty = 2) + 
  geom_text(aes(x=3, label="x = 3",y=0.33), colour="black", 
  angle=90, vjust = 1.2, text=element_text(size=11))



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#########################################  8. CUSTOMISE THE FIRST WAVE  #######################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


### Solution to the task on customising the first wave ###

# Produce diagnostics for all the emulators trained in wave 1
vd <- validation_diagnostics(ems_wave1, validation = validation, targets = targets, plt=TRUE)

# Define a vector indicating the factor by which each sigma should be multiplied
inflations <- c(1.5,1.5,1,1,1,1,1.5,1.5,2,2,1.5,2)

# Loop through the vector `inflations` and modify the sigma values of emulators
for (i in 1:length(ems_wave1)) {
ems_wave1[[i]] <- ems_wave1[[i]]$mult_sigma(inflations[[i]])
}

# Produce diagnostics for all the emulators trained in wave 1, with the modified sigmas
vd <- validation_diagnostics(ems_wave1, validation = validation, targets = targets, plt=TRUE)

# Generate 180 new parameter sets using the modified emulators and plot them
new_points <- generate_new_design(ems_wave1, 180, targets, verbose = TRUE)
plot_wrap(new_points, ranges)

### End of the solution ###



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###############################################  9. SECOND WAVE  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


###  Solution to the task on training and customising new emulators ###

# Start by evaluating the function `get_results` on `new_points` 
new_initial_results <- setNames(data.frame(t(apply(new_points, 1, get_results,
                                            c(25, 40, 100, 200, 300, 350), c('I', 'R')))), names(targets))

# Bind by columns `new_points` to the model output `new_initial_results` to create the data.frame `wave1`
wave1 <- cbind(new_points, new_initial_results)

# Split `wave1` into training and validation sets
new_t_sample <- sample(1:nrow(wave1), 90)
new_training <- wave1[new_t_sample,]
new_validation <- wave1[-new_t_sample,]

# Train wave 2 emulators using `emulator_from_data`, passing the new ranges to it
ems_wave2 <- emulator_from_data(new_training, names(targets), ranges, check.ranges = TRUE, c_lengths= rep(0.55,length(targets)))
# Produce diagnostics for all wave 2 emulators

vd <- validation_diagnostics(ems_wave2, validation = new_validation, targets = targets, plt=TRUE)
# Define a vector indicating the factor by which each sigma should be multiplied
inflations <- c(2,4,2,2,2,1,3,2,2,2,2,2)
for (i in 1:length(ems_wave2)) {
ems_wave2[[i]] <- ems_wave2[[i]]$mult_sigma(inflations[[i]])
}

# Produce diagnostics for the modified wave 2 emulators 
vd <- validation_diagnostics(ems_wave2, validation =  new_validation, targets = targets, plt=TRUE)

# Generate 180 new parameter sets and plot them
new_new_points <- generate_new_design(c(ems_wave2, ems_wave1), 180, targets, verbose=TRUE)
plot_wrap(new_new_points, ranges)

### End of the solution ###


# Examine the adjusted R^2 squared of the regression hypersurface for all wave 2 emulators
R_squared_new <- list()
for (i in 1:length(ems_wave2)) {
  R_squared_new[[i]] <- summary(ems_wave2[[i]]$model)$adj.r.squared
}
names(R_squared_new) <- names(ems_wave2)
unlist(R_squared_new)

# Train new wave 1 emulators, setting `quadratic=FALSE` to assume a linear regression term  
ems_wave1_linear <- emulator_from_data(training, names(targets), 
                                        ranges, quadratic=FALSE, c_lengths= rep(0.55,length(targets)))
# Examine the adjusted R^2 squared of the linear regression hypersurface for these new emulators
R_squared_linear <- list()
 for (i in 1:length(ems_wave1_linear)) {
   R_squared_linear[[i]] <- summary(ems_wave1_linear[[i]]$model)$adj.r.squared
 }
 names(R_squared_linear) <- names(ems_wave1_linear)
 unlist(R_squared_linear)

# Plot the emulator variance for the number of infectious people at t=200 in the (beta1,gamma)-plane
emulator_plot(ems_wave1_linear$I200, plot_type = 'var', 
               params = c('beta1', 'gamma'))

# Plot the emulator implausibility for the number of infectious people at t=200 in the (beta1,gamma)-plane
emulator_plot(ems_wave1_linear$I200, plot_type = 'imp', targets = targets, 
              params = c('beta1', 'gamma'), cb=TRUE)

# Increase theta by a factor of three in the linear emulator for the number of infectious people at t=200 
ems_wave1_linear$I200 <- ems_wave1_linear$I20$set_hyperparams(
              list(theta=ems_wave1_linear$I200$corr$hyper_p$theta *3))
# Plot the variance of the modified linear emulator for the number of infectious people at t=200 in the 
# (beta1,gamma)-plane
emulator_plot(ems_wave1_linear$I200, plot_type = 'var', 
              params = c('beta1', 'gamma'))

# Plot the implausibility of the modified linear emulator for the number of infectious people at t=200 in 
# the (beta1,gamma)-plane
emulator_plot(ems_wave1_linear$I200, plot_type = 'imp', targets = targets, 
              params = c('beta1', 'gamma'), cb=TRUE)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#############################  10. VISUALISATIONS OF NON-IMPLAUSIBLE SPACE BY WAVE  ###########################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Show the distribution of the non-implausible space before the wave 1, at the end of wave 1 and at the end of
# wave 2 using the function `wave_points`
wave_points(list(initial_points, new_points, new_new_points), input_names = names(ranges), p_size = 1)

# Create a dataframe `wave2` binding the parameters sets generated at the end of wave 2 with the corresponding 
# model outputs
new_new_initial_results <- setNames(data.frame(t(apply(new_new_points, 1, 
                                    get_results, c(25, 40, 100, 200, 300, 350), 
                                    c('I', 'R')))), names(targets))
wave2 <- cbind(new_new_points, new_new_initial_results)

# Assess how much better parameter sets at later waves perform compared to the original `initial_points` through 
#`simulator_plot`
all_points <- list(wave0, wave1, wave2)
simulator_plot(all_points, targets)


### Solution to the task on `simulator_plot` function ###

# Set `normalize=TRUE` to rescale and uniformize the targets bounds in the plot
simulator_plot(all_points, targets, normalize = TRUE)

# Set `logscale=TRUE` to plot log-scaled target bounds
simulator_plot(all_points, targets, logscale = TRUE)

### End of the solution ###

# For each combination of two outputs, show the output values for non-implausible parameter sets at each wave.
wave_values(all_points, targets, l_wod = 1, p_size = 1)
