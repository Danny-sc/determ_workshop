library(hmer)
library(deSolve)
library(ggplot2)
library(reshape2)
library(purrr)
library(tidyverse)
set.seed(123)
ode_results <- function(parms, end_time = 365*2) {
  forcer = matrix(c(0, parms['beta1'], 100, parms['beta2'], 180, parms['beta2'], 270, parms['beta3']), ncol = 2, byrow = TRUE)
  force_func = approxfun(x = forcer[,1], y = forcer[,2], method = "linear", rule = 2)
  des = function(time, state, parms) {
    with(as.list(c(state, parms)), {
      dS <- b*(S+E+I+R)-force_func(time)*I*S/(S+E+I+R) + omega*R - mu*S
      dE <- force_func(time)*I*S/(S+E+I+R) - sigma*E - mu*E
      dI <- sigma*E - alpha*I - gamma*I - mu*I
      dR <- gamma*I - omega*R - mu*R
      return(list(c(dS, dE, dI, dR)))
    })
  }
  yini = c(S = 900, E = 100, I = 0, R = 0)
  times = seq(0, end_time, by = 1)
  out = deSolve::ode(yini, times, des, parms)
  return(out)
}
get_results <- function(params, times, outputs) {
  t_max <- max(times)
  all_res <- ode_results(params, t_max)
  actual_res <- all_res[all_res[,'time'] %in% times, c('time', outputs)]
  shaped <- reshape2::melt(actual_res[,outputs])
  return(setNames(shaped$value, paste0(shaped$Var2, actual_res[,'time'], sep = "")))
}

ranges = list(
  b = c(1e-5, 1e-4), # birth rate
  mu = c(1e-5, 1e-4), # rate of death from other causes
  beta1 = c(0.2, 0.3), # infection rate at time t=0
  beta2 = c(0.1, 0.2), # infection rates at time t=100
  beta3 = c(0.3, 0.5), # infection rates at time t=270
  sigma = c(0.07, 0.21), # rate of becoming infectious after infection
  alpha = c(0.01, 0.025), # rate of death from the disease
  gamma = c(0.05, 0.08), # recovery rate
  omega = c(0.002, 0.004) # rate at which immunity is lost following recovery
)

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

chosen_params <- list(b = 1/(76*365), mu = 1/(76*365), beta1 = 0.214, 
                      beta2 = 0.107, beta3 = 0.428, sigma = 1/7, alpha = 1/50, gamma = 1/14, omega = 1/365)

initial_LHS <- lhs::randomLHS(180, 9)

initial_points <- setNames(data.frame(t(apply(initial_LHS, 1, 
                                              function(x) x*unlist(lapply(ranges, function(x) x[2]-x[1])) + 
                                                unlist(lapply(ranges, function(x) x[1]))))), names(ranges))

initial_results <- setNames(data.frame(t(apply(initial_points, 1, get_results, 
                                               c(25, 40, 100, 200, 300, 350), c('I', 'R')))), names(targets))

wave0 <- cbind(initial_points, initial_results)

t_sample <- sample(1:nrow(wave0), 90)
training <- wave0[t_sample,]
validation <- wave0[-t_sample,]

ems_wave1 <- emulator_from_data(training, names(targets), ranges, 
                                c_lengths= rep(0.55,length(targets)))

op <- par(mfrow = c(1, 3))
cd <- comparison_diag(ems_wave1$R200, validation = validation, targets = targets)
ce <- classification_diag(ems_wave1$R200, validation = validation, targets =targets)
se <- standard_errors(ems_wave1$R200, validation = validation, targets = targets)
par(op)


sigmadoubled_emulator <- ems_wave1$R200$mult_sigma(1.1)
op <- par(mfrow = c(1, 3))
cd <- comparison_diag(sigmadoubled_emulator, validation = validation, targets = targets)
ce <- classification_diag(sigmadoubled_emulator, validation = validation, targets = targets)
se <- standard_errors(sigmadoubled_emulator, validation = validation, targets = targets)
par(op)

vd <- validation_diagnostics(ems_wave1, validation = validation, targets = targets)

inflations <- c(1.5,1.5,2,1.5,2,1.5,2,1.5,1,1.3,2,2)
for (i in 1:length(ems_wave1)) {
  ems_wave1[[i]] <- ems_wave1[[i]]$mult_sigma(inflations[[i]])
}

vd <- validation_diagnostics(ems_wave1, validation = validation, targets = targets)

new_points <- generate_new_runs(ems_wave1, 180, targets, verbose = TRUE)
plot_wrap(new_points, ranges)

min_val <- list()
max_val <- list()
new_ranges <- list()
for (i in 1:length(ranges)) {
  par <- names(ranges)[[i]]
  min_val[[par]] <- max(min(new_points[,par])-0.05*diff(range(new_points[,par])), 
                        ranges[[par]][1])
  max_val[[par]] <- min(max(new_points[,par])+0.05*diff(range(new_points[,par])),
                        ranges[[par]][2])
  new_ranges[[par]] <- c(min_val[[par]], max_val[[par]])
}

new_initial_results <- setNames(data.frame(t(apply(new_points, 1, get_results, 
                                                   c(25, 40, 100, 200, 300, 350), c('I', 'R')))), names(targets))

wave1 <- cbind(new_points, new_initial_results)


new_t_sample <- sample(1:nrow(wave1), 90)
new_training <- wave1[new_t_sample,]
new_validation <- wave1[-new_t_sample,]

ems_wave2 <- emulator_from_data(new_training, names(targets), new_ranges, c_lengths= rep(0.55,length(targets)))

vd <- validation_diagnostics(ems_wave2, validation = new_validation, targets = targets)

inflations <- c(4,6,4,4,4,2,4,4,4,4,4,4)
for (i in 1:length(ems_wave2)) {
  ems_wave2[[i]] <- ems_wave2[[i]]$mult_sigma(inflations[[i]])
}
vd <- validation_diagnostics(ems_wave2, validation =  new_validation, targets = targets)

new_new_points <- generate_new_runs(c(ems_wave2, ems_wave1), 180, targets, verbose=TRUE)
plot_wrap(new_new_points, ranges)

R_squared_new <- list()
for (i in 1:length(ems_wave2)) {
  R_squared_new[[i]] <- summary(ems_wave2[[i]]$model)$adj.r.squared
}
names(R_squared_new) <- names(ems_wave2)
unlist(R_squared_new)

ems_wave1_linear <- emulator_from_data(training, names(targets), 
                                       ranges, quadratic=FALSE, c_lengths= rep(0.55,length(targets)))
R_squared_linear <- list()
for (i in 1:length(ems_wave1_linear)) {
  R_squared_linear[[i]] <- summary(ems_wave1_linear[[i]]$model)$adj.r.squared
}
names(R_squared_linear) <- names(ems_wave1_linear)
unlist(R_squared_linear)

emulator_plot(ems_wave1_linear$I200, plot_type = 'var', 
              params = c('beta1', 'gamma'))

emulator_plot(ems_wave1_linear$I200, plot_type = 'imp', targets = targets, 
              params = c('beta1', 'gamma'), cb=TRUE)

ems_wave1_linear$I200 <- ems_wave1_linear$I20$set_hyperparams(
  list(theta=ems_wave1_linear$I200$corr$hyper_p$theta *3))
emulator_plot(ems_wave1_linear$I200, plot_type = 'var', 
              params = c('beta1', 'gamma'))

emulator_plot(ems_wave1_linear$I200, plot_type = 'imp', targets = targets, 
              params = c('beta1', 'gamma'), cb=TRUE)

wave_points(list(initial_points, new_points, new_new_points), input_names = names(ranges))

new_new_initial_results <- setNames(data.frame(t(apply(new_new_points, 1, 
                                                       get_results, c(25, 40, 100, 200, 300, 350), 
                                                       c('I', 'R')))), names(targets))
wave2 <- cbind(new_new_points, new_new_initial_results)

all_points <- list(wave0, wave1, wave2)
simulator_plot(all_points, targets)

simulator_plot(all_points, targets, normalize = TRUE)
simulator_plot(all_points, targets, logscale = TRUE)

