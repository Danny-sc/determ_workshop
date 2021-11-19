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
  b = c(1e-5, 1e-4),
  mu = c(1e-5, 1e-4),
  beta1 = c(0.2, 0.3),
  beta2 = c(0.1, 0.2),
  beta3 = c(0.3, 0.5),
  sigma = c(0.07, 0.21),
  alpha = c(0.01, 0.025),
  gamma = c(0.05, 0.08),
  omega = c(0.002, 0.004)
)

targets <- list(
  I25 = list(val = 115.88, sigma = 3.48),
  I40 = list(val = 137.84, sigma = 4.14),
  I100 = list(val = 26.34, sigma = 0.795),
  I200 = list(val = 0.68, sigma = 0.015),
  I300 = list(val = 29.55, sigma = 0.885),
  I350 = list(val = 68.89, sigma = 2.07),
  R25 = list(val = 125.12, sigma = 3.75),
  R40 = list(val = 256.80, sigma = 7.71),
  R100 = list(val = 538.99, sigma = 16.17),
  R200 = list(val = 444.23, sigma = 13.32),
  R300 = list(val = 371.08, sigma = 11.13),
  R350 = list(val = 549.42, sigma = 16.485)
)

test_grid <- expand.grid(beta1 = seq(ranges$beta1[1], ranges$beta1[2], length.out = 20), 
                         +                          gamma = seq(ranges$gamma[1], ranges$gamma[2], length.out = 20))
for (i in names(ranges)[!names(ranges) %in% c('beta1', 'gamma')]) test_grid[,i] <- rep(mean(ranges[[i]]), 400)

R200300_results <- setNames(data.frame(t(apply(test_grid, 1, get_results, c(200,300), c('R')))), names(targets)[10:11])

ggplot(data = test_grid, aes(x = beta1, y = gamma)) +geom_contour_filled(aes(z = R200300_results[,1]), bins=20)
