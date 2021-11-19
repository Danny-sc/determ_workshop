--- 
title: "Workshop 1: calibrating a deterministic model"
author: "Danny Scarponi, Andy Iskauskas"
site: bookdown::bookdown_site
output:
    bookdown::pdf_book:
        includes:
            in_header: header.tex
    bookdown::gitbook:
        config:
            sharing: null
        css: 'style.css'
        highlight: tango
        includes:
            in_header: _toggle.html
        keep_md: TRUE
linkcolor: blue
documentclass: book
link-citations: yes
description: "An interactive introduction to the HoMER package"
---







# Introduction to the model {#intro}

This tutorial offers an interactive introduction to the main functionality of the 
[HoMER](https://github.com/Tandethsquire/hmer) package, using a synthetic example of an epidemiological model. In this workshop, you will be invited to carry out a series of tasks (see "Task" boxes) which will enhance your understanding of the package and its tools. Thanks to these activities, you will learn to calibrate deterministic models using history matching and model emulation,  and to use your judgement to customise the process. Although self-contained, this tutorial should be considered as a natural continuation of Tutorial 1 (hyperlink), which gives a more general overview of the history matching with emulation process for deterministic models, and shows how to perform it using [HoMER](https://github.com/Tandethsquire/hmer). Following Workshop 1, you may also want to read Tutorial 2 and Workshop 2, where we demonstrate how to calibrate stochastic models. 

The model that we chose for demonstration purposes is a deterministic  <span class="abbr" title="A model consisting of four compartments 

- $S$: Susceptible individuals,
- $E$: Exposed individuals (i.e. people that are infected but not infectious yet), 
- $I$: Infectious individuals,  
- $R$: Recovered individuals, 

and four possible transitions

- $S \rightarrow E$, when a susceptible individual becomes infected, 
- $E \rightarrow I$, when an infected individual becomes infectious,
- $I \rightarrow R$, when an infectious individual recovers,
- $R \rightarrow S$, when a recovered individual becomes susceptible again.

SEIRS models are used to study those infectious diseases that do not confer permanent immunity."><abbr title="A model consisting of four compartments 

- S: Susceptible individuals,
- E: Exposed individuals (i.e. people that are infected but not infectious yet), 
- I: Infectious individuals,  
- R: Recovered individuals, 

and four possible transitions

- S to E, when a susceptible individual becomes infected, 
- E to I, when an infected individual becomes infectious,
- I to R, when an infectious individual recovers,
- R to S, when a recovered individual becomes susceptible again.

SEIRS models are suitable to study those infectious diseases that have an incubation period and do not confer permanent immunity.">
SEIRS</abbr></span>
model, described by the following differential equations:
\begin{align}
\frac{dS}{dt} &= b N - \frac{\beta(t)IS}{N} + \omega R -\mu S  \\ 
\frac{dE}{dt} &= \frac{\beta(t)IS}{N} - \sigma E - \mu E \\ 
\frac{dI}{dt} &= \sigma E - \gamma I - (\mu + \alpha) I \\ 
\frac{dR}{dt} &= \gamma I - \omega R - \mu R
\end{align}
where $N$ is the total population, varying over time, and the parameters are as follows:

- $b$ is the birth rate, 

- $\mu$ is the  rate of death from other causes, 

- $\beta(t)$ is the infection rate between each infectious and susceptible individual, 

- $\sigma$ is the rate of becoming infectious after infection, 

- $\alpha$ is the rate of death from the disease, 

- $\gamma$ is the recovery rate and  

- $\omega$ is the rate at which immunity is lost following recovery. 

<div class="figure" style="text-align: center">
<img src="SEIRSdiagram.png" alt="SEIRS Diagram"  />
<p class="caption">(\#fig:unnamed-chunk-2)SEIRS Diagram</p>
</div>

The rate of infection $\beta(t)$ is set to be a simple linear function interpolating between points, where the points in question are $\beta(0)=\beta_1$, $\beta(100) = \beta(180) = \beta_2$, $\beta(270) = \beta_3$ and where $\beta_2 < \beta_1 < \beta_3$. This choice was made to represent an infection rate that initially drops due to external (social) measures and then raises when a more infectious variant appears. Here $t$ is taken to measure days. Below we show a graph of the infection rate over time when $\beta_1=0.3, \beta_2=0.1$ and $\beta_3=0.4$:

<div class="figure" style="text-align: center">
<img src="infection_rate.png" alt="Infection rate graph" width="60%" />
<p class="caption">(\#fig:unnamed-chunk-3)Infection rate graph</p>
</div>

In order to obtain the solution of the differential equations for a given set of parameters, we will use a helper function, `ode_results`. The function assumes an initial population of 900 susceptible individuals, 100 exposed individuals, and no infectious or recovered individuals. Below we use `ode_results` with an example set of parameters and plot the model output over time.


```r
example_params <- c(
  b = 1/(60*365),
  mu = 1/(76*365),
  beta1 = 0.2, beta2 = 0.1, beta3 = 0.3,
  sigma = 0.13,
  alpha = 0.01,
  gamma = 0.08,
  omega = 0.003
)
solution <- ode_results(example_params)
plot(solution)
```

<img src="_main_files/figure-html/unnamed-chunk-4-1.png" style="display: block; margin: auto;" />

<div class="panel panel-default"><div class="panel-heading"> Task </div><div class="panel-body"> 
Familiarise yourself with the model. Investigate how the plots change as you change the values of the parameters.

<span class="abbr" title=""><abbr title="Copy the code below, modify the value of (some) parameters and run it.

```r
example_params <- c(
  b = 1/(60*365),
  mu = 1/(76*365),
  beta1 = 0.2, beta2 = 0.1, beta3 = 0.3,
  sigma = 0.13,
  alpha = 0.01,
  gamma = 0.08,
  omega = 0.003
)
solution <- ode_results(example_params)
plot(solution)
```
">
R tip</abbr></span> </div></div>



<button id="displayTextunnamed-chunk-6" onclick="javascript:toggle('unnamed-chunk-6');">Show: Solution</button>

<div id="toggleTextunnamed-chunk-6" style="display: none"><div class="panel panel-default"><div class="panel-heading panel-heading1"> Solution </div><div class="panel-body">
Let us see what happens when a higher force of infection is considered:

```r
higher_foi_params <- c(
  b = 1/(60*365),
  mu = 1/(76*365),
  beta1 = 0.3, beta2 = 0.1, beta3 = 0.5,
  sigma = 0.13,
  alpha = 0.01,
  gamma = 0.08,
  omega = 0.003
)
higher_foi_solution <- ode_results(higher_foi_params)
plot(solution, higher_foi_solution)
```

<img src="_main_files/figure-html/unnamed-chunk-63-1.png" style="display: block; margin: auto;" />

Here the black line shows the model output when it is run using the original parameters and the red dotted line when it is run using a higher force of infection. As expected, the number of susceptible individuals decreases, while the size of outbreaks increases.

Let us now also increase the rate of becoming infectious following infection $\sigma$:

```r
higher_sigma_params <- c(
  b = 1/(60*365),
  mu = 1/(76*365),
  beta1 = 0.3, beta2 = 0.1, beta3 = 0.5,
  sigma = 0.21,
  alpha = 0.01,
  gamma = 0.08,
  omega = 0.003
)
higher_sigma_solution <- ode_results(higher_sigma_params)
plot(higher_foi_solution,higher_sigma_solution)
```

<img src="_main_files/figure-html/unnamed-chunk-64-1.png" style="display: block; margin: auto;" />

Here the black line is the model output when the model is run with the previous parameter set, and the red dotted line is the model output when we also increase sigma. We observe a decrease in the number of exposed individuals. Again, this is in agreement with our expectation: a higher rate of becoming infectious means that people leave the exposed compartmental to enter the infectious compartment faster than before.


Finally, what happens when a lower value of the recovery rate $\gamma$ is used?

```r
smaller_recovery_params <- c(
  b = 1/(60*365),
  mu = 1/(76*365),
  beta1 = 0.3, beta2 = 0.1, beta3 = 0.5,
  sigma = 0.21,
  alpha = 0.01,
  gamma = 0.05,
  omega = 0.003
)
smaller_recovery_solution <- ode_results(smaller_recovery_params)
plot(higher_sigma_solution, smaller_recovery_solution)
```

<img src="_main_files/figure-html/unnamed-chunk-65-1.png" style="display: block; margin: auto;" />

Here the black line is the model output when the model is run with the previous parameter set, and the red dotted line is the model output when we also decrease the recovery rate. Again, as one expects, this causes the number of susceptible individuals to decrease and the number of infectious individuals to increase, at least during the first peak.</div></div></div>

# Defining 'wave0' data

Before we tackle the emulation, we need to define some objects. First of all, let us set the parameter ranges:


```r
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
```

We then turn to the targets we will match: the number of infectious individuals $I$ and the number of recovered individuals $R$ at times $t=25, 40, 100, 200, 200, 350$. For each of these outputs, we define a pair (val, sigma), where ‘val’ represents the mean value of the output and ‘sigma’ represents its standard deviation:


```r
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
```
  
The 'sigmas' in our `targets` list represent the uncertainty we have about the observations. Note that in general we can also choose to include model uncertainty in the 'sigmas', to reflect how accurate we think our model is. 
  
<infobutton id="displayTextunnamed-chunk-9" onclick="javascript:toggle('unnamed-chunk-9');">Show: More on how targets were set</infobutton>

<div id="toggleTextunnamed-chunk-9" style="display: none"><div class="panel panel-default"><div class="panel-body">

Since our model is synthetic, we couldn't rely on observed data to define our targets. Instead, we chose the parameter set 

```r
chosen_params <- list(b = 1/(76*365), mu = 1/(76*365), beta1 = 0.214, 
                      beta2 = 0.107, beta3 = 0.428, sigma = 1/7, alpha = 1/50, gamma = 1/14, omega = 1/365)
```

ran the model with it and used the relevant model outputs as the ‘val’ in `targets`. The ‘sigma’ components were chosen to be $5\%$ of the corresponding ‘val’.</div></div></div>

Finally we need a set of `wave0` data to start. When using your own model, you can create the dataframe 'wave0' following the steps below:
  
- build a named dataframe `initial_points` containing a space filling set of points, which can be generated using a Latin Hypercube Design or another sampling method of your choice. A rule of thumb is to select at least $10p$ parameter sets, where $p$ is the number of parameters in the model. The columns of `initial_points` should be named exactly as in the list `ranges`;

- run your model on the parameter sets in `initial_points` and define a dataframe `initial_results` in the following way: the nth row of `initial_results` should contain the model outputs corresponding to the chosen targets for the parameter set in the nth row of `initial_points`. The columns of `initial_results` should be named exactly as in the list `targets`;

- bind `initial_points` and `initial_results` by their columns to form the dataframe `wave0`.

For this workshop, we generate parameter sets using a [Latin Hypercube](https://en.wikipedia.org/wiki/Latin_hypercube_sampling) design (see fig. \@ref(fig:figlhs) for a Latin hypercube example in two dimensions). 

<div class="figure" style="text-align: center">
<img src="LHS.png" alt="An example of Latin hypercube with 10 points in two dimensions: there is only one sample point in each row and each column." width="60%" />
<p class="caption">(\#fig:figlhs)An example of Latin hypercube with 10 points in two dimensions: there is only one sample point in each row and each column.</p>
</div>

Through the function `randomLHS` in the package `lhs` we create a hypercube design with 180 parameter sets for our nine inputs, 20 times the dimensionality of the parameter space. To access the external function `randomLHS`, we use the syntax `package::function()`:


```r
initial_LHS <- lhs::randomLHS(180, 9)
```

Note that in `initial_LHS` each parameter is distributed on $[0,1]$. This is not exactly what we need, since each parameter has a different range. We therefore re-scale each component in `initial_LHS` multiplying it by the difference between the upper and lower bounds of the range of the  corresponding parameter and then we add the lower bound for that parameter. In this way we obtain `initial_points`, which contains parameter values in the correct ranges. 


```r
initial_points <- setNames(data.frame(t(apply(initial_LHS, 1, 
                  function(x) x*unlist(lapply(ranges, function(x) x[2]-x[1])) + 
                  unlist(lapply(ranges, function(x) x[1]))))), names(ranges))
```

We then run the model for the parameter sets in `initial_points` through the `get_results` function. This is a helper function that acts as `ode_results`, but has the additional feature of allowing us to decide which outputs and times should be returned: in our case we need the values of $I$ and $R$ at $t=25, 40, 100, 200, 300, 350$.


```r
initial_results <- setNames(data.frame(t(apply(initial_points, 1, get_results, 
                   c(25, 40, 100, 200, 300, 350), c('I', 'R')))), names(targets))
```

Finally, we bind the parameter sets `initial_points` to the model runs `initial_results` and save everything in the data.frame `wave0`: 


```r
wave0 <- cbind(initial_points, initial_results)
```

Note that when using your own model, you can create the 'wave0' dataframe however you wish. The info box below will help you build `wave0` when working on your own model. 


# Emulators {#constr}

Let us start by splitting `wave0` in two parts: the training set, on which we will train the emulators, and a validation set, which will be used to do diagnostics of the emulators. 


```r
t_sample <- sample(1:nrow(wave0), 90)
training <- wave0[t_sample,]
validation <- wave0[-t_sample,]
```


## A brief recap on emulators
Before building emulators, let us quickly remind ourselves what an emulator is and how it is structured. Note that a more detailed discussion about the structure of an emulator can be found in Tutorial 2 (Section 3, and Appendix A and B). 

An [emulator](https://en.wikipedia.org/wiki/Emulator) is a way of representing our 
<span class="abbr" title=""><abbr title="In Bayesian statistics, probability expresses a degree of belief in an event. Such belief can be based either on prior knowledge or on personal beliefs about the event. 
">beliefs</abbr></span>
about the behaviour of a function whose output is unknown. In this workshop what is unknown is the behaviour of our SEIRS model and we will emulate each of the model outputs separately. Given a training dataset, i.e. a set of model runs, we can train an emulator and use it to get expectation and variance for a model output at any parameter set, without the need to run the model at the chosen set. We think of the expectation as the prediction provided by the emulator at the chosen parameter set, and we think of the variance as the uncertainty associated to that prediction. 

The general structure of a univariate emulator is as follows:
$$f(x) = g(x)^T \xi + u(x),$$
where $g(x)^T \xi$ is a regression term and $u(x)$ is a Gaussian process with mean zero. The role of the regression term is to mimic the global behaviour of the model output, while $u(x)$ represents localised deviations of the output from this global behaviour near to $x$.  

The regression term is specified by:

- a vector of functions of the parameters $g(x)$ which determine the shape and complexity of the regression 
<span class="abbr" title=""><abbr title="A hypersurface is a mathematical object that generalizes the concept of surface from the three-dimensional space to hyperspaces, i.e. spaces of dimension higher than three.
">hypersurface</abbr></span>
 we fit to the training data. For example, if $x$ is one dimensional, i.e. we have just one parameter, setting $g(x)=(1,x)$  corresponds to fitting a straight line to the training data. Similarly, setting $g(x)=(1,x,x^2)$ corresponds to fitting a parabola to the training data. Fig \@ref(fig:regresresid) shows a one-dimensional example with a quadratic global behaviour: the model output to emulate is in black, while the best fitting parabola is in red.

- a vector of regression coefficients $\xi$. In the one dimensional case for example, if we set $g(x)=(1,x)$, then $\xi=(\xi_0,\xi_1)$, where $\xi_0$ is the $y$-intercept and $\xi_1$ is the gradient of the straight line fitted to the training data. 

<div class="figure" style="text-align: center">
<img src="regres_resid_plot.png" alt="Regression term and residuals in one dimensional example" width="100%" />
<p class="caption">(\#fig:regresresid)Regression term and residuals in one dimensional example</p>
</div>

In general, and especially when dealing with complex models, we cannot expect the regression hypersurface to perfectly explain the behaviour of the output. For this reason it is necessary to account for the local deviations of the output from the regression hypersurface. These local deviations, also referred to as residuals, are shown in blue in Fig \@ref(fig:regresresid). When the parameter space is one-dimensional, they indicate how far the regression term is from the model output at each point. Since residuals are unknown, we treat them as random variables: for each parameter $x$, we then have a random variable $u(x)$ representing the residual at $x$. In the [HoMER](https://github.com/Tandethsquire/hmer) package we assume this collection of random variables $u(x)$ to be a [Gaussian process](https://en.wikipedia.org/wiki/Gaussian_process), with mean zero. Informally this means the following:

- for each parameter set $x$, we consider the residual $u(x)$ as a normally distributed random variable with mean zero. Note that the mean is assumed to be zero, since, even if we expect to see local deviations, we do not expect the output to be systematically above (or below) the regression hypersurface;

- given any pair of parameter sets $(x,x')$, the pair $(u(x),u(x'))$ is a multivariate normal variable, with mean $(0,0)$.

To fully describe the Gaussian process $u(x)$ we need to define the covariance structure, i.e. we need to say how correlated the residuals at $x$ and $x'$ are, for any pair $(x,x')$. A common covariance structure used in Gaussian processes is given by 

$$\text{Cov}(u(x), u(x'))= \sigma^2  c(x,x^{\prime}) $$
where  $c$ is the square-exponential correlation function

$$c(x,x^{\prime}) :=  \exp\left(\frac{-\sum\limits_{i}(x_{i}-x_{i}^{\prime})^2}{\theta^2}\right)$$

where $x_i$ is the ith-component of the parameter set $x.$ This covariance structure is the default option in the [HoMER](https://github.com/Tandethsquire/hmer) package, even though other structures are available.

Let us look at the various terms in this covariance structure:

- $\sigma^2$ is the **emulator variance**, i.e the variance of $u(x)$, for all parameter sets $x$. The value of $\sigma$ reflects how far from the regression hypersurface the model output can be. The larger the value of $\sigma$, the farthest the model output can be from the regression hypersurface. In particular, larger values of $\sigma$ correspond to more uncertain emulators. For example, Fig. \@ref(fig:regresresid) was generated with a $\sigma$ of $0.3$. A higher $\sigma$ of $1$, would create wider residuals, as in the plot below:

<div class="figure" style="text-align: center">
<img src="regres_resid_plot_highsigma.png" alt="Regression term and residuals in one dimensional example, with higher $\sigma$" width="100%" />
<p class="caption">(\#fig:unnamed-chunk-15)Regression term and residuals in one dimensional example, with higher $\sigma$</p>
</div>

- $\theta$ is the **correlation length** of the Gaussian process. For a given pair $(x,x')$, the larger $\theta$ is, the larger is the covariance between $u(x)$ and $u(x')$. This means that the size of $\theta$ determines how close two parameter sets must be in order for the corresponding residuals to be non-negligibly correlated.  Informally, we can think of $\theta$ in the following way: if the distance of two parameters sets is no more than $\theta$, then their residuals will be well correlated. In particular, a larger $\theta$ results in a smoother (less wiggly) emulator. In the one dimensional example in Fig. \@ref(fig:regresresid), a $\theta$ of $1$ was used. A value of $\theta$ equal to $0.4$ would result in less smooth residuals:

<div class="figure" style="text-align: center">
<img src="regres_resid_plot_lowtheta.png" alt="Regression term and residuals in one dimensional example, with lower $\theta$" width="100%" />
<p class="caption">(\#fig:unnamed-chunk-16)Regression term and residuals in one dimensional example, with lower $\theta$</p>
</div>

Choosing values for $\sigma$ and $\theta$ corresponds to making a judgment about how far we expect the output to be from the regression hypersurface ($\sigma$) and about its smoothness ($\theta$). While the [HoMER](https://github.com/Tandethsquire/hmer) package, and in particular the function `emulator_from_data`, selects values of $\sigma$ and $\theta$ for us based on the provided training data, we will see in this workshop how we can intervene to customise the choice of these hyperparameters and the benefits that this operation brings.

## Training emulators

We are now ready to train the emulators using the `emulator_from_data` function, which just needs the training set, the names of the targets we want to emulate and the ranges of the parameters:


```r
ems_wave1 <- emulator_from_data(training, names(targets), ranges)
```

In `ems_wave1` we have information about all emulators. Let us take a look at the emulator of the number of recovered individuals at time $t=200$:


```r
ems_wave1$R200
```

```
## Parameters and ranges:  b: c(0, 1e-04): mu: c(0, 1e-04): beta1: c(0.2, 0.3): beta2: c(0.1, 0.2): beta3: c(0.3, 0.5): sigma: c(0.07, 0.21): alpha: c(0.01, 0.025): gamma: c(0.05, 0.08): omega: c(0.002, 0.004) 
## Specifications: 
## 	 Basis functions:  (Intercept); b; mu; beta1; beta2; sigma; alpha; gamma; omega; I(beta1^2); I(beta2^2); I(sigma^2); I(alpha^2); I(gamma^2); I(omega^2); mu:sigma; beta1:beta2; beta1:alpha; beta1:gamma; beta1:omega; beta2:sigma; beta2:gamma; sigma:alpha; sigma:gamma; sigma:omega; alpha:gamma; alpha:omega; gamma:omega 
## 	 Active variables b; mu; beta1; beta2; sigma; alpha; gamma; omega 
## 	 Regression Surface Expectation:  496.2728; 1.9876; -4.6397; 8.2317; 25.0217; -15.1378; -57.9129; -4.992; -57.2594; -2.5115; -12.4326; -4.9669; 1.0659; -8.4424; 6.5641; 2.4632; 1.5658; 5.5449; 11.0464; 9.2804; 1.8736; 9.3984; -7.801; -5.4154; -6.0825; 4.3294; -1.9138; 5.5718 
## 	 Regression surface Variance (eigenvalues):  0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0 
## Correlation Structure: 
## Bayes-adjusted emulator - prior specifications listed. 
## 	 Variance (Representative):  2.040089 
## 	 Expectation:  0 
## 	 Correlation type: exp_sq 
## 	 Hyperparameters:  theta: 0.5009 
## 	 Nugget term: 0 
## Mixed covariance:  0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```

The print statement provides an overview of the emulator specifications, which refer to the global part, and correlation structure, which refers to the local part:

- Active variables: these are the variables that have the most explanatory power for the chosen output. In our case  all variables but $\beta_3$ are active.

- Basis Functions: these are the functions composing the vector $g(x)$. Note that, since by default `emulator_from_data` uses quadratic regression for the global part of the emulator, the list of basis functions contains not only the active variables but also products of them.  
 
- First and second order specifications for $\xi$ and $u(x)$. Note that  by default `emulator_from_data` assumes that the regression surface is known and its coefficients are fixed. This explains why Regression Surface Variance and Mixed Covariance (which shows the covariance of $\xi$ and $u(x)$) are both zero. The term Variance refers to $\sigma^2$ in $u(x)$. 

We can also plot the emulators to see how they represent the output space: the `emulator_plot` function does this for emulator expectation (default option), variance, standard deviation, and implausibility.
The emulator expectation plots show the structure of the regression surface, which is at most quadratic in its parameters, through a 2D slice of the input space. 


```r
emulator_plot(ems_wave1$R200, params = c('beta1', 'gamma'))
```

<img src="_main_files/figure-html/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />

Here for each pair $(\bar \beta_1,\bar \gamma)$ the plot shows the expected value produced by the emulator `ems_wave1$R200` at the parameter set having $\beta_1=\bar \beta_1$, $\gamma=\bar \gamma$ and all other parameters equal to their mid-range value (the ranges of parameters are those that were passed to `emulator_from_data` to train `ems_wave1`). Note that we chose to display $\beta_1$ and $\gamma$, but any other pair can be selected. For consistency, we will use $\beta_1$ and $\gamma$ throughout this workshop.

<div class="panel panel-default"><div class="panel-heading"> Task </div><div class="panel-body"> 
The argument `fixed_vals` allows us to change where to slice the parameters that are not shown in the plots. Try changing this argument: what happens if we set all hidden parameters equal to the values in `chosen_parameters`? What do you expect the emulator expectation to be in this case? 
  
<span class="abbr" title=""><abbr title="To make the unshown parameters equal to the values in `chosen_params`, set `fixed_vals` to `chosen_params[!names(chosen_params) %in% c('beta1', 'gamma')]`">
R tip</abbr></span>  </div></div>

<button id="displayTextunnamed-chunk-21" onclick="javascript:toggle('unnamed-chunk-21');">Show: Solution</button>

<div id="toggleTextunnamed-chunk-21" style="display: none"><div class="panel panel-default"><div class="panel-heading panel-heading1"> Solution </div><div class="panel-body">
When slicing at the values of `chosen_parameters` we would expect the emulator expectation to be somewhere close to the target ($444.23$) when $\beta_1$ and $\gamma$ are equal to their values in `chosen_params`, i.e. $0.214$ and $1/14$. Let us plot the emulator expectation with the `emulator_plot` function and add the point $(\beta_1=0.214,\gamma=1/14)$ with `geom_point`:
  

```r
emulator_plot(ems_wave1$R200, params = c('beta1', 'gamma'), 
              fixed_vals = chosen_params[!names(chosen_params) %in% c('beta1', 'gamma')])+ 
  geom_point(aes(x=0.214, y=1/14), size=3)
```

<img src="_main_files/figure-html/unnamed-chunk-67-1.png" style="display: block; margin: auto;" />

The plot is in agreement with our reasoning: when $\beta_1$ and $\gamma$ are equal to $0.214$ and $1/14$, the emulator expectation is around in the interval $[443,446]$ (see black point in the box). It is worth noting that when using a real, more complex model, we would not necessarily expect the emulated output at wave zero to be as close to the model output as it is for the synthetic model we chose for this workshop.</div></div></div>



<div class="panel panel-default"><div class="panel-heading"> Task </div><div class="panel-body"> 
Is $\beta_3$ active for all emulators? Why? 
  
<span class="abbr" title=""><abbr title="To show what variables are active for an emulator 'em' you can access the parameter `active_vars` of the emulator, typing `em$active_vars`.">
R tip</abbr></span>  </div></div>


<button id="displayTextunnamed-chunk-23" onclick="javascript:toggle('unnamed-chunk-23');">Show: Solution</button>

<div id="toggleTextunnamed-chunk-23" style="display: none"><div class="panel panel-default"><div class="panel-heading panel-heading1"> Solution </div><div class="panel-body">
Let us take a look at the emulators. To show what variables are active for an emulator 'em' we can simply type `em$active_vars`. For example:
  

```r
ems_wave1$R200$active_vars
```

```
## [1]  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE
```

shows what variables are active for the $R200$ emulator. Since $\beta_3$ is the fifth parameter (see for example `ranges`), then we can type

```r
ems_wave1$R200$active_vars[5]
```

```
## [1] FALSE
```

to just focus on $\beta_3$. To look at the role of $\beta_3$ in all emulators at once, we create a logical vector looping through `ems_wave1`:

```r
beta3_role <- logical()
for (i in 1:length(ems_wave1)) beta3_role[i] <- ems_wave1[[i]]$active_vars[5]
names(beta3_role) <- names(ems_wave1)
beta3_role
```

```
##   I25   I40  I100  I200  I300  I350   R25   R40  R100  R200  R300  R350 
## FALSE FALSE FALSE FALSE  TRUE  TRUE FALSE  TRUE FALSE FALSE  TRUE  TRUE
```

Here we see that $\beta_3$ tends to be more active at later times. This is in fact what we would expect: the later infection/recovery rates  don't have an impact on early time outputs.</div></div></div>

Looking at what variables are active for different emulators is often an instructive exercise. The code below produces a plot that shows all dependencies at once.


```r
plot_actives(ems_wave1)
```

<img src="_main_files/figure-html/unnamed-chunk-24-1.png" style="display: block; margin: auto;" />

From this table, we can immediately see that $\mu$ is inactive for most outputs, while $\beta_1$, $\beta_2$, $\sigma$, $\alpha$, $\gamma$ are active for most outputs. We also notice again that $\beta_3$ tends to be active for outputs at later times and inactive for outputs at earlier times. 

As mentioned above, `emulator_plot` can also plot the variance of a given emulator: 


```r
emulator_plot(ems_wave1$R200, plot_type = 'var', params = c('beta1', 'gamma'))
```

<img src="_main_files/figure-html/unnamed-chunk-25-1.png" style="display: block; margin: auto;" />

This plot shows the presence of a training point (purple-blue area on the right) close to the chosen slice of the input space. As discussed above, by default `emulator_plot` fixes all non-shown parameters to their mid-range, but different slices can be explored, through the argument `fixed_vals`. The purple-blue area indicates that the variance is low when we are close to the training point, which is in accordance with our expectation. 

Now that we have taken a look at the emulator expectation and the emulator variance, we might want to compare the relative contributions of the global and the residual terms to the overall emulator expectation. This can be done simply by examining the adjusted $R^2$ of the regression hypersurface:


```r
summary(ems_wave1$R200$model)$adj.r.squared
```

```
## [1] 0.9988402
```

We see that we have a very high value of the adjusted $R^2$. This means that the regression term explains most of the behaviour of the output $R200$. In particular, the residuals contribute little to the emulator predictions. This is not surprising, considering that we are working with a relatively simple SEIRS model. When dealing with more complex models, the regression term may be able to explain the model output less well. In such cases the residuals play a more important role.


# Implausibility

In this section we focus on implausibility and its role in the history matching process. Once emulators are built, we want to use them to systematically explore the input space. For any chosen parameter set, the emulator provides us with an approximation of the corresponding model output. This value is what we need to assess the implausibility of the parameter set in question.

For a given model output and a given target, the implausibility measures the difference between the emulator output and the target, taking into account all sources of uncertainty. For a parameter set $x$, the general form for the implausibility $I(x)$ is

$$I(x) = \frac{|f(x)-z|}{\sqrt{V_0 + V_c(x)+V_s+V_m}},$$

where $f(x)$ is the emulator output, $z$ the target, and the terms in the denominator refer to various forms of uncertainty. In particular

- $V_0$ is the variance associated with the observation uncertainty (i.e. uncertainty in estimates from observed data);
- $V_c(x)$ refers to the uncertainty one introduces when using the emulator output instead of the model output itself. Note that this term depends on $x$, since the emulator is more/less certain about its predictions based on how close/far $x$ is from training parameter sets;
- $V_s$ is the ensemble variability and represents the stochastic nature of the model (this term is not present in this tutorial, since the model is deterministic);
- $V_m$ is the model discrepancy, accounting for possible mismatches between the model and reality.

Since in this case study we want to emulate our model, without reference to a real-life analogue, the model represents the reality perfectly. For this reason we have $V_m=0$. Similarly we have $V_s=0$, since our model is deterministic. The observation uncertainty $V_0$ is represented by the 'sigma' values in the `targets` list, while $V_c$ is the emulator variance, which we discussed in the previous section.

A very large value of $I(x)$ means that we can be confident that the parameter set $x$ does not provide a good match to the observed data, even factoring in the additional uncertainty that comes with the use of emulators. When $I(x)$ is small, it could mean that the emulator output is very close 
to the model output or it could mean that the uncertainty in the denominator of $I(x)$ is large. In the former case, the emulator retains the parameter set, since it is likely to give a good fit to the observation for that output. In the latter case, the emulator does not have enough information to rule the parameter set out and therefore keeps it to explore it further in the next wave.

An important aspect to consider is the choice of cut-off for the implausibility measure. A rule of thumb follows [Pukelsheim's $3\sigma$ rule](https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule), a very general result which states that for any continuous unimodal distribution $95\%$ of the probability lies within $3$ sigma of the mean, regardless of asymmetry (or skewness etc). Following this rule, we set the implausibility threshold to be $3$: this means that a parameter $x$ is classified as non-implausible only if its implausibility is below $3$. 

For a given emulator, we can plot the implausibility through the function `emulator_plot` by setting `plot_type='imp'`. Note that we also set `cb=TRUE` to ensure that the produced plots are colour blind friendly:


```r
emulator_plot(ems_wave1$R200, plot_type = 'imp', 
              targets = targets, params = c('beta1', 'gamma'), cb=TRUE)
```

<img src="_main_files/figure-html/unnamed-chunk-27-1.png" style="display: block; margin: auto;" />

This is a 2D slice through the input space: for a chosen pair $(\bar\beta_1,\bar\gamma)$, the plot shows the implausibility of the parameter set having $\beta_1=\bar \beta_1$, $\gamma=\bar \gamma$ and all other parameters set to their mid-range value. Parameter sets with implausibility more than $3$ are highly unlikely to give a good fit and will be discarded when forming the parameters sets for the next wave.

Given multiple emulators, we can visualise the implausibility of several emulators at once:


```r
emulator_plot(ems_wave1, plot_type = 'imp', 
              targets = targets, params = c('beta1', 'gamma'), cb=TRUE)
```

<img src="_main_files/figure-html/unnamed-chunk-28-1.png" style="display: block; margin: auto;" />

This plot is useful to get an overall idea of which emulators have higher/lower implausibility, but how do we measure overall implausibility? We want a single measure for the implausibility at a given parameter set, but for each emulator we obtain an individual value for $I$. The simplest way to combine them is to consider maximum implausibility at each parameter set:
$$I_M(x) = \max_{i=1,\dots,N}I_{i}(x),$$ where $I_i(x)$ is the implausibility at $x$ coming from the $i$th emulator.
Note that Pukelsheim's rule applies for each emulator separately, but when we combine several emulators' implausibilities together a threshold of $3$ might be overly restrictive. For this reason, for large collections of emulators, it may be useful to replace the maximum implausibility with the second- or third-maximum implausibility. This also provides robustness to the failure of one or two of the emulators. 

<div class="panel panel-default"><div class="panel-heading"> Task </div><div class="panel-body"> 
Explore the functionalities of `emulator_plot` and produce a variety of implausibility plots. Here are a few suggestions: set `plot_type` to 'imp' to get implausibility plots or to 'nimp' to display the maximum implausibility plot; use the argument `nth` to obtain the second- or third- maximum implausibility plot; select a subset of all targets to pass to `emulator_plot`; change the value of the argument `fixed_vals` to decide where to slice the parameters that are not shown in the plots. </div></div>

<button id="displayTextunnamed-chunk-30" onclick="javascript:toggle('unnamed-chunk-30');">Show: Solution</button>

<div id="toggleTextunnamed-chunk-30" style="display: none"><div class="panel panel-default"><div class="panel-heading panel-heading1"> Solution </div><div class="panel-body">
Let us start by visualising the maximum implausibility passing all emulators to `emulator_plot` and setting `plot_type='nimp'`:

```r
emulator_plot(ems_wave1, plot_type = 'nimp', targets = targets, params = c('beta1', 'gamma'), cb=T)
```

<img src="_main_files/figure-html/unnamed-chunk-71-1.png" style="display: block; margin: auto;" />

This plot shows very high values of the implausibility for most points in the box. During the first few waves of history matching, one can consider second-maximum implausibility, rather than maximum implausibility. This means that instead of requiring the implausibility measure to be under the chosen threshold for all outputs, we allow (at most) one of them to be over it. This approach, which may result in less space cut out during the first few waves, has the advantage of being more conservative, reducing the risk that parts of the input space may be incorrectly cut. The more strict maximum implausibility measure can then be adopted in later waves, when the space to search is considerably smaller than the original input space, and the emulators will be less uncertain. To work with second-maximum implausibility we simply add nth=2 to the previous function call:

```r
emulator_plot(ems_wave1, plot_type = 'nimp', targets = targets, params = c('beta1', 'gamma'), cb=T, nth=2)
```

<img src="_main_files/figure-html/unnamed-chunk-72-1.png" style="display: block; margin: auto;" />

One of the advantages of history matching and emulation is that we are not obliged to emulate all outputs at each wave. This flexibility comes in handy, for example, when the emulator for a given output does not perform well at a certain wave: we can simply exclude that output and emulate it at a later wave. Another common situation where it may be useful to select a subset of emulators is when we have early outputs and late outputs, as in this workshop. It is often the case that later outputs have greater variability compared to earlier outputs, since they have more time to diverge. As a consequence, including emulators for later outputs in the first few waves may not be particularly convenient: it would both increase the number of calculations to make (since we would train more emulators), and would probably contribute to a lesser extent to the reduction of the parameter space.

We can focus on early times outputs (up to $t=200$), and produce implausibility plots for them:
  

```r
restricted_ems <- ems_wave1[c(1,2,3,4,7,8,9,10)]
emulator_plot(restricted_ems, plot_type = 'imp', targets = targets, params = c('beta1', 'gamma'), cb=T)
```

<img src="_main_files/figure-html/unnamed-chunk-73-1.png" style="display: block; margin: auto;" />

Finally let us set the unshown parameters to be as in `chosen_params`:

```r
emulator_plot(restricted_ems, plot_type = 'nimp', targets = targets[c(1,2,3,4,7,8,9,10)], 
              params = c('beta1', 'gamma'), 
              fixed_vals = chosen_params[!names(chosen_params) %in% c('beta1', 'gamma')], cb=T)+
 geom_point(aes(x=0.214, y=1/14), size=3)
```

<img src="_main_files/figure-html/unnamed-chunk-74-1.png" style="display: block; margin: auto;" />

The plot shows what we expected: when $\beta_1$ and $\gamma$ are equal to their values in `chosen_params`, i.e. $0.214$ and $1/14$, the implausibility measure is well below the threshold $3$ (cf. black point in the box). Note that when working with real models, one usually cannot check if the implausibility is low around fitting parameters, simply because these are not known. However, if one happens to have first hand fitted the model and has therefore a set of fitting parameters, then the above check can be performed.
</div></div></div>

# Emulator diagnostics

For a given set of emulators, we want to assess how accurately they reflect the model outputs over the input space. For a given validation set, we can ask the following questions:

- Within uncertainties, does the emulator output accurately represent the equivalent model output?

- Does the emulator adequately classify parameter sets as implausible or non-implausible?

- What are the standardised errors of the emulator outputs in light of the model outputs?

The function `validation_diagnostics` provides us with three diagnostics, addressing the three questions above.


```r
vd <- validation_diagnostics(ems_wave1$R200, validation = validation, targets = targets)
```

<img src="_main_files/figure-html/unnamed-chunk-32-1.png" style="display: block; margin: auto;" />

The first plot shows the emulator outputs plotted against the model outputs. In particular, the emulator expectation is plotted against the model output for each validation point, providing the dots in the graph. The emulator uncertainty at each validation point is shown in the form of a vertical interval that goes from $3\sigma$ below to $3\sigma$ above the emulator expectation, where $\sigma$ is the emulator variance at the considered point. The uncertainty interval can be expressed by the formula: $E[f(x)]\pm 3 \sqrt{Var(f(x)}$. An 'ideal' emulator would exactly reproduce the model results: this behaviour is represented by the green line $f(x)=E[f(x)]$ (this is a diagonal line, visible here only in the bottom left and top right corners). Any parameter set whose emulated prediction lies more than $3\sigma$ away from the model output is highlighted in red. Note that we do not need to have no red points for the test to be passed: since we are plotting $3\sigma$ bounds, statistically speaking it is ok to have up to $5\%$ of validation points in red (see [Pukelsheim's $3\sigma$ rule](https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule)). Apart from the number of points failing the diagnostic, it is also worth looking at whether the points that fail the diagnostic do so systematically. For example: are they all overestimates/underestimates of the model output?

The second column compares the emulator implausibility to the equivalent model implausibility (i.e. the implausibility calculated replacing the emulator output with the model output). There are three cases to consider:

- The emulator and model both classify a set as implausible or non-implausible (bottom-left and top-right quadrants). This is fine. Both are giving the same classification for the parameter set. 

- The emulator classifies a set as non-implausible, while the model rules it out (top-left quadrant): this is also fine. The emulator should not be expected to shrink the parameter space as much as the model does, at least not on a single wave. Parameter sets classified in this way will survive this wave, but may be removed on subsequent waves as the emulators grow more accurate on a reduced parameter space.

- The emulator rules out a set, but the model does not (bottom-right quadrant): these are the problem sets, suggesting that the emulator is ruling out parts of the parameter space that it should not be ruling out.

As for the first test, we should be alarmed only if we spot a systematic problem, with $5\%$ or more of the points in the bottom-right quadrant. Note, however, that it is always up to the user to decide how serious a misclassification is. For instance, a possible check is to identify points that are incorrectly ruled out by one emulator, and see if they would be considered non-implausible by all other emulators. If they are, then we should think about changing the misclassifying emulator.

Finally, the third column gives the standardised errors of the emulator outputs in light of the model output: for each validation point, the difference between the emulator output and the model output is calculated, and then divided by the standard deviation $\sigma$ of the emulator at the point. The general rule is that we want our standardised errors to be somewhat normally distributed around $0$, with $95\%$ of the probability mass between $-2$ and $2$. When looking at the standard errors plot, we should ask ourselves at least the following questions:

- Is more than $5\%$ of the probability mass outside the interval $[-2,2]$? If the answer is yes, this means that, even factoring in all the uncertainties in the emulator and in the observed data, the emulator output is too often far from the model output. 

- Is $95\%$ of the probability mass concentrated in a considerably smaller interval than $[-2,2]$ (say, for example, $[-0.5,0.5]$)? For this to happen, the emulator uncertainty must be quite large. In such case the emulator, being extremely cautious, will cut out a small part of the parameter space and we will end up needing many more waves of history matching than are necessary, or, even worse, we just won't be able to reduce the non-implausible parameter space.

- Is the histogram skewing significantly in one direction or the other? If this is the case, the emulator tends to either overestimate or underestimate the model output.

The diagnostic plot on the left shows very few points in red: this is a good result! The plot in the middle has only one red point, which is promising too. Finally, in the third diagnostic if we consider all standardised errors below $-2$ and above $2$, we clearly get more than $5\%$ of all errors: this is not ideal, since we would like the standardised errors to have $95\%$ of the probability mass within $\pm 2$.

A way of improving the performance of an emulator is by changing the variance $\sigma^2$ in the Gaussian process $u$:
$$\sigma^2 \left[(1-\delta) c(x,x^{\prime}) + \delta I_{\{x=x^\prime\}}\right].$$
The lower the value of $\sigma$, the more 'certain' the emulator will be. This means that when an emulator is a little too overconfident (as in our case above), we can try increasing $\sigma$. Below we train a new emulator setting $\sigma$ to be twice as much as its default value, through the method `mult_sigma`:



```r
sigmadoubled_emulator <- ems_wave1$R200$mult_sigma(2)
vd <- validation_diagnostics(sigmadoubled_emulator, 
                             validation = validation, targets = targets)
```

<img src="_main_files/figure-html/unnamed-chunk-34-1.png" style="display: block; margin: auto;" />

A higher value of $\sigma$ has therefore allowed us to build a more conservative emulator that performs better than before.


<div class="panel panel-default"><div class="panel-heading"> Task </div><div class="panel-body"> 
Explore different values of $\sigma$. What happens for very small/large values of $\sigma$?
  
<span class="abbr" title=""><abbr title="If `em` is an emulator, you can change its sigma by a factor x through the following line of code: `ems$mult_sigma(x)`
">
R tip</abbr></span> </div></div>

<button id="displayTextunnamed-chunk-36" onclick="javascript:toggle('unnamed-chunk-36');">Show: Solution</button>

<div id="toggleTextunnamed-chunk-36" style="display: none"><div class="panel panel-default"><div class="panel-heading panel-heading1"> Solution </div><div class="panel-body">

Let us set $\sigma$ to be ten times smaller than its default value:
  

```r
tinysigma_emulator <- ems_wave1$R200$mult_sigma(0.1)
vd <- validation_diagnostics(tinysigma_emulator, validation = validation, targets = targets)
```

<img src="_main_files/figure-html/unnamed-chunk-76-1.png" style="display: block; margin: auto;" />

In this case we built a very overconfident emulator. This is shown by the very small uncertainty intervals in the first column: as a consequence many points are in red. Similarly, if we look at the third column we notice that the standardised errors are extremely large, well beyond the value of $2$. 

Let us now set $\sigma$ to be ten times larger than its default value: 
  

```r
hugesigma_emulator <- ems_wave1$R200$mult_sigma(10)
vd <- validation_diagnostics(hugesigma_emulator, validation = validation, targets = targets)
```

<img src="_main_files/figure-html/unnamed-chunk-78-1.png" style="display: block; margin: auto;" />

With this choice of $\sigma$, we see that our emulator is extremely cautious. If we look at the plot in the middle, we see that now a lot of points in the validation set have an implausibility less or equal to $3$. This implies that this emulator will reduce the input space slowly. As explained above, having consistent very small standardised errors is not positive: it implies that, even though we trained a regression hypersurface in order to catch the global behaviour of the output, the sigma is so large that the emulator is being dominated by the correlation structure. This means at best that we will have to do many more waves of history matching than are necessary, and at worst that our emulators won’t be able to reduce the non-implausible parameter space.

The above exploration highlights the importance of finding a value of $\sigma$ that produces an emulator which on one hand is not overconfident and on the other is able to quickly reduce the input space. Note that there is not a universal rule to navigate this tradeoff: the role of the scientist's judgment is fundamental.  
</div></div></div>

# Proposing new points

The function `generate_new_runs` is designed to generate new sets of parameters; its default behaviour is as follows.

STEP 1. If  prior parameter sets are provided, step 1 is skipped, otherwise a set is generated using a [Latin Hypercube Design](https://en.wikipedia.org/wiki/Latin_hypercube_sampling), rejecting implausible parameter sets. Figure \@ref(fig:lhsamp) shows an example of LH sampling (not for our model), with implausible parameter sets in red and non-implausible ones in blue:

<div class="figure" style="text-align: center">
<img src="LHSamp.png" alt="LH sampling" width="75%" />
<p class="caption">(\#fig:lshamp)LH sampling</p>
</div>

STEP 2. Pairs of parameter sets are selected at random and more sets are sampled from lines connecting them, with particular importance given to those that are close to the non-implausible boundary. Figure \@ref(fig:linesamp) shows line sampling on top of the previously shown LH sampling: non-implausible parameter sets provided by LH sampling are in grey, parameter sets proposed by line sampling are in red if implausbile and in blue if non-implausible.

<div class="figure" style="text-align: center">
<img src="LineSamp.png" alt="Line sampling" width="75%" />
<p class="caption">(\#fig:linesamp)Line sampling</p>
</div>

STEP 3. Using these as seeding points, more parameter sets are generated using [importance sampling](https://en.wikipedia.org/wiki/Importance_sampling) to attempt to fully cover the non-implausible region. Figure \@ref(fig:impsamp) has non-implausible parameter sets provided by LH and line sampling in grey and non-implausibile parameter sets found by importance sampling in blue. 

<div class="figure" style="text-align: center">
<img src="ImpSamp.png" alt="Importance sampling" width="75%" />
<p class="caption">(\#fig:impsamp)Importance sampling</p>
</div>

The combination of the three steps above brings to the following final design of non-implausible parameter sets:

<div class="figure" style="text-align: center">
<img src="FinalSamp.png" alt="Overall design of new non-implausible parameter sets"  />
<p class="caption">(\#fig:unnamed-chunk-37)Overall design of new non-implausible parameter sets</p>
</div>

Let us generate $180$ new sets of parameters, using the emulators for the time up to $t=200$. Note that the generation of new points might require a few minutes.


```r
new_points_restricted <- generate_new_runs(restricted_ems, 180, targets, verbose=TRUE)
```

```
## [1] "Performing Latin Hypercube sampling..."
## [1] "Proposing at implausibility I=3.5"
## [1] "47 points generated from LHS at I=3.5"
## [1] "Performing line sampling..."
## [1] "Line sampling generated 9 more points."
## [1] "Performing importance sampling..."
## [1] "Importance sampling generated 207 more points."
## [1] "Selecting final points using maximin criterion..."
## [1] "128 points generated from LHS at I=3"
## [1] "Performing line sampling..."
## [1] "Line sampling generated 13 more points."
## [1] "Performing importance sampling..."
## [1] "Importance sampling generated 47 more points."
## [1] "Selecting final points using maximin criterion..."
## [1] "Resampling: Performing line sampling..."
## [1] "Line sampling generated 18 more points."
## [1] "Resampling: Performing importance sampling..."
## [1] "Importance sampling generated 78 more points."
## [1] "Selecting final points using maximin criterion..."
```

We now plot `new_points_restricted` through `plot_wrap`. Note that we pass `ranges` too to `plot_wrap` to ensure that 
the plot shows the entire range for each parameter: this allows us to see how the new set of parameters compares with respect to the original input space.


```r
plot_wrap(new_points_restricted, ranges)
```

<img src="_main_files/figure-html/unnamed-chunk-39-1.png" style="display: block; margin: auto;" />
     
By looking at the plot we can learn a lot about the non-implausible space. For example, it seems clear that low values of $\gamma$ cannot produce a match (cf. fourth column). We can also deduce relationships between parameters: $\beta_1$ and $\sigma$ are an example of negatively-correlated parameters. If $\beta_1$ is large then $\sigma$ needs to be small, and vice versa.

<div class="panel panel-default"><div class="panel-heading"> Task </div><div class="panel-body"> 
What do you think would happen if we tried generating new parameter sets using all emulators, instead of those relating to times up to $t=200$ only? Try changing the code and seeing what happens. </div></div>

<button id="displayTextunnamed-chunk-41" onclick="javascript:toggle('unnamed-chunk-41');">Show: Solution</button>

<div id="toggleTextunnamed-chunk-41" style="display: none"><div class="panel panel-default"><div class="panel-heading panel-heading1"> Solution </div><div class="panel-body">
We expect that using all emulators (instead of just a subset of them) will make the non-implausible space smaller than before, since more conditions will have to be met.  Let us check if our intuition corresponds to what the plots show:
  

```r
new_points <- generate_new_runs(ems_wave1, 180, targets, verbose = TRUE)
```

```
## [1] "Performing Latin Hypercube sampling..."
## [1] "Proposing at implausibility I=4"
## [1] "52 points generated from LHS at I=4"
## [1] "Performing line sampling..."
## [1] "Line sampling generated 24 more points."
## [1] "Performing importance sampling..."
## [1] "Importance sampling generated 216 more points."
## [1] "Selecting final points using maximin criterion..."
## [1] "64 points generated from LHS at I=3"
## [1] "Performing line sampling..."
## [1] "Line sampling generated 17 more points."
## [1] "Performing importance sampling..."
## [1] "Importance sampling generated 233 more points."
## [1] "Selecting final points using maximin criterion..."
## [1] "Resampling: Performing line sampling..."
## [1] "Line sampling generated 24 more points."
## [1] "Resampling: Performing importance sampling..."
## [1] "Importance sampling generated 161 more points."
## [1] "Selecting final points using maximin criterion..."
```

```r
plot_wrap(new_points, ranges)
```

<img src="_main_files/figure-html/fignewpts-1.png" style="display: block; margin: auto;" />
The answer is yes: points now look slightly less spread than before, i.e. a higher part of the input space has been cut-out.   </div></div></div>

In order to quantify how much of the input space is removed by a given set of emulators, we can use the function `space_removed`, which takes a list of emulators and a list of targets we want to match to. The output is a plot that shows the percentage of space that is removed by the emulators as a function of the implausibility cut-off. Note that here we also set the argument `ppd`, which determines the number of points per input dimension to sample at. In this workshop we have $9$ parameters, so a `ppd` of $3$ means that `space_removed` will test the emulators on $3^9$ sets of parameters.


```r
space_removed(ems_wave1, targets, ppd=3) + geom_vline(xintercept = 3, lty = 2) + 
  geom_text(aes(x=3, label="x = 3",y=0.33), colour="black", 
  angle=90, vjust = 1.2, text=element_text(size=11))
```

```
## Warning: Ignoring unknown parameters: text
```

<img src="_main_files/figure-html/unnamed-chunk-42-1.png" style="display: block; margin: auto;" />

By default the plot shows the percentage of space that is deemed implausible both when the observational errors are exactly the ones in `targets` and when the observational errors are $80\%$ (resp. $90\%$, $110\%$ and $120\%$)  of the values in `targets`. Here we see that with an implausibility cut-off of $3$, the percentage of space removed is around $99\%$. If we use a cut-off of $5\%$ instead, then we would discard around $93\%$ (resp. $90\%$) of the space with the observational errors provided by `targets` (resp. with observational errors equal to $120\%$ of the values in `targets`). As expected, larger observational errors correspond to lower percentages of space removed.

# Customise the first wave

<div class="panel panel-default"><div class="panel-heading"> Task </div><div class="panel-body"> 
Now that we have learnt how to customise the various steps of the process, try to improve the performance of the first wave of emulation and history matching. First, have a look at the emulator diagnostics and see if you can improve the performance of the emulators. Then generate new points using your improved emulators, and compare them to those shown in the task at the end of last section.  </div></div>

<button id="displayTextunnamed-chunk-44" onclick="javascript:toggle('unnamed-chunk-44');">Show: Solution</button>

<div id="toggleTextunnamed-chunk-44" style="display: none"><div class="panel panel-default"><div class="panel-heading panel-heading1"> Solution </div><div class="panel-body">
First of all let us take a look at the diagnostics of all the emulators trained in wave 1:
  

```r
vd <- validation_diagnostics(ems_wave1, validation = validation, targets = targets)
```

<img src="_main_files/figure-html/unnamed-chunk-79-1.png" style="display: block; margin: auto;" /><img src="_main_files/figure-html/unnamed-chunk-79-2.png" style="display: block; margin: auto;" /><img src="_main_files/figure-html/unnamed-chunk-79-3.png" style="display: block; margin: auto;" /><img src="_main_files/figure-html/unnamed-chunk-79-4.png" style="display: block; margin: auto;" />

All emulators seem to have more than $5\%$ of the probability mass for the standardised errors outside $[-2,2]$. We would then benefit from having slightly more conservative emulators, which we can obtain by increasing $\sigma$. After some trial and error, we chose the following values of sigma for our emulators:

```r
inflations <- c(4,2,4,2,2,2,2,4,4,2,2,4)
for (i in 1:length(ems_wave1)) {
  ems_wave1[[i]] <- ems_wave1[[i]]$mult_sigma(inflations[[i]])
}
```

```r
vd <- validation_diagnostics(ems_wave1, validation = validation, targets = targets)
```

<img src="_main_files/figure-html/unnamed-chunk-81-1.png" style="display: block; margin: auto;" /><img src="_main_files/figure-html/unnamed-chunk-81-2.png" style="display: block; margin: auto;" /><img src="_main_files/figure-html/unnamed-chunk-81-3.png" style="display: block; margin: auto;" /><img src="_main_files/figure-html/unnamed-chunk-81-4.png" style="display: block; margin: auto;" />

The diagnostics look good now. We are then ready to generate new points with our improved emulators:

```r
new_points <- generate_new_runs(ems_wave1, 180, targets, verbose = TRUE)
```

```
## [1] "Performing Latin Hypercube sampling..."
## [1] "75 points generated from LHS at I=3"
## [1] "Performing line sampling..."
## [1] "Line sampling generated 19 more points."
## [1] "Performing importance sampling..."
## [1] "Importance sampling generated 194 more points."
## [1] "Selecting final points using maximin criterion..."
## [1] "Resampling: Performing line sampling..."
## [1] "Line sampling generated 16 more points."
## [1] "Resampling: Performing importance sampling..."
## [1] "Importance sampling generated 87 more points."
## [1] "Selecting final points using maximin criterion..."
```

```r
plot_wrap(new_points, ranges)
```

<img src="_main_files/figure-html/unnamed-chunk-82-1.png" style="display: block; margin: auto;" />

If we compare the parameter sets we just generated with those generated using non-customised emulators, we note that the space has now been reduced less than before. This happened because our customisation helped us to build more conservative emulators, decreasing the risk of rejecting good parts of the input space. Building emulators carefully ensures that we end up with a set of points that are truly representative of the set of all points that fit the data (rather than a subset of it).</div></div></div>

# Second wave

To perform a second wave of history matching and emulation we follow the same procedure as in the previous sections, with two caveats. We start by forming a dataframe `wave1` using parameters sets in `new_points`, as we did with `wave0`, i.e. we evaluate the function `get_results` on `new_points` and then bind the obtained outputs to `new_points`. Half of `wave1` should be used as the training set for the new emulators, and the other half as the validation set to evaluate the new emulators' performance. Note that when dealing with computationally expensive models, using the same number of points for the training and validation sets may not feasible. If $p$ is the number of parameters, a good rule of thumb is to build a training set with at least $10p$ points, and a validation set with at least $p$ points. 

Now note that parameter sets in `new_points` tend to lie in a small region inside the original input space, since `new_points` contains only non-implausible points, according to the first wave emulators. The first caveat is then that it is preferable to train the new emulators only on the non-implausible region found in wave one. To do this we define new ranges for the parameters:


```r
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
```

The list `new_ranges` contains lower and upper bounds for each parameter. The upper bound for a given parameter is determined in the following way: 

- Among all points in `new_points`, the maximum value of the parameter is identified.

- $5\%$ of the original range of the parameter is added to the maximum value found in the previous point. This step enlarges slightly the new range, to make sure that we are including all non-implausible points.

- The minimum between the value found above and the upper bound in the original `ranges` list is selected: this ensure that we do not end up with a new upper bound which is larger than the original one.

A similar calculation was adopted to determine the new lower bounds of parameters.

Since wave two emulators will be trained only on the non-implausible space from wave one, their implausibility cannot be assessed everywhere in the original input space. For this reason, the second caveat is that when generating new parameter sets at the end of wave two, we must consider implausibility across both wave one and wave two emulators, i.e. we need to pass emulators from both waves to the function `generate_new_points`. To do this we simply make the first argument of `generate_new_points` equal to the vector `c(ems_wave2, ems_wave1)`, where `ems_wave2` are the wave two emulators. Note that `generate_new_points` picks the ranges of parameters from the first emulator in the list. For this reason, it is important to put the second wave emulators first (in `c(ems_wave2, ems_wave1)`), which have smaller ranges. 

In the task below, you can have a go at wave 2 of the emulation and history matching process yourself.

<div class="panel panel-default"><div class="panel-heading"> Task </div><div class="panel-body"> 
Using `new_points` and `new_ranges`, train new emulators. Customise them and generate new parameter sets. </div></div>

<button id="displayTextunnamed-chunk-47" onclick="javascript:toggle('unnamed-chunk-47');">Show: Solution</button>

<div id="toggleTextunnamed-chunk-47" style="display: none"><div class="panel panel-default"><div class="panel-heading panel-heading1"> Solution </div><div class="panel-body">
We start by evaluating the function `get_results` on `new_points`. In other words, we run the model using the parameter sets we generated at the end of wave 1:

```r
new_initial_results <- setNames(data.frame(t(apply(new_points, 1, get_results, 
                                               c(25, 40, 100, 200, 300, 350), c('I', 'R')))), names(targets))
```

and then binding `new_points` to the model runs `new_initial_results` to create the data.frame `wave1`, containing the input parameter sets and model outputs:

```r
wave1 <- cbind(new_points, new_initial_results)
```

We split `wave1` into training and validation sets
 

```r
new_t_sample <- sample(1:nrow(wave1), 90)
new_training <- wave1[new_t_sample,]
new_validation <- wave1[-new_t_sample,]
```
 
and train wave two emulators on the space defined by `new_ranges`: 

```r
ems_wave2 <- emulator_from_data(new_training, names(targets), new_ranges)
```
Let us check their diagnostics:
  

```r
vd <- validation_diagnostics(ems_wave2, validation = new_validation, targets = targets)
```

<img src="_main_files/figure-html/unnamed-chunk-87-1.png" style="display: block; margin: auto;" /><img src="_main_files/figure-html/unnamed-chunk-87-2.png" style="display: block; margin: auto;" /><img src="_main_files/figure-html/unnamed-chunk-87-3.png" style="display: block; margin: auto;" /><img src="_main_files/figure-html/unnamed-chunk-87-4.png" style="display: block; margin: auto;" />
 
All emulators fail at least one of the three dagnostics and most of them fail the first and third diagnostics quite badly. Let us modify the sigmas in order to build more conservative emulators that pass the three diagnostics. After some trial and error, we chose the following values of sigma for our emulators:

```r
inflations <- c(4,6,4,4,4,2,4,4,4,4,4,4)
for (i in 1:length(ems_wave2)) {
ems_wave2[[i]] <- ems_wave2[[i]]$mult_sigma(inflations[[i]])
}
vd <- validation_diagnostics(ems_wave2, validation =  new_validation, targets = targets)
```

<img src="_main_files/figure-html/unnamed-chunk-88-1.png" style="display: block; margin: auto;" /><img src="_main_files/figure-html/unnamed-chunk-88-2.png" style="display: block; margin: auto;" /><img src="_main_files/figure-html/unnamed-chunk-88-3.png" style="display: block; margin: auto;" /><img src="_main_files/figure-html/unnamed-chunk-88-4.png" style="display: block; margin: auto;" />

The diagnostics look ok now. Let us try to generate new parameter sets using all emulators build so far:
  

```r
new_new_points <- generate_new_runs(c(ems_wave2, ems_wave1), 180, targets, verbose=TRUE)
```

```
## [1] "Performing Latin Hypercube sampling..."
## [1] "Proposing at implausibility I=4"
## [1] "52 points generated from LHS at I=4"
## [1] "Performing line sampling..."
## [1] "Line sampling generated 21 more points."
## [1] "Performing importance sampling..."
## [1] "Importance sampling generated 223 more points."
## [1] "Selecting final points using maximin criterion..."
## [1] "68 points generated from LHS at I=3"
## [1] "Performing line sampling..."
## [1] "Line sampling generated 28 more points."
## [1] "Performing importance sampling..."
## [1] "Importance sampling generated 90 more points."
## [1] "Selecting final points using maximin criterion..."
## [1] "Resampling: Performing line sampling..."
## [1] "Line sampling generated 23 more points."
## [1] "Resampling: Performing importance sampling..."
## [1] "Importance sampling generated 75 more points."
## [1] "Selecting final points using maximin criterion..."
```

```r
plot_wrap(new_new_points, ranges)
```

<img src="_main_files/figure-html/unnamed-chunk-89-1.png" style="display: block; margin: auto;" />

This worked well: the new non-implausible region is clearly smaller than the one we had at the end of wave one. In the next section we will show how to make visualisations to direcly compare the non-implausible space at the end of different waves of the process.
</div></div></div>

As we did for the wave 1 emulators, let us check the values of the adjusted $R^2$ for the new emulators:


```r
R_squared_new <- list()
for (i in 1:length(ems_wave2)) {
  R_squared_new[[i]] <- summary(ems_wave2[[i]]$model)$adj.r.squared
}
names(R_squared_new) <- names(ems_wave2)
unlist(R_squared_new)
```

```
##       I25       I40      I100      I200      I300      I350       R25       R40 
## 0.9998347 0.9993282 0.9988811 0.9940101 0.9950220 0.9310200 0.9999403 0.9997465 
##      R100      R200      R300      R350 
## 0.9992135 0.9996098 0.9844951 0.9942361
```

All $R^2$ values are very high, meaning that the regression term is contributing far more than the residuals. As all of the emulators we have seen so far have had high $R^2$ values, we have not discussed the customisation of $\theta$. We now want to briefly comment on what happens when instead the $R^2$ are lower and the residuals play a more substantial role. In such cases, the extent to which residuals at different parameter sets are correlated is a key ingredient in the training of emulators, since it determines how informative the model outputs at training parameter sets are. For example, if residuals are highly correlated even at parameter sets that are far away from each other, then knowing the model output at a given parameter set gives us information about a wide region around it. This results in rather confident emulators, which cut a lot of space out. If instead residuals are correlated only for parameter sets that are close to each other, then knowing the model output at a given parameter set gives us information about a small region around it. This creates more uncertain emulators, which can rule out a lot of input parameter space. It is then clear that when we don't have very high $R^2$ values, we can use $\theta$ to increase or decrease the amount  of space cut out by emulators. 

In practice, the default value of $\theta$ is chosen very carefully by the HoMER package, and most users calibrating deterministic models will not have to vary the value of $\theta$ for most of their emulators. If, however, you find that the non-implausible space is shrinking very slowly, particularly in later waves (see section 9 for details of how to check this), then the value of $\theta$ may be too conservative. If this occurs, then you can increase the $\theta$ of the emulators to increase the rate at which space is reduced. You should only do this if you are confident that your outputs are smooth enough to justify the choice of $\theta$ however, or you risk the emulators incorrectly excluding space when model fits could be found. We discuss the choice of $\theta$ further in workshop on calibrating stochastic models.

To see all this in practice, we train new wave one emulators, assuming a linear regression term by setting `quadratic=FALSE`: 


```r
 ems_wave1_linear <- emulator_from_data(training, names(targets), 
                                        ranges, quadratic=FALSE)
 R_squared_linear <- list()
 for (i in 1:length(ems_wave1_linear)) {
   R_squared_linear[[i]] <- summary(ems_wave1_linear[[i]]$model)$adj.r.squared
 }
 names(R_squared_linear) <- names(ems_wave1_linear)
 unlist(R_squared_linear)
```

```
##       I25       I40      I100      I200      I300      I350       R25       R40 
## 0.9697263 0.9491857 0.9007032 0.7157199 0.8191644 0.1823831 0.9827103 0.9809724 
##      R100      R200      R300      R350 
## 0.9310675 0.9645519 0.8112647 0.8390023
```
 
By forcing the regression hypersurface to be linear, we obtain emulators where the global term is not sufficient to explain the model output on its own. As a rough guide, $R^2$ values of below $0.9$ indicate that the residuals are playing an important role. Let us see what happens if we plot the variance and the implausibility for the linear $I200$ emulator before and after increasing its $\theta$ by a factor of $3$:
   

```r
 emulator_plot(ems_wave1_linear$I200, plot_type = 'var', 
               params = c('beta1', 'gamma'))
```

<img src="_main_files/figure-html/unnamed-chunk-50-1.png" style="display: block; margin: auto;" />


```r
emulator_plot(ems_wave1_linear$I200, plot_type = 'imp', targets = targets, 
              params = c('beta1', 'gamma'), cb=TRUE)
```

<img src="_main_files/figure-html/unnamed-chunk-51-1.png" style="display: block; margin: auto;" />
 

```r
ems_wave1_linear$I200 <- ems_wave1_linear$I20$set_hyperparams(
              list(theta=ems_wave1_linear$I200$corr$hyper_p$theta *3))
emulator_plot(ems_wave1_linear$I200, plot_type = 'var', 
              params = c('beta1', 'gamma'))
```

<img src="_main_files/figure-html/unnamed-chunk-52-1.png" style="display: block; margin: auto;" />


```r
emulator_plot(ems_wave1_linear$I200, plot_type = 'imp', targets = targets, 
              params = c('beta1', 'gamma'), cb=TRUE)
```

<img src="_main_files/figure-html/unnamed-chunk-53-1.png" style="display: block; margin: auto;" />

First of all, the blue-purple area in the variance plot becomes larger after $\theta$ increased: this shows that a higher $\theta$ results in the model output at training points influencing a wider region around itself. Second, we see that a higher $\theta$ 
causes the implausibility measure to have higher values: as a consequence, more space will be ruled out as implausible.

# Multi-wave visualisations

In this last section we present two visualisations that can be used to compare the non-implausible space identified at different waves of the process.

The first visualisation, obtained through the function `wave_points`, shows the distribution of the non-implausible space for the waves of interest. For example, let us plot the distribution of parameter sets at the beginning, at the end of wave one and at the end of wave two:


```r
wave_points(list(initial_points, new_points, new_new_points), input_names = names(ranges))
```

<img src="_main_files/figure-html/unnamed-chunk-54-1.png" style="display: block; margin: auto;" />

Here `initial_points` are in yellow, `new_points` are in green and `new_new_points` are in purple. The plots in the main diagonal show the distribution of each parameter singularly: we can easily see that the distributions tend to become more and more narrow wave after wave. In the off-diagonal boxes we have plots for all possible pairs of parameters. Again, the non-implausible region identified at the end of each wave clearly becomes smaller and smaller.

The second visualisation allows us to assess how much better parameter sets at later waves perform compared to the original `initial_points`. Let us first create the dataframe `wave2`:



```r
new_new_initial_results <- setNames(data.frame(t(apply(new_new_points, 1, 
                                    get_results, c(25, 40, 100, 200, 300, 350), 
                                    c('I', 'R')))), names(targets))
wave2 <- cbind(new_new_points, new_new_initial_results)
```

We now produce the plots using the function `simulator_plot`:


```r
all_points <- list(wave0, wave1, wave2)
simulator_plot(all_points, targets)
```

<img src="_main_files/figure-html/unnamed-chunk-56-1.png" style="display: block; margin: auto;" />

We can see that, compared to the space-filling random parameter sets used to train the first emulators, the new parameter sets are in much closer agreement with our targets. While the `initial_points` did not match any of the targets, we can see that at the end of wave two, we have several targets matched (I25, I40, R25, R40, R100, R200). Subsequent waves, trained on the new parameter sets, will be more confident in the new non-implausible region: this will allow them to refine the region and increase the number of targets met.

<div class="panel panel-default"><div class="panel-heading"> Task </div><div class="panel-body"> 
In the plot above, some targets are easier to read than others: this is due to the targets having quite different values and ranges. To help with this, `simulator_plot` has the argument `normalize`, which can be set to `TRUE` to rescale the target bounds in the plot. Similarly, the argument `logscale` can be used to plot log-scaled target bounds. Explore these options and get visualisations that are easier to interpret. </div></div>

<button id="displayTextunnamed-chunk-58" onclick="javascript:toggle('unnamed-chunk-58');">Show: Solution</button>

<div id="toggleTextunnamed-chunk-58" style="display: none"><div class="panel panel-default"><div class="panel-heading panel-heading1"> Solution </div><div class="panel-body">

```r
simulator_plot(all_points, targets, normalize = TRUE)
```

<img src="_main_files/figure-html/unnamed-chunk-90-1.png" style="display: block; margin: auto;" />

```r
simulator_plot(all_points, targets, logscale = TRUE)
```

<img src="_main_files/figure-html/unnamed-chunk-90-2.png" style="display: block; margin: auto;" />
These normalised/logscaled plots allow us to investigate better targets such as I200: it is now clear that this is not matched yet, even at the end of wave two. </div></div></div>



<!--
# TJ material



```
rmd_files:
    html: ['index.Rmd', 'ch1.Rmd']
    latex: ['index.Rmd', 'ch1.Rmd', 'ch_appendix.Rmd']
output_dir: "docs"
```


-->

<!--chapter:end:index.Rmd-->
