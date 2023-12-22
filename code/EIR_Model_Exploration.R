# This script shows the dynamics for the temperature/rainfall relatiomships in Ukuwaba et al

#### parameter values ####

T = 30
# eggs laid per gonotrophic cycle
epsilon = 50
# adult mosquito mortality (no intervention)
mu_m = -log(.98)

# egg-adult sensitivity to surface moisture
a.R = .5
# standardizd anomaly of cumulative weekly rainfall- NO DEFAULT, IS DATA
c.RD = 1
# mean anomaly level at which 0.5 survival is expected
b.R = 0 


# epidemiological parameters
PHM = 0.125
PMH = 0.5
X = .01


#### functions for model ####

# length of the gonotrophic cycle
GP = function(T){
  1 / (.017 * T - .165)
}

# eggs produced over the lifespan
B = function(T, epsilon, mu_m){
  epsilon / (exp(GP(T) * mu_m) - 1)
}


# egg-to-adult survivorship due to rainfall
PEAr = function(R, a.R, c.RD, b.R){
  1 / (1 + exp(-a.R*(c.RD-b.R)))
}

# water temperature relative to air temperature
Tw = function(T){
  return(T + 2)
}

# rate of immature survival
PEAt = function(T){
  TW = Tw(T)
  return(-.00924*(TW^2) + .453 * TW - 4.77)
}

# duration of egg-adult stage
tau_EA = function(T){
  TW = Tw(T)
  return(((.000111*TW*(TW - 14.7)) * sqrt(34 - TW))^(-1))
}

# mosquito biting rate
a = function(T){
  0.017*T - .165
}

# Sporogony
n = function(T){
  (.0001122*T*(T - 15.384) * sqrt(35-T))^(-1)
}

L = function(T, R, epsilon, mu_m){
  B(T, epsilon, mu_m) * PEAt(T) / tau_EA(T)
}

# mosquito density
m = function(T, R, epsilon, mu_m){
  L(T, R, epsilon, mu_m) / mu_m
}


# primary function
EIR = function(T, X, PHM, R, epsilon, mu_m){
  numerator = m(T, R, epsilon, mu_m) * (a(T)) ^ 2 * PHM * X * exp(-mu_m * n(T))
  denominator = mu_m + a(T) * PHM *X
  return(numerator / denominator)
}

# convert to FOI
mu_SE = function(T, X, PHM, R, epsilon, mu_m, PMH){
  1 - exp(-PMH * EIR(T, X, PHM, R, epsilon, mu_m))
}

mu_SE(T, X, PHM, R, epsilon, mu_m, PMH)
