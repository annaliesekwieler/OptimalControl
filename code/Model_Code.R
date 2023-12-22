# This script codes the ODE model 

#### load packages ####

library(deSolve)

#### functions for mosquito model ####

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

#### basic model (no weather) ####

SEI_model = function(time, y, params){
  with(as.list(c(y, params)), {
    
    FOI = temps[time]
    
    # rate of change
    dS = -FOI*S - deltaT*S + mut + rc*Ic + rs*Is
    dE = FOI*S - tau*E - deltaT*E
    dIc = -rc*Ic + c * tau * E - deltaT*Ic
    dIs = -rs * Is + (1-c) * tau * E - deltaT*Is
    
    return(list(c(dS, dE, dIc, dIs)))
    
  })
  
}


params = c(
  rc = 1/11,
  rs = 1/240,
  tau = 1/10,
  c = .1,
  FOI = .01,
  deltaT = .0001,
  mut = .0001
)

state = c(S = .9, E = .01, Ic = .05, Is = .04)

times = seq(1, 1000, by =1)

out = data.frame(ode(y = state, times = times, func = SEI_model, parms = params))
plot(out$time, out$S, type="l", xlim=c(0, 1000), ylim=c(0,1), lwd=2, xlab = "Time (Days)", 
     ylab="Proportion of Population")
lines(out$time, out$E, col="red", lwd=2)
lines(out$time, out$Ic, col="green", lwd=2)
lines(out$time, out$Is, col="blue", lwd=2)
legend(400, .5, c("S", "E", "Ic", "It"), col=c("black", "red", "green", "blue"), lwd=2)


#### adding in time-varying FOI ####

# dummy vector of mean temperatures for the day
temps = c(rep(.01, 500), rep(.02, 500))

# parameters sans FOI
params = c(
  rc = 1/11,
  rs = 1/240,
  tau = 1/10,
  c = .1,
  deltaT = .0001,
  mut = .0001
)

out = data.frame(ode(y = state, times = times, func = SEI_model, parms = params))
plot(out$time, out$S, type="l", xlim=c(0, 1000), ylim=c(0,1), lwd=2, xlab = "Time (Days)", 
     ylab="Proportion of Population")
lines(out$time, out$E, col="red", lwd=2)
lines(out$time, out$Ic, col="green", lwd=2)
lines(out$time, out$Is, col="blue", lwd=2)
legend(400, .5, c("S", "E", "Ic", "Is"), col=c("black", "red", "green", "blue"), lwd=2)

