
#### functions for mosquito model ####


# length of the gonotrophic cycle
GP = function(temp){
  1 / (.017 * temp - .165)
}

# eggs produced over the lifespan
B = function(temp, epsilon, mu_m){
  epsilon / (exp(GP(temp) * mu_m) - 1)
}


# egg-to-adult survivorship due to rainfall
PEAr = function(R, a.R, c.RD, b.R){
  1 / (1 + exp(-a.R*(c.RD-b.R)))
}

# water temperature relative to air temperature
Tw = function(temp){
  return(temp + 2)
}

# rate of immature survival
# temp and cRD are data
PEAt_PEAr = function(temp, c.RD, aR, bR){
  TW = Tw(temp)
  PEAt = -.00924*(TW^2) + .453 * TW - 4.77
  
  PEAr = 1 / (1 + exp(-aR * (c.RD - bR)))
  
  return(PEAt * PEAr)
  
}

# duration of egg-adult stage
tau_EA = function(temp){
  TW = Tw(temp)
  TW[TW > 34] = 34
  return(((.000111*TW*(TW - 14.7)) * sqrt(34 - TW))^(-1))
}

# mosquito biting rate
a = function(temp){
  0.017*temp - .165
}

# Sporogony
n = function(temp){
  (.0001122*temp*(temp - 15.384) * sqrt(35-temp))^(-1)
}

L = function(temp, epsilon, mu_m, c.RD, aR, bR){
  B(temp, epsilon, mu_m) * PEAt_PEAr(temp, c.RD, aR, bR) / tau_EA(temp)
}

# mosquito density
m = function(temp, epsilon, mu_m, c.RD, aR, bR){
  L(temp, epsilon, mu_m, c.RD, aR, bR) / mu_m
}

# multiply these by epsilon to get true value
B_no_epsilon = function(temp,  mu_m){
  1 / (exp(GP(temp) * mu_m) - 1)
}

L_no_epsilon = function(temp, mu_m, c.RD, aR, bR){
  B_no_epsilon(temp,  mu_m) * PEAt_PEAr(temp, c.RD, aR, bR) / tau_EA(temp)
}

m_no_epsilon = function(temp, mu_m, c.RD, aR, bR){
  L_no_epsilon(temp, mu_m, c.RD, aR, bR) / mu_m
}


# primary function
EIR = function(temp, X, PHM, epsilon, mu_m, c.RD, aR, bR){
  numerator = m(temp, epsilon, mu_m, c.RD, aR, bR) * (a(temp)) ^ 2 * PHM * X * exp(-mu_m * n(temp))
  denominator = mu_m + a(temp) * PHM *X
  return(numerator / denominator)
}

# convert to FOI
mu_SE = function(temp, X, PHM, epsilon, mu_m, c.RD, aR, bR, PMH){
  1 - exp(-PMH * EIR(temp, X, PHM, epsilon, mu_m, c.RD, aR, bR))
}

#### set default parameter values #####

# eggs laid per gonotrophic cycle
# epsilon = 50
# adult mosquito mortality (no intervention)
mu_m = -log(.98)

# epidemiological parameters

# transmission probabilities
PHM = 0.125
# PMH = 0.5

# incubation rate in humans
tau = 1/10


# cumulative rainfall parameters
aR = 0.5

# recovery rates
# rc = 1/5
# rs = 1 / 100


#### function that optimizes bR  for a particular raster ####

# function of the temperature, rainfall, and proportion of cases treated
# returns the bR that is optimal for that raster

# requires all other parameters to be specified globally

SEI_model = odin::odin({
  
  # Derivatives 
  
  X = Ic + Is
  
  M_true = M * epsilon
  
  FOI = 1 - 2.71828182845904523536 ^ (-PMH * (M_true * (A ^ 2) * PHM * X *
                                                (2.71828182845904523536 ^ (-mu_m * N)))/
                                        (mu_m + A * PHM * X))
  
  
  
  # prevalence is current IC + IS
  
  
  deriv(S) <- -FOI*S - deltaT*S + mut + rc*Ic + rs*Is
  deriv(E) <- FOI*S - tau*E - deltaT*E
  deriv(Ic) <- -rc*Ic + c * tau * E - deltaT*Ic
  deriv(Is) <- -rs * Is + (1-c) * tau * E - deltaT*Is
  
  # Initial conditions
  initial(S) <- S_init
  initial(E) <- E_init
  initial(Ic) <- Ic_init
  initial(Is) <- Is_init
  
  M = interpolate(tvec_spline, m_in, 'spline')
  A = interpolate(tvec_spline, a_in, 'spline')
  N = interpolate(tvec_spline, n_in, 'spline')
  
  tvec_spline[] = user()
  m_in[] = user()
  a_in[] = user()
  n_in[] = user()
  
  dim(tvec_spline) = 365
  dim(m_in) = 365
  dim(a_in) = 365
  dim(n_in) = 365
  
  ## parameters
  deltaT = user()
  mut = user()
  rc = user()
  rs = user()
  tau = user()
  c = user()
  PMH = user()
  mu_m = user()
  PHM = user()
  epsilon = user()
  
  # Initial conditions
  S_init = user()
  E_init = user()
  Ic_init = user()
  Is_init = user()
  
})


calibration_optimizer = function(temps, rainfall, prop_cases_treated, prevalence_2019, incidence_2020,
                        prevalence_2020){
  
  temps[temps > 33.7] = 33.7
  temps[temps > 15.4] = 15.4
  
  predicted_prevalence = function(bR, rs, PMH, rc, tau, epsilon){
    
    
    m_in_vec = m_no_epsilon(temps,  mu_m, c.RD = rainfall, aR, bR)
    
    a_in_vec = a(temps)
    
    n_in_vec = a(temps)
    
    t = seq(1, 365, by =1/8)
    
    model = SEI_model$new(S_init = 1 - .01 - prevalence_2019,
                          E_init=.01,
                          Ic_init = prevalence_2019 * prop_cases_treated,
                          Is_init = prevalence_2019 * (1-prop_cases_treated),
                          deltaT = population_data_kenya_2020$death_rate / 365,
                          mut = population_data_kenya_2020$birth_rate / 365,
                          rc = rc,
                          rs = rs,
                          tau = tau,
                          epsilon = epsilon,
                          c = prop_cases_treated,
                          PHM = PHM,
                          PMH = PMH,
                          mu_m = mu_m,
                          m_in = m_in_vec,
                          a_in = a_in_vec,
                          n_in = n_in_vec,
                          tvec_spline = seq(1, 365))
    
    out = model$run(t)
    
    out=as.data.frame(out)
    
    return(mean(out$Ic + out$Is))
    
  }
  
  predicted_incidence = function(bR, rs, PMH, rc, tau, epsilon){
  
    
    m_in_vec = m_no_epsilon(temps, mu_m, c.RD = rainfall, aR, bR)
    
    a_in_vec = a(temps)
    
    n_in_vec = a(temps)
    
    t = seq(1, 365, by =1/8)
    
    model = SEI_model$new(S_init = 1 - .01 - prevalence_2019,
                          E_init=.01,
                          Ic_init = prevalence_2019 * prop_cases_treated,
                          Is_init = prevalence_2019 * (1-prop_cases_treated),
                          deltaT = population_data_kenya_2020$death_rate / 365,
                          mut = population_data_kenya_2020$birth_rate / 365,
                          rc = rc,
                          rs = rs,
                          tau = tau,
                          epsilon = epsilon,
                          c = prop_cases_treated,
                          PHM = PHM,
                          PMH = PMH,
                          mu_m = mu_m,
                          m_in = m_in_vec,
                          a_in = a_in_vec,
                          n_in = n_in_vec,
                          tvec_spline = seq(1, 365))
    
    out = model$run(t)
    
    out=as.data.frame(out)
    
    prev = out$Ic + out$Is
    
    m_vector = rep(m_in_vec, each = 8)[1:2913]
    a_vector = rep(a_in_vec, each = 8)[1:2913]
    n_vector = rep(n_in_vec, each = 8)[1:2913]
    
    return(mean(out$S * (1 - 2.71828182845904523536 ^ (-PMH * (m_vector * (a_vector ^ 2) * PHM * prev *
                                                                (2.71828182845904523536 ^
                                                                   (-mu_m * n_vector)))/
                                                        (mu_m + a_vector * PHM * prev)))))
    
  }
  
  
  optimizing_function = function(parameters){
    bR = parameters[1]
    rs = parameters[2]
    PMH = parameters[3]
    rc = parameters[4]
    tau = parameters[5]
    epsilon = parameters[6]
    ((predicted_prevalence(bR, rs, PMH, rc, tau, epsilon) - prevalence_2020) ^ 2)  + 
      (abs(predicted_incidence(bR, rs, PMH, rc, tau, epsilon) - incidence_2020))
  }
  
  optimized_result = optim(par = c(0, .1, .5, .1, .3, 50), optimizing_function, method="L-BFGS-B",
                           lower=c(-10000, 1/60, 0, 1/2000, 1/60, 1/100), upper = c(10000, 1/7, 1, 1/2, 1/5, 1000))
  
  return(list(predicted_prevalence(optimized_result$par[1], optimized_result$par[2],
                                   optimized_result$par[3], optimized_result$par[4],
                                   optimized_result$par[5], optimized_result$par[6]),
              predicted_incidence(optimized_result$par[1], optimized_result$par[2],
                                  optimized_result$par[3], optimized_result$par[4],
                                  optimized_result$par[5], optimized_result$par[6]),
              list(optimized_result$par)))
  
}

# obj = bR_optimizer(temps = temps, rainfall = rainfall,
#                    prop_cases_treated = kenya_2020$prop_cases[1],
#              prevalence_2019 = kenya_2019_prevalence$prevalence[[1]],
#              prevalence_2020 =  kenya_2020$prevalence[1])
# 
# 
# 

