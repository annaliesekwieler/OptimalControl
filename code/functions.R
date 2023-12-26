
library(zoo)

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
PMH = 0.5

# incubation rate in humans
tau = 1/10
rc = 1/5


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
  
  PMH = .5
  PHM = 0.125
  rc = .2
  tau = 1/10
  
  A_modified = A * (1-C) + A * A_sr * C
  mu_modified = mu_m * (1-C) + mu_m * G_sr * C
  
  # exp
  FOI = 1 - 2.72 ^ (-PMH * (M * (A_modified ^ 2) * PHM * I *
                                                (2.71828182845904523536 ^ (-mu_modified * N)))/
                                        (mu_modified + A_modified * PHM * I))
  
  
  
  # prevalence is current IC + IS
  
  
  deriv(S) <- -FOI*S - deltaT*S + mut*S + r*I
  deriv(E) <- FOI*S - tau*E - deltaT*E
  deriv(I) <- -r*I + tau * E - deltaT*I
  deriv(C) <- U * (1 - C / Cmax) - rho * C
  
  
  # Initial conditions
  initial(S) <- S_init
  initial(E) <- E_init
  initial(I) <- I_init
  initial(C) <- C_init
  
  M = interpolate(tvec_spline, m_in, 'spline')
  A = interpolate(tvec_spline, a_in, 'spline')
  N = interpolate(tvec_spline, n_in, 'spline')
  U = interpolate(tvec_spline, u_in, 'spline')
  
  tvec_spline[] = user()
  m_in[] = user()
  a_in[] = user()
  n_in[] = user()
  u_in[] = user()
  
  dim(tvec_spline) = 1825
  dim(m_in) = 1825
  dim(a_in) = 1825
  dim(n_in) = 1825
  dim(u_in) = 1825
  
  ## parameters
  deltaT = user()
  mut = user()
  r = user()
  c = user()
  mu_m = user()
  epsilon = user()
  G_sr = user()
  A_sr = user()
  Cmax = user()
  rho = user()
  
  # Initial conditions
  S_init = user()
  E_init = user()
  I_init = user()
  C_init = user()
  
})

epsilon = 50


PHM = .125
PMH = 0.5

rs = 1/224

relative_infectiousness = .1


calibration_optimizer = function(temps, rainfall, prop_cases_treated, prevalence_2019, incidence_2020,
                        prevalence_2020, default_bR = 0, x, y){
  
  temps[temps > 31.7] = 31.7
  temps[temps < 13.4] = 13.4
  
  temps = c(rollmean(temps, k = 50)[268:316],
                  rollmean(temps, k = 50))
  
  rainfall = c(rollmean(rainfall, k = 40)[288:326],
            rollmean(rainfall, k = 40))
  
  temps[temps > 31.7] = 31.7
  temps[temps < 13.4] = 13.4
  
  predicted_prevalence = function(bR){
    
    m_in_vec_one_year = rep(m(temps,  epsilon, mu_m, c.RD = rainfall, aR, bR)) 
    
    m_in_vec = rep(m_in_vec_one_year, times = 5)
    
    a_in_vec_one_year = a(temps)
    a_in_vec = rep(a_in_vec_one_year, times = 5)
    
    n_in_vec_one_year = n(temps)
    n_in_vec = rep(n_in_vec_one_year, times = 5)
    
    t = seq(1, 365 * 5, by =1/8)

    model = SEI_model$new(S_init = 1 - .01 - prevalence_2019,
                          E_init=.01,
                          I_init = prevalence_2019,
                          C_init = 0,
                          deltaT = population_data_kenya_2020$death_rate / 365,
                          mut = population_data_kenya_2020$birth_rate / 365,
                          rs = rs,
                          epsilon = epsilon,
                          relative_infectiousness = relative_infectiousness,
                          c = prop_cases_treated,
                          mu_m = mu_m,
                          m_in = m_in_vec,
                          a_in = a_in_vec,
                          n_in = n_in_vec,
                          u_in = rep(0, 1825),
                          tvec_spline = seq(1, 365*5),
                          A_sr = 1,
                          G_sr = 1, 
                          rho = 1/180,
                          Cmax = .8)
    
    out = model$run(t)
    
    out=as.data.frame(out)[11675:nrow(out),]
    
    return(mean(out$Ic + out$Is))
    
  }
  
  predicted_incidence = function(bR){
  
    
    m_in_vec_one_year = rep(m(temps,  epsilon, mu_m, c.RD = rainfall, aR, bR)) 
    
    m_in_vec = rep(m_in_vec_one_year, times = 5)
    
    a_in_vec_one_year = a(temps)
    a_in_vec = rep(a_in_vec_one_year, times = 5)
    
    n_in_vec_one_year = n(temps)
    n_in_vec = rep(n_in_vec_one_year, times = 5)
    
    t = seq(1, 365 * 5, by =1/8)

    
    
    model = SEI_model$new(S_init = 1 - .01 - prevalence_2019,
                          E_init=.01,
                          Ic_init = prevalence_2019 * prop_cases_treated,
                          Is_init = prevalence_2019 * (1-prop_cases_treated),
                          deltaT = population_data_kenya_2020$death_rate / 365,
                          mut = population_data_kenya_2020$birth_rate / 365,
                          rs = rs,
                          epsilon = epsilon,
                          c = prop_cases_treated,
                          relative_infectiousness = 0.5,
                          mu_m = mu_m,
                          m_in = m_in_vec,
                          a_in = a_in_vec,
                          n_in = n_in_vec,
                          tvec_spline = seq(1, 365*5))
    
    out = model$run(t)
    
    out=as.data.frame(out)[11674:nrow(out),]
    
    prev = out$Ic + out$Is
    
    m_vector = rep(m_in_vec_one_year, each = 8)
    a_vector = rep(a_in_vec_one_year, each = 8)
    n_vector = rep(n_in_vec_one_year, each = 8)
    
    eir = (m_vector * (a_vector ^ 2) * PHM * prev * exp(-mu_m*n_vector)) / 
      (mu_m + a_vector * PHM * prev)
    
    foi = 1 - exp(-PMH * eir)
    
    return(mean(out$S * foi))
    
  }
  
  
  optimizing_function = function(parameters){
    bR = parameters[1]
  #  rs = parameters[2]
     abs(predicted_prevalence(bR) - prevalence_2020)  #  + 
   #   ((predicted_incidence(bR) - incidence_2020)^2)
  }
  
  optimized_result = optim(par = c( 
  default_bR
  ), optimizing_function, method="L-BFGS-B",
                            lower=c(
                              -10000
                              ), upper = c(
                                 10000
                                ),
                           control=list(factr=1e-32, maxit=100000, pgtol=c(1e-16)))
  
  return(list(predicted_prevalence(optimized_result$par[1]),
              predicted_incidence(optimized_result$par[1]),
              list(optimized_result$par), optimized_result$convergence, x, y))
  
}


#### Predicted prevalence and incidence with effects of SRs ####

# load SR effects
load("../data/posterior_samples.RData")
mean(result$a_mult)
mean(result$g_mult)

