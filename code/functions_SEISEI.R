
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

# adult mosquito mortality
g = function(temp){
  1 - (-.000828 * (temp^2) + .0367 * temp + .522)
}

#### set default parameter values #####

# eggs laid per gonotrophic cycle
epsilon = 50

# epidemiological parameters

# transmission probabilities
b = 0.125
c = 0.5

# cumulative rainfall parameters
aR = 0.5

#### function that optimizes bR  for a particular raster ####

# function of the temperature, rainfall, and proportion of cases treated
# returns the bR that is optimal for that raster

# requires all other parameters to be specified globally


SEI_model = odin::odin({
  
  # Derivatives
  
  A_modified = A * (1-C) + A * .75 * C
  mu_modified = G * (1-C) + G * .9 * C

  # prevalence is current IC + IS
  
  
  deriv(Sh) <- -A_modified * c * Im * Sh + r*Ih
  deriv(Eh) <- A_modified * c * Im * Sh - tau * Eh
  deriv(Ih) <- tau * Eh - r * Ih
  
  deriv(Sm) = M * mu_modified - A_modified * b * Ih * Sm - mu_modified * Sm
  deriv(Em) =  A_modified * b * Ih * Sm - (1/N)*Em - mu_modified * Em
  deriv(Im) = (1/N)*Em - mu_modified * Im
  
  deriv(C) = U * (1 - C / Cmax) - rho * C
  
  
  # Initial conditions
  initial(Sh) <- Sh_init
  initial(Eh) <- Eh_init
  initial(Ih) <- Ih_init
  
  initial(Sm) <- Sm_init
  initial(Em) <- Em_init
  initial(Im) <- Im_init
  
  initial(C) <- C_init
  
  M = interpolate(tvec_spline, m_in, 'spline')
  A = interpolate(tvec_spline, a_in, 'spline')
  N = interpolate(tvec_spline, n_in, 'spline')
  G = interpolate(tvec_spline, g_in, 'spline')
  U = interpolate(tvec_spline, u_in, 'spline')
  
  tvec_spline[] = user()
  m_in[] = user()
  a_in[] = user()
  n_in[] = user()
  g_in[] = user()
  u_in[] = user()
  
  dim(tvec_spline) = 1825
  dim(m_in) = 1825
  dim(a_in) = 1825
  dim(n_in) = 1825
  dim(g_in) = 1825
  dim(u_in) = 1825
  
  ## parameters
  r = user()
  Cmax = user()
  rho = user()
  tau = user()
  b = user()
  c = user()
  phi_a = user()
  phi_g = user()
  
  # Initial conditions
  Sh_init = user()
  Eh_init = user()
  Ih_init = user()
  Sm_init = user()
  Em_init = user()
  Im_init = user()
  C_init = user()
  
})

run_SEISEI = function(temps, rainfall, prop_cases_treated, prevalence_2019,
                      incidence_2020, prevalence_2020,
                      u_in = rep(0, 1825), rho = 1/180, Cmax = .8,
                      b = .125, c = 0.5, tau = 1/10,
                      phi_a = 0.75, phi_g = 0.1,
                      # start times to get result (default is last year)
                      start = 11675, bR = 0) {
  
  # use rolling mean on temperature
  temps = c(rollmean(temps, k = 50)[268:316],
            rollmean(temps, k = 50))
  rainfall = c(rollmean(rainfall, k = 40)[288:326],
               rollmean(rainfall, k = 40))
  # replace values that are above or below the thresholds
  temps[temps > 31.7] = 31.7
  temps[temps < 15.39] = 15.39
  
  # calculate m, g, a, and n
  
  r = (1-prop_cases_treated) / ( (prevalence_2020 / incidence_2020) - 
                                   (1 / 7))

  g_in_vec = rep(g(temps), times = 5)
  m_in_vec = m(temps,  epsilon = 50, mu_m = g_in_vec, c.RD = rainfall, aR, bR)
  a_in_vec = rep(a(temps), times = 5)
  n_in_vec = rep(n(temps), times = 5)
 
  
  
  t = seq(1, 365 * 5, by =1/8)
  
  model = SEI_model$new(Sh_init = 1 - .01 - prevalence_2019,
                        Eh_init=.01,
                        Ih_init = prevalence_2019,
                        Sm_init = 1 - .01 - prevalence_2019,
                        Em_init=.01,
                        Im_init = prevalence_2019,
                        C_init = 0,
                        r = r,
                        g_in = g_in_vec,
                        m_in = m_in_vec,
                        a_in = a_in_vec,
                        n_in = n_in_vec,
                        u_in = rep(0, 1825),
                        tvec_spline = seq(1, 365*5),
                        rho = rho,
                        Cmax = Cmax,
                        b = b,
                        c = c,
                        tau = tau,
                        phi_a = phi_a,
                        phi_g = phi_g)
  
  out = model$run(t)
  out=as.data.frame(out)[start:nrow(out),]
  
  return(out)
  
}

calibration_optimizer = function(temps, rainfall, prop_cases_treated, prevalence_2019, incidence_2020,
                                 prevalence_2020, default_bR = 0, x, y){
  
  predicted_prevalence = function(bR){
    
    model_res = run_SEISEI(temps = temps, rainfall = rainfall,
                           prop_cases_treated = prop_cases_treated,
                           prevalence_2019,
                          incidence_2020, prevalence_2020, bR = bR)
  
    return(mean(model_res$Ih))
    
  }
  
  
  optimizing_function = function(parameters){
    bR = parameters[1]
    abs(predicted_prevalence(bR) - prevalence_2020)  
  }
  
  optimized_result = optim(par = c( 
    20
  ), optimizing_function, method="L-BFGS-B",
  lower=c(
    -10000
  ), upper = c(
    10000
  ),
  control=list(factr=1e-32, maxit=100000, pgtol=c(1e-16)))
  
  return(list(predicted_prevalence(optimized_result$par[1]),
              list(optimized_result$par), optimized_result$convergence, x, y))
  
}
