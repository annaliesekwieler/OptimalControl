# load data
load("../data/kenya_data.RData")

# load climate data
load("../data/kenya_climate_data/kenya_climate_data.RData")

# load population data
load("../data/population_data_kenya_2020.RData")

# take out rasters that weren't able to be calibrated
to_keep = which(!is.na(kenya_data$bR))
kenya_data = kenya_data[to_keep,]
kenya_temp = kenya_temp[to_keep]
kenya_rainfall = kenya_rainfall[to_keep]

# source functions with model code
source("functions_SEISEI.R")

# take out the value for a single raster
prevalence_2019 = kenya_data[1, "prevalence_2019"]
prevalence_2020 = kenya_data[1, "prevalence_2020"]

temps = kenya_temp[[1]]
rainfall = kenya_rainfall[[1]]

temps[temps > 31.7] = 31.7
temps[temps < 13.4] = 13.4

temps = c(rollmean(temps, k = 50)[268:316],
          rollmean(temps, k = 50))

rainfall = c(rollmean(rainfall, k = 40)[288:326],
             rollmean(rainfall, k = 40))

temps[temps > 31.7] = 31.7
temps[temps < 13.4] = 13.4

#### calibration ####

# time period for simulation
# start after burn in period, then run
# see if we need a burn in here or is relatively stable
tvec_spline = seq(0, 365 * 6,by=1/12)
tvec = seq(0, 365 * 6,by=1/365)

# set values of epidemiological and entomological parameters
r = 1 / 157
mu_m = mu_m
rho = 1/14
Cmax = .8
PHM = .125
PMH = 0.5
tau = 1/10
aR = 0.5
bR = kenya_data$bR[1]


# vectors of a, n, and m. For a and m, under no coverage and coverage
m_n_in = rep(m(temps,  epsilon = 50, mu_m, c.RD = rainfall, aR, bR), times = 5) 
m_c_in = rep(m(temps,  epsilon = 50, mu_m * .9, c.RD = rainfall, aR, bR), times = 5) 
a_n_in = rep(a(temps), times = 5)
a_c_in = rep(a(temps)*.75, times = 5)
n_in = n(temps)

# initial control function
u.fun = function(t){
  rep(0,length(t))
}

# make up cost per case
cost_per_clinical_case = 300

parms = c()
# set min and max application rate
parms['umin'] = 0
parms['umax'] = 1/30

# set values of parameters for min and max application rate
umin = rep(parms['umin'],length(tvec))
umax = rep(parms['umax'],length(tvec))


# final conditions of adjoint variables
stateAdjoint = c(
  lambdaSh = 0,
  lambdaEh = 0,
  lambdaIh = 0,
  lambdaSm = 0,
  lambdaEm = 0,
  lambdaIm = 0,
  lambdaC = 0
  )

# tuning parameters for forward-backward sweep
maxit = 5000
tol = 1e-4
omega = 0.99

# track the objective function value and optimal control across iterations
J = rep(1e5,maxit)
rel.error = rep(0, maxit)
u.proposed = rep(NA,length(tvec))
u.outfbsm = rep(NA,length(tvec))
# umat.proposed = matrix(NA,maxit,length(tvec))
# umat.outfbsm = matrix(NA,maxit,length(tvec))



# iterate through the forward-backward sweep for a certain number of iterations
iteration = 1
u.proposed = rep(0, length(tvec))
u.outfbsm = rep(0, length(tvec))
# umat.proposed[iteration,] = rep(0, length(tvec))
# umat.outfbsm[iteration,] = rep(0, length(tvec))
iteration = 2

# start with a single raster

start.time = Sys.time()
raster = 1
repeat{
  # solve for state variables forward in time
  mod = SEI_model$new(
    S_init =1 - kenya_data$prevalence_2020[raster] - .01,
    E_init=.01,
    Ic_init = kenya_data$prevalence_2020[raster] * kenya_data$prop_cases[raster],
    Is_init = kenya_data$prevalence_2020[raster] * (1 - kenya_data$prop_cases[raster]),
    deltaT = population_data_kenya_2020$death_rate / 365,
    mut = population_data_kenya_2020$birth_rate / 365,
    rs = 1/224,
    epsilon = 50,
    relative_infectiousness = .1,
    c = kenya_data$prop_cases[raster],
    mu_m = mu_m * ,
    m_in = m_in_vec,
    a_in = a_in_vec,
    n_in = n_in_vec,
    tvec_spline = seq(1, 365*5))
  outState = mod$run(tvec)
  # outState = ode(
  #   y = c(states.2025/sum(states.2025),C=0),
  #   times = tvec, func = model_control, parms = parms, method = 'lsoda')
  
  # get functions approximating state variable solutions
  S0fun = approxfun(tvec, outState[,'S0'], rule = 2)
  I1fun = approxfun(tvec, outState[,'I1'], rule = 2)
  R1fun = approxfun(tvec, outState[,'R1'], rule = 2)
  S1fun = approxfun(tvec, outState[,'S1'], rule = 2)
  I2fun = approxfun(tvec, outState[,'I2'], rule = 2)
  R2fun = approxfun(tvec, outState[,'R2'], rule = 2)
  S2fun = approxfun(tvec, outState[,'S2'], rule = 2)
  I3fun = approxfun(tvec, outState[,'I3'], rule = 2)
  R3fun = approxfun(tvec, outState[,'R3'], rule = 2)
  S3fun = approxfun(tvec, outState[,'S3'], rule = 2)
  I4fun = approxfun(tvec, outState[,'I4'], rule = 2)
  R4fun = approxfun(tvec, outState[,'R4'], rule = 2)
  Cfun  = approxfun(tvec, outState[,'C'],  rule = 2)
  ufun  = approxfun(tvec, mod$contents()$u_in, rule = 2)
  
  # calculate objective function
  I1int = integrate(
    function(t){I1fun(t)},min(outState[,'t']),
    max(outState[,'t']),stop.on.error=F)$value
  I2int = integrate(
    function(t){I2fun(t)},min(outState[,'t']),
    max(outState[,'t']),stop.on.error=F)$value
  I3int = integrate(
    function(t){I3fun(t)},min(outState[,'t']),
    max(outState[,'t']),stop.on.error=F)$value
  I4int = integrate(
    function(t){I4fun(t)},min(outState[,'t']),
    max(outState[,'t']),stop.on.error=F)$value
  uint = integrate(
    function(t){ufun(t)},min(outState[,'t']),
    max(outState[,'t']),stop.on.error=F)$value
  J[iteration] =
    parms['w_I1'] * I1int * parms['gamma'] / 365 +
    parms['w_I2'] * I2int * parms['gamma'] / 365 +
    parms['w_I3'] * I3int * parms['gamma'] / 365 +
    parms['w_I4'] * I4int * parms['gamma'] / 365 +
    parms['w_C'] * uint
  
  # solve for adjoint variables backward in time
  mod = model_adjoint_bkwd$new(
    beta_mag = exp(beta.opt$minimum),
    w_I1 = parms['w_I1'],
    w_I2 = parms['w_I2'],
    w_I3 = parms['w_I3'],
    w_I4 = parms['w_I4'],
    tvec_spline = tvec_spline,
    tvec = tvec,
    beta_c_in = rev(beta_c_in),
    beta_n_in = rev(beta_n_in),
    delta_0_in = rev(delta.0.fun(tvec_spline)),
    delta_1_in = rev(delta.1.fun(tvec_spline)),
    delta_2_in = rev(delta.2.fun(tvec_spline)),
    delta_3_in = rev(delta.3.fun(tvec_spline)),
    delta_4_in = rev(delta.4.fun(tvec_spline)),
    S0_in = rev(S0fun(tvec)),
    I1_in = rev(I1fun(tvec)),
    R1_in = rev(R1fun(tvec)),
    S1_in = rev(S1fun(tvec)),
    I2_in = rev(I2fun(tvec)),
    R2_in = rev(R2fun(tvec)),
    S2_in = rev(S2fun(tvec)),
    I3_in = rev(I3fun(tvec)),
    R3_in = rev(R3fun(tvec)),
    S3_in = rev(S3fun(tvec)),
    I4_in = rev(I4fun(tvec)),
    R4_in = rev(R4fun(tvec)),
    C_in  = rev(Cfun(tvec)),
    u_in  = rev(ufun(tvec)),
    lambdaS0_init = stateAdjoint[1],
    lambdaI1_init = stateAdjoint[2],
    lambdaR1_init = stateAdjoint[3],
    lambdaS1_init = stateAdjoint[4],
    lambdaI2_init = stateAdjoint[5],
    lambdaR2_init = stateAdjoint[6],
    lambdaS2_init = stateAdjoint[7],
    lambdaI3_init = stateAdjoint[8],
    lambdaR3_init = stateAdjoint[9],
    lambdaS3_init = stateAdjoint[10],
    lambdaI4_init = stateAdjoint[11],
    lambdaR4_init = stateAdjoint[12],
    lambdaC_init  = stateAdjoint[13])
  outAdjoint = mod$run(tvec)
  
  # beta_n_in = spline(tvec_spline,beta_n_in,method='natural')
  # beta_n_fun = approxfun(beta_n_in$x,beta_n_in$y)
  # beta_c_in = spline(tvec_spline,beta_c_in,method='natural')
  # beta_c_fun = approxfun(beta_c_in$x,beta_c_in$y)
  # u.fun = ufun
  # outAdjoint = ode(
  #   y = stateAdjoint, times = rev(tvec), func = adjointEqns, parms = parms, method = 'lsoda')
  
  # get functions approximating adjoint variable solutions
  lambdaS0fun = approxfun(rev(tvec), outAdjoint[,'lambdaS0'], rule = 2)
  lambdaI1fun = approxfun(rev(tvec), outAdjoint[,'lambdaI1'], rule = 2)  
  lambdaR1fun = approxfun(rev(tvec), outAdjoint[,'lambdaR1'], rule = 2)
  lambdaS1fun = approxfun(rev(tvec), outAdjoint[,'lambdaS1'], rule = 2)
  lambdaI2fun = approxfun(rev(tvec), outAdjoint[,'lambdaI2'], rule = 2)
  lambdaR2fun = approxfun(rev(tvec), outAdjoint[,'lambdaR2'], rule = 2)
  lambdaS2fun = approxfun(rev(tvec), outAdjoint[,'lambdaS2'], rule = 2)
  lambdaI3fun = approxfun(rev(tvec), outAdjoint[,'lambdaI3'], rule = 2)
  lambdaR3fun = approxfun(rev(tvec), outAdjoint[,'lambdaR3'], rule = 2)
  lambdaS3fun = approxfun(rev(tvec), outAdjoint[,'lambdaS3'], rule = 2)
  lambdaI4fun = approxfun(rev(tvec), outAdjoint[,'lambdaI4'], rule = 2)
  lambdaR4fun = approxfun(rev(tvec), outAdjoint[,'lambdaR4'], rule = 2)
  lambdaCfun  = approxfun(rev(tvec), outAdjoint[,'lambdaC'],  rule = 2)
  
  # get function approximating control that maximizes Hamiltonian
  psi = dHdu(tvec, parms)
  sign.psi = sign(unlist(psi))
  temp.u = ifelse(sign.psi==-1,parms['umax'],parms['umin'])
  u.fun = approxfun(tvec, temp.u, rule=2)
  
  # umat.outfbsm[iteration,] = temp.u
  u.outfbsm.old = u.outfbsm
  u.outfbsm = temp.u
  
  if((iteration == maxit)){
    # umat.proposed[iteration,] = temp.u
    u.proposed = temp.u
  } else {
    # umat.proposed[iteration,] = omega * umat.proposed[iteration-1,] + (1-omega) * temp.u
    u.proposed = omega * u.proposed + (1-omega) * temp.u
    u.fun = approxfun(tvec, u.proposed, rule=2)
  }
  
  rel.error[iteration] =
    sum(abs(u.outfbsm - u.outfbsm.old)) /
    max(1e-6,sum(u.outfbsm))
  # print(rel.error[iteration])
  
  # exit once the control stops changing sufficiently
  if(iteration == maxit){
    break
  }
  
  # update iteration
  iteration = iteration + 1
}

# report how long it took
end.time = Sys.time()
end.time - start.time

# simulate with the optimal control and calculated associated J
u.opt = u.fun = approxfun(tvec, u.outfbsm, rule=2)
# u.opt = u.fun = approxfun(tvec, umat.outfbsm[iteration-1,], rule=2)
mod = model_control_cuminc$new(
  beta_mag = exp(beta.opt$minimum),
  tvec_spline = tvec_spline,
  tvec = tvec,
  beta_c_in = beta_c_in,
  beta_n_in = beta_n_in,
  delta_0_in = delta.0.fun(tvec_spline),
  delta_1_in = delta.1.fun(tvec_spline),
  delta_2_in = delta.2.fun(tvec_spline),
  delta_3_in = delta.3.fun(tvec_spline),
  delta_4_in = delta.4.fun(tvec_spline),
  mu_in = mu.fun(tvec_spline),
  u_in = u.fun(tvec),
  S0_init = states.2025[1],
  I1_init = states.2025[2],
  R1_init = states.2025[3],
  S1_init = states.2025[4],
  I2_init = states.2025[5],
  R2_init = states.2025[6],
  S2_init = states.2025[7],
  I3_init = states.2025[8],
  R3_init = states.2025[9],
  S3_init = states.2025[10],
  I4_init = states.2025[11],
  R4_init = states.2025[12])
outState.control = mod$run(tvec)
# outState.control = ode(
#   y = c(states.2025/sum(states.2025),C=0,cum.inc=0),
#   times = tvec, func = model_control.cum.inc, parms = parms, method = 'lsoda')
S0fun = approxfun(tvec, outState.control[,'S0'], rule = 2)
I1fun = approxfun(tvec, outState.control[,'I1'], rule = 2)
R1fun = approxfun(tvec, outState.control[,'R1'], rule = 2)
S1fun = approxfun(tvec, outState.control[,'S1'], rule = 2)
I2fun = approxfun(tvec, outState.control[,'I2'], rule = 2)
R2fun = approxfun(tvec, outState.control[,'R2'], rule = 2)
S2fun = approxfun(tvec, outState.control[,'S2'], rule = 2)
I3fun = approxfun(tvec, outState.control[,'I3'], rule = 2)
R3fun = approxfun(tvec, outState.control[,'R3'], rule = 2)
S3fun = approxfun(tvec, outState.control[,'S3'], rule = 2)
I4fun = approxfun(tvec, outState.control[,'I4'], rule = 2)
R4fun = approxfun(tvec, outState.control[,'R4'], rule = 2)
Cfun  = approxfun(tvec, outState.control[,'C'],  rule = 2)
I1int = integrate(
  function(t){I1fun(t)},min(outState.control[,'t']),
  max(outState.control[,'t']),stop.on.error=F)$value
I2int = integrate(
  function(t){I2fun(t)},min(outState.control[,'t']),
  max(outState.control[,'t']),stop.on.error=F)$value
I3int = integrate(
  function(t){I3fun(t)},min(outState.control[,'t']),
  max(outState.control[,'t']),stop.on.error=F)$value
I4int = integrate(
  function(t){I4fun(t)},min(outState.control[,'t']),
  max(outState.control[,'t']),stop.on.error=F)$value
uint = integrate(
  function(t){u.fun(t)},min(outState.control[,'t']),
  max(outState.control[,'t']),stop.on.error=F)$value
J.control = unname(
  parms['w_I1'] * I1int * parms['gamma'] / 365 +
    parms['w_I2'] * I2int * parms['gamma'] / 365 +
    parms['w_I3'] * I3int * parms['gamma'] / 365 +
    parms['w_I4'] * I4int * parms['gamma'] / 365 +
    parms['w_C'] * uint)
J.control.control = unname(parms['w_C'] * uint)
J.control.dengue = J.control - J.control.control

# simulate without the optimal control
u.fun = function(t){rep(0,length(t))}
mod = model_control_cuminc$new(
  beta_mag = exp(beta.opt$minimum),
  tvec_spline = tvec_spline,
  tvec = tvec,
  beta_c_in = beta_c_in,
  beta_n_in = beta_n_in,
  delta_0_in = delta.0.fun(tvec_spline),
  delta_1_in = delta.1.fun(tvec_spline),
  delta_2_in = delta.2.fun(tvec_spline),
  delta_3_in = delta.3.fun(tvec_spline),
  delta_4_in = delta.4.fun(tvec_spline),
  mu_in = mu.fun(tvec_spline),
  u_in = u.fun(tvec),
  S0_init = states.2025[1],
  I1_init = states.2025[2],
  R1_init = states.2025[3],
  S1_init = states.2025[4],
  I2_init = states.2025[5],
  R2_init = states.2025[6],
  S2_init = states.2025[7],
  I3_init = states.2025[8],
  R3_init = states.2025[9],
  S3_init = states.2025[10],
  I4_init = states.2025[11],
  R4_init = states.2025[12])
outState.baseline = mod$run(tvec)
# outState.baseline = ode(
#   y = c(states.2025/sum(states.2025),C=0,cum.inc=0),
#   times = tvec, func = model_control.cum.inc, parms = parms, method = 'lsoda')
S0fun = approxfun(tvec, outState.baseline[,'S0'], rule = 2)
I1fun = approxfun(tvec, outState.baseline[,'I1'], rule = 2)
R1fun = approxfun(tvec, outState.baseline[,'R1'], rule = 2)
S1fun = approxfun(tvec, outState.baseline[,'S1'], rule = 2)
I2fun = approxfun(tvec, outState.baseline[,'I2'], rule = 2)
R2fun = approxfun(tvec, outState.baseline[,'R2'], rule = 2)
S2fun = approxfun(tvec, outState.baseline[,'S2'], rule = 2)
I3fun = approxfun(tvec, outState.baseline[,'I3'], rule = 2)
R3fun = approxfun(tvec, outState.baseline[,'R3'], rule = 2)
S3fun = approxfun(tvec, outState.baseline[,'S3'], rule = 2)
I4fun = approxfun(tvec, outState.baseline[,'I4'], rule = 2)
R4fun = approxfun(tvec, outState.baseline[,'R4'], rule = 2)
Cfun  = approxfun(tvec, outState.baseline[,'C'],  rule = 2)
I1int = integrate(
  function(t){I1fun(t)},min(outState.baseline[,'t']),
  max(outState.baseline[,'t']),
  stop.on.error=F)$value
I2int = integrate(
  function(t){I2fun(t)},min(outState.baseline[,'t']),
  max(outState.baseline[,'t']),
  stop.on.error=F)$value
I3int = integrate(
  function(t){I3fun(t)},min(outState.baseline[,'t']),
  max(outState.baseline[,'t']),
  stop.on.error=F)$value
I4int = integrate(
  function(t){I4fun(t)},min(outState.baseline[,'t']),
  max(outState.baseline[,'t']),
  stop.on.error=F)$value
uint = integrate(
  function(t){u.fun(t)},min(outState.baseline[,'t']),
  max(outState.baseline[,'t']),
  stop.on.error=F)$value
J.baseline = unname(
  parms['w_I1'] * I1int * parms['gamma'] / 365 +
    parms['w_I2'] * I2int * parms['gamma'] / 365 +
    parms['w_I3'] * I3int * parms['gamma'] / 365 +
    parms['w_I4'] * I4int * parms['gamma'] / 365 +
    parms['w_C'] * uint)
J.baseline.control = unname(parms['w_C'] * uint)
J.baseline.dengue = J.baseline - J.baseline.control

# statistics about control
control.yearly = aggregate(u.opt(tvec)>0,by=list(floor(tvec)),FUN=mean)[1:10,2]
coverage.yearly = aggregate(outState.control[,'C'],by=list(floor(tvec)),FUN=mean)[1:10,2]
infections.averted.absolute = 
  max(outState.baseline[,'cuminc'],na.rm=T) -
  max(outState.control[,'cuminc'],na.rm=T)
infections.averted.pctdiff = 
  (max(outState.baseline[,'cuminc'],na.rm=T) -
     max(outState.control[,'cuminc'],na.rm=T)) /
  max(outState.baseline[,'cuminc'],na.rm=T)
money.spent = J.control.control - J.baseline.control
money.saved = J.baseline.dengue - J.control.dengue

# save to file
toSave = c('rel.error','J','J.control','J.baseline',
           'J.control.dengue','J.control.control',
           'J.baseline.dengue','J.baseline.control',
           'infections.averted.absolute',
           'infections.averted.pctdiff',
           'control.yearly','control.propn',
           'money.spent','money.saved',
           'parms','states.2025','tvec',
           'outState.control','outState.baseline','u.opt')
save(list=toSave,file=paste('optimalcontrol/optimalcontrol_',cc.ind,'.RData',sep=''))