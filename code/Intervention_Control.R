library(dplyr)
library(fitdistrplus)
library(mvtnorm)
library(BayesianTools)
library(postpack)
library(coda)
library(doParallel)


#### load functions and trial data ####

source("functions.R")

human.location <- read.csv("HumanLocation.csv")
cluster.info <- read.csv("cluster.info.csv")
coverage = read.csv("coverage.csv")

psi = mean(coverage$coverage)



# intervention data
incidence.intervention <- read.csv("intervention.epi.csv")
hlc.intervention <- read.csv("intervention.hlc.csv")
parity.intervention = read.csv("intervention.parity.csv")
prevalence_intervention = read.csv("prevalence_intervention.csv")


# remove one subject that had 0 followup time
incidence.intervention = incidence.intervention %>%
  filter(followup.time > 0)


# mean inside and outside proportions
prop_time_inside = mean(human.location$prop.inside)
prop_time_outside = mean(human.location$prop.outside)


#### sample sizes ####


ss_parity_intervention_control = parity.intervention %>%
  filter(Treatment.Allocation == "U") %>%
  summarize(number.tested = sum(number.tested)) %>%
  pull(number.tested)
ss_parity_intervention_treatment = parity.intervention %>%
  filter(Treatment.Allocation == "T") %>%
  summarize(number.tested = sum(number.tested)) %>%
  pull(number.tested)

# HLC sample size
ss_hlc_intervention_control = nrow(hlc.intervention %>%
                                     filter(Treatment.Allocation == "U"))
ss_hlc_intervention_treatment = nrow(hlc.intervention %>%
                                       filter(Treatment.Allocation == "T"))

# incidence- get the followup times for each person

ss_incidence_intervention_control= incidence.intervention %>%
  filter(Treatment.Allocation == "U") %>%
  pull(followup.time)
ss_incidence_intervention_treatment = incidence.intervention %>%
  filter(Treatment.Allocation == "T") %>%
  pull(followup.time)

# prevalence 

ss_prevalence_int_control = prevalence_intervention %>%
  filter(Treatment.Allocation == "U") %>%
  summarize(total_tests = sum(total_tests)) %>%
  pull(total_tests)
ss_prevalence_int_treatment = prevalence_intervention %>%
  filter(Treatment.Allocation == "T") %>%
  summarize(total_tests = sum(total_tests)) %>%
  pull(total_tests)



#### function that performs inference with a given data set ####


#### function that performs inference with a given data set ####


# unpack the data

parity.intervention.control.data = sum(parity.intervention$number.parous[
  parity.intervention$Treatment.Allocation == "U"])
parity.intervention.treatment.data = sum(parity.intervention$number.parous[
  parity.intervention$Treatment.Allocation == "T"])


prev.intervention.control.data = sum(prevalence_intervention$positive_tests[
  prevalence_intervention$Treatment.Allocation == "U"])
prev.intervention.treatment.data = sum(prevalence_intervention$positive_tests[
  prevalence_intervention$Treatment.Allocation == "T"])


hlc.inside.intervention.control.data = hlc.intervention$inside[
  hlc.intervention$Treatment.Allocation == "U"]
hlc.outside.intervention.control.data = hlc.intervention$outside[
  hlc.intervention$Treatment.Allocation == "U"]
hlc.inside.intervention.treatment.data = hlc.intervention$inside[
  hlc.intervention$Treatment.Allocation == "T"]
hlc.outside.intervention.treatment.data = hlc.intervention$outside[
  hlc.intervention$Treatment.Allocation == "T"]



infections.intervention.control.data = incidence.intervention$total.infections[
  incidence.intervention$Treatment.Allocation == "U"]
infections.intervention.treatment.data = incidence.intervention$total.infections[
  incidence.intervention$Treatment.Allocation == "T"]

#### define likelihood function ####
likelihood = function(params){
  
  
  # unpack the parameters and transform them back to their natural scale
  params=as.numeric(params)
  
  a_in_u= exp(params[1])
  a_out= exp(params[2])
  gu = sigmoid(params[3])
  lambda = exp(params[4])
  phi = exp(params[5])
  prop_indoor_u = sigmoid(params[6])
  b = sigmoid(params[7])
  c = sigmoid(params[8])
  r=exp(params[9])
  n = exp(params[10])
  q = sigmoid(params[11])
  
  
  # parity for control clusters- 0 coverage 
  parity.rate.control = parity(a_in_u, a_out, a_mult = 0, gu, prop_indoor_u,
                               g_mult = 0, p_mult = 1, psi = 0, q)
  
  # if the parity is 0, add a very small value to it to avoid an error
  if(parity.rate.control == 0){
    parity.rate.control = parity.rate.control + 1e-20
  }
  
  # Human biting rates indoor and outdoor for HLC- equations 11 and 12. Note psi = 0 
  H_in_u = HBR_iu(a_in_u, gu, lambda ,
                  prop_indoor_u, 
                  g_mult = 1, p_mult = 1, psi = 0)
  H_out_u = HBR_ou(a_out, gu, lambda ,
                   prop_indoor_u, 
                   g_mult = 1, p_mult = 0, psi = 0)
  
  if(H_in_u == 0){
    H_in_u  = H_in_u  + 1e-20
  }
  
  if(H_out_u == 0){
    H_out_u  = H_out_u  + 1e-20
  }
  
  # likelihood for parity
  parity.baseline.ll = sum(dbinom(x = parity.intervention.control.data,
                                  size = ss_parity_intervention_control,
                                  prob = parity.rate.control, log = T))
  
  # since HLC was only conducted for 50 min / hr, multiply the bite number by 5/6
  indoor.hlc.baseline.ll = sum(dnbinom(hlc.inside.intervention.control.data,
                                       mu = H_in_u * (5/6),
                                       size=phi, log = T))
  
  outdoor.hlc.baseline.ll = sum(dnbinom(hlc.outside.intervention.control.data,
                                        mu = H_out_u * (5/6), size=phi,
                                        log = T))
  
  # calculate the prevalence, again noting 0 coverage
  prevalence.control = prevalence_untreated(a_in_u, a_out, gu, lambda ,
                                            prop_indoor_u, b, c,
                                            r, n, a_mult = 0,
                                            g_mult = 0, p_mult = 1, psi=0, Qi = prop_time_inside,
                                            Qo = prop_time_outside)
  
  if(is.nan(prevalence.control)){
    prevalence.control = 1e-20
  }
  
  # find incidence
  incidence.control = incidence(prevalence.control, r)
  
  if(incidence.control == 0){
    incidence.control = 1e-20
  }
  
  if(prevalence.control <= 0){
    prevalence.control = 1e-20
  }
  
  inc.baseline.ll = sum(dpois(x = infections.intervention.control.data, 
                              lambda = (incidence.control +1e-10) *
                                ss_incidence_intervention_control / 365, log=T))
  
  prev.baseline.ll = sum(dbinom(x = prev.intervention.control.data,size=ss_prevalence_int_control,
                                prob=prevalence.control, log=T))
  
  
  return(
    (parity.baseline.ll + 
       inc.baseline.ll + 
       prev.baseline.ll + 
       indoor.hlc.baseline.ll + outdoor.hlc.baseline.ll
    )
  )
}


#### Run Bayesian tools ####

# prior
load("MultivariatePrior.RData")

# prior density
prior_density = function(params){
  dens_mv = dmvnorm(params, mu_mvn,
                    sigma_mvn, log = T)
  return(dens_mv)
}

# sampler function
prior_sampler = function(n=1){
  sample_mv = rmvnorm(n, mu_mvn, sigma_mvn)
  return(cbind(sample_mv[1], sample_mv[2],sample_mv[3], sample_mv[4],
               sample_mv[5], sample_mv[6],sample_mv[7], sample_mv[8], sample_mv[9], 
               sample_mv[10], sample_mv[11]))
}

# combine into Bayesian setup
bayesianSetup = createBayesianSetup(likelihood, prior_density, prior_sampler)

# initialize the convergence
convergence = F

#run for an initial 1000 iterations
output= runMCMC(bayesianSetup, settings=list(iterations = 1000),
                sampler="DEzs")
total_iterations_so_far = 1000


# until it converges, keep adding iterations
while(!convergence){
  
  # add 10000 iterations each time and keep track
  output= runMCMC(output, settings=list(iterations=10000))
  total_iterations_so_far = total_iterations_so_far + 10000
  
  # determine convergence
  gelman_output = gelmanDiagnostics(output)
  
  # once it's reached the convergence threshold, save the samples
  if(max(gelman_output$psrf[,2])< 1.1){
    
    convergence = T
    print("We have converged")
    
    # get a sample of posterior
    result = getSample(output, start = 1)
    result = tail(result, 100000)
    result = as.data.frame(result)
    
    # convert to natural units
    result_exp = exp(result)
    result_sigmoid = sigmoid(result)
    result = cbind(result_exp[,1:2], result_sigmoid[,3], result_exp[,4:5], result_sigmoid[,6:8],
                   result_exp[,9:10], result_sigmoid[,11])
    # give parameter names
    names(result) =c("a_in_u", "a_out", "gu", "lambda", "phi", "prop_indoor_u","b", "c", "r", "n",
                     "q")
    
    # save the result
    save(result, file="Intervention_Control.RData")
  } else{
    print(paste("Did not converge on", total_iterations_so_far, "iterations"))
  }
}
