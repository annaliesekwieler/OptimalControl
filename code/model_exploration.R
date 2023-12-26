# functions that run epi and mosquito models
source("functions_SEISEI.R")

# set the temperature, rainfall, and input data 
load("../data/kenya_data.RData")
load("../data/kenya_climate_data/kenya_climate_data.RData")
i=1
temps = kenya_temp[[i]]$daily_temp
temps[temps > 31.7] = 31.7
temps[temps < 15.39] = 15.39
temps = c(rollmean(temps, k = 50)[268:316],
          rollmean(temps, k = 50))
rainfall = c(rollmean(rainfall, k = 40)[288:326],
             rollmean(rainfall, k = 40))
temps[temps > 31.7] = 31.7
temps[temps < 15.39] = 15.39
rainfall = kenya_rainfall[[i]]
prevalence_2019 = kenya_data[i,"prevalence_2019"]
incidence_2020 = kenya_data[i,"incidence_2020"]/365
prevalence_2020 = kenya_data[i, "prevalence_2020"]

m_in_vec_one_year = rep(m(temps,  epsilon = 50, mu_m, c.RD = rainfall, aR, bR)) 

m_in_vec = rep(m_in_vec_one_year, times = 5)

a_in_vec_one_year = a(temps)
a_in_vec = rep(a_in_vec_one_year, times = 5)

n_in_vec_one_year = n(temps)
n_in_vec = rep(n_in_vec_one_year, times = 5)

t = seq(1, 365 * 5, by =1/8)

model = SEI_model$new(Sh_init = 1 - .01 - prevalence_2019,
                      Eh_init=.01,
                      Ih_init = prevalence_2019,
                      Sm_init = 1 - .01 - prevalence_2019,
                      Em_init=.01,
                      Im_init = prevalence_2019,
                      C_init = 0,
                      r = 1 / 157,
                      mu_m = mu_m,
                      m_in = m_in_vec,
                      a_in = a_in_vec,
                      n_in = n_in_vec,
                      u_in = rep(0, 1825),
                      tvec_spline = seq(1, 365*5),
                      rho = 1/14,
                      Cmax = .8,
                      PHM = .125,
                      PMH = 0.5,
                      tau = 1/10)

out = model$run(t)

out=as.data.frame(out)[11675:nrow(out),]