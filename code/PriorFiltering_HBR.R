#### intro ####

# This script takes the results from PriorFiltering_CRC (which are 
# sets of parameters generated according to marginal priors and 
# fall within a reasonable incidence and prevalence range)

# This script converts these samples into the log/inverse sigmoid scale


#### load CRC results ####

load("../results/PriorFiltering_CRC_HBR.RData")

new_prior_df = new_prior_df[complete.cases(new_prior_df),]

# add in qu
new_prior_df$q = rbeta(nrow(new_prior_df), 1, 1)

source("functions.R")

# convert to log and inverse sigmoid scales
prior_df_log = log(new_prior_df)
prior_df_inv = inverse_sigmoid(new_prior_df)

# glue together columns 
prior_df = cbind(prior_df_log[,1:2], prior_df_inv[,3],
                 prior_df_log[,4:5], prior_df_inv[,6:8],
                 prior_df_log[,9:10], prior_df_inv[,11])
names(prior_df) = names(new_prior_df)

# get the mvn distribution
mu_mvn = colMeans(prior_df)
sigma_mvn = cov(prior_df)

save(mu_mvn, sigma_mvn, file = "../results/MultivariatePrior_HBR.RData")



#### Prior for running on only intervention or on baseline and intervention together ####

result = new_prior_df 


result$a_mult = runif(nrow(result), .5, 1.5)
result$g_mult = runif(nrow(result), .5, 1.5)
result$p_mult = runif(nrow(result), .5, 1.5)

result_log = log(result)
result_inv_sig = inverse_sigmoid(result)


result_transformed= cbind(result_log[,1:2], result_inv_sig[,3],
                          result_log[,4:5], result_inv_sig[,6:8],
                          result_log[,9:10], result_inv_sig[11],
                          result_log[,12:14])

mu_mvn_baseline = colMeans(result_transformed)
sigma_mvn_baseline = var(result_transformed)
save(mu_mvn_baseline, sigma_mvn_baseline, file="../results/All_Prior_HBR.RData")


#### separate indoor and outdoor ####

new_prior_df = new_prior_df %>%
  rename(phi_indoor = phi)

new_prior_df$phi_outdoor = rgamma(nrow(new_prior_df), 1, 1)

new_prior_df_log = log(new_prior_df)
new_prior_df_inv_sig = inverse_sigmoid(new_prior_df)


new_prior_df_transformed= cbind(new_prior_df_log[,1:2], new_prior_df_inv_sig[,3],
                                new_prior_df_log[,4:5], new_prior_df_inv_sig[,6:8],
                                new_prior_df_log[,9:10], new_prior_df_inv_sig[11],
                                new_prior_df_log[,12])

mu_mvn_baseline = colMeans(new_prior_df_transformed)
sigma_mvn_baseline = var(new_prior_df_transformed)
save(mu_mvn_baseline, sigma_mvn_baseline, file="../results/MultivariatePrior_IndoorOutdoor.RData")
