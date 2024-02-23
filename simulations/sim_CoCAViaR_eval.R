library(tidyverse)
library(dplyr)
library(tidyr)
library(reshape2)
library(SystemicRisk)


# Load simulation results
res_df_MC <- readRDS(file = "simulations/data/sim_CoCAViaR.rds")

# Set options
model_choice <- "CoCAViaR_SAV_diag"
diag_entries <- which(MTS::Vech(diag(6))==1)


# Estimated Covariance Matrix and standard deviations (diagonal covariance entries only) as data frames
res_MC_cov <- res_df_MC %>%
  filter(model==model_choice, type %in% c("cov_est_asy")) %>%
  group_by(theta_index, model, n, alpha, beta, risk_measure, init_method, bandwidth_factor, type) %>%
  summarize(theta_cov_mean = mean(value, na.rm = TRUE),
            theta_cov_median = median(value, na.rm = TRUE))

res_MC_sd <- res_MC_cov %>%
  group_by(model, n, alpha, beta, risk_measure, init_method, bandwidth_factor) %>%
  filter(theta_index %in% diag_entries) %>%
  mutate(theta_index = 1:length(theta_index)) %>% # Reset theta index after selecting diagonal entries only
  group_by(theta_index) %>%
  mutate(sd_asy_mean = sqrt(theta_cov_mean),
         sd_asy_median = sqrt(theta_cov_median))

# True parameter values
res_MC_true <- res_df_MC %>%
  filter(model==model_choice, type=="true_value") %>%
  group_by(theta_index, model, n, alpha, beta, risk_measure, init_method, bandwidth_factor, type) %>%
  summarize(theta_true = mean(value, na.rm = TRUE))

# Extract the mean parameter estimates
res_MC_est <- res_df_MC %>%
  dplyr::filter(model==model_choice, type=="param_est") %>%
  group_by(theta_index, model, n, alpha, beta, risk_measure, init_method, bandwidth_factor, type) %>%
  summarize(theta_mean = mean(value),
            theta_median = median(value),
            theta_sd_emp = sd(value)) %>%
  ungroup()


# Join data frames together
res_MC <- res_MC_sd %>% ungroup %>% dplyr::select(-type) %>%
  left_join(res_MC_est %>% ungroup %>% dplyr::select(-c(type,bandwidth_factor)),
            by=c("theta_index", "model", "n", "alpha", "beta", "risk_measure", "init_method")) %>%
  left_join(res_MC_true %>% ungroup %>% dplyr::select(-c(type,bandwidth_factor)),
            by=c("theta_index", "model", "n", "alpha", "beta", "risk_measure", "init_method")) %>%
  dplyr::select(theta_index, model, n, beta, alpha, risk_measure, init_method, bandwidth_factor, theta_true, theta_mean, theta_median, theta_sd_emp, sd_asy_mean, sd_asy_median)


## Confidence interval coverage
res_MC_CI <- as_tibble(res_df_MC) %>%
  filter(model==model_choice, type %in% c("cov_est_asy")) %>%
  group_by(i_MC, model, n, alpha, beta, risk_measure, init_method, bandwidth_factor, type) %>%
  filter(theta_index %in% diag_entries) %>%
  mutate(theta_index=1:length(diag_entries)) %>% # This line should be robustified!!!
  reshape2::dcast(i_MC + theta_index + model + n + alpha + beta + risk_measure + init_method + bandwidth_factor ~ type) %>%
  left_join(res_df_MC %>% dplyr::filter(type=="param_est") %>% dplyr::select(-c(type,bandwidth_factor)) %>% rename(param_est = value),
            by=c("i_MC", "theta_index", "model", "n", "alpha", "beta", "risk_measure", "init_method")) %>%
  mutate(CI_asy_lower = param_est - qnorm(0.95)*sqrt(cov_est_asy),
         CI_asy_upper = param_est + qnorm(0.95)*sqrt(cov_est_asy)) %>%
  left_join(res_df_MC %>% filter(type=="true_value") %>% dplyr::select(-c(type,bandwidth_factor)) %>% rename(true_value = value),
            by=c("i_MC", "theta_index", "model", "n", "alpha", "beta", "risk_measure", "init_method")) %>%
  as_tibble()


# Compute coverage rates
res_MC_CIcoverage <- res_MC_CI %>%
  group_by(theta_index, model, n, alpha, beta, risk_measure, init_method, bandwidth_factor) %>%
  summarize(CI_asy_coverage = mean( (true_value >= CI_asy_lower) & (true_value <= CI_asy_upper)),
            CI_asy_length = median(CI_asy_upper - CI_asy_lower))


# Look at some results
res_MC %>% filter(n==4000)  %>% arrange(alpha)
res_MC_CIcoverage %>% filter(n==4000)  %>% arrange(alpha)





##########
##   Automatically print to LaTeX tables
##########

# Joint table displaying both, consistency and covariance estimation results
res_MC %>%
  filter(bandwidth_factor==1, init_method=="omega") %>%
  dplyr::mutate(bias_mean=(theta_mean-theta_true),
                bias_median=(theta_median-theta_true)) %>%
  dplyr::select(theta_index, beta, alpha, n, bias_mean, bias_median, theta_sd_emp, sd_asy_median) %>%
  full_join(res_MC_CIcoverage %>%
              filter(bandwidth_factor==1, init_method=="omega") %>%
              dplyr::select(c(theta_index, n, alpha, beta, CI_asy_coverage)),
            by=c("theta_index", "beta", "alpha", "n")) %>%
  mutate(theta_type=ifelse(theta_index<=3, "VaR", "CoVaR"),
         theta_index_partial=ifelse(theta_index<=3, theta_index, theta_index-3))  %>%   # This transforms the CoVaR parameters with number (3,4) to (1,2), etc
  tidyr::pivot_wider(id_cols = c("beta", "alpha", "n", "theta_type"),
                     names_from = theta_index_partial,
                     values_from = c(bias_mean, bias_median, theta_sd_emp, sd_asy_median, CI_asy_coverage)) %>%
  arrange(desc(theta_type), beta, alpha, n) %>%
  dplyr::filter(n<=4000) %>%
  mutate(empty1=NA, empty2=NA, empty3=NA) %>%
  dplyr::select(beta, n, empty1,
                bias_mean_1, bias_median_1, theta_sd_emp_1, sd_asy_median_1, CI_asy_coverage_1, empty2,
                bias_mean_2, bias_median_2, theta_sd_emp_2, sd_asy_median_2, CI_asy_coverage_2, empty3,
                bias_mean_3, bias_median_3, theta_sd_emp_3, sd_asy_median_3, CI_asy_coverage_3) %>%
  xtable::xtable(digits=c(0, 2, 0, 0,   c(4,4,3,3,2), 0,   c(4,4,3,3,2), 0,   c(4,4,3,3,2))) %>%
  print(file="simulations/output/CoCAViaR_AllResults_RR.txt", include.rownames=FALSE, booktabs=TRUE)





# Table analyzing the bandwidth choice
res_MC %>%
  filter(init_method=="omega") %>%
  dplyr::mutate(bias_mean=(theta_mean-theta_true),
                bias_median=(theta_median-theta_true)) %>%
  dplyr::select(theta_index, beta, alpha, n, bandwidth_factor, theta_sd_emp, sd_asy_median) %>%
  full_join(res_MC_CIcoverage %>%
              filter(init_method=="omega") %>%
              dplyr::select(c(theta_index, n, alpha, beta, bandwidth_factor, CI_asy_coverage)),
            by=c("theta_index", "beta", "alpha", "n", "bandwidth_factor")) %>%
  mutate(theta_type=ifelse(theta_index<=3, "VaR", "CoVaR"),
         theta_index_partial=ifelse(theta_index<=3, theta_index, theta_index-3))  %>%   # This transforms the CoVaR parameters with number (3,4) to (1,2), etc
  tidyr::pivot_wider(id_cols = c("beta", "alpha", "n", "theta_type", "bandwidth_factor"),
                     names_from = theta_index_partial,
                     values_from = CI_asy_coverage,
                     names_prefix = "param_") %>%
  tidyr::pivot_wider(id_cols = c("beta", "alpha", "n", "theta_type"),
                     names_from = bandwidth_factor,
                     values_from = c(param_1, param_2, param_3),
                     names_prefix="CI_") %>%
  arrange(desc(theta_type), beta, alpha, n) %>%
  mutate(empty1=NA, empty2=NA, empty3=NA) %>%
  dplyr::select(beta, n, empty1,
                param_1_CI_0.25, param_1_CI_0.5, param_1_CI_1, param_1_CI_2, param_1_CI_4, empty2,
                param_2_CI_0.25, param_2_CI_0.5, param_2_CI_1, param_2_CI_2, param_2_CI_4, empty3,
                param_3_CI_0.25, param_3_CI_0.5, param_3_CI_1, param_3_CI_2, param_3_CI_4) %>%
  xtable::xtable(digits=c(0, 2, 0, 0,   c(2,2,2,2,2), 0,   c(2,2,2,2,2), 0,   c(2,2,2,2,2))) %>%
  print(file="simulations/output/CoCAViaR_Bandwidth_RR.txt", include.rownames=FALSE, booktabs=TRUE)



# Table considering initialization with omega or at 0
res_MC %>%
  filter(bandwidth_factor==1) %>%
  dplyr::mutate(bias_mean=(theta_mean-theta_true),
                bias_median=(theta_median-theta_true)) %>%
  dplyr::select(theta_index, beta, alpha, n, init_method, bias_mean, theta_sd_emp) %>%
  mutate(theta_type=ifelse(theta_index<=3, "VaR", "CoVaR"),
         theta_index_partial=ifelse(theta_index<=3, theta_index, theta_index-3))  %>%   # This transforms the CoVaR parameters with number (3,4) to (1,2), etc
  tidyr::pivot_wider(id_cols = c("beta", "alpha", "n", "theta_type", "init_method"),
                     names_from = theta_index_partial,
                     values_from = c(bias_mean, theta_sd_emp)) %>%
  tidyr::pivot_wider(id_cols = c("beta", "alpha", "n", "theta_type"),
                     names_from = init_method,
                     values_from = c(bias_mean_1, bias_mean_2, bias_mean_3, theta_sd_emp_1, theta_sd_emp_2, theta_sd_emp_3),
                     names_prefix="") %>%
  arrange(desc(theta_type), beta, alpha, n) %>%
  mutate(empty1=NA, empty2=NA, empty3=NA) %>%
  dplyr::select(beta, n, empty1,
                bias_mean_1_omega, bias_mean_1_1, bias_mean_2_omega, bias_mean_2_1, bias_mean_3_omega, bias_mean_3_1, empty2,
                theta_sd_emp_1_omega, theta_sd_emp_1_1, theta_sd_emp_2_omega, theta_sd_emp_2_1, theta_sd_emp_3_omega, theta_sd_emp_3_1) %>%
  xtable::xtable(digits=c(0, 2, 0, 0,   c(3,3,3,3,3,3), 0,   c(3,3,3,3,3,3))) %>%
  print(file="simulations/output/CoCAViaR_Initialization_RR.txt", include.rownames=FALSE, booktabs=TRUE)

