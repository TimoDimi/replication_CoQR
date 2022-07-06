library(dplyr)
library(tidyr)
library(reshape2)
library(CoQR)


# Load simulated data
res_df_MC <- readRDS(file = "simulations/data/sim_crossDGP.rds")


# Set options
model_choice <- "CoQR_6p"
q1 <- 3
diag_entries <- which(MTS::Vech(diag(2*q1))==1)

# Estimated Covariance Matrix and standard deviations (diagonal covariance entries only) as data frames
res_MC_cov <- res_df_MC %>%
  filter(model==model_choice, type %in% c("cov_est_asy")) %>%
  group_by(theta_index, model, n, alpha, beta, SRM, type) %>%
  summarize(theta_cov_mean = mean(value, na.rm = TRUE),
            theta_cov_median = median(value, na.rm = TRUE))

res_MC_sd <- res_MC_cov %>%
  filter(theta_index %in% diag_entries) %>%
  group_by(model, n, alpha, beta, SRM) %>%
  mutate(theta_index = 1:length(theta_index)) %>% # Reset theta index after selecting diagonal entries only
  group_by(theta_index) %>%
  mutate(sd_asy_mean = sqrt(theta_cov_mean),
         sd_asy_median = sqrt(theta_cov_median))

# True parameter values
res_MC_true <- res_df_MC %>%
  filter(model==model_choice, type=="true_value") %>%
  group_by(theta_index, model, n, alpha, beta, SRM) %>%
  summarize(theta_true = mean(value,na.rm = TRUE))

# Extract the mean parameter estimates
res_MC_est <- res_df_MC %>%
  dplyr::filter(model==model_choice, type=="param_est") %>%
  group_by(theta_index, model, n, alpha, beta, SRM) %>%
  summarize(theta_mean = mean(value),
            theta_median = median(value),
            theta_sd_emp = sd(value)) %>%
  ungroup()

# Join data frames together
res_MC <- res_MC_est %>%
  left_join(res_MC_true,
            by=c("theta_index", "model", "n", "alpha", "beta", "SRM")) %>%
  left_join(res_MC_sd,
            by=c("theta_index", "model", "n", "alpha", "beta", "SRM")) %>%
  dplyr::select(theta_index, model, SRM, beta, alpha, n, theta_true, theta_mean, theta_median, theta_sd_emp, sd_asy_mean, sd_asy_median)



## Confidence interval coverage
res_MC_CI <- res_df_MC %>%
  filter(model==model_choice, type %in% c("cov_est_asy")) %>%
  group_by(i_MC, model, n, alpha, beta, SRM, type) %>%
  filter(theta_index %in% diag_entries) %>%
  mutate(theta_index=1:length(diag_entries)) %>% # This line is dangerous!
  reshape2::dcast(i_MC + theta_index + model + n + alpha + beta + SRM ~ type) %>%
  left_join(res_df_MC %>% dplyr::filter(type=="param_est") %>% select(-type) %>% rename(param_est = value),
            by=c("i_MC", "theta_index", "model", "n", "alpha", "beta", "SRM")) %>%
  mutate(CI_asy_lower = param_est - qnorm(0.975)*sqrt(cov_est_asy),
         CI_asy_upper = param_est + qnorm(0.975)*sqrt(cov_est_asy)) %>%
  left_join(res_df_MC %>% filter(type=="true_value") %>% select(-type) %>% rename(true_value = value),
            by=c("i_MC", "theta_index", "model", "n", "alpha", "beta", "SRM")) %>%
  as_tibble()


# Compute CI coverage
res_MC_CIcoverage <- res_MC_CI %>%
  group_by(theta_index, model, n, alpha, beta, SRM) %>%
  summarize(CI_asy_coverage = mean( (true_value >= CI_asy_lower) & (true_value <= CI_asy_upper)),
            CI_asy_length = median(CI_asy_upper - CI_asy_lower))



# Look at specific results
res_MC %>% filter(n==4000) %>% arrange(alpha)
res_MC_CIcoverage %>% filter(n==4000) %>% arrange(alpha)




##########
##   Automatically print results to LaTeX tables
##########



# Joint table displaying both, consistency and covariance estimation results
res_MC %>%
  dplyr::mutate(bias_mean=(theta_mean-theta_true),
                bias_median=(theta_median-theta_true)) %>%
  dplyr::select(theta_index, beta, alpha, n, bias_mean, bias_median, theta_sd_emp, sd_asy_median) %>%
  full_join(res_MC_CIcoverage %>% dplyr::select(c(theta_index, n, alpha, beta, CI_asy_coverage)),
            by=c("theta_index", "beta", "alpha", "n")) %>%
  mutate(theta_type=ifelse(theta_index<=q1, "VaR", "CoVaR"),
         theta_index_partial=ifelse(theta_index<=q1, theta_index, theta_index-q1))  %>%   # This appears weird, transform CoVaR parameters with number (3,4) to (1,2)...
  tidyr::pivot_wider(id_cols = c("beta", "alpha", "n", "theta_type"),
                     names_from = theta_index_partial,
                     values_from = c(bias_mean, bias_median, theta_sd_emp, sd_asy_median, CI_asy_coverage)) %>%
  arrange(desc(theta_type), beta, alpha, n) %>%
  mutate(empty1=NA, empty2=NA, empty3=NA) %>%
  dplyr::select(beta, n, empty1, bias_mean_1, bias_median_1, theta_sd_emp_1, sd_asy_median_1, CI_asy_coverage_1,
                empty2, bias_mean_2, bias_median_2, theta_sd_emp_2, sd_asy_median_2, CI_asy_coverage_2,
                empty3, bias_mean_3, bias_median_3, theta_sd_emp_3, sd_asy_median_3, CI_asy_coverage_3) %>%
  xtable::xtable(digits=c(0, 2, 0, 0, c(4,4,3,3,2), 0, c(4,4,3,3,2), 0, c(4,4,3,3,2) )) %>%
  print(file="simulations/output/CoQR_AllResults.txt", include.rownames=FALSE, booktabs=TRUE)








# Separate tables for consistency and inference results

#
# # Save consistency results to a LaTeX table
# res_MC %>%
#   dplyr::mutate(bias=(theta_mean-theta_true),
#                 median_bias=(theta_median-theta_true)) %>%
#   dplyr::select(theta_index, beta, alpha, n, bias, median_bias) %>%
#   tidyr::pivot_wider(id_cols = c("beta", "alpha", "n"),
#                      names_from = theta_index,
#                      values_from = c(bias, median_bias)) %>%
#   arrange(beta, alpha, n) %>%
#   select(n, bias_1, median_bias_1, bias_2, median_bias_2, bias_3, median_bias_3, bias_4, median_bias_4) %>%
#   xtable::xtable(digits=c(0, 0, rep(4,8))) %>%
#   print(file="simulations/output/CoQR_consistency.txt", include.rownames=FALSE, booktabs=TRUE)
#
#
# # Save covariance results to a LaTeX table
# res_MC %>%
#   dplyr::select(c(theta_index, n, alpha, beta, theta_sd_emp, sd_asy_median)) %>%
#   full_join(res_MC_CIcoverage %>% dplyr::select(c(theta_index, n, alpha, beta, CI_asy_coverage)),
#             by=c("theta_index", "beta", "alpha", "n")) %>%
#   tidyr::pivot_wider(id_cols = c("beta", "alpha", "n"),
#                      names_from = theta_index,
#                      values_from = c(theta_sd_emp, sd_asy_median, CI_asy_coverage)) %>%
#   arrange(beta, alpha, n) %>%
#   select(n, theta_sd_emp_1, sd_asy_median_1, CI_asy_coverage_1,
#          theta_sd_emp_2, sd_asy_median_2, CI_asy_coverage_2,
#          theta_sd_emp_3, sd_asy_median_3, CI_asy_coverage_3,
#          theta_sd_emp_4, sd_asy_median_4, CI_asy_coverage_4,
#          theta_sd_emp_5, sd_asy_median_5, CI_asy_coverage_5,
#          theta_sd_emp_6, sd_asy_median_6, CI_asy_coverage_6) %>%
#   xtable::xtable(digits=c(0, 0, rep(3,18))) %>%
#   print(file="simulations/output/CoQR_AsyCovEst.txt", include.rownames=FALSE, booktabs=TRUE)





