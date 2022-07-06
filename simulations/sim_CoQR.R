library(MASS)
library(dplyr)
library(tibble)
library(ggplot2)
library(doParallel)
library(MTS)
library(tsibble)
library(CoQR)

# Load CoQR DGP function
source("R/CoQR_DGP.R")
source("R/GARCH_utils.R")


# Preliminary settings
M <- 5000
optim_replications <- c(1,5)
n_set <- c(500, 1000, 2000, 4000)
prob_level_list <- list(c(beta=0.9, alpha=0.9),
                        c(beta=0.95, alpha=0.95))

SRM_set <- c("CoVaR")
model_set <- c("CoQR_4p", "CoQR_6p")
gamma_list <- list(c(1, 1.5, 0, 0.25, 0.5, 0),
                   c(1, 1.5, 2, 0.25, 1, 0))

sigma1 <- 2
sigma2 <- 1
rho <- 0.5
Sigma2 <- matrix(c(sigma1^2, rho*sigma1*sigma2, rho*sigma1*sigma2, sigma2^2), nrow=2) # Covariance matrix
Sigma <- t(chol(Sigma2)) # Square-root of Covariance matrix
phi <- c(0.5, 0.8)
nu <- 8


# Cluster Settings
core.max <- 90
cl <- makeCluster(min(parallel::detectCores()-1, M, core.max) )
registerDoParallel(cl)
start.time <- Sys.time()
res_df_MC <- foreach(
  i_MC = 1:M,
  .combine=rbind,
  .packages=c("dplyr", "tibble", "MTS", "MASS", "mvtnorm", "cubature", "rmgarch", "lubridate", "abind", "tsibble", "CoQR"),
  .errorhandling="pass"
)%dopar%{
  set.seed(i_MC) # set seed for reproducibility
  res_df <- tibble()

  # Load CoQR DGP function
  source("R/CoQR_DGP.R")
  source("R/GARCH_utils.R")

  for (model in model_set){
    for (n in n_set){
      for (prob_level in prob_level_list){
        beta <- prob_level[[1]]
        alpha <- prob_level[[2]]

        # Simulate data
        gamma <- switch(model,
                        CoQR_4p = {gamma_list[[1]]},
                        CoQR_6p = {gamma_list[[2]]})

        sim_help <- sim_CoQR_DGP(n=n, gamma=gamma, Sigma=Sigma, nu=nu, phi=phi, beta=beta, alpha=alpha)
        theta_true <- sim_help$theta_true

        # Transform data to tsibble
        data <- sim_help$data %>%
          mutate(Date=1:n()) %>%
          as_tsibble(index=Date)

        for (SRM in SRM_set){
          # Parameter estimate
          if (model=="CoQR_4p"){
            theta0 <- theta_true[-c(3,6)]
            est_obj <- CoQR(data = data %>% dplyr::select(-z3), model="joint_linear", SRM=SRM, beta=beta, alpha=alpha, theta0=theta0, optim_replications=optim_replications)
          } else if (model=="CoQR_6p"){
            theta0 <- theta_true
            est_obj <- CoQR(data=data, model="joint_linear", SRM=SRM, beta=beta, alpha=alpha, theta0=theta0, optim_replications=optim_replications)
          } else {
            stop("enter a correct model name!")
          }

          df_theta_est <- data.frame(i_MC = i_MC,
                                     model=model,
                                     n=n,
                                     alpha=alpha,
                                     beta=beta,
                                     SRM = SRM,
                                     type = "param_est",
                                     theta_index = 1:length(est_obj$theta),
                                     value = est_obj$theta)

          df_theta_true <- data.frame(i_MC = i_MC,
                                      model=model,
                                      n=n,
                                      alpha=alpha,
                                      beta=beta,
                                      SRM = SRM,
                                      type = "true_value",
                                      theta_index = 1:length(theta0),
                                      value = theta0)


          # Covariance estimation
          sum_obj_asy <- summary(est_obj)

          df_theta_cov_asy <-  data.frame(i_MC = i_MC,
                                          model=model,
                                          n=n,
                                          alpha=alpha,
                                          beta=beta,
                                          SRM = SRM,
                                          type = "cov_est_asy",
                                          theta_index = 1:length(MTS::Vech(sum_obj_asy$cov)),
                                          value = MTS::Vech(sum_obj_asy$cov))

          res_df <- rbind(res_df,
                          df_theta_est,
                          df_theta_true,
                          df_theta_cov_asy)
        }
      }
    }
  }
  res_df
}
stopCluster(cl)
end.time <- Sys.time()
(run.time <- end.time-start.time)

head(res_df_MC)
saveRDS(res_df_MC, file = "simulations/data/sim_crossDGP.rds")


