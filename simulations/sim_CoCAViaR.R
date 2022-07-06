library(MASS)
library(dplyr)
library(tibble)
library(ggplot2)
library(doParallel)
library(mvtnorm)
library(cubature)
library(MTS)
library(rmgarch)
library(CoQR)

# Load GARCH DGP function
source("R/GARCH_utils.R")

# setwd("/home/dimitriadis_520k/Documents/Code/CoQR/")

M <- 5000
optim_replications <- c(1,5)
n_set <- c(500, 1000, 2000, 4000)
prob_level_list <- list(c(beta=0.9, alpha=0.9), c(beta=0.95, alpha=0.95))
model_set <- c("CoCAViaR_SAV_diag")
SRM_set <- c("CoVaR")

# Cluster Settings
core.max <- 90
cl <- makeCluster(min(parallel::detectCores()-1, M, core.max) )
registerDoParallel(cl)
start.time <- Sys.time()
res_df_MC <- foreach(
  i_MC = 1:M,
  .combine=rbind,
  .packages=c("dplyr", "tibble", "MTS", "MASS", "mvtnorm", "cubature","rmgarch","lubridate","abind","tsibble", "CoQR"),
  .errorhandling="pass"
)%dopar%{

  # Load GARCH DGP function
  source("R/GARCH_utils.R")

  # Set seed for reproducibility
  set.seed(i_MC)
  res_df <- tibble()

  for (model in model_set){
    for (n in n_set){
      for (prob_level in prob_level_list){
        beta <- prob_level[1]
        alpha <- prob_level[2]

        # Simulate data from a multivariate GARCH model; switch parameters depending on the model
        nu <- 8
        rho <- 0.5
        if (model == "CoCAViaR_SAV_diag"){
          GARCH.coef.1 <- c(0.04, 0.15, 0, 0.75, 0)
          GARCH.coef.2 <- c(0.02, 0, 0.1, 0, 0.8)
          mGARCH_sim <- CCC.GARCH.sim(GARCH.coef.1, GARCH.coef.2, alpha, beta, nu, rho, n)
          theta0_CoVaR <- c(mGARCH_sim$VaR.par[c(1,2,4)], mGARCH_sim$CoVaR.par[c(1,3,5)])
          theta0_MES <- c(mGARCH_sim$VaR.par[c(1,2,4)], mGARCH_sim$MES.par[c(1,3,5)])
        } else if (model == "CoCAViaR_SAV_fullA"){
          GARCH.coef.1 <- c(0.04, 0.08, 0.07, 0.8, 0)
          GARCH.coef.2 <- c(0.02, 0.05, 0.1, 0, 0.75)
          mGARCH_sim <- CCC.GARCH.sim(GARCH.coef.1, GARCH.coef.2, alpha, beta, nu, rho, n)
          theta0_CoVaR <- c(mGARCH_sim$VaR.par[-5], mGARCH_sim$CoVaR.par[-4])
          theta0_MES <- c(mGARCH_sim$VaR.par[-5], mGARCH_sim$MES.par[-4])
        }


        # Transform data to tsibble
        data <- tsibble(x = mGARCH_sim$X,
                        y = mGARCH_sim$Y,
                        Date=1:length(mGARCH_sim$X),
                        index=Date)

        for (SRM in SRM_set){
          # True parameters
          if (SRM == "CoVaR"){
            theta0 <- theta0_CoVaR
          } else {
            theta0 <- theta0_MES
          }

          est_obj <- CoQR(data=data,
                          model=model, SRM=SRM, beta=beta, alpha=alpha,
                          theta0=theta0, optim_replications=optim_replications)

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
                                      theta_index = 1:length(est_obj$theta),
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
saveRDS(res_df_MC, file = "simulations/data/sim_CoCAViaR.rds")

