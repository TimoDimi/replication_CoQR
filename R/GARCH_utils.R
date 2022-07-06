# This file contains utilities for the simulation of CoVaR models from a CCC GARCH model
# and functions to compute CoVaR (and MES) for a bivariate t-distribution

library(rmgarch)
library(mvtnorm)
library(cubature)



#============================================================
# Simulate CCC-GARCH(1,1) process with t-innovations
# ===========================================================

# Inputs:
# GARCH.coef.1 = GARCH parameters of X
# GARCH.coef.2 = GARCH parameters of Y
# alpha        = confidence level of CoVaR
# beta         = confidence level of VaR
# nu           = degrees of freedom of innovations
# rho          = correlation of innovations
# n            = sample size
# ntrans       = length of burn-in period


CCC.GARCH.sim <- function(GARCH.coef.1, GARCH.coef.2, alpha, beta, nu, rho, n, ntrans=1000){
  m   <- n + ntrans

  # 0. Calculate "t"-innovations eps
  sigma <- (nu - 2) / nu * matrix(c(1, rho, rho, 1), nrow=2)
  eps   <- rmvt(m, sigma, df = nu)  # correlated t_{nu}-variables with unit variance and correlation rho

  # 0.1 Calculate VaR, ES and MES of eps
  sd.wt <- sqrt(nu / (nu-2) )

  VaR   <- qt(p = beta, df = nu) / sd.wt
  MES   <- MES.true.t(nu, rho, beta)
  CoVaR <- CoVaR.true.t(alpha, beta, nu, matrix(c(1, rho, rho, 1), nrow=2))


  # 1. CCC-GARCH(1,1) process
  X <- Y <- sigma.X <- sigma.Y <- numeric(m)  # initializes vectors of X's and sigma^2's to be m long
  X[1] <- Y[1] <- rnorm(1)                    # (stochastic) start value
  sigma.X[1] <- sigma.Y[1] <- rnorm(1)^2      # (stochastic) start value
  for(t in 2 : m){
    sigma.X[t] <- GARCH.coef.1[1] + GARCH.coef.1[2] * abs( X[t-1] ) + GARCH.coef.1[3] * abs( Y[t-1] ) + GARCH.coef.1[4] * sigma.X[t-1] + GARCH.coef.1[5] * sigma.Y[t-1]
    sigma.Y[t] <- GARCH.coef.2[1] + GARCH.coef.2[2] * abs( X[t-1] ) + GARCH.coef.2[3] * abs( Y[t-1] ) + GARCH.coef.2[4] * sigma.X[t-1] + GARCH.coef.2[5] * sigma.Y[t-1]

    X[t]        <- sigma.X[t] * eps[t, 1]    # Simulate GARCH(1,1) process
    Y[t]        <- sigma.Y[t] * eps[t, 2]    # Simulate GARCH(1,1) process
  }

  # 2. Calcuate true conditional VaR, CoVaR and MES
  sigma.np1 <- GARCH.coef.2[1] + GARCH.coef.2[2] * abs( X[m] ) + GARCH.coef.2[3] * abs( Y[m] ) + GARCH.coef.2[4] * sigma.X[m] + GARCH.coef.2[5] * sigma.Y[m]  # true conditional variance at t=n+1

  VaR.true   <- sigma.np1 * VaR
  CoVaR.true <- sigma.np1 * CoVaR
  MES.true   <- sigma.np1 * MES

  # 3. Compute parameters in VaR / CoVaR / MES equations in (4)-(6) of the main paper
  VaR.par   <- c( GARCH.coef.1[1:3] * VaR,   GARCH.coef.1[4], GARCH.coef.1[5] * VaR / CoVaR )
  CoVaR.par <- c( GARCH.coef.2[1:3] * CoVaR, GARCH.coef.2[4] * CoVaR / VaR, GARCH.coef.2[5] )
  MES.par   <- c( GARCH.coef.2[1:3] * MES,   GARCH.coef.2[4] * MES / VaR,   GARCH.coef.2[5] )

  return( list(X = X[(ntrans + 1) : m], Y = Y[(ntrans + 1) : m],
               VaR.true = VaR.true, CoVaR.true = CoVaR.true, MES.true = MES.true,
               VaR.par  = VaR.par,  CoVaR.par  = CoVaR.par,  MES.par  = MES.par ) )
}





########
# TEST #
########

# GARCH.coef.1 <- c(0.01, 0.1, 0.02, 0.85, 0.05)
# GARCH.coef.2 <- c(0.01, 0.02, 0.1, 0.05, 0.85)
# alpha=0.95; beta=0.95; nu=5; rho=0.3; n=2000
#
# data <- CCC.GARCH.sim(GARCH.coef.1, GARCH.coef.2, alpha, beta, nu, rho, n)
#
# par( mfrow=c(2,1) )
# plot(data$X)
# plot(data$Y)





# ===========================================================
# Function for CoVaR calculation for bivariate, zero-mean t-distribution
# ===========================================================

# Inputs:
# alpha = confidence level of CoVaR
# beta  = confidence level of VaR
# nu    = degrees of freedom
# H     = variance-covariance matrix

CoVaR.true.t <- function(alpha, beta, nu, H){
  root.CoVaR.t <- function(x){
    sigma <-  (nu - 2) / nu * H
    sd.wt <- sqrt(nu / (nu-2) )
    VaR   <- qt(p = beta, df = nu) / sd.wt * sqrt( H[1,1] )
    prob  <- pmvt(lower=c(VaR, x), upper=c(Inf, Inf), df=nu, sigma=sigma)

    return( prob - (1-alpha) * (1-beta) )
  }
  uniroot(root.CoVaR.t, interval = c(0, 10), extendInt = "yes")$root
}


# ===========================================================
# Function for MES calculation for bivariate t-distribution
# ===========================================================

# Inputs:
# nu   = degrees of freedom
# rho  = correlation
# beta = risk level

MES.true.t <- function(nu, rho, beta){
  if(rho==0){
    MES <- 0
  }
  else{
    sd.wt <- sqrt(nu / (nu-2) )
    f.MES <- function(x, nu, rho) { x[2] * dmvt(x, sigma = (nu - 2) / nu * matrix(c(1, rho, rho, 1), nrow=2), df = nu, log = FALSE) } # "x" is vector
    VaR   <- qt(p = beta, df = nu) / sd.wt
    MES   <- 1/(1 - beta) * adaptIntegrate(f.MES, lowerLimit = c(VaR, -Inf), upperLimit = c(Inf, Inf), nu, rho)$integral
  }
  return( MES )
}




