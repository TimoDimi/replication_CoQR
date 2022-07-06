library(dplyr)
library(ggplot2)
library(reshape2)
library(lubridate)
library(CoQR)


####################################################################################
###
###     IMPORTANT:
###     The following code requires the file "data_appl_AB16.rds" which we cannot make
###     available publicly as we do not have the licence to publish the raw data
###
####################################################################################


set.seed(2022)
data.set <- readRDS(file = "application/data/data_appl_AB16.rds") %>%
  dplyr::filter(Date <= "2021-12-31")

# Lag all covariates by one day and transform to tsibble
  summarize(Date=Date,
            x=JPM.loss,
            y=SPF.loss,
            Intercept=1,
            Spread = lag(Spread),
            ChangeSpread = lag(ChangeSpread),
            TEDSpread = lag(TEDSpread),
            SP500.ret = lag(SP500.ret),
            RV.data = lag(RV.data)) %>%
  na.omit() %>%
  as_tsibble(index=Date)


CoQR_obj <- CoQR(data=data_lagged,
     model="joint_linear",
     SRM="CoVaR",
     beta=0.95,
     alpha=0.95,
     optim_replications=c(5,20))

summary(CoQR_obj)
# VaR Coefficients:
#   Estimate  Std. Error t value  Pr(>|t|)
#   Intercept     0.01259995  0.00310781  4.0543 5.099e-05 ***
#   Spread        0.00247625  0.00077196  3.2078 0.0013455 **
#   TEDSpread     0.02658856  0.00653465  4.0689 4.792e-05 ***
#   ChangeSpread -0.03113515  0.02747426 -1.1332 0.2571604
#   SP500.ret    -0.15070443  0.06931602 -2.1742 0.0297360 *
#   RV.data       0.71508069  0.19592008  3.6499 0.0002649 ***
#
# CoVaR Coefficients:
#   Estimate Std. Error t value  Pr(>|t|)
#   Intercept     0.0053722  0.0076382  0.7033  0.481884
#   Spread        0.0023487  0.0026643  0.8815  0.378076
#   TEDSpread     0.0631499  0.0093250  6.7721 1.404e-11 ***
#   ChangeSpread  0.1923377  0.0924236  2.0810  0.037477 *
#   SP500.ret    -0.1371597  0.4910525 -0.2793  0.780012
#   RV.data       3.0606485  1.1471962  2.6679  0.007655 **
#
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


plot(CoQR_obj)


