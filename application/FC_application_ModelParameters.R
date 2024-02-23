
library(dplyr)
library(ggplot2)
library(tidyr)
library(rmgarch)
library(lubridate)
library(tsibble)
library(SystemicRisk)

source("R/SystemicDCCroll.R")
source("R/GARCH_utils.R")


####################################################################################
###
###     IMPORTANT:
###     The following code requires the file "data_Assets.rds" which we cannot make
###     available publicly as we do not have the licence to publish the raw data
###
####################################################################################


# Set options
alpha <- 0.95
beta <- 0.95
length_IS <- 3000
CoCAViaR_model_set <- c("CoCAViaR_SAV_diag", "CoCAViaR_SAV_fullA", "CoCAViaR_SAV_full",
                        "CoCAViaR_AS_pos", "CoCAViaR_AS_signs", "CoCAViaR_AS_mixed")

# Load data and define window length
data_Assets_full <- readRDS(file = "application/data/data_Assets.rds")
data_Assets <- data_Assets_full %>%
  dplyr::filter(Date <= "2021-12-31")

Date_max <- data_Assets$Date %>% unique() %>% .[length_IS]


# Fit models
CoCAViaR_modelfits <- list()
for (jj in 1:length(CoCAViaR_model_set)){
  set.seed(1)
  CoCAViaR_model <- CoCAViaR_model_set[jj]

  Mest_obj <- SRM(data = data_Assets %>%
                    dplyr::filter(Asset %in% c("JPM", "SP500"), Date <= Date_max) %>%
                    tidyr::pivot_wider(names_from = "Asset", values_from = "NegReturn") %>%
                    rename(x=JPM, y=SP500) %>%
                    as_tsibble(index=Date),
                  model=CoCAViaR_model, risk_measure="CoVaR", beta=beta, alpha=alpha)
  CoCAViaR_modelfits[[paste(CoCAViaR_model)]] <- summary(Mest_obj)
}


CoCAViaR_modelfits$CoCAViaR_SAV_diag

saveRDS(CoCAViaR_modelfits, file = "application/data/CoCAViaR_ModelParameters_RR.rds")


#  Plot in-sample model fits
CoCAViaR_SAV_fullA_fit <- SRM(data = data_Assets %>%
                                dplyr::filter(Asset %in% c("JPM", "SP500"), Date <= Date_max) %>%
                                tidyr::pivot_wider(names_from = "Asset", values_from = "NegReturn") %>%
                                rename(x=JPM, y=SP500) %>%
                                as_tsibble(index=Date),
                              model="CoCAViaR_SAV_full", risk_measure="CoVaR", beta=beta, alpha=alpha)

plot(CoCAViaR_SAV_fullA_fit, facet_names=c("JPM / VaR","SP500 / CoVaR"))
ggsave("application/output/appl_ISplot_RR.pdf", width=20, height=15, units="cm")



#  Plot out-of-sample forecasts
CoCAViaR_SAV_fullA_roll <- SRMroll(data = data_Assets %>%
                                 dplyr::filter(Asset %in% c("JPM", "SP500")) %>%
                                 tidyr::pivot_wider(names_from = "Asset", values_from = "NegReturn") %>%
                                 rename(x=JPM, y=SP500) %>%
                                 as_tsibble(index=Date),
                               model="CoCAViaR_SAV_full",
                               length_IS=3000, refit_freq=100,
                               risk_measure="CoVaR", beta=beta, alpha=alpha)

plot(CoCAViaR_SAV_fullA_roll, facet_names=c("JPM / VaR","SP500 / CoVaR"))
ggsave("application/output/appl_OOSplot_RR.pdf", width=20, height=15, units="cm")

