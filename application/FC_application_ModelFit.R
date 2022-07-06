library(dplyr)
library(ggplot2)
library(reshape2)
library(rmgarch)
library(mvtnorm)
library(cubature)
library(lubridate)
library(doParallel)
library(tsibble)
library(CoQR)

# Load functions for the DCC models
source("R/GARCH_utils.R")
source("R/SystemicDCCroll.R")

####################################################################################
###
###     IMPORTANT:
###     The following code requires the file "data_Assets.rds" which we cannot make
###     available publicly as we do not have the licence to publish the raw data
###
####################################################################################


###    This file should be ran on a cluster. It takes approximately 6-12 hours per
###    iteration in the foreach loop


# 0. Set options
alpha <- 0.95
beta <- 0.95
length_IS <- 3000
refit_freq <- 100
SRM <- "CoVaR"

# CoCAViaR models
CoCAViaR_model_set <- c("CoCAViaR_SAV_diag", "CoCAViaR_SAV_fullA", "CoCAViaR_SAV_full",
                        "CoCAViaR_AS_pos", "CoCAViaR_AS_signs", "CoCAViaR_AS_mixed")


# DCC GARCH specifications
GARCH11_norm_spec <- ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean=FALSE),
                                variance.model = list(garchOrder = c(1,1), model = "sGARCH"),
                                distribution.model = "norm")
GARCH11_std_spec <- ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean=FALSE),
                               variance.model = list(garchOrder = c(1,1), model = "sGARCH"),
                               distribution.model = "std")
gjrGARCH11_std_spec <- ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean=FALSE),
                                  variance.model = list(garchOrder = c(1,1), model = "gjrGARCH"),
                                  distribution.model = "std")

DCC_GARCH11_norm_spec <- dccspec(uspec = multispec(replicate(2, GARCH11_norm_spec)), dccOrder = c(1,1), distribution = "mvnorm")
DCC_GARCH11_std_spec <- dccspec(uspec = multispec(replicate(2, GARCH11_std_spec)), dccOrder = c(1,1), distribution = "mvt")
DCC_gjrGARCH11_std_spec <- dccspec(uspec = multispec(replicate(2, gjrGARCH11_std_spec)), dccOrder = c(1,1), distribution = "mvt")

DCC_spec_list <- list("DCC-GARCH-n" = DCC_GARCH11_norm_spec,
                      "DCC-GARCH-t" = DCC_GARCH11_std_spec,
                      "DCC-gjrGARCH-t" = DCC_gjrGARCH11_std_spec)

H_sqrt_set <- c("Cholesky", "symmetric")

# Load assets
data_Assets <- readRDS(file = "application/data/data_Assets.rds") %>%
  dplyr::filter(Date <= "2021-12-31")

asset_names <- c("BAC", "C", "GS", "JPM", "SP500", "SPF")

# All (non-identical) asset combinations
AssetCombinations <- rbind(t(combn(asset_names, 2)),
                           t(combn(rev(asset_names), 2)))

# Only combinations with SP500 or SPF
i_SP500 <- which(as.logical(rowSums(AssetCombinations == "SP500")))
i_SPF <- which(as.logical(rowSums(AssetCombinations == "SPF")))
AssetCombinations_subset <- AssetCombinations[unique(c(i_SP500, i_SPF)),]
M_AssetPairs <- nrow(AssetCombinations_subset)



# Cluster Settings
core.max <- 40
cl <- makeCluster(min(parallel::detectCores()-1, M_AssetPairs, core.max) )
registerDoParallel(cl)
start_time <- Sys.time()
FCs_application <- foreach(
  i_AssetPair = 1:M_AssetPairs,
  .combine=rbind,
  .packages=c("dplyr", "ggplot2", "tibble", "MTS", "reshape2", "rmgarch", "mvtnorm", "cubature", "lubridate", "CoQR", "tsibble"),
  .errorhandling="remove"
)%dopar%{
  # Start parallel loop

  # Load DCC functions
  source("R/GARCH_utils.R")
  source("R/SystemicDCCroll.R")

  # Select data
  AssetPair <- AssetCombinations_subset[i_AssetPair,]

  data_AssetPair <- data_Assets %>%
    dplyr::filter(Asset %in% AssetPair) %>%
    tidyr::pivot_wider(names_from = "Asset", values_from = "NegReturn") %>%
    dplyr::rename(x=AssetPair[1], y=AssetPair[2]) %>%
    tsibble::as_tsibble(index=Date) %>%
    dplyr::select(Date, x,y)

  # Dynamic CoQR forecasts for CoVaR
  FC_tbl <- tibble()
  for (CoCAViaR_model in CoCAViaR_model_set){
    CoQRroll_obj_tmp <- CoQRroll(data=data_AssetPair,
                                 length_IS=length_IS, refit_freq=refit_freq,
                                 model=CoCAViaR_model, SRM=SRM, beta=beta, alpha=alpha)

    FC_tbl <- bind_rows(FC_tbl,
                        CoQRroll_obj_tmp$FC_df %>% mutate(SRM=SRM, model=CoCAViaR_model))
  }

  # DCC GARCH forecasts (for CoVaR and MES)
  for (i_DCC_spec in 1:length(DCC_spec_list)){
    for (H_sqrt_method in H_sqrt_set){

      SystemicDCCroll_obj <- SystemicDCCroll(data=data_AssetPair,
                                             DCC_spec=DCC_spec_list[[i_DCC_spec]],
                                             length_IS=length_IS, refit_freq=refit_freq, beta=beta, alpha=alpha,
                                             DCC_FC_dist=c("np"), H_sqrt=H_sqrt_method)

      # Only consider CoVaR forecasts
      FC_tbl <- bind_rows(FC_tbl,
                          SystemicDCCroll_obj$SRM_FC %>%
                            dplyr::select(-MES) %>%
                            mutate(SRM="CoVaR", model=paste0(names(DCC_spec_list)[i_DCC_spec],"-", H_sqrt_method)))
    }
  }

  FC_tbl %>% mutate(assetX=AssetPair[1], assetY=AssetPair[2])
}
stopCluster(cl)
end_time <- Sys.time()
(run_time <- end_time-start_time)


head(FCs_application)
saveRDS(FCs_application, file = "application/data/application_FCs.rds")


