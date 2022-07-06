library(dplyr)
library(ggplot2)
library(reshape2)
library(CoQR)

####################################################################################
###
###     IMPORTANT:
###     The following code requires the file "application_FCs.rds" which we cannot make
###     available publicly as we do not have the licence to publish the raw data
###
####################################################################################


# Load data:
FCs_application <- readRDS(file = "application/data/application_FCs.rds")

# Rename models
FCs_application <- FCs_application %>%
  mutate(model=recode(model,
                      "CoCAViaR_SAV_fullA" = "CoCAViaR-SAV-fullA",
                      "CoCAViaR_SAV_diag" = "CoCAViaR-SAV-diag",
                      "CoCAViaR_SAV_full" = "CoCAViaR-SAV-full",
                      "CoCAViaR_AS_signs" = "CoCAViaR-AS-signs",
                      "CoCAViaR_AS_mixed" = "CoCAViaR-AS-mixed",
                      "CoCAViaR_AS_pos" = "CoCAViaR-AS-pos",
                      "DCC-GARCH-n-Cholesky" = "DCC-n-Chol",
                      "DCC-GARCH-n-symmetric" = "DCC-n-sym",
                      "DCC-GARCH-t-Cholesky" = "DCC-t-Chol",
                      "DCC-GARCH-t-symmetric" = "DCC-t-sym",
                      "DCC-gjrGARCH-t-Cholesky" = "DCC-gjr-t-Chol",
                      "DCC-gjrGARCH-t-symmetric" = "DCC-gjr-t-sym"))
cov_est_method <- "HAC"

# Forecast Evaluation:
# compute average losses, ranks, and violations (hits) per asset and model
FC_eval_application <- FCs_application %>%
  group_by(assetX, assetY, model, SRM) %>%
  summarize(loss_VaR = mean(loss_VaR(x=x, VaR=VaR, beta=beta)),
            loss_CoVaR = mean(loss_CoVaR(x=x, y=y, VaR=VaR, CoVaR=CoVaR, alpha=alpha)),
            loss_CoVaR_std = sum(loss_CoVaR(x=x, y=y, VaR=VaR, CoVaR=CoVaR, alpha=alpha)) / sum((x>=VaR)),
            VaR_violations = mean(x >= VaR),
            CoVaR_violations = sum((x >= VaR & y >= CoVaR), na.rm=T)/sum((x >= VaR), na.rm=T)) %>%
  filter(SRM %in% c("CoVaR", "both")) %>%
  ungroup() %>%
  group_by(assetX, assetY, SRM) %>%
  mutate(rk_loss_VaR = rank(loss_VaR),
         rk_loss_CoVaR = rank(loss_CoVaR))


# One-half sided tests of Fissler and Hoga (2021, arXiv)
# First create "VAR_baseline" and "CoVaR_baseline" columns with repeated entries (throughout different models)
# Then run the "SystemicRiskFCeval" function for each model
baseline_model <-  "DCC-t-Chol"
FCs_baseline <- left_join(FCs_application,
                          FCs_application %>% dplyr::filter(model==baseline_model),
                          by=c("Date", "Date_index", "x", "y", "SRM", "assetX", "assetY"),
                          suffix = c("", "_baseline")) %>%
  filter(model!=baseline_model)

FC_eval_inference <- FCs_baseline %>%
  group_by(assetX, assetY, SRM, model) %>%
  summarize(baseline_model=unique(model_baseline),
            zone = SystemicRiskFCeval(x=x, y=y,
                                      VaR1=VaR_baseline,
                                      VaR2=VaR,
                                      CoVaR1=CoVaR_baseline,
                                      CoVaR2=CoVaR,
                                      SRM="CoVaR", beta=0.95, alpha=0.95,
                                      sided="onehalf", cov_method=cov_est_method,
                                      sig_level=0.1)$zone,
            pval = SystemicRiskFCeval(x=x, y=y,
                                      VaR1=VaR_baseline,
                                      VaR2=VaR,
                                      CoVaR1=CoVaR_baseline,
                                      CoVaR2=CoVaR,
                                      SRM="CoVaR", beta=0.95, alpha=0.95,
                                      sided="onehalf", cov_method=cov_est_method,
                                      sig_level=0.1)$pval %>% as.numeric())

# Join the two data frames
FC_eval_joint <- left_join(FC_eval_application,
          FC_eval_inference,
          by=c("assetX", "assetY", "SRM", "model")) %>%
  arrange(assetX, assetY, loss_CoVaR) %>%
  ungroup()


# Specify the order of the models
model_order <- c("CoCAViaR-SAV-fullA",
                 "CoCAViaR-SAV-diag",
                 "CoCAViaR-SAV-full",
                 "CoCAViaR-AS-signs",
                 "CoCAViaR-AS-mixed",
                 "CoCAViaR-AS-pos",
                 "DCC-n-Chol",
                 "DCC-n-sym",
                 "DCC-t-Chol",
                 "DCC-t-sym",
                 "DCC-gjr-t-Chol",
                 "DCC-gjr-t-sym")



# Print the results for Y_t = SP500 to a LaTeX table
FC_eval_joint %>%
  dplyr::filter(assetX %in% c("BAC", "C", "GS", "JPM", "SPF"),
                assetY=="SP500") %>%
  dplyr::mutate(loss_VaR = 10*loss_VaR, loss_CoVaR=1000*loss_CoVaR,
                VaR_violations=100*VaR_violations, CoVaR_violations=100*CoVaR_violations,
                empty1=NA, empty2=NA, empty3=NA) %>%
  dplyr::select(assetX, model, empty1, loss_VaR, rk_loss_VaR, VaR_violations, empty2, loss_CoVaR, rk_loss_CoVaR, CoVaR_violations, empty3, zone, pval) %>%
  group_by(assetX) %>%
  # arrange(match(model, model_order), .by_group = TRUE) %>%
  arrange(loss_CoVaR, loss_VaR, .by_group = TRUE) %>%
  xtable::xtable(digits=c(rep(0,4), 3,0,1, 0, 3,0,1, 0, 0, 2)) %>%
  print(file="application/output/FC_application_Banks_SP500.txt", include.rownames=FALSE, booktabs=TRUE)



# Print the results for Y_t = SP Financials (SPF) to a LaTeX table
FC_eval_joint %>%
  dplyr::filter(assetX %in% c("BAC", "C", "GS", "JPM"),
                assetY=="SPF") %>%
  dplyr::mutate(loss_VaR = 10*loss_VaR, loss_CoVaR=1000*loss_CoVaR,
                VaR_violations=100*VaR_violations, CoVaR_violations=100*CoVaR_violations,
                empty1=NA, empty2=NA, empty3=NA) %>%
  dplyr::select(assetX, model, empty1, loss_VaR, rk_loss_VaR, VaR_violations, empty2, loss_CoVaR, rk_loss_CoVaR, CoVaR_violations, empty3, zone, pval) %>%
  group_by(assetX) %>%
  arrange(loss_CoVaR, loss_VaR, .by_group = TRUE) %>%
  # arrange(match(model,model_order), .by_group = TRUE) %>%
  xtable::xtable(digits=c(rep(0,4), 3,0,1, 0, 3,0,1, 0, 0, 2)) %>%
  print(file="application/output/FC_application_Banks_SPF.txt", include.rownames=FALSE, booktabs=TRUE)





# Exemplary plots:
# JPM as systematically most relevant bank
FCs_plot <- FCs_baseline %>% filter(assetX=="JPM", assetY=="SP500", model=="CoCAViaR-SAV-fullA")
SysEvalObj <- with(FCs_plot, SystemicRiskFCeval(x=x, y=y,
                                 VaR1=VaR_baseline,
                                 VaR2=VaR,
                                 CoVaR1=CoVaR_baseline,
                                 CoVaR2=CoVaR,
                                 SRM="CoVaR", beta=0.95, alpha=0.95,
                                 sided="onehalf", cov_method=cov_est_method,
                                 sig_level=0.1))
p <- autoplot(SysEvalObj) +
  xlab("Average VaR Loss Difference") +
  ylab("Average CoVaR Loss Difference") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))

p

ggsave("application/output/FC_eval_JPMSP500_SAV.pdf", p, width=15, height=12, units="cm")

# Financial Index
FCs_plot <- FCs_baseline %>% filter(assetX=="JPM", assetY=="SP500", model=="CoCAViaR-AS-mixed")
SysEvalObj <- with(FCs_plot, SystemicRiskFCeval(x=x, y=y,
                                                VaR1=VaR_baseline,
                                                VaR2=VaR,
                                                CoVaR1=CoVaR_baseline,
                                                CoVaR2=CoVaR,
                                                SRM="CoVaR", beta=0.95, alpha=0.95,
                                                sided="onehalf", cov_method=cov_est_method,
                                                sig_level=0.1))


p <- autoplot(SysEvalObj) +
  xlab("Average VaR Loss Difference") +
  ylab("Average CoVaR Loss Difference") +
  theme(axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16))
p

ggsave("application/output/FC_eval_JPMSP500_AS.pdf", p, width=15, height=12, units="cm")




