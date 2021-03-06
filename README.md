
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Replication Material for the Paper “Dynamic Co-Quantile Regression”

<!-- badges: start -->
<!-- badges: end -->

The code in this repository generates the results of the simulations and
the empirical applications of the working paper Dimitriadis, T. and
Hoga, Y. (2022), Dyanmic Co-Quantile Regressions available on
[arXiv](https://arxiv.org/abs/2206.14275).

## Data Availability

We unfortunately do not have the licence to make the data publicly
available. The code of the empirical applications does of course not
work without these files. We still publish the code to make it available
for replication. It can be tested with publicly available or simulated
data. We will make sure to submit the original data files as
supplementary files to a journal in the peer-review process such that
referees can replicate our results. The simulations of course work
without any additional files.

## Simulations

The files for the simulations are available in the folder ‘simulations’.
The file ‘sim_CoQR.R’ generates the parameter estimates for the
predictive CoQR in Section 3.1, which are saved under
‘simulations/data’. It is recommended to run the file on a cluster on
20+ kernels. The file ‘sim_CoQR_eval.R’ evaluates these results and
generates the output in the ‘simulations/output’ folder.

Equivalently, the files ‘sim_CoCAViaR.R’ and ‘sim_CoCAViaR_eval.R’
generate the parameter estimates and forecast evaluation results for the
CoCAViaR models of Section 3.2. It is recommended to run the file on a
cluster on 20+ kernels.

## Applications

### 4.1 Predictive CoVaR Regression

The results of Section 4.1 are generated by the file
‘application/AB16reg_application.R’. Note again that it requires the
publicly unavailable file ‘data_appl_AB16.rds’.

### 4.2 CoVaR Forecasting with CoCAViaR Models

The predictions of Section 4.2 are generated by the file
‘application/FC_application_ModelFit.R’, which requires the publicly
unavailable financial returns in the file ‘data_Assets.rds’. The file
‘application/FC_application_FCEvaluation.R’ evaluates these predictions
and generates the output in the folder ‘application/output’. Finally,
the file ‘application/FC_application_ModelParameters.R’ generates the
parameter estimates for Table 5.
