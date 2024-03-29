
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Replication Material for the Paper “Dynamic CoVaR Modeling”

<!-- badges: start -->
<!-- badges: end -->

The code in this repository generates the results of the simulations and
the empirical applications of the working paper Dimitriadis, T. and
Hoga, Y. (2024), Dynamic CoVaR Modeling available on
[arXiv](https://arxiv.org/abs/2206.14275).

## Data Availability

We unfortunately do not have the licence to make the data publicly
available. The code of the empirical applications does of course not
work without these files. We still publish the code to make it available
for replication. It can be tested with publicly available or simulated
data. The simulations of course work without any additional files.

## Simulations

The files for the simulations are available in the folder ‘simulations’.
The file ‘sim_CoCAViaR.R’ generates the parameter estimates for the
CoCAViaR models in Section 3, which are saved under ‘simulations/data’.
It is recommended to run the file on a cluster on 10+ kernels. The file
‘sim_CoCAViaR_eval.R’ evaluates these results and generates the output
in the ‘simulations/output’ folder.

## Application: CoVaR Forecasting with CoCAViaR Models

The predictions of Section 4 are generated by the file
‘application/FC_application_ModelFit.R’, which requires the publicly
unavailable financial returns in the file ‘data_Assets.rds’. The file
‘application/FC_application_FCEvaluation.R’ evaluates these predictions
and generates the output in the folder ‘application/output’. Finally,
the file ‘application/FC_application_ModelParameters.R’ generates the
parameter estimates for Table 3.
