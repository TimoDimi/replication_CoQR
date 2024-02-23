###  This file runs the respective files generating the results for the Simulations
###  and the Applications

### Possibly install the SystemicRisk (and devtools) packages
# install.packages("devtools")
# devtools::install_github("TimoDimi/SystemicRisk")


### Section 3: Dynamic CoCAViaR Simulations
source("simulations/sim_CoCAViaR.R")
source("simulations/sim_CoCAViaR_eval.R")


### Section 4: Forecasting application with dynamic CoCAViaR models

# Important note: The following files only run with the publicly unavailable
# data files mentioned in the README file.

source("application/FC_application_ModelFit.R")
source("application/FC_application_FCEvaluation.R")
source("application/FC_application_ModelParameters.R")



