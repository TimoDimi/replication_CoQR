###  This file runs the respective files generating the results for the Simulations
###  and the Applications


### Simulations

# Section 3.1: Predictive CoQR simulations
source("simulations/sim_CoQR.R")
source("simulations/sim_CoQR_eval.R")

# Section 3.2: Dynamic CoCAViaR simulations
source("simulations/sim_CoCAViaR.R")
source("simulations/sim_CoCAViaR_eval.R")


### Applications

# Important note: The following files only run with the publicly unavailable
# data files mentioned in the README file.

# Section 4.1: Predictive CoQR application
source("application/CS_application.R")

# Section 4.2: Forecasting application with dynamic CoCAViaR models
source("application/FC_application_ModelParameters.R")
source("application/FC_application_FCEvaluation.R")
source("application/FC_application_ModelParameters.R")



