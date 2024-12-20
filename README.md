# NMA_IPD_RMST
R code for the article "Network meta-analysis of time-to-event endpoints with individual participant data using restricted mean survival time regression"

`simulation results.R`: This file contains R codes for generating simulation data and results (saved in `sim_out.RData` and `sim_res.RData`)

`figures for simulation.R`: This file contains R codes for generating all figures for the simulation study

`tables for simulation.R`: This file contains R codes for generating all tables for the simulation study

source files:
1) `source_dat_sim_aft.R`: source code for genearting data and true values for RMST estimands, using the method from Wang and Schaubel, Lifetime Data Analysis, 2018; Zhong and Schaubel, Biometrics, 2022
2) `source_adj_two_stage.R`: source code for proposed two-stage RMST IPD-NMA model
3) `source_adj_one_stage.R`: source code for proposed one-stage RMST IPD-NMA model
4) `source_NPF_two_stage.R`: source code for existing non-parametric frequentist model
5) `source_NPB_two_stage.R`: source code for existing non-parametric Bayesian model


