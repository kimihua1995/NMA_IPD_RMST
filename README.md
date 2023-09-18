# NMA_IPD_RMST
R code for the article "Network meta-analysis of time-to-event endpoints with individual participant data using restricted mean survival time regression"

simulation.R: This file contains all R codes for simulation study and figures generation

source files:
1) source_dat_sim_aft.R: source code for genearting data and true values for RMST estimands, using the method from Wang and Schaubel, Lifetime Data Analysis, 2018; Zhong and Schaubel, Biometrics, 2022
2) source_adj_two_stage.R: source code for proposed two-stage RMST IPD-NMA model
3) source_adj_one_stage.R: source code for proposed one-stage RMST IPD-NMA model
4) source_NPF_two_stage.R: source code for existing non-parametric frequentist model
5) source_NPB_two_stage.R: source code for existing non-parametric Bayesian model
6) source_plot.R: source code for creating the plots
