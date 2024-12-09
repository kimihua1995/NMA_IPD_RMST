# This R script generates all tables from the simulation study, including:
# Web Tables 3 to 9 in the Supplementary Material

library(Matrix)
library(writexl)
library(zoo)
library(MASS)
library(lme4)
library(nlme)
library(survival)
library(survminer)
library(R2jags)
library(mvmeta)
library(ggplot2)
library(gridExtra)
library(xtable)


# input source files
source("source_dat_sim_aft.R")
source("source_adj_two_stage.R")
source("source_adj_one_stage.R")
source("source_NPF_two_stage.R")
source("source_NPB_two_stage.R")



# load cleaned-up results from "simulation results.R"
load("sim_res.RData")



# function for generating tables
sim.table <- function(sim,subgroup){
  two_X0 <- sim$two[,1:18]
  two_X1 <- sim$two[,19:36]
  two_tau <- sim$two[,37:42]
  one_X0 <- sim$one[,1:18]
  one_X1 <- sim$one[,19:36]
  one_tau <- sim$one[,37:42]
  NPF_X0 <- sim$NPF[,1:18]
  NPF_X1 <- sim$NPF[,19:36]
  NPB_X0 <- sim$NPB[,1:18]
  NPB_X1 <- sim$NPB[,19:36]
  tab_X0 <- rbind(two_X0, one_X0, NPF_X0, NPB_X0)
  tab_X1 <- rbind(two_X1, one_X1, NPF_X1, NPB_X1)
  tab_X0 <- cbind(c(subgroup,rep("",11)),
                  c("Two","","","One","","","NPF","","","NPB","",""),
                  rep(c(0.1,0.3,0.5),4),tab_X0)
  tab_X1 <- cbind(c(subgroup,rep("",11)),
                  c("Two","","","One","","","NPF","","","NPB","",""),
                  rep(c(0.1,0.3,0.5),4),tab_X1)
  colnames(tab_X0) <- colnames(tab_X1) <- c("","Method","tau",rep(c("TrtA","TrtB","TrtC"),6))
  #rownames(tab_X0) <- rownames(tab_X1) <- NULL
  tab_tau <- rbind(two_tau, one_tau)
  tab_tau <- cbind(c(subgroup,rep("",5)),c("Two","","","One","",""),rep(c(0.1,0.3,0.5),2), tab_tau)
  colnames(tab_tau) <- c("","Method","tau","a1","a2","a3","b1","b2","b3")
  #rownames(tab_tau) <- NULL
  return(list(tab_X0 = tab_X0, tab_X1 = tab_X1, tab_tau = tab_tau))
}




# Scenario 1
sim_tab1_1 <- sim.table(sim1_1, "n=200")
sim_tab1_2 <- sim.table(sim1_2, "n=500")
sim_tab1_3 <- sim.table(sim1_3, "n=1000")

# Scenario 2
sim_tab2_1 <- sim.table(sim2_1, "nt=10")
sim_tab2_2 <- sim.table(sim2_2, "nt=20")
sim_tab2_3 <- sim.table(sim2_3, "nt=30")


# Scenario 3
sim_tab3_1 <- sim.table(sim3_1, "Network 1")
sim_tab3_2 <- sim.table(sim3_2, "Network 2")
sim_tab3_3 <- sim.table(sim3_3, "Network 3")





# Web Table 3: Scenario 1 with subgroup x=0
print(xtable(sim_tab1_1$tab_X0, type = "latex"), include.rownames=FALSE)
print(xtable(sim_tab1_2$tab_X0, type = "latex"), include.rownames=FALSE)
print(xtable(sim_tab1_3$tab_X0, type = "latex"), include.rownames=FALSE)


# Web Table 4: Scenario 2 with subgroup x=0
print(xtable(sim_tab2_1$tab_X0, type = "latex"), include.rownames=FALSE)
print(xtable(sim_tab2_2$tab_X0, type = "latex"), include.rownames=FALSE)
print(xtable(sim_tab2_3$tab_X0, type = "latex"), include.rownames=FALSE)



# Web Table 5: Scenario 3 with subgroup x=0
print(xtable(sim_tab3_1$tab_X0, type = "latex"), include.rownames=FALSE)
print(xtable(sim_tab3_2$tab_X0, type = "latex"), include.rownames=FALSE)
print(xtable(sim_tab3_3$tab_X0, type = "latex"), include.rownames=FALSE)


# Web Table 6: Scenario 1 with subgroup x=1
print(xtable(sim_tab1_1$tab_X1, type = "latex"), include.rownames=FALSE)
print(xtable(sim_tab1_2$tab_X1, type = "latex"), include.rownames=FALSE)
print(xtable(sim_tab1_3$tab_X1, type = "latex"), include.rownames=FALSE)



# Web Table 7: Scenario 2 with subgroup x=1
print(xtable(sim_tab2_1$tab_X1, type = "latex"), include.rownames=FALSE)
print(xtable(sim_tab2_2$tab_X1, type = "latex"), include.rownames=FALSE)
print(xtable(sim_tab2_3$tab_X1, type = "latex"), include.rownames=FALSE)


# Web Table 8: Scenario 3 with subgroup x=1
print(xtable(sim_tab3_1$tab_X1, type = "latex"), include.rownames=FALSE)
print(xtable(sim_tab3_2$tab_X1, type = "latex"), include.rownames=FALSE)
print(xtable(sim_tab3_3$tab_X1, type = "latex"), include.rownames=FALSE)



# Web Table 9: between-study heterogeneity
sim_tab_tau <- rbind(sim_tab1_1$tab_tau,sim_tab1_2$tab_tau,sim_tab1_3$tab_tau,
                     sim_tab2_1$tab_tau,sim_tab2_2$tab_tau,sim_tab2_3$tab_tau,
                     sim_tab3_1$tab_tau,sim_tab3_2$tab_tau,sim_tab3_3$tab_tau)
print(xtable(sim_tab_tau, type = "latex"), include.rownames=FALSE)

