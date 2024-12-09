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


# function for the simulation study
sim.output <- function(seed, S, nt, n1, n2, alphas, a.list, betas, px, sigma,
                       tau, t_trt, Omega, df, method = c(1:4)){
  # method (all arm-based)
  # 1: Adjusted Two-Stage
  # 2: Adjusted One-Stage
  # 3: Non-Parametric Frequentist Two-Stage
  # 4: Non-Parametric Bayesian Two-Stage
  output_two <- output_one <- output_NPF <- output_NPB <- list()
  set.seed(seed)
  ntrt <- length(alphas)
  for (i in 1:length(a.list)){
    data.list <- list()
    for (s in 1:S){
      set.seed(s)
      data.list[[s]] <- dat.sim.aft(nt,n1,n2,alphas,a.list[i],betas,px,sigma,t_trt)
    }
    if (1 %in% method){
      output_two[[i]] <- sim.adj.two(S,data.list,ntrt,nt,tau)
    }else{output_two[[i]] <- NULL}
    if (2 %in% method){
      output_one[[i]] <- sim.adj.one(S,data.list,ntrt,nt,tau)
    }else{output_one[[i]] <- NULL}
    if (3 %in% method){
      output_NPF[[i]] <- sim.NPF(S,data.list,ntrt,nt,tau)
    }else{output_NPF[[i]] <- NULL}
    if (4 %in% method){
      output_NPB[[i]] <- sim.NPB(S,data.list,ntrt,nt,tau,Omega,df)
    }else{output_NPB[[i]] <- NULL}
    
  }
  
  
  return(list(two = output_two, one = output_one, 
              NPF = output_NPF, NPB = output_NPB))
}




# original results
alphas = c(0.5,1.5,1)
betas = c(0.3,0.5,0.7)
sigma = c(1,1.5,2)
a.list=c(.1,.3,.5)
tau=4
seed=12321
df=5
Omega=diag(1,3)
px=0.5
S=1000

## Senario 1
nt = 20
t_trt = rep(list(1:3),nt)
out1_1 <- sim.output(seed,S,nt,n1=200,n2=200,alphas,a.list,betas,px,sigma,tau,
                     t_trt,Omega,df,1:4)
out1_2 <- sim.output(seed,S,nt,n1=500,n2=500,alphas,a.list,betas,px,sigma,tau,
                     t_trt,Omega,df,1:4)
out1_3 <- sim.output(seed,S,nt,n1=1000,n2=1000,alphas,a.list,betas,px,sigma,tau,
                     t_trt,Omega,df,1:4)


## Senario 2
n1 = n2 = 500
out2_1 <- sim.output(seed,S,nt=10,n1,n2,alphas,a.list,betas,px,sigma,tau,
                     t_trt=rep(list(1:3),10),Omega,df,1:4)
out2_2 <- sim.output(seed,S,nt=20,n1,n2,alphas,a.list,betas,px,sigma,tau,
                     t_trt=rep(list(1:3),20),Omega,df,1:4)
out2_3 <- sim.output(seed,S,nt=30,n1,n2,alphas,a.list,betas,px,sigma,tau,
                     t_trt=rep(list(1:3),30),Omega,df,1:4)

## Senario 3
nt = 20
n1 = 300
n2 = 700
out3_1 <- sim.output(seed,S,nt,n1,n2,alphas,a.list,betas,px,sigma,tau,
                     t_trt=c(rep(list(1:3),14), rep(list(c(1,2)),2), rep(list(c(1,3)),2), rep(list(c(2,3)),2)),
                     Omega,df,1:4)
out3_2 <- sim.output(seed,S,nt,n1,n2,alphas,a.list,betas,px,sigma,tau,
                     t_trt=c(rep(list(1:3),8), rep(list(c(1,2)),4), rep(list(c(1,3)),4), rep(list(c(2,3)),4)),
                     Omega,df,1:4)
out3_3 <- sim.output(seed,S,nt,n1,n2,alphas,a.list,betas,px,sigma,tau,
                     t_trt=c(rep(list(c(1,2)),7), rep(list(c(1,3)),7), rep(list(c(2,3)),6)),
                     Omega,df,1:4)



save(out1_1,out1_2,out1_3,out2_1,out2_2,out2_3,
     out3_1,out3_2,out3_3,
     file = "sim_out.RData")
#load("sim_out.RData")




# clean-up the original results from the sim.output
sim.res <- function(output, ntrt, fit0, method){
  output_two <- output$two
  output_one <- output$one
  output_NPF <- output$NPF
  output_NPB <- output$NPB
  
  
  if (1 %in% method){
    res_two <- res.output.adj.two(output_two,ntrt,fit0)
  }else{res_two <- NULL}
  if (2 %in% method){
    res_one <- res.output.adj.one(output_one,ntrt,fit0)
  }else{res_one <- NULL}
  if (3 %in% method){
    res_NPF <- res.output.NPF(output_NPF,ntrt,fit0)
  }else{res_NPF <- NULL}
  if (4 %in% method){
    res_NPB <- res.output.NPB(output_NPB,ntrt,fit0)
  }else{res_NPB <- NULL}
  
  
  return(list(two = res_two,
              one = res_one, 
              NPF = res_NPF, 
              NPB = res_NPB))
}



fit0 <- true.coef(alphas,betas,px,sigma,tau)  # get true values of the parameters
ntrt=3
sim1_1 <- sim.res(out1_1, ntrt, fit0, 1:4)
sim1_2 <- sim.res(out1_2, ntrt, fit0, 1:4)
sim1_3 <- sim.res(out1_3, ntrt, fit0, 1:4)
sim2_1 <- sim.res(out2_1, ntrt, fit0, 1:4)
sim2_2 <- sim.res(out2_2, ntrt, fit0, 1:4)
sim2_3 <- sim.res(out2_3, ntrt, fit0, 1:4)
sim3_1 <- sim.res(out3_1, ntrt, fit0, 1:4)
sim3_2 <- sim.res(out3_2, ntrt, fit0, 1:4)
sim3_3 <- sim.res(out3_3, ntrt, fit0, 1:4)


save(sim1_1,sim1_2,sim1_3,sim2_1,sim2_2,sim2_3,sim3_1,sim3_2,sim3_3,
     file = "sim_res.RData")

#load("sim_res.RData")


