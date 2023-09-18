library(MASS)
library(lme4)
library(nlme)
library(survival)
library(R2jags)
library(mvmeta)
library(ggplot2)
library(gridExtra)
source("source_dat_sim_aft.R")
source("source_adj_two_stage.R")
source("source_adj_one_stage.R")
source("source_NPF_two_stage.R")
source("source_NPB_two_stage.R")
source("source_plot.R")




sim.output <- function(seed, S, nt, n1, n2, alphas, a.list, betas, px, 
                       sigma, tau, t_trt, Omega, df, method = c(1:4)){
  # method (all arm-based)
  # 1: Adjusted Two-Stage
  # 2: Adjusted One-Stage
  # 3: Non-Parametric Frequentist Two-Stage
  # 4: Non-Parametric Bayesian Two-Stage
  output_two <- output_one <- output_NPF <- output_NPB <- list()
  set.seed(seed)
  for (i in 1:length(a.list)){
    data.list <- list()
    for (s in 1:S){
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





out <- sim.output(seed = seed, S = S,
                  nt = nt, n1 = n1, n2 = n2,
                  alphas = alphas, a.list = a.list,
                  betas = betas, px = px, sigma = sigma,
                  tau = tau, t_trt = t_trt, 
                  Omega = Omega, df = df,
                  method = 1:4)





######################
# Simulation Results
######################

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

# Senario 1
nt = 20
t_trt = rep(list(1:3),nt)
out1_1 <- sim.output(seed,S,nt,n1=200,n2=200,alphas,a.list,betas,px,sigma,tau,
                     t_trt,Omega,df,1:4)
out1_2 <- sim.output(seed,S,nt,n1=500,n2=500,alphas,a.list,betas,px,sigma,tau,
                     t_trt,Omega,df,1:4)
out1_3 <- sim.output(seed,S,nt,n1=1000,n2=1000,alphas,a.list,betas,px,sigma,tau,
                     t_trt,Omega,df,1:4)


# Senario 2
n1 = n2 = 500
out2_1 <- sim.output(seed,S,nt=10,n1,n2,alphas,a.list,betas,px,sigma,tau,
                     t_trt=rep(list(1:3),10),Omega,df,1:4)
out2_2 <- sim.output(seed,S,nt=20,n1,n2,alphas,a.list,betas,px,sigma,tau,
                     t_trt=rep(list(1:3),20),Omega,df,1:4)
out2_3 <- sim.output(seed,S,nt=30,n1,n2,alphas,a.list,betas,px,sigma,tau,
                     t_trt=rep(list(1:3),30),Omega,df,1:4)

# Senario 3
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










#############################
# Clean-up the results
#############################
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


fit0 <- true.coef(alphas,betas,px,sigma,tau)
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















######################
# Creating the plot data
######################
tab1_1 <- combine_result(sim1_1, "n=200")
tab1_2 <- combine_result(sim1_2, "n=500")
tab1_3 <- combine_result(sim1_3, "n=1000")
tab1 <- rbind(tab1_1, tab1_2, tab1_3)
tab1$scenario <- factor(tab1$scenario,
                        levels = c("n=200","n=500","n=1000"),
                        labels = c("n=200","n=500","n=1000"))
tab2_1 <- combine_result(sim2_1, "nt=10")
tab2_2 <- combine_result(sim2_2, "nt=20")
tab2_3 <- combine_result(sim2_3, "nt=30")
tab2 <- rbind(tab2_1, tab2_2, tab2_3)
tab2$scenario <- factor(tab2$scenario,
                        levels = c("nt=10","nt=20","nt=30"),
                        labels = c("nt=10","nt=20","nt=30"))
tab3_1 <- combine_result(sim3_1, "Network1")
tab3_2 <- combine_result(sim3_2, "Network2")
tab3_3 <- combine_result(sim3_3, "Network3")
tab3 <- rbind(tab3_1, tab3_2, tab3_3)
tab3$scenario <- factor(tab3$scenario,
                        levels = c("Network1","Network2","Network3"),
                        labels = c("Network1","Network2","Network3"))


#######################
# Making the plot
#######################
tab1$rmst0_CP[tab1$rmst0_CP < 50] <- 50; tab1$rmst1_CP[tab1$rmst1_CP < 50] <- 50
tab2$rmst0_CP[tab2$rmst0_CP < 50] <- 50; tab2$rmst1_CP[tab2$rmst1_CP < 50] <- 50
tab3$rmst0_CP[tab3$rmst0_CP < 50] <- 50; tab3$rmst1_CP[tab3$rmst1_CP < 50] <- 50
plot_rmst1 <- plot_rmst(tab1, "Scenario 1")
plot_rmst2 <- plot_rmst(tab2, "Scenario 2")
plot_rmst3 <- plot_rmst(tab3, "Scenario 3")




pdf("Figure1.pdf", width = 10, height = 10)
grid.arrange(plot_rmst1$p_bias0+theme(legend.position="none"),
             plot_rmst2$p_bias0+theme(legend.position="none"),
             plot_rmst3$p_bias0+theme(legend.position="none"),
             get_legend(plot_rmst1$p_bias0),
             ncol = 2,heights = c(5,5),widths = c(5,5))
dev.off()


pdf("Figure2.pdf", width = 10, height = 10)
grid.arrange(plot_rmst1$p_mse0+theme(legend.position="none"),
             plot_rmst2$p_mse0+theme(legend.position="none"),
             plot_rmst3$p_mse0+theme(legend.position="none"),
             get_legend(plot_rmst1$p_mse0),
             ncol = 2,heights = c(5,5),widths = c(5,5))
dev.off()


pdf("Figure3.pdf", width = 10, height = 10)
grid.arrange(plot_rmst1$p_cp0+theme(legend.position="none"),
             plot_rmst2$p_cp0+theme(legend.position="none"),
             plot_rmst3$p_cp0+theme(legend.position="none"),
             get_legend(plot_rmst1$p_cp0),
             ncol = 2,heights = c(5,5),widths = c(5,5))
dev.off()




pdf("eFigure2.pdf", width = 10, height = 10)
grid.arrange(plot_rmst1$p_bias1+theme(legend.position="none"),
             plot_rmst2$p_bias1+theme(legend.position="none"),
             plot_rmst3$p_bias1+theme(legend.position="none"),
             get_legend(plot_rmst1$p_bias1),
             ncol = 2,heights = c(5,5),widths = c(5,5))
dev.off()


pdf("eFigure3.pdf", width = 10, height = 10)
grid.arrange(plot_rmst1$p_mse1+theme(legend.position="none"),
             plot_rmst2$p_mse1+theme(legend.position="none"),
             plot_rmst3$p_mse1+theme(legend.position="none"),
             get_legend(plot_rmst1$p_mse1),
             ncol = 2,heights = c(5,5),widths = c(5,5))
dev.off()


pdf("eFigure4.pdf", width = 10, height = 10)
grid.arrange(plot_rmst1$p_cp1+theme(legend.position="none"),
             plot_rmst2$p_cp1+theme(legend.position="none"),
             plot_rmst3$p_cp1+theme(legend.position="none"),
             get_legend(plot_rmst1$p_cp1),
             ncol = 2,heights = c(5,5),widths = c(5,5))
dev.off()







