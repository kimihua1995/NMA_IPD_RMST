library(Matrix)
library(writexl)
library(zoo)
library(MASS)
library(lme4)
library(nlme)
library(survival)
library(R2jags)
library(mvmeta)
library(ggplot2)

source("source_dat_sim.R")
source("source_adj_two.R")
source("source_adj_one.R")
source("source_NPF.R")
source("source_NPB.R")



sim.output <- function(seed, nt, n1, n2, alphas, betas, px, a.list, lambda, gamma, 
                       tcen, tau, t_trt, Omega, df, method = c(1:4)){
  # method (all arm-based)
  # 1: Covariate-Adjusted Two-Stage
  # 2: Covariate-Adjusted One-Stage
  # 3: Non-Parametric Frequentist
  # 4: Non-Parametric Bayesian
  output_glmm <- output_mvmeta <- output_bayesian_arm <- 
    output_mvmeta2 <- list()
  for (i in 1:length(a.list)){
    if (1 %in% method){
      output_mvmeta2[[i]] <- sim.RMST.mvmeta.adj(seed, nt, n1, n2, alphas, betas, px,
                                                 a.list[i], lambda, gamma, tcen, t_trt, tau)
    }else{output_mvmeta2[[i]] <- NULL}
    if (2 %in% method){
      output_glmm[[i]] <- sim.RMST.GLMM.adj(seed, nt, n1, n2, alphas, betas, px,
                                            a.list[i],lambda, gamma, tcen, t_trt, tau)
    }else{output_glmm[[i]] <- NULL}
    if (3 %in% method){
      output_mvmeta[[i]] <- sim.RMST.mvmeta(seed, nt, n1, n2, alphas, betas, px,
                                            a.list[i],lambda, gamma, tcen, t_trt, tau)
    }else{output_mvmeta[[i]] <- NULL}
    if (4 %in% method){
      output_bayesian_arm[[i]] <- sim.RMST.bayesian.arm(seed, nt, n1, n2, alphas, betas, px,
                                          a.list[i], lambda, gamma, tcen, t_trt, tau, Omega, df)
    }else{output_bayesian_arm[[i]] <- NULL}
    
  }
  

  return(list(mvmeta2 = output_mvmeta2, glmm = output_glmm, 
              mvmeta = output_mvmeta, bayesian_arm = output_bayesian_arm))
}



# parameters
a.list=c(.1,.4,.7)
alphas=c(0.3,-1,-0.5)
lambda=c(0.3,0.4,0.2)
gamma=c(1,0.5,1.2)
tcen=10
tau=5
seed=12321
df <- 5
Omega <- diag(1,3)
betas <- c(-0.3,-0.3,-0.3)
px <- 0.5


# Senario 1
nt = 20
t_trt = rep(list(1:3),nt)
out1_1 <- sim.output(seed, nt, n1=200, n2=200, alphas, betas, px, a.list, lambda, gamma, tcen, tau, t_trt, Omega, df, 1:4)
out1_2 <- sim.output(seed, nt, n1=500, n2=500, alphas, betas, px, a.list, lambda, gamma, tcen, tau, t_trt, Omega, df, 1:4)
out1_3 <- sim.output(seed, nt, n1=1000, n2=1000, alphas, betas, px, a.list, lambda, gamma, tcen, tau, t_trt, Omega, df, 1:4)




# Senario 2
n1 = n2 = 500
out2_1 <- sim.output(seed, nt=10, n1, n2, alphas, betas, px, a.list, lambda, gamma, tcen, tau, 
                      t_trt = rep(list(1:3),10), Omega, df, 1:4)
out2_2 <- sim.output(seed, nt=20, n1, n2, alphas, betas, px, a.list, lambda, gamma, tcen, tau, 
                      t_trt = rep(list(1:3),20), Omega, df, 1:4)
out2_3 <- sim.output(seed, nt=30, n1, n2, alphas, betas, px, a.list, lambda, gamma, tcen, tau, 
                      t_trt = rep(list(1:3),30), Omega, df, 1:4)



# Senario 3
nt = 20
n1 = 300
n2 = 700
out3_1 <- sim.output(seed, nt, n1, n2, alphas, betas, px, a.list, lambda, gamma, tcen, tau, 
                      t_trt = c(rep(list(1:3),14), rep(list(c(1,2)),2), rep(list(c(1,3)),2), rep(list(c(2,3)),2)),
                     Omega, df, 1:4)
out3_2 <- sim.output(seed, nt, n1, n2, alphas, betas, px, a.list, lambda, gamma, tcen, tau, 
                      t_trt = c(rep(list(1:3),8), rep(list(c(1,2)),4), rep(list(c(1,3)),4), rep(list(c(2,3)),4)),
                      Omega, df, 1:4)
out3_3 <- sim.output(seed, nt, n1, n2, alphas, betas, px, a.list, lambda, gamma, tcen, tau, 
                      t_trt = c(rep(list(c(1,2)),7), rep(list(c(1,3)),7), rep(list(c(2,3)),6)),
                      Omega, df, 1:4)




# output results from the sim.output
sim.res <- function(output, alphas, betas, a.list, gamma, lambda, tau, method){
  output_mvmeta2 <- output$mvmeta2
  output_glmm <- output$glmm
  output_mvmeta <- output$mvmeta
  output_bayesian_arm <- output$bayesian_arm
  
  
  if (1 %in% method){
    res_mvmeta2 <- res.output.RMST.adj(output_mvmeta2,alphas,betas,a.list,gamma,lambda,tau)
  }else{res_mvmeta2 <- NULL}
  if (2 %in% method){
    res_glmm <- res.output.RMST.adj(output_glmm,alphas,betas,a.list,gamma,lambda,tau)
  }else{res_glmm <- NULL}
  if (3 %in% method){
    res_mvmeta <- res.output.RMST.adj(output_mvmeta,alphas,betas,a.list,gamma,lambda,tau)
  }else{res_mvmeta <- NULL}
  if (4 %in% method){
    res_bayesian_arm <- res.output.RMST.adj(output_bayesian_arm,alphas,betas,a.list,gamma,lambda,tau)
  }else{res_bayesian_arm <- NULL}
  
  
  return(list(mvmeta2 = res_mvmeta2,
              glmm = res_glmm, 
              mvmeta = res_mvmeta, 
              bayesian_arm = res_bayesian_arm))
}



sim1_1 <- sim.res(out1_1, alphas, betas, a.list, gamma, lambda, tau, 1:4)
sim1_2 <- sim.res(out1_2, alphas, betas, a.list, gamma, lambda, tau, 1:4)
sim1_3 <- sim.res(out1_3, alphas, betas, a.list, gamma, lambda, tau, 1:4)
sim2_1 <- sim.res(out2_1, alphas, betas, a.list, gamma, lambda, tau, 1:4)
sim2_2 <- sim.res(out2_2, alphas, betas, a.list, gamma, lambda, tau, 1:4)
sim2_3 <- sim.res(out2_3, alphas, betas, a.list, gamma, lambda, tau, 1:4)
sim3_1 <- sim.res(out3_1, alphas, betas, a.list, gamma, lambda, tau, 1:4)
sim3_2 <- sim.res(out3_2, alphas, betas, a.list, gamma, lambda, tau, 1:4)
sim3_3 <- sim.res(out3_3, alphas, betas, a.list, gamma, lambda, tau, 1:4)












##########################################################################################
##########################################################################################
# plot for RMST
combine_result <- function(sim,scenario){
  comb_tab <- NULL
  method <- names(sim)
  for (i in 1:length(sim)){
    if (!is.null(sim[[i]])){
      temp <- data.frame(alpha_bias = c(sim[[i]][,"alpha1 bias"],sim[[i]][,"alpha2 bias"],sim[[i]][,"alpha3 bias"]),
                         beta_bias = c(sim[[i]][,"beta1 bias"],sim[[i]][,"beta2 bias"],sim[[i]][,"beta3 bias"]),
                         alpha_CP = c(sim[[i]][,"alpha1 CP"],sim[[i]][,"alpha2 CP"],sim[[i]][,"alpha3 CP"]),
                         beta_CP = c(sim[[i]][,"beta1 CP"],sim[[i]][,"beta2 CP"],sim[[i]][,"beta3 CP"]),
                         alpha_MSE = c(sim[[i]][,"alpha1 MSE"],sim[[i]][,"alpha2 MSE"],sim[[i]][,"alpha3 MSE"]),
                         beta_MSE = c(sim[[i]][,"beta1 MSE"],sim[[i]][,"beta2 MSE"],sim[[i]][,"beta3 MSE"]),
                         alpha_tau = c(sim[[i]][,"alpha1 tau"],sim[[i]][,"alpha2 tau"],sim[[i]][,"alpha3 tau"]),
                         beta_tau = c(sim[[i]][,"beta1 tau"],sim[[i]][,"beta2 tau"],sim[[i]][,"beta3 tau"]),
                         trt = rep(c("A","B","C"), each = length(a.list)),
                         a = rep(a.list,3))
      temp$method <- method[i]
      comb_tab <- rbind(comb_tab, temp)
    }
  }
  comb_tab$method <- factor(comb_tab$method,
                            levels = c("mvmeta2","glmm","mvmeta","bayesian_arm"),
                            labels = c("mvmeta2","glmm","mvmeta","bayesian_arm"))
  comb_tab$a <- factor(comb_tab$a, levels = c(0.1,0.4,0.7),
                       labels = c("low","median","high"))
  comb_tab$scenario <- scenario
  return(comb_tab)
}

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

plot_rmst <- function(tab){
  p_bias1 <- ggplot(tab, aes(x = a, y = abs(alpha_bias), color = method)) +
    #geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.3)+
    geom_point(aes(shape = method),size = 2) +
    coord_cartesian(ylim = c(0,0.3)) +
    facet_grid(cols = vars(scenario), rows = vars(trt)) +
    #facet_wrap(trt ~ scenario, ncol = 2, nrow = 3) +
    theme_bw() +
    theme(legend.position="bottom",
          legend.box.margin = margin(-10,-10,-10,-10)) +
    scale_color_manual(name = "Method",
                       values = c("darkorange2","seagreen","darkslateblue","mediumorchid2"),
                       labels = c("E-Two-Stage","One-Stage","F-Two-Stage","B-Two-Stage"))+
    scale_shape_manual(name = "Method",
                       values = c(16,17,15,18),
                       labels = c("E-Two-Stage","One-Stage","F-Two-Stage","B-Two-Stage"))+
    scale_y_continuous(breaks = seq(0,0.3,by=0.05))+
    labs(y = "Absolute Bias of RMST (alpha)",
         x = "Between-Trial Heterogeneity")
  
  p_bias2 <- ggplot(tab, aes(x = a, y = abs(beta_bias), color = method)) +
    #geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.3)+
    geom_point(aes(shape = method),size = 2) +
    coord_cartesian(ylim = c(0,0.05)) +
    facet_grid(cols = vars(scenario), rows = vars(trt)) +
    #facet_wrap(trt ~ scenario, ncol = 2, nrow = 3) +
    theme_bw() +
    theme(legend.position="bottom",
          legend.box.margin = margin(-10,-10,-10,-10)) +
    scale_color_manual(name = "Method",
                       values = c("darkorange2","seagreen","darkslateblue","mediumorchid2"),
                       labels = c("E-Two-Stage","One-Stage","F-Two-Stage","B-Two-Stage"))+
    scale_shape_manual(name = "Method",
                       values = c(16,17,15,18),
                       labels = c("E-Two-Stage","One-Stage","F-Two-Stage","B-Two-Stage"))+
    scale_y_continuous(breaks = seq(0,0.05,by=0.01))+
    labs(y = "Absolute Bias of RMST (beta)",
         x = "Between-Trial Heterogeneity")
  
  p_cp1 <- ggplot(tab, aes(x = a, y = alpha_CP, color = method)) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red1", size = 0.3)+
    geom_point(aes(shape = method),size = 2) +
    coord_cartesian(ylim = c(50,100)) +
    facet_grid(cols = vars(scenario), rows = vars(trt)) +
    #facet_wrap(trt ~ scenario, ncol = 2, nrow = 3) +
    theme_bw() +
    theme(legend.position="bottom",
          legend.box.margin = margin(-10,-10,-10,-10)) +
    scale_color_manual(name = "Method",
                       values = c("darkorange2","seagreen","darkslateblue","mediumorchid2"),
                       labels = c("E-Two-Stage","One-Stage","F-Two-Stage","B-Two-Stage"))+
    scale_shape_manual(name = "Method",
                       values = c(16,17,15,18),
                       labels = c("E-Two-Stage","One-Stage","F-Two-Stage","B-Two-Stage"))+
    scale_y_continuous(breaks = seq(50,100,by=10))+
    labs(y = "CP(%) of RMST (alpha)",
         x = "Between-Trial Heterogeneity")
  
  
  p_cp2 <- ggplot(tab, aes(x = a, y = beta_CP, color = method)) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red1", size = 0.3)+
    geom_point(aes(shape = method),size = 2) +
    coord_cartesian(ylim = c(50,100)) +
    facet_grid(cols = vars(scenario), rows = vars(trt)) +
    #facet_wrap(trt ~ scenario, ncol = 2, nrow = 3) +
    theme_bw() +
    theme(legend.position="bottom",
          legend.box.margin = margin(-10,-10,-10,-10)) +
    scale_color_manual(name = "Method",
                       values = c("darkorange2","seagreen","darkslateblue","mediumorchid2"),
                       labels = c("E-Two-Stage","One-Stage","F-Two-Stage","B-Two-Stage"))+
    scale_shape_manual(name = "Method",
                       values = c(16,17,15,18),
                       labels = c("E-Two-Stage","One-Stage","F-Two-Stage","B-Two-Stage"))+
    scale_y_continuous(breaks = seq(50,100,by=10))+
    labs(y = "CP(%) of RMST (beta)",
         x = "Between-Trial Heterogeneity")
  
  p_mse1 <- ggplot(tab, aes(x = a, y = alpha_MSE, color = method)) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red1", size = 0.3)+
    geom_point(aes(shape = method),size = 2) +
    coord_cartesian(ylim = c(0,0.1)) +
    facet_grid(cols = vars(scenario), rows = vars(trt)) +
    #facet_wrap(trt ~ scenario, ncol = 2, nrow = 3) +
    theme_bw() +
    theme(legend.position="bottom",
          legend.box.margin = margin(-10,-10,-10,-10)) +
    scale_color_manual(name = "Method",
                       values = c("darkorange2","seagreen","darkslateblue","mediumorchid2"),
                       labels = c("E-Two-Stage","One-Stage","F-Two-Stage","B-Two-Stage"))+
    scale_shape_manual(name = "Method",
                       values = c(16,17,15,18),
                       labels = c("E-Two-Stage","One-Stage","F-Two-Stage","B-Two-Stage"))+
    scale_y_continuous(breaks = seq(0,0.1,by=0.02))+
    labs(y = "MSE of RMST (alpha)",
         x = "Between-Trial Heterogeneity")
  
  p_mse2 <- ggplot(tab, aes(x = a, y = beta_MSE, color = method)) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red1", size = 0.3)+
    geom_point(aes(shape = method),size = 2) +
    coord_cartesian(ylim = c(0,0.1)) +
    facet_grid(cols = vars(scenario), rows = vars(trt)) +
    #facet_wrap(trt ~ scenario, ncol = 2, nrow = 3) +
    theme_bw() +
    theme(legend.position="bottom",
          legend.box.margin = margin(-10,-10,-10,-10)) +
    scale_color_manual(name = "Method",
                       values = c("darkorange2","seagreen","darkslateblue","mediumorchid2"),
                       labels = c("E-Two-Stage","One-Stage","F-Two-Stage","B-Two-Stage"))+
    scale_shape_manual(name = "Method",
                       values = c(16,17,15,18),
                       labels = c("E-Two-Stage","One-Stage","F-Two-Stage","B-Two-Stage"))+
    scale_y_continuous(breaks = seq(0,0.1,by=0.02))+
    labs(y = "MSE of RMST (beta)",
         x = "Between-Trial Heterogeneity")
  
  p_tau1 <- ggplot(tab, aes(x = a, y = alpha_tau, color = method)) +
    #geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.3)+
    geom_point(aes(shape = method),size = 2) +
    coord_cartesian(ylim = c(0,1)) +
    facet_grid(cols = vars(scenario), rows = vars(trt)) +
    #facet_wrap(trt ~ scenario, ncol = 2, nrow = 3) +
    theme_bw() +
    theme(legend.position="bottom",
          legend.box.margin = margin(-10,-10,-10,-10)) +
    scale_color_manual(name = "Method",
                       values = c("darkorange2","seagreen","darkslateblue","mediumorchid2"),
                       labels = c("E-Two-Stage","One-Stage","F-Two-Stage","B-Two-Stage"))+
    scale_shape_manual(name = "Method",
                       values = c(16,17,15,18),
                       labels = c("E-Two-Stage","One-Stage","F-Two-Stage","B-Two-Stage"))+
    scale_y_continuous(breaks = seq(0,1,by=0.2))+
    labs(y = "Tau (alpha)",
         x = "Between-Trial Heterogeneity")
  
  
  p_tau2 <- ggplot(tab, aes(x = a, y = beta_tau, color = method)) +
    #geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.3)+
    geom_point(aes(shape = method),size = 2) +
    coord_cartesian(ylim = c(0,2)) +
    facet_grid(cols = vars(scenario), rows = vars(trt)) +
    #facet_wrap(trt ~ scenario, ncol = 2, nrow = 3) +
    theme_bw() +
    theme(legend.position="bottom",
          legend.box.margin = margin(-10,-10,-10,-10)) +
    scale_color_manual(name = "Method",
                       values = c("darkorange2","seagreen","darkslateblue","mediumorchid2"),
                       labels = c("E-Two-Stage","One-Stage","F-Two-Stage","B-Two-Stage"))+
    scale_shape_manual(name = "Method",
                       values = c(16,17,15,18),
                       labels = c("E-Two-Stage","One-Stage","F-Two-Stage","B-Two-Stage"))+
    scale_y_continuous(breaks = seq(0,2,by=0.4))+
    labs(y = "Tau (beta)",
         x = "Between-Trial Heterogeneity")
  
  return(list(p_bias1 = p_bias1, p_bias2 = p_bias2, 
              p_cp1 = p_cp1, p_cp2 = p_cp2, 
              p_mse1 = p_mse1, p_mse2 = p_mse2, 
              p_tau1 = p_tau1, p_tau2 = p_tau2))
}




tab1$alpha_CP[tab1$alpha_CP < 50] <- 50
tab2$alpha_CP[tab2$alpha_CP < 50] <- 50
tab2$alpha_CP[tab3$alpha_CP < 50] <- 50
plot_rmst1 <- plot_rmst(tab1)
plot_rmst2 <- plot_rmst(tab2)
plot_rmst3 <- plot_rmst(tab3)


























