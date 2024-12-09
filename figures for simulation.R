# This R script generates all figures from the simulation study, including:
# Figures 1 to 3 in the main manuscript, 
# and Web Figures 1 to 4 in the Supplementary Material


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



# generate plot data
combine_result <- function(sim,scenario){
  comb_tab <- NULL
  method <- names(sim)
  for (i in 1:length(sim)){
    if (!is.null(sim[[i]])){
      temp <- data.frame(rmst0_bias = c(sim[[i]][,"rmst0_1 bias"],sim[[i]][,"rmst0_2 bias"],sim[[i]][,"rmst0_3 bias"]),
                         rmst1_bias = c(sim[[i]][,"rmst1_1 bias"],sim[[i]][,"rmst1_2 bias"],sim[[i]][,"rmst1_3 bias"]),
                         rmst0_CP = c(sim[[i]][,"rmst0_1 CP"],sim[[i]][,"rmst0_2 CP"],sim[[i]][,"rmst0_3 CP"]),
                         rmst1_CP = c(sim[[i]][,"rmst1_1 CP"],sim[[i]][,"rmst1_2 CP"],sim[[i]][,"rmst1_3 CP"]),
                         rmst0_MSE = c(sim[[i]][,"rmst0_1 MSE"],sim[[i]][,"rmst0_2 MSE"],sim[[i]][,"rmst0_3 MSE"]),
                         rmst1_MSE = c(sim[[i]][,"rmst1_1 MSE"],sim[[i]][,"rmst1_2 MSE"],sim[[i]][,"rmst1_3 MSE"]),
                         #alpha_tau = c(sim[[i]][,"alpha1 tau"],sim[[i]][,"alpha2 tau"],sim[[i]][,"alpha3 tau"]),
                         #beta_tau = c(sim[[i]][,"beta1 tau"],sim[[i]][,"beta2 tau"],sim[[i]][,"beta3 tau"]),
                         trt = rep(c("Trt A","Trt B","Trt C"), each = length(c(0.1,0.4,0.7))),
                         a = rep(c(0.1,0.3,0.5),3))
      temp$method <- method[i]
      comb_tab <- rbind(comb_tab, temp)
    }
  }
  comb_tab$method <- factor(comb_tab$method,
                            levels = c("two","one","NPF","NPB"),
                            labels = c("two","one","NPF","NPB"))
  comb_tab$a <- factor(comb_tab$a, levels = c(0.1,0.3,0.5),
                       labels = c("low","moderate","high"))
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

plot_rmst <- function(tab, title){
  p_bias0 <- ggplot(tab, aes(x = a, y = rmst0_bias, color = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red1", size = 0.3)+
    geom_point(aes(shape = method),size = 2, position = position_dodge(width = 0.6)) +
    coord_cartesian(ylim = c(-0.15,0.05)) +
    facet_grid(cols = vars(scenario), rows = vars(trt)) +
    #facet_wrap(trt ~ scenario, ncol = 2, nrow = 3) +
    theme_bw() +
    theme(legend.position="right",
          legend.box.margin = margin(-10,-10,-10,-10)) +
    scale_color_manual(name = "Method",
                       values = c("darkorange2","seagreen","darkslateblue","mediumorchid2"),
                       labels = c("Proposed Two-Stage","Proposed One-Stage",
                                  "NPF Two-Stage",
                                  "NPB Two-Stage"))+
    scale_shape_manual(name = "Method",
                       values = c(16,17,15,18),
                       labels = c("Proposed Two-Stage","Proposed One-Stage",
                                  "NPF Two-Stage",
                                  "NPB Two-Stage"))+
    scale_y_continuous(breaks = c(-0.15,-0.1,-0.05,0,0.05)) +
    labs(y = "Bias of log-RMST for x=0",
         x = "Between-Study Heterogeneity",
         title = title)
  
  p_bias1 <- ggplot(tab, aes(x = a, y = rmst1_bias, color = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red1", size = 0.3)+
    geom_point(aes(shape = method),size = 2, position = position_dodge(width = 0.6)) +
    coord_cartesian(ylim = c(-0.15,0.05)) +
    facet_grid(cols = vars(scenario), rows = vars(trt)) +
    #facet_wrap(trt ~ scenario, ncol = 2, nrow = 3) +
    theme_bw() +
    theme(legend.position="right",
          legend.box.margin = margin(-10,-10,-10,-10)) +
    scale_color_manual(name = "Method",
                       values = c("darkorange2","seagreen","darkslateblue","mediumorchid2"),
                       labels = c("Adjusted Two-Stage","Adjusted One-Stage",
                                  "Non-Parametric Frequentist Two-Stage",
                                  "Non-Parametric Bayesian Two-Stage"))+
    scale_shape_manual(name = "Method",
                       values = c(16,17,15,18),
                       labels = c("Adjusted Two-Stage","Adjusted One-Stage",
                                  "Non-Parametric Frequentist Two-Stage",
                                  "Non-Parametric Bayesian Two-Stage"))+
    scale_y_continuous(breaks = c(-0.15,-0.1,-0.05,0,0.05))+
    labs(y = "Bias of log-RMST for x=1",
         x = "Between-Study Heterogeneity",
         title = title)
  
  
  p_mse0 <- ggplot(tab, aes(x = a, y = rmst0_MSE, color = method)) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red1", size = 0.3)+
    geom_point(aes(shape = method),size = 2, position = position_dodge(width = 0.6)) +
    coord_cartesian(ylim = c(0,0.03)) +
    facet_grid(cols = vars(scenario), rows = vars(trt)) +
    #facet_wrap(trt ~ scenario, ncol = 2, nrow = 3) +
    theme_bw() +
    theme(legend.position="right",
          legend.box.margin = margin(-10,-10,-10,-10)) +
    scale_color_manual(name = "Method",
                       values = c("darkorange2","seagreen","darkslateblue","mediumorchid2"),
                       labels = c("Proposed Two-Stage","Proposed One-Stage",
                                  "NPF Two-Stage",
                                  "NPB Two-Stage"))+
    scale_shape_manual(name = "Method",
                       values = c(16,17,15,18),
                       labels = c("Proposed Two-Stage","Proposed One-Stage",
                                  "NPF Two-Stage",
                                  "NPB Two-Stage"))+
    scale_y_continuous(breaks = seq(0,0.03,by=0.005))+
    labs(y = "MSE of log-RMST for x=0",
         x = "Between-Study Heterogeneity",
         title = title)
  
  
  p_mse1 <- ggplot(tab, aes(x = a, y = rmst1_MSE, color = method)) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red1", size = 0.3)+
    geom_point(aes(shape = method),size = 2, position = position_dodge(width = 0.6)) +
    coord_cartesian(ylim = c(0,0.03)) +
    facet_grid(cols = vars(scenario), rows = vars(trt)) +
    #facet_wrap(trt ~ scenario, ncol = 2, nrow = 3) +
    theme_bw() +
    theme(legend.position="right",
          legend.box.margin = margin(-10,-10,-10,-10)) +
    scale_color_manual(name = "Method",
                       values = c("darkorange2","seagreen","darkslateblue","mediumorchid2"),
                       labels = c("Adjusted Two-Stage","Adjusted One-Stage",
                                  "Non-Parametric Frequentist Two-Stage",
                                  "Non-Parametric Bayesian Two-Stage"))+
    scale_shape_manual(name = "Method",
                       values = c(16,17,15,18),
                       labels = c("Adjusted Two-Stage","Adjusted One-Stage",
                                  "Non-Parametric Frequentist Two-Stage",
                                  "Non-Parametric Bayesian Two-Stage"))+
    scale_y_continuous(breaks = seq(0,0.03,by=0.005))+
    labs(y = "MSE of log-RMST for x=1",
         x = "Between-Study Heterogeneity",
         title = title)
  
  
  p_cp0 <- ggplot(tab, aes(x = a, y = rmst0_CP, color = method)) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red1", size = 0.3)+
    geom_point(aes(shape = method),size = 2, position = position_dodge(width = 0.6)) +
    coord_cartesian(ylim = c(50,100)) +
    facet_grid(cols = vars(scenario), rows = vars(trt)) +
    #facet_wrap(trt ~ scenario, ncol = 2, nrow = 3) +
    theme_bw() +
    theme(legend.position="right",
          legend.box.margin = margin(-10,-10,-10,-10)) +
    scale_color_manual(name = "Method",
                       values = c("darkorange2","seagreen","darkslateblue","mediumorchid2"),
                       labels = c("Proposed Two-Stage","Proposed One-Stage",
                                  "NPF Two-Stage",
                                  "NPB Two-Stage"))+
    scale_shape_manual(name = "Method",
                       values = c(16,17,15,18),
                       labels = c("Proposed Two-Stage","Proposed One-Stage",
                                  "NPF Two-Stage",
                                  "NPB Two-Stage"))+
    scale_y_continuous(breaks = seq(50,100,by=10))+
    labs(y = "Coverage Probability(%) of log-RMST for x=0",
         x = "Between-Study Heterogeneity",
         title = title)
  
  
  p_cp1 <- ggplot(tab, aes(x = a, y = rmst1_CP, color = method)) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red1", size = 0.3)+
    geom_point(aes(shape = method),size = 2, position = position_dodge(width = 0.6)) +
    coord_cartesian(ylim = c(50,100)) +
    facet_grid(cols = vars(scenario), rows = vars(trt)) +
    #facet_wrap(trt ~ scenario, ncol = 2, nrow = 3) +
    theme_bw() +
    theme(legend.position="right",
          legend.box.margin = margin(-10,-10,-10,-10)) +
    scale_color_manual(name = "Method",
                       values = c("darkorange2","seagreen","darkslateblue","mediumorchid2"),
                       labels = c("Adjusted Two-Stage","Adjusted One-Stage",
                                  "Non-Parametric Frequentist Two-Stage",
                                  "Non-Parametric Bayesian Two-Stage"))+
    scale_shape_manual(name = "Method",
                       values = c(16,17,15,18),
                       labels = c("Adjusted Two-Stage","Adjusted One-Stage",
                                  "Non-Parametric Frequentist Two-Stage",
                                  "Non-Parametric Bayesian Two-Stage"))+
    scale_y_continuous(breaks = seq(50,100,by=10))+
    labs(y = "Coverage Probability(%) of log-RMST for x=1",
         x = "Between-Study Heterogeneity",
         title = title)
  
  
  return(list(p_bias0 = p_bias0, p_bias1 = p_bias1, 
              p_cp0 = p_cp0, p_cp1 = p_cp1, 
              p_mse0 = p_mse0, p_mse1 = p_mse1))
}




tab1$rmst0_CP[tab1$rmst0_CP < 50] <- 50; tab1$rmst1_CP[tab1$rmst1_CP < 50] <- 50
tab2$rmst0_CP[tab2$rmst0_CP < 50] <- 50; tab2$rmst1_CP[tab2$rmst1_CP < 50] <- 50
tab3$rmst0_CP[tab3$rmst0_CP < 50] <- 50; tab3$rmst1_CP[tab3$rmst1_CP < 50] <- 50
plot_rmst1 <- plot_rmst(tab1, "Scenario 1")
plot_rmst2 <- plot_rmst(tab2, "Scenario 2")
plot_rmst3 <- plot_rmst(tab3, "Scenario 3")




# a function to get legend from the ggplot
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}




#######################################################
# Figures 1 to 3, and Web Figures 2 to 4:
#######################################################
pdf("Figure1.pdf", width = 10, height = 10)
grid.arrange(plot_rmst1$p_bias0+theme(legend.position="none"),
             plot_rmst2$p_bias0+theme(legend.position="none"),
             plot_rmst3$p_bias0+theme(legend.position="none"),
             get_legend(plot_rmst1$p_bias1),
             ncol = 2,heights = c(5,5),widths = c(5,5))
dev.off()


pdf("Figure2.pdf", width = 10, height = 10)
grid.arrange(plot_rmst1$p_mse0+theme(legend.position="none"),
             plot_rmst2$p_mse0+theme(legend.position="none"),
             plot_rmst3$p_mse0+theme(legend.position="none"),
             get_legend(plot_rmst1$p_mse1),
             ncol = 2,heights = c(5,5),widths = c(5,5))
dev.off()


pdf("Figure3.pdf", width = 10, height = 10)
grid.arrange(plot_rmst1$p_cp0+theme(legend.position="none"),
             plot_rmst2$p_cp0+theme(legend.position="none"),
             plot_rmst3$p_cp0+theme(legend.position="none"),
             get_legend(plot_rmst1$p_cp1),
             ncol = 2,heights = c(5,5),widths = c(5,5))
dev.off()




pdf("AP_Figure2.pdf", width = 10, height = 10)
grid.arrange(plot_rmst1$p_bias1+theme(legend.position="none"),
             plot_rmst2$p_bias1+theme(legend.position="none"),
             plot_rmst3$p_bias1+theme(legend.position="none"),
             get_legend(plot_rmst1$p_bias1),
             ncol = 2,heights = c(5,5),widths = c(5,5))
dev.off()


pdf("AP_Figure3.pdf", width = 10, height = 10)
grid.arrange(plot_rmst1$p_mse1+theme(legend.position="none"),
             plot_rmst2$p_mse1+theme(legend.position="none"),
             plot_rmst3$p_mse1+theme(legend.position="none"),
             get_legend(plot_rmst1$p_mse1),
             ncol = 2,heights = c(5,5),widths = c(5,5))
dev.off()


pdf("AP_Figure4.pdf", width = 10, height = 10)
grid.arrange(plot_rmst1$p_cp1+theme(legend.position="none"),
             plot_rmst2$p_cp1+theme(legend.position="none"),
             plot_rmst3$p_cp1+theme(legend.position="none"),
             get_legend(plot_rmst1$p_cp1),
             ncol = 2,heights = c(5,5),widths = c(5,5))
dev.off()







#######################################################
# Web Figure 1
#######################################################
set.seed(12321)
dat <- dat.sim.aft(nt=1,n1=10000,n2=10000,
                   alphas,a=0,betas,px,sigma,t_trt=list(1:3))
dat1 <- dat[dat$X == 1,]
dat0 <- dat[dat$X == 0,]
fit1 <- survfit(Surv(time,status) ~ trt,data = dat1)
fit0 <- survfit(Surv(time,status) ~ trt,data = dat0)
curve1 <- ggsurvplot(fit1,data = dat1,size = 1.5,
                     censor.size = 1,conf.int = F,
                     pval = F, #pval.coord = c(0,0.15),
                     xlim = c(0,5),
                     break.time.by = 1,
                     #xscale = 365.25/12,
                     xlab = "t",
                     ylab = "S(t)",
                     #ylim = c(0,0.12),
                     #break.y.by = 0.02,
                     title = "",
                     surv.plot.height = 5,
                     # risk table
                     risk.table = "nrisk_cumevents",
                     risk.table.fontsize = 5,
                     legend.labs = c("Trt A","Trt B","Trt C"),
                     palette = c("#E7B800", "#2E9FDF","red"),
                     cumevents = F,
                     legend.title = "Survival Curves for X=1",
                     tables.height = 0.2,
                     cumevents.height = 0.2,
                     fontsize = 2.0)
curve0 <- ggsurvplot(fit0,data = dat0,size = 1.5,
                     censor.size = 1,conf.int = F,
                     pval = F, #pval.coord = c(0,0.15),
                     xlim = c(0,5),
                     break.time.by = 1,
                     #xscale = 365.25/12,
                     xlab = "t",
                     ylab = "S(t)",
                     #ylim = c(0,0.12),
                     #break.y.by = 0.02,
                     title = "",
                     surv.plot.height = 5,
                     # risk table
                     risk.table = "nrisk_cumevents",
                     risk.table.fontsize = 5,
                     legend.labs = c("Trt A","Trt B","Trt C"),
                     palette = c("#E7B800", "#2E9FDF","red"),
                     cumevents = F,
                     legend.title = "Survival Curves for X=0",
                     tables.height = 0.2,
                     cumevents.height = 0.2,
                     fontsize = 2.0)

pdf("AP_Figure1.pdf", height = 7, width = 7)
grid.arrange(curve0$plot + geom_vline(xintercept=4, linetype="dashed", 
                                      color = "black", size=1),
             curve1$plot + geom_vline(xintercept=4, linetype="dashed", 
                                      color = "black", size=1),
             ncol = 1,heights = c(3.5,3.5),widths = 7)
dev.off()