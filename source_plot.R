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


plot_rmst <- function(tab, title){
  p_bias0 <- ggplot(tab, aes(x = a, y = rmst0_bias, color = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red1", size = 0.3)+
    geom_point(aes(shape = method),size = 2) +
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
    scale_y_continuous(breaks = c(-0.15,-0.1,-0.05,0,0.05)) +
    labs(y = "Bias of log-RMST for x=0",
         x = "Between-Study Heterogeneity",
         title = title)
  
  p_bias1 <- ggplot(tab, aes(x = a, y = rmst1_bias, color = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red1", size = 0.3)+
    geom_point(aes(shape = method),size = 2) +
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
    geom_point(aes(shape = method),size = 2) +
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
    labs(y = "MSE of log-RMST for x=0",
         x = "Between-Study Heterogeneity",
         title = title)
  
  
  p_mse1 <- ggplot(tab, aes(x = a, y = rmst1_MSE, color = method)) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red1", size = 0.3)+
    geom_point(aes(shape = method),size = 2) +
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
    geom_point(aes(shape = method),size = 2) +
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
    labs(y = "Coverage Probability(%) of log-RMST for x=0",
         x = "Between-Study Heterogeneity",
         title = title)
  
  
  p_cp1 <- ggplot(tab, aes(x = a, y = rmst1_CP, color = method)) +
    geom_hline(yintercept = 95, linetype = "dashed", color = "red1", size = 0.3)+
    geom_point(aes(shape = method),size = 2) +
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




get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
