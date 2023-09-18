###################################################################
data.adj.two <- function(data,tau){
  nt <- length(table(data$trial))
  ntrt <- length(table(data$trt))
  weights <- cen.weight(data,ntrt=ntrt,nt=nt,t=tau)
  data_w <- merge(data, weights, by = "id")
  data_w$Y <- pmin(data_w$time,tau)
  nt <- length(table(data_w$trial))
  rmst_m <- S <- NULL
  for (j in 1:nt){
    dataj <- data_w[data_w$trial==j,]
    fit <- glm(Y ~ -1 + trt1 + trt2 + trt3 + trt1:X + trt2:X + trt3:X, 
              family = quasipoisson(), weights = weight, data = dataj)
    rmst_m <- rbind(rmst_m, coef(fit))
    vcov_fit <- vcov(fit)
    S <- rbind(S, as.vector(vcov_fit[lower.tri(vcov_fit,diag = T)]))
  }
  return(list(rmst_m = rmst_m, S = S))
}



sim.adj.two <- function(S,data.list,ntrt,nt,tau){
  beta.list <- cov_beta.list <- sigmab.list <- NULL
  beta0.list <- beta1.list <- cov_beta0.list <- cov_beta1.list <- NULL
  
  for (s in 1:S){
    data <- data.list[[s]]
    
    mvmeta_df <- data.adj.two(data, tau)
    skip_to_next <- FALSE
    tryCatch(fit <- mvmeta(mvmeta_df$rmst_m, S = mvmeta_df$S, 
                           method = "reml", bscov = "unstr", control = list(maxiter=1000)),
             error = function(e){skip_to_next <<- TRUE})
    if (skip_to_next) {next} 
    
    
    
    beta <- fit$coefficients
    vcov_beta <- fit$vcov
    M <- matrix(c(1,0,0,1,0,0,
                  0,1,0,0,1,0,
                  0,0,1,0,0,1), nrow = 3, ncol = 6, byrow = T)
    beta.list <- rbind(beta.list, beta)
    cov_beta.list <- rbind(cov_beta.list, sqrt(diag(vcov_beta)))
    sigmab.list <- rbind(sigmab.list, sqrt(diag(fit$Psi)))
    beta0.list <- rbind(beta0.list, beta[1:ntrt])
    beta1.list <- rbind(beta1.list, beta[1:ntrt] + beta[(ntrt+1):(2*ntrt)])
    cov_beta0.list <- rbind(cov_beta0.list, sqrt(diag(vcov_beta)[1:ntrt]))
    cov_beta1.list <- rbind(cov_beta1.list, sqrt(diag(M %*% vcov_beta %*% t(M))))
    
    if (!s%%500) print(s)
    #print(s)
  }
  
  
  return(list(beta = beta.list,
              cov_beta = cov_beta.list,
              sigmab = sigmab.list,
              beta0 = beta0.list,
              beta1 = beta1.list,
              cov_beta0 = cov_beta0.list,
              cov_beta1 = cov_beta1.list))
}






res.output.adj.two <- function(output,ntrt,fit0){
  tab <- NULL

  #true_rmst0 <- exp(fit0$coefficients[1:ntrt])
  #true_rmst1 <- exp(fit0$coefficients[1:ntrt] + fit0$coefficients[(ntrt+1):(2*ntrt)])
  true_rmst0 <- fit0$coefficients[1:ntrt]
  true_rmst1 <- fit0$coefficients[1:ntrt] + fit0$coefficients[(ntrt+1):(2*ntrt)]
  
  for (i in 1:length(output)){
    res <- output[[i]]
    
    #rmst0 = colMeans(exp(res$beta0))
    #rmst1 = colMeans(exp(res$beta1))
    rmst0 = colMeans(res$beta0)
    rmst1 = colMeans(res$beta1)
    rmst0_bias = round(rmst0 - true_rmst0,3)
    rmst1_bias = round(rmst1 - true_rmst1,3)
    rmst0_cov = res$cov_beta0
    rmst1_cov = res$cov_beta1
    rmst0_sd = colMeans(rmst0_cov)
    rmst1_sd = colMeans(rmst1_cov)
    rmst0_se = apply(res$beta0,2,sd)
    rmst1_se = apply(res$beta1,2,sd)
    
    #rmst0_l = exp(res$beta0 - 1.96*rmst0_cov)
    #rmst1_l = exp(res$beta1 - 1.96*rmst1_cov)
    #rmst0_h = exp(res$beta0 + 1.96*rmst0_cov)
    #rmst1_h = exp(res$beta1 + 1.96*rmst1_cov)
    rmst0_l = res$beta0 - 1.96*rmst0_cov
    rmst1_l = res$beta1 - 1.96*rmst1_cov
    rmst0_h = res$beta0 + 1.96*rmst0_cov
    rmst1_h = res$beta1 + 1.96*rmst1_cov
    
    cover.p <- function(beta_l,beta_h,true_beta){
      cover <- (beta_l <= c(true_beta))*(beta_h >= c(true_beta))
      #cp <- sum(cover == 1,na.rm = T)/length(cover)
      cp <- mean(cover, na.rm = T)
      round(cp*100,1)
    }
    
    CP0 <- CP1 <- NULL
    for (j in 1:ntrt){
      CP0 <- c(CP0, cover.p(rmst0_l[,j], rmst0_h[,j], true_rmst0[j]))
      CP1 <- c(CP1, cover.p(rmst1_l[,j], rmst1_h[,j], true_rmst1[j]))
    }
    
    
    #MSE0 <- round(rowMeans((t(exp(res$beta0)) - true_rmst0)^2),3)
    #MSE1 <- round(rowMeans((t(exp(res$beta1)) - true_rmst1)^2),3)
    MSE0 <- round(rowMeans((t(res$beta0) - true_rmst0)^2),3)
    MSE1 <- round(rowMeans((t(res$beta1) - true_rmst1)^2),3)
    
    tab <- rbind(tab, c(round(rmst0,3), rmst0_bias, 
                        round(rmst0_sd,3), round(rmst0_se,3),
                        CP0, MSE0,
                        round(rmst1,3), rmst1_bias, 
                        round(rmst1_sd,3), round(rmst1_se,3),
                        CP1, MSE1,
                        round(colMeans(res$sigmab),3)))
  }
  
  colnames(tab) <- c(paste0("rmst0_",1:ntrt), paste0("rmst0_",1:ntrt," bias"),
                     paste0("rmst0_",1:ntrt," sd"), paste0("rmst0_",1:ntrt," se"),
                     paste0("rmst0_",1:ntrt," CP"), paste0("rmst0_",1:ntrt," MSE"),
                     paste0("rmst1_",1:ntrt), paste0("rmst1_",1:ntrt," bias"),
                     paste0("rmst1_",1:ntrt," sd"), paste0("rmst1_",1:ntrt," se"),
                     paste0("rmst1_",1:ntrt," CP"), paste0("rmst1_",1:ntrt," MSE"),
                     paste0("trt",1:ntrt," tau"), paste0("trt",1:ntrt,"_X tau"))
  return(as.data.frame(tab))
}
