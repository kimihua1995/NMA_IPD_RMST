my.rmst=function(time, status, tau){
  ft <- survfit(Surv(time, status)~1)
  idx <- ft$time<=tau
  
  wk.time <- sort(c(ft$time[idx],min(tau,ft$time)))
  wk.surv <- ft$surv[idx]
  wk.n.risk  <- ft$n.risk[idx]
  wk.n.event <- ft$n.event[idx]
  
  time.diff <- diff(c(0, wk.time))
  areas <- time.diff * c(1, wk.surv)
  rmst <-  sum(areas)
  
  wk.var <- ifelse((wk.n.risk-wk.n.event)==0, 0,
                   wk.n.event /(wk.n.risk *(wk.n.risk - wk.n.event)))
  wk.var <- c(wk.var,0)
  rmst.var <- sum( cumsum(rev(areas[-1]))^2 * rev(wk.var)[-1])
  rmst.se  <- sqrt(rmst.var)
  
  return(list(rmst = rmst, se = rmst.se))
}


#########################################################
ad.rmst <- function(dat,tau){
  nt <- length(table(dat$trial))
  ntrt <- length(table(dat$trt))
  
  rmst_m <- rmst_sd <- matrix(NA, nrow = nt, ncol = ntrt)
  for (j in 1:nt){
    datj <- dat[dat$trial == j,]
    for (i in 1:ntrt){
      datij <- datj[datj$trt == i,]
      if (nrow(datij) > 0){
        res <- my.rmst(datij$time, datij$status, tau)
        rmst_m[j,i] <- res$rmst
        rmst_sd[j,i] <- res$se
      }
    }
  }
  
  return(list(rmst = rmst_m, se = rmst_sd))
}


###################################################################
dat.NPF <- function(dat, tau){
  res <- ad.rmst(dat,tau)
  rmst_m <- res$rmst
  rmst_sd <- res$se
  nt <- length(table(dat$trial))
  ntrt <- length(table(dat$trt))
  
  theta <- as.vector(t(rmst_m))
  S <- matrix(NA, nrow = nt, ncol = ntrt*(ntrt+1)/2)
  for (j in 1:nt){
    s=1
    for (i in 1:ntrt){
      for (l in i:ntrt){
        S[j,s] <- ifelse(i==l,rmst_sd[j,i]*rmst_sd[j,l],0)
        s <- s + 1
      }
    }
  }
  return(list(rmst_m = rmst_m, theta = theta, S = S))
}


sim.NPF <- function(S,data.list,ntrt,nt,tau){
  beta.list <- cov_beta.list <- sigmab.list <- NULL
  beta0.list <- beta1.list <- cov_beta0.list <- cov_beta1.list <- NULL

  for (s in 1:S){
    data <- data.list[[s]]
    
    mvmeta_df_0 <- dat.NPF(data[data$X==0,], tau)
    mvmeta_df_1 <- dat.NPF(data[data$X==1,], tau)
    skip_to_next <- FALSE
    tryCatch(fit0 <- mvmeta(mvmeta_df_0$rmst_m, S = mvmeta_df_0$S, 
                           method = "reml", bscov = "unstr", control = list(maxiter=1000)),
             error = function(e){skip_to_next <<- TRUE})
    tryCatch(fit1 <- mvmeta(mvmeta_df_1$rmst_m, S = mvmeta_df_1$S, 
                            method = "reml", bscov = "unstr", control = list(maxiter=1000)),
             error = function(e){skip_to_next <<- TRUE})
    if (skip_to_next) {next} 
    
    
    
    beta0 <- fit0$coefficients
    beta1 <- fit1$coefficients
    cov0 <- sqrt(diag(fit0$vcov))
    cov1 <- sqrt(diag(fit1$vcov))
    sigmab0 <- sqrt(diag(fit0$Psi))
    sigmab1 <- sqrt(diag(fit1$Psi))
    
    beta <- c(beta0, beta1 - beta0)
    cov_beta <- c(cov0, sqrt(cov0^2 + cov1^2))
    sigmab <- c(sigmab0, sigmab1)
    
    beta.list <- rbind(beta.list, beta)
    cov_beta.list <- rbind(cov_beta.list, cov_beta)
    sigmab.list <- rbind(sigmab.list, sigmab)
    beta0.list <- rbind(beta0.list, beta0)
    beta1.list <- rbind(beta1.list, beta1)
    cov_beta0.list <- rbind(cov_beta0.list, cov0)
    cov_beta1.list <- rbind(cov_beta1.list, cov1)
    
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




res.output.NPF <- function(output,ntrt,fit0){
  tab <- NULL
  
  #true_rmst0 <- exp(fit0$coefficients[1:ntrt])
  #true_rmst1 <- exp(fit0$coefficients[1:ntrt] + fit0$coefficients[(ntrt+1):(2*ntrt)])
  true_rmst0 <- fit0$coefficients[1:ntrt]
  true_rmst1 <- fit0$coefficients[1:ntrt] + fit0$coefficients[(ntrt+1):(2*ntrt)]
  
  for (i in 1:length(output)){
    res <- output[[i]]
    
    rmst0 = colMeans(log(res$beta0))
    rmst1 = colMeans(log(res$beta1))
    #rmst0 = colMeans(res$beta0)
    #rmst1 = colMeans(res$beta1)
    rmst0_bias = round(rmst0 - true_rmst0,3)
    rmst1_bias = round(rmst1 - true_rmst1,3)
    rmst0_cov = res$cov_beta0
    rmst1_cov = res$cov_beta1
    rmst0_sd = colMeans(rmst0_cov)
    rmst1_sd = colMeans(rmst1_cov)
    rmst0_se = apply(res$beta0,2,sd)
    rmst1_se = apply(res$beta1,2,sd)
    
    rmst0_l = log(res$beta0 - 1.96*rmst0_cov)
    rmst1_l = log(res$beta1 - 1.96*rmst1_cov)
    rmst0_h = log(res$beta0 + 1.96*rmst0_cov)
    rmst1_h = log(res$beta1 + 1.96*rmst1_cov)
    #rmst0_l = res$beta0 - 1.96*rmst0_cov
    #rmst1_l = res$beta1 - 1.96*rmst1_cov
    #rmst0_h = res$beta0 + 1.96*rmst0_cov
    #rmst1_h = res$beta1 + 1.96*rmst1_cov
    
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
    
    MSE0 <- round(rowMeans((t(log(res$beta0)) - true_rmst0)^2),3)
    MSE1 <- round(rowMeans((t(log(res$beta1)) - true_rmst1)^2),3)
    #MSE0 <- round(rowMeans((t(res$beta0) - true_rmst0)^2),3)
    #MSE1 <- round(rowMeans((t(res$beta1) - true_rmst1)^2),3)
    
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
                     paste0("rmst0_",1:ntrt," tau"), paste0("rmst1_",1:ntrt," tau"))
  return(as.data.frame(tab))
}
