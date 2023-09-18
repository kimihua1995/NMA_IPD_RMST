# calculate rmst and its sd
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
# generate the trail by arm matrix for rmst and sd of the data 
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




dat.NPB <- function(dat,tau,Omega,df){
  rmst <- ad.rmst(dat,tau)
  rmst_m <- c(t(rmst$rmst))
  rmst_sd <- c(t(rmst$se))
  nt <- nrow(rmst$rmst)
  ntrt <- ncol(rmst$se)
  
  trt <- rep(1:ntrt,nt); trt <- trt[!is.na(rmst_m)]
  trial <- rep(1:nt, each=ntrt); trial <- trial[!is.na(rmst_m)]
  rmst_m <- as.vector(na.omit(rmst_m))
  rmst_sd <- as.vector(na.omit(rmst_sd))
  return(list("y" = rmst_m, "s" = trial, "t" = trt,
              "sd" = rmst_sd, "NT" = ntrt, "NS" = nt, 
              "N" = length(rmst_m), "zero.AB" = rep(0,ntrt),
              "Omega" = Omega, "df" = df))
}


rmst_ad_model <- function(){
  for (i in 1:N) {
    y[i] ~ dnorm(mean[s[i],t[i]], prec[i])
    prec[i] <- 1/pow(sd[i],2)
    mean[s[i],t[i]] <- mu[t[i]] + v[s[i],t[i]]
  }

  
  for (j in 1:NS) { v[j, 1:NT] ~ dmnorm(zero.AB[1:NT], invR[1:NT, 1:NT]) }
  
  invR[1:NT, 1:NT] ~ dwish(Omega[ , ],df)
  R[1:NT, 1:NT] <- inverse(invR[1:NT, 1:NT])
  
  
  for (k in 1:NT) { 
    mu[k] ~ dnorm(0, 0.001) 
    tau[k] <- sqrt(R[k,k])
  }
}




sim.NPB <- function(S,data.list,ntrt,nt,tau,Omega,df){
  #beta.list <- beta.m.list <- beta_sd.list <- beta_l.list <- beta_h.list <- NULL
  #sigmab.list <- sigmab.m.list <- NULL
  beta.list <- cov_beta.list <- sigmab.list <- NULL
  beta0.list <- beta1.list <- cov_beta0.list <- cov_beta1.list <- NULL

  for (s in 1:S){
    data <- data.list[[s]]
  
    dat_model0 <- dat.NPB(data[data$X==0,],tau,Omega,df)
    dat_model1 <- dat.NPB(data[data$X==1,],tau,Omega,df)
    par_model <- c("mu","tau")
    skip_to_next <- FALSE
    tryCatch(fit_jag0 <- jags(data = dat_model0, inits = NULL, 
                             parameters.to.save = par_model,
                             n.iter=5000, n.burnin=1000, n.chains=2, n.thin=1,
                             DIC = TRUE, jags.seed=12321,
                             model.file=rmst_ad_model,
                             quiet = TRUE, progress.bar = "none"),
             error = function(e){skip_to_next <<- TRUE})
    tryCatch(fit_jag1 <- jags(data = dat_model1, inits = NULL, 
                              parameters.to.save = par_model,
                              n.iter=5000, n.burnin=1000, n.chains=2, n.thin=1,
                              DIC = TRUE, jags.seed=12321,
                              model.file=rmst_ad_model,
                              quiet = TRUE, progress.bar = "none"),
             error = function(e){skip_to_next <<- TRUE})
    if (skip_to_next) {next} 
    
    res0 <- round(fit_jag0$BUGSoutput$summary[-1,],3)
    res1 <- round(fit_jag1$BUGSoutput$summary[-1,],3)
    
    
    beta0 <- res0[1:ntrt,"mean"]
    beta1 <- res1[1:ntrt,"mean"]
    cov0 <- res0[1:ntrt,"sd"]
    cov1 <- res1[1:ntrt,"sd"]
    sigmab0 <- res0[(ntrt+1):(2*ntrt),"mean"]
    sigmab1 <- res1[(ntrt+1):(2*ntrt),"mean"]
    
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









res.output.NPB <- function(output,ntrt,fit0){
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



