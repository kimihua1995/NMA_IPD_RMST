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




arm.dat <- function(dat,tau,Omega,df){
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




sim.RMST.bayesian.arm <- function(seed,nt,n1,n2,alphas,betas,px,a,lambda,gamma,tcen,t_trt,tau,Omega,df){
  #beta.list <- beta.m.list <- beta_sd.list <- beta_l.list <- beta_h.list <- NULL
  #sigmab.list <- sigmab.m.list <- NULL
  beta.list <- cov_beta.list <- sigmab.list <- NULL
  beta0.list <- beta1.list <- cov_beta0.list <- cov_beta1.list <- NULL
  ntrt <- length(alphas)
  set.seed(seed)
  for (s in 1:1000){
    data <- my.meta.survial.weibull.sim.adj(nt = nt, n1 = n1, n2 = n2, alphas = alphas, 
                                            betas = betas, px = px, a = a,
                                            lambda = lambda, gamma = gamma, tcen = tcen, t_trt = t_trt)
    
    dat_model0 <- arm.dat(data[data$X==0,],tau,Omega,df)
    dat_model1 <- arm.dat(data[data$X==1,],tau,Omega,df)
    par_model <- c("mu","tau")
    skip_to_next <- FALSE
    tryCatch(fit_jag0 <- jags(data = dat_model0, inits = NULL, 
                             parameters.to.save = par_model,
                             n.iter=5000, n.burnin=1000, n.chains=2, n.thin=1,
                             DIC = TRUE, jags.seed=seed,
                             model.file=rmst_ad_model,
                             quiet = TRUE, progress.bar = "none"),
             error = function(e){skip_to_next <<- TRUE})
    tryCatch(fit_jag1 <- jags(data = dat_model1, inits = NULL, 
                              parameters.to.save = par_model,
                              n.iter=5000, n.burnin=1000, n.chains=2, n.thin=1,
                              DIC = TRUE, jags.seed=seed,
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
    sigmab <- c(sigmab0, sqrt(sigmab0^2 + sigmab1^2))
    
    beta.list <- rbind(beta.list, beta)
    cov_beta.list <- rbind(cov_beta.list, cov_beta)
    sigmab.list <- rbind(sigmab.list, sigmab)
    beta0.list <- rbind(beta0.list, beta0)
    beta1.list <- rbind(beta1.list, beta1)
    cov_beta0.list <- rbind(cov_beta0.list, cov0)
    cov_beta1.list <- rbind(cov_beta1.list, cov1)
    
  }
  
  
  return(list(beta = beta.list,
              cov_beta = cov_beta.list,
              sigmab = sigmab.list,
              beta0 = beta0.list,
              beta1 = beta1.list,
              cov_beta0 = cov_beta0.list,
              cov_beta1 = cov_beta1.list))
}





res.output.RMST.bayesian.arm <- function(output,alphas,a.list,gamma,lambda,tau){
  tab <- NULL
  ntrt <- length(alphas)
  for (i in 1:length(output)){
    res <- output[[i]]
    
    true_beta <- true_sigmab <- NULL
    for (j in 1:length(alphas)){
      true_beta <- c(true_beta, weibull.true.rmst(tau,alpha = alphas[j], a = a.list[i],
                                                  gamma = gamma[j], lambda = lambda[j], link = "identity"))
      true_sigmab <- c(true_sigmab, weibull.true.rmst.sd(tau,alpha = alphas[j], a = a.list[i],
                                                         gamma = gamma[j], lambda = lambda[j], link = "identity"))
    }
    
    true_rmstd <- c()
    for (k in 1:(ntrt-1)){
      true_rmstd[k] <- true_beta[k+1] - true_beta[1]
    }
    
    cover.p <- function(beta_l,beta_h,true_beta){
      cover <- (beta_l <= c(true_beta))*(beta_h >= c(true_beta))
      cp <- sum(cover == 1,na.rm = T)/length(cover)
      round(cp*100,1)
    }
    
    beta_bias <- round(colMeans(res$beta) - true_beta, 3)
    beta_sd = round(colMeans(res$beta_sd),3)
    beta_se = round(apply(res$beta,2,sd),3)
    
    beta_CP <- NULL
    for (j in 1:ntrt){
      beta_CP <- c(beta_CP, cover.p(res$beta_l[,j],res$beta_h[,j],true_beta[j]))
    }
    
    #MSE <- round(beta_se^2 + beta_bias^2,3)
    MSE <- round(rowMeans((t(res$beta) - true_beta)^2),3)
    
    best_arm <- apply(res$beta, 1, which.max)
    true_best <- which.max(true_beta)
    p_rank <- round(mean(best_arm == true_best)*100,2)
    
    
    
    rmstd_CP <- NULL
    rmstd_bias <- round(colMeans(res$d) - true_rmstd, 3)
    for (j in 1:(ntrt-1)){
      rmstd_CP <- c(rmstd_CP, cover.p(res$d_l[,j], res$d_h[,j], true_rmstd[j]))
    }
    rmstd_sd <- round(colMeans(res$d_sd),3)
    rmstd_se <- round(apply(res$d,2,sd),3)
    
    sigmab = colMeans(res$sigmab)
    sigmab_bias = round(sigmab - true_sigmab,3)
    #sigmab_bias = round((sigmab - true_sigmab)/true_sigmab*100,1)
    
    tab <- rbind(tab, c(round(colMeans(res$beta),3), beta_bias, beta_sd, beta_se, beta_CP, MSE,
                        p_rank, round(sigmab,3), sigmab_bias,
                        round(colMeans(res$d),3), rmstd_bias, rmstd_sd, rmstd_se,rmstd_CP))
  }
  
  colnames(tab) <- c(paste0("beta",1:ntrt),paste0("beta",1:ntrt," bias"),
                     paste0("sd",1:ntrt), paste0("se",1:ntrt),
                     paste0("beta",1:ntrt," CP"), paste0("beta",1:ntrt," MSE"),"rank",
                     paste0("tau",1:ntrt),paste0("tau",1:ntrt," bias"),
                     paste0("rmstd",2:ntrt),paste0("rmstd",2:ntrt," bias"),
                     paste0("rmstd",2:ntrt," sd"), paste0("rmstd",2:ntrt, " se"),
                     paste0("rmstd",2:ntrt," CP"))
  return(as.data.frame(tab))
}







