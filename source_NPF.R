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
mvmeta.dat <- function(dat, tau){
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


sim.RMST.mvmeta <- function(seed,nt,n1,n2,alphas,betas,px,a,lambda,gamma,tcen,t_trt,tau){
  beta.list <- cov_beta.list <- sigmab.list <- NULL
  beta0.list <- beta1.list <- cov_beta0.list <- cov_beta1.list <- NULL
  ntrt <- length(alphas)
  set.seed(seed)
  for (s in 1:1000){
    data <- my.meta.survial.weibull.sim.adj(nt = nt, n1 = n1, n2 = n2, alphas = alphas, 
                                        betas = betas, px = px, a = a,
                                        lambda = lambda, gamma = gamma, tcen = tcen, t_trt = t_trt)
    mvmeta_df_0 <- mvmeta.dat(data[data$X==0,], tau)
    mvmeta_df_1 <- mvmeta.dat(data[data$X==1,], tau)
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




