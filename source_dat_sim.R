# simulate data for non-PH hazards
# different baseline hazard functions for each arm

my.meta.survial.weibull.sim.adj <- function(nt,n1,n2,alphas,betas,px,a,lambda,gamma,tcen,t_trt){
  # nt: number of trials
  # n: number of patients in each trial
  # alphas: coefficient in the weibull arm-based PH model
  # a: heterogeneity of alphas
  # t_trt: included treatments in each trial (list of vectors)
  #set.seed(seed)
  ntrt <- length(alphas)
  alphaj <- matrix(NA, nrow = nt, ncol = ntrt)
  for (k in 1:ntrt){
    alphaj[,k] <- rnorm(nt, mean = alphas[k], sd = a)
  }
  n <- runif(nt,n1,n2)
  dat_comb <- list()
  
  for (j in 1:nt){
    n_trt_j <- length(t_trt[[j]])
    prob <- rep(0,ntrt)
    prob[t_trt[[j]]] <- 1
    trt0 <- sample(x = 1:ntrt, n[j], replace = T, prob = prob)
    trt <- matrix(0,nrow = n[j],ncol = ntrt)
    for (i in 1:n[j]){
      trt[i,trt0[i]] <- 1
    }
    u <- runif(n[j],0.00001,1)
    Xj <- rbinom(n[j],1,px)
    timee <- (-log(u)*exp(-trt %*% alphaj[j,] - trt %*% betas * Xj)/(trt %*% lambda))^(1/(trt %*% gamma))
    timec <- runif(n[j],1,tcen)
    time <- pmin(timee,timec)
    status <- as.numeric(timee <= timec)
    dat <- as.data.frame(cbind(time,status,trt,Xj))
    names(dat) <- c("time","status",paste0("trt",1:ntrt),"X")
    dat$trt <- trt0
    dat$study_id <- 1:n[j]
    dat$trial <- j
    dat_comb[[j]] <- dat
  }
  
  data <- do.call("rbind", dat_comb)
  data$trial <- factor(data$trial)
  data$id <- 1:nrow(data)
  
  return(data)
}







weibull.true.rmst.adj <- function(t,alpha,a,gamma,lambda,beta,x){
  require(zoo)
  mu_beta <- function(y){
    dnorm(y,alpha,a) * integrate(function(z){exp(-lambda*(z^gamma)*exp(beta*x+y))},lower=0,upper=t,rel.tol = 1e-10)$value
  }

  y.seq <- seq(alpha - 5, alpha + 5, 0.001)
  return(sum(0.001*rollapply(sapply(y.seq,mu_beta),2,sum)/2))
}






res.output.RMST.adj <- function(output,alphas,betas,a.list,gamma,lambda,tau){
  tab <- NULL
  ntrt <- length(alphas)
  for (i in 1:length(output)){
    res <- output[[i]]

    true_rmst1 <- true_rmst0 <- NULL
    for (j in 1:ntrt){
      true_rmst0 <- c(true_rmst0, weibull.true.rmst.adj(tau,alpha = alphas[j], a = a.list[i],
                                                          gamma = gamma[j], lambda = lambda[j], beta = betas[j], x = 0))
      true_rmst1 <- c(true_rmst1, weibull.true.rmst.adj(tau,alpha = alphas[j], a = a.list[i],
                                                          gamma = gamma[j], lambda = lambda[j], beta = betas[j], x = 1))
    }
    
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
    
    
    MSE0 <- round(rowMeans((t(res$beta0) - true_rmst0)^2),3)
    MSE1 <- round(rowMeans((t(res$beta1) - true_rmst1)^2),3)
    
    tab <- rbind(tab, c(round(rmst0,3), rmst0_bias, 
                        round(rmst0_sd,3), round(rmst0_se,3),
                        CP0, MSE0,
                        round(rmst0,3), rmst0_bias, 
                        round(rmst0_sd,3), round(rmst0_se,3),
                        CP0, MSE0,
                        round(colMeans(res$sigmab),3)))
  }
  
  colnames(tab) <- c(paste0("rmst0_",1:ntrt), paste0("rmst0_",1:ntrt," bias"),
                     paste0("rmst0_",1:ntrt," sd"), paste0("rmst0_",1:ntrt," se"),
                     paste0("rmst0_",1:ntrt," CP"), paste0("rmst0_",1:ntrt," MSE"),
                     paste0("rmst1_",1:ntrt), paste0("rmst1_",1:ntrt," bias"),
                     paste0("rmst1_",1:ntrt," sd"), paste0("rmst1_",1:ntrt," se"),
                     paste0("rmst1_",1:ntrt," CP"), paste0("rmst1_",1:ntrt," MSE"),
                     paste0("alpha",1:ntrt," tau"),paste0("beta",1:ntrt," tau"))
  return(as.data.frame(tab))
}
