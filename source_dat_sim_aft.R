dat.sim.aft <- function(nt,n1,n2,alphas,a,betas,px,sigma,t_trt){
  ntrt <- length(alphas)
  alphaj <- matrix(NA, nrow = nt, ncol = ntrt)
  for (k in 1:ntrt){
    alphaj[,k] <- rnorm(nt, mean = alphas[k], sd = a)
  }
  n <- runif(nt,n1,n2)
  dat_comb <- list()
  for (j in 1:nt){
    prob <- rep(0,ntrt)
    prob[t_trt[[j]]] <- 1
    trt0 <- sample(x = 1:ntrt, n[j], replace = T, prob = prob)
    trt <- matrix(0,nrow = n[j],ncol = ntrt)
    for (i in 1:n[j]){
      trt[i,trt0[i]] <- 1
    }
    
    Xj <- rbinom(n[j],1,px)
    D <- matrix(NA,nrow = n[j],ncol = ntrt)
    #W <- revd(n[j])
    W <- rnorm(n[j])
    for (k in 1:ntrt){
      D[,k] <- exp(alphaj[j,k] + betas[k]*Xj + sigma[k]*W)
    }
    timee <- apply(D*trt, 1, sum)
    timec <- rexp(n[j],0.15)
    time <- pmin(timee, timec)
    status <- as.numeric(timee <= timec)
    
    dat <- as.data.frame(cbind(time, status, timee, timec, trt0, trt, Xj))
    colnames(dat) <- c("time","status","timee","timec",
                       "trt","trt1","trt2","trt3","X")
    dat$trial <- j
    dat_comb[[j]] <- dat
  }
  data <- do.call("rbind", dat_comb)
  data$trial <- factor(data$trial)
  data$id <- 1:nrow(data)
  return(data)
}



true.coef <- function(alphas,betas,px,sigma,tau){
  set.seed(12321)
  n <- 10000000
  ntrt <- length(alphas)
  prob <- rep(1,ntrt)
  trt0 <- sample(x = 1:ntrt, n, replace = T, prob = prob)
  trt <- matrix(0,nrow = n,ncol = ntrt)
  for (i in 1:n){
    trt[i,trt0[i]] <- 1
  }
  X <- rbinom(n,1,px)
  D <- matrix(NA,nrow = n,ncol = ntrt)
  #W <- revd(n)
  W <- rnorm(n)
  for (k in 1:ntrt){
    D[,k] <- exp(alphas[k] + betas[k]*X + sigma[k]*W)
  }
  timee <- apply(D*trt, 1, sum)
  
  data_true <- as.data.frame(cbind(timee, trt0, trt, X))
  colnames(data_true) <- c("D","trt","trt1","trt2","trt3","X")
  data_true$Y <- pmin(data_true$D, tau)
  fit0 <- glm(Y ~ -1 + trt1 + trt2 + trt3 + trt1:X + trt2:X + trt3:X,
              data = data_true, family = quasipoisson())
  return(fit0)
}






res.output.RMST.adj <- function(output,alphas,betas,gamma,lambda,
                                tau,fit0){
  tab <- NULL
  ntrt <- length(alphas)
  
  true_rmst0 <- exp(fit0[1:ntrt])
  true_rmst1 <- exp(fit0[1:ntrt] + fit0[(ntrt+1):(2*ntrt)])
  
  for (i in 1:length(output)){
    res <- output[[i]]
    
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

