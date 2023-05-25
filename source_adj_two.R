###################################################################
mvmeta.dat.adj <- function(data,tau){
  nt <- length(table(data$trial))
  ntrt <- length(table(data$trt))
  weights <- cen.weight(data,ntrt=ntrt,nt=nt,t=tau)
  data_w <- merge(data, weights, by = "id")
  data_w$Y <- pmin(data_w$time,tau)
  nt <- length(table(data_w$trial))
  rmst_m <- S <- NULL
  for (j in 1:nt){
    dataj <- data_w[data_w$trial==j,]
    fit <- lm(Y ~ -1 + trt1 + trt2 + trt3 + trt1:X + trt2:X + trt3:X, 
              weights = weight, data = dataj)
    rmst_m <- rbind(rmst_m, coef(fit))
    vcov_fit <- vcov(fit)
    S <- rbind(S, as.vector(vcov_fit[lower.tri(vcov_fit,diag = T)]))
  }
  return(list(rmst_m = rmst_m, S = S))
}



sim.RMST.mvmeta.adj <- function(seed,nt,n1,n2,alphas,betas,px,a,lambda,gamma,tcen,t_trt,tau){
  beta.list <- cov_beta.list <- sigmab.list <- NULL
  beta0.list <- beta1.list <- cov_beta0.list <- cov_beta1.list <- NULL
  ntrt <- length(alphas)
  set.seed(seed)
  for (s in 1:1000){
    data <- my.meta.survial.weibull.sim.adj(nt = nt, n1 = n1, n2 = n2, alphas = alphas, 
                                            betas = betas, px = px, a = a, 
                                            lambda = lambda, gamma = gamma, tcen = tcen, t_trt = t_trt)
    mvmeta_df <- mvmeta.dat.adj(data, tau)
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
    
  }
  
  
  return(list(beta = beta.list,
              cov_beta = cov_beta.list,
              sigmab = sigmab.list,
              beta0 = beta0.list,
              beta1 = beta1.list,
              cov_beta0 = cov_beta0.list,
              cov_beta1 = cov_beta1.list))
}



