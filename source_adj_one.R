func_surv <- function(y, d, id0){
  #--input--
  #y=time
  #d=status
  
  #--
  id=order(y)
  y=y[id]
  d=d[id]
  
  #--
  t_idx = unique(c(0,y))
  ny = length(y)
  
  #--
  Y = N = C = S = H = D = E = rep(0,length(t_idx))
  
  #i=1
  Y[1] = ny
  N[1] = 0
  C[1] = 0
  S[1] = 1
  H[1] = 0
  D[1] = 0
  E[1] = 0
  
  #i>=2
  for(i in 2:length(t_idx)){
    Y[i] = Y[i-1] - N[i-1] - C[i-1]
    N[i] = ifelse(sum(y==t_idx[i] & d==1)>0, sum(y==t_idx[i] & d==1), 0)
    C[i] = ifelse(sum(y==t_idx[i] & d==0)>0, sum(y==t_idx[i] & d==0), 0)
    
    if(Y[i]<0){Y[i] = 0}
    
    S[i] = ifelse(Y[i]==0, S[i-1], S[i-1]*(1-(N[i]/Y[i])))
    H[i] = ifelse(Y[i]*(Y[i]-N[i])==0, 0, N[i]/(Y[i]*(Y[i]-N[i])))
    
    if(S[i]<0){S[i] = 0}
    
    D[i] = sum(H[2:i])
    E[i] = sqrt((S[i]**2)*D[i])
    
    if(is.na(S[i])){S[i] = 0}
    if(is.na(E[i])){E[i] = 0}
  }
  
  #--output--
  out           = as.data.frame(cbind(t_idx, Y, N, C, S, E))
  colnames(out) = c("t_idx", "n_risk", "n_event", "n_censor", "surv", "se")
  
  #--to match the output of survfit--
  out2 = out[t_idx!=0,]
  
  #--
  Z2 = list()
  Z2$out      = out2
  Z2$t_idx    = out2[,"t_idx"]
  Z2$n_risk   = out2[,"n_risk"]
  Z2$n_event  = out2[,"n_event"]
  Z2$n_censor = out2[,"n_censor"]
  Z2$surv     = out2[,"surv"]
  Z2$se       = out2[,"se"]
  Z2$id       = id0[id]
  return(Z2)
}





cen.weight <- function(data,ntrt,nt,t){
  weight <- id <- NULL
  for (j in 1:nt){
    data.j <- data[data$trial == j,]
    y <- pmin(data.j$time, t)
    data.j$status[y == t] <- 1
    for (k in 1:ntrt){
      yk <- y[data.j$trt == k]
      dk <- data.j$status[data.j$trt == k]
      idk <- data.j$id[data.j$trt == k]
      if (length(yk) > 0){
        fitck <- func_surv(yk, 1-dk, idk)
        weightk <- dk[order(yk)]/rep(pmax(fitck$surv,0.001), table(yk[order(yk)]))
        idk_w <- fitck$id
        weight <- c(weight, weightk)
        id <- c(id, idk_w)
      }
    }
  }
  return(data.frame(weight = weight,
                    id = id))
}



RMST.GLMM.adj <- function(data,ntrt,nt,tau){
  #require(MASS)
  #require(lme4)
  #require(nmle)
  weights <- cen.weight(data,ntrt=ntrt,nt=nt,t=tau)
  data_w <- merge(data, weights, by = "id")
  data_w$Y <- pmin(data_w$time,tau)
  
  
  fit_PQL <- MASS::glmmPQL(fixed = Y ~ -1 + trt1 + trt2 + trt3 +
                           trt1:X + trt2:X + trt3:X, 
                           random = list(trial = ~ -1 + trt1 + trt2 + trt3 +
                                trt1:X + trt2:X + trt3:X),
                           family = quasi(link="identity",variance="constant"), 
                           data = data_w[data_w$weight != 0,], weights = weight, 
                           control = list(opt = "optim"), niter = 100, verbose = F)

  
  
  rand_effect <- VarCorr(fit_PQL)
  beta <- fit_PQL$coefficients$fixed
  vcov_beta <- fit_PQL$varFix
  M <- matrix(c(1,0,0,1,0,0,
                0,1,0,0,1,0,
                0,0,1,0,0,1), nrow = 3, ncol = 6, byrow = T)
  return(list(beta = beta,
              cov_beta = sqrt(diag(vcov_beta)),
              sigmab = as.numeric(rand_effect[1:(2*ntrt),"StdDev"]),
              beta0 = beta[1:ntrt],
              beta1 = beta[1:ntrt] + beta[(ntrt+1):(2*ntrt)],
              cov_beta0 = sqrt(diag(vcov_beta)[1:ntrt]),
              cov_beta1 = sqrt(diag(M %*% vcov_beta %*% t(M)))
              )
         )
}



sim.RMST.GLMM.adj <- function(seed,nt,n1,n2,alphas,betas,px,a,lambda,gamma,tcen,t_trt,tau){
  beta.list <- cov_beta.list <- sigmab.list <- NULL
  beta0.list <- beta1.list <- cov_beta0.list <- cov_beta1.list <- NULL
  ntrt <- length(alphas)
  set.seed(seed)
  for (s in 1:1000){
    data <- my.meta.survial.weibull.sim.adj(nt = nt, n1 = n1, n2 = n2, alphas = alphas, 
                                            betas = betas, px = px, a = a, 
                            lambda = lambda, gamma = gamma, tcen = tcen, t_trt = t_trt)
    
    skip_to_next <- FALSE
    tryCatch(res <- RMST.GLMM.adj(data,ntrt=ntrt,nt=nt,tau=tau),
             error = function(e){skip_to_next <<- TRUE})
    if (skip_to_next) {next} 
    beta.list <- rbind(beta.list, res$beta)
    cov_beta.list <- rbind(cov_beta.list, res$cov_beta)
    sigmab.list <- rbind(sigmab.list, res$sigmab)
    beta0.list <- rbind(beta0.list, res$beta0)
    beta1.list <- rbind(beta1.list, res$beta1)
    cov_beta0.list <- rbind(cov_beta0.list, res$cov_beta0)
    cov_beta1.list <- rbind(cov_beta1.list, res$cov_beta1)
    
  
  }
  
  
  return(list(beta = beta.list,
              cov_beta = cov_beta.list,
              sigmab = sigmab.list,
              beta0 = beta0.list,
              beta1 = beta1.list,
              cov_beta0 = cov_beta0.list,
              cov_beta1 = cov_beta1.list))
}






