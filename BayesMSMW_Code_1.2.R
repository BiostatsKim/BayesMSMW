library(expm)
library(MASS)
library(truncnorm)
library(dplyr)
library(abind)
# This function is for multi-way analysis (assuming rank R structure)
# Let X be n x p x d, let y be length n, let XTest be N x p x d, YTest be length N, let TrueVal be p x d
BayesMSMW <- function(X.list,Y.list,N_sim=11000,burn_in=1000,Alpha=1,Beta=NULL,multi.source="yes",rank=2,outcome="binary",zsd=1){ #If ms=1 then multi-source analysis is performed, else, single-source
  # Progress bar
  pb = txtProgressBar(min = 0, max = N_sim, initial = 0,style = 3) 
  # Derive dimensions
  ## Vector of source sizes
  Ps <- as.vector(data.matrix(as.data.frame(lapply(X.list, dim)))[2,])
  ## Total number of covariates
  P <- sum(Ps)
  ## Sample size
  N <- as.vector(data.matrix(as.data.frame(lapply(X.list, dim)))[1,])[1]
  ## number of ways (eg time points)
  d <- as.vector(data.matrix(as.data.frame(lapply(X.list, dim)))[3,])[1]
  ## outcome
  y <- unlist(Y.list)
  # Set m= number of sources
  M <- if(is.vector(Ps)==TRUE){length(Ps)}else{1}
  ## concatenate datasets from list into single object X
  Xmat <- X.list[[1]]
  if(length(X.list)>1){
    for(m in 2:length(X.list)){Xmat <- abind(Xmat,X.list[[m]],along = 2)}
    X <- Xmat
  }
  else{X <- Xmat}
  
  # Check that sources are dimensionally linked
  check1 <- var(data.matrix(as.data.frame(lapply(X.list, dim)))[1,])
  check2 <- var(data.matrix(as.data.frame(lapply(X.list, dim)))[3,])
  check3 <- var(data.matrix(as.data.frame(lapply(Y.list, length)))[1,])
  check4 <- if(check1!=0){FALSE}else if(check2!=0){FALSE}else if(check3!=0){FALSE}else{TRUE}

  if(outcome=="binary"){
    N1  <- sum(y)         # Number of successes
    N0  <- N - N1         # Number of failures
  } else{
    N1  <- length(which(y>0))         # Number of successes
    N0  <- N - N1         # Number of failures
  }
  #true_B <-TrueVal      # True values of parameters
  R <- rank

  
  # Initialize objects/parameters
  tau2 <- rep(1,M)
  Var <- list()
  Ww <- rep(1,P)
  cPss <- cumsum(Ps)
  Pss <- c(0,cPss) + 1
  
  
  # Conjugate prior on the coefficients \beta ~ N(beta_0, Q_0)
  w_0 <- rep(0, P)      # w corresponds to covariates (p)
  v_0 <- rep(0, d)      # v corresponds to time points (d)
  
  # Initialize parameters
  tau2_v <- c(1)        # variance for v
  y_var <- c(1)         # variance for y
  w <- rep(0, P)        # initialize w
  v <- rep(0, d)        # initialize v
  z <- rep(0, N)        # initialize z (used in Gibbs sampler)
  
  # Matrices storing samples of the parameter
  w_list <- list()
  tauv_chain <- matrix(0, nrow = N_sim, ncol = 1)
  v_list <- list()
  prod_chain <- matrix(0, nrow = N_sim, ncol = (P*d)) 
  y_var_chain <- matrix(0, nrow = N_sim, ncol = 1)
  
  vecX <- matrix(X,nrow = N, ncol = P*d)
  Ws <- c(rep(0,P*R))
  Ws_0 <- c(rep(0,P*R))
  Vs <- c(rep(0,d*R))
  Vs_0 <- c(rep(0,d*R))
  Xa <- matrix(X,nrow = N, ncol = P*R)
  Xb <- matrix(X,nrow = N, ncol = d*R)
  probtrain_chain <- matrix(0, nrow = N_sim, ncol = N) 
  probs_train <- Y
  Xar <- matrix(nrow= P,ncol = rank)
  Xbr <- matrix(nrow= d,ncol = rank)
  Xr1 <- matrix(nrow = P,ncol = rank)
  Vs_r <- matrix(Vs,ncol = R)
  Ws_r <- matrix(Ws,ncol = R)
  
  
  
  
  
  if(multi.source == "yes"){  # Gibbs sampling algorithm, Multi-Source Multi-Way
    
    tau2_chain <- matrix(0, nrow = N_sim, ncol = M)  
    
    for (t in 2:N_sim) {
      
      # Create aX = n x rP
      for(a in 1:N){
        for(r in 1:R){
          Xar[,r] <- X[a,,] %*% Vs_r[,r]
        }
        Xa[a,] <- as.vector(Xar)
      }
      
      # Update Mean of z
      mu_z <- Xa %*% Ws
      # Draw latent variable z from its full conditional: z | W, y, X        
      if(outcome=="binary"){
        z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = zsd, a = -Inf, b = 0)
        z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = zsd, a = 0, b = Inf)
      } else{
        z <- Y
      }
      
      
      # Compute posterior variance of W
      Q_0 <- diag( rep( rep(tau2,times=Ps), R ) )
      #Q_0 <- diag(rep(c(rep(tau2_1, P1), rep(tau2_2, P2)),R),R*P)
      prec_0 <- solve(Q_0)
      yv <- if(zsd==1){diag(zsd,R*P)}else{diag(1/y_var,R*P)}
      Var <- chol2inv(chol(prec_0 + (yv)%*%crossprod(Xa, Xa)))
      
      # Alt method to help with convergence
      s <- 10^-10
      Var_alt <- chol2inv(chol(s*(prec_0 + (yv)%*%crossprod(Xa, Xa))))
      Mean <- Var %*% (prec_0 %*% Ws_0 + (yv)%*%crossprod(Xa, z)) 
      alt <- function(){
        c(mvrnorm(1, Mean, Var_alt * s))
      }
      Ws <- tryCatch( {c(mvrnorm(1, Mean, Var))}, #<---- Var * 10^-10
                      error = function(e){
                        alt()
                      }
      )
      
      # Posterior variance, prior = IG(Alpha,1/sqrt(P))
      Ws_r <- matrix(Ws,ncol = R)
      if(is.null(Beta)==TRUE){
      for(m in 1:M){tau2[m]<-1/rgamma(1,Alpha + ((Ps[m]*R)/2)  ,1/sqrt(Ps[m]*R) + (1/2)*(sum((Ww[seq(from=Pss[m],to=cPss[m])])^2)))}
      }else{
      for(m in 1:M){tau2[m]<-1/rgamma(1,Alpha + ((Ps[m]*R)/2)  ,Beta + (1/2)*(sum((Ww[seq(from=Pss[m],to=cPss[m])])^2)))}
      }

      
      # Create Xb = n x rd
      for(a in 1:N){
        for(r in 1:R){
          Xbr[,r] <- t(Ws_r[,r]) %*% X[a,,]
        }
        Xb[a,] <- as.vector(Xbr)
      }
      
      mu_z <- Xb %*% Vs
      # Draw latent variable z from its full conditional: z | V, y, X        
      if(outcome=="binary"){
        z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = zsd, a = -Inf, b = 0)
        z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = zsd, a = 0, b = Inf)
      } else{
        z <- Y
      }
      
      # Compute posterior variance of V
      Q_0 <- diag(tau2_v, d*R)
      prec_0 <- solve(Q_0)
      yv <- if(zsd==1){diag(zsd,R*d)}else{diag(1/y_var,R*d)}
      Var <- chol2inv(chol(prec_0 + (yv)%*%crossprod(Xb, Xb))) ### crossprod == X'X
      
      # Alt method to help with convergence issues
      s <- 10^-10
      Var_alt <- chol2inv(chol(s*(prec_0 + (yv)%*%crossprod(Xb, Xb))))
      # Compute posterior mean of V
      Mean <- Var %*% (prec_0 %*% Vs_0 + (yv)%*%crossprod(Xb, z)) 
      alt <- function(){
        c(mvrnorm(1, Mean, Var_alt * s))
      }
      Vs <- tryCatch( {c(mvrnorm(1, Mean, Var))}, #<---- Var * 10^-10
                      error = function(e){
                        alt()
                      }
      )
      
      # v Variance fixed at 1
      tau2_v <- 1
      Vs_r <- matrix(Vs,ncol = R)
      
      
      # Store the draws
      iter <- paste("t:",t,sep='')
      w_list[[iter]] <- Ws
      tauv_chain[t, ] <- tau2_v
      tau2_chain[t, ] <- tau2
      v_list[[iter]] <- Vs
      prod_chain[t, ] <- as.vector(Ws_r %*% t(Vs_r))
      for(n in 1:N){
        for(r in 1:R){
          Xr1[,r] <-  X[n,,] %*% Vs_r[,r]
        }
        probs_train[n] <- if(outcome=="binary"){pnorm(as.vector(Ws_r) %*% as.vector(Xr1))}else{as.vector(Ws_r) %*% as.vector(Xr1)}
      }
      
      probtrain_chain[t, ] <- probs_train
      
      if(zsd==1){
        y_var <- zsd
      }else{
        b_param <- sum((Y - probs_train)^2)
        y_var <- 1/rgamma(1,0.001 + N/2, 0.001 + b_param)
      }
      y_var_chain[t, ] <- y_var
      setTxtProgressBar(pb,t)
    }
    
  } else{

    tau2_chain <- matrix(0, nrow = N_sim, ncol = 1)  
    
    for (t in 2:N_sim) {
      
      # Create aX = n x rP
      for(a in 1:N){
        for(r in 1:R){
          Xar[,r] <- X[a,,] %*% Vs_r[,r]
        }
        Xa[a,] <- as.vector(Xar)
      }
      
      # Update Mean of z
      mu_z <- Xa %*% Ws
      # Draw latent variable z from its full conditional: z | W, y, X        
      if(outcome=="binary"){
        z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = zsd, a = -Inf, b = 0)
        z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = zsd, a = 0, b = Inf)
      } else{
        z <- Y
      }
      
      # Compute posterior variance of W
      Q_0 <- diag(c(rep(tau2_1, R*P)),R*P)
      prec_0 <- solve(Q_0)
      yv <- diag(1/y_var,R*P)
      Var <- chol2inv(chol(prec_0 + (yv)%*%crossprod(Xa, Xa)))
      
      # Alt method to help with convergence
      s <- 10^-10
      Var_alt <- chol2inv(chol(s*(prec_0 + (yv)%*%crossprod(Xa, Xa))))
      Mean <- Var %*% (prec_0 %*% Ws_0 + (yv)%*%crossprod(Xa, z)) 
      alt <- function(){
        c(mvrnorm(1, Mean, Var_alt * s))
      }
      Ws <- tryCatch( {c(mvrnorm(1, Mean, Var))}, #<---- Var * 10^-10
                      error = function(e){
                        alt()
                      }
      )
      
      # Posterior variance, prior = IG(Alpha,1/sqrt(P))
      if(is.null(Beta)==TRUE){
      tau2 <- 1/rgamma(1,Alpha + ((P*R)/2) ,1/sqrt(P*R) + (1/2)*(sum((Ws)^2)))
      }else{
      tau2 <- 1/rgamma(1,Alpha + ((P*R)/2) ,Beta + (1/2)*(sum((Ws)^2)))
      }
      
      Ws_r <- matrix(Ws,ncol = R)
      
      # Create Xb = n x rd
      for(a in 1:N){
        for(r in 1:R){
          Xbr[,r] <- t(Ws_r[,r]) %*% X[a,,]
        }
        Xb[a,] <- as.vector(Xbr)
      }
      
      mu_z <- Xb %*% Vs
      # Draw latent variable z from its full conditional: z | V, y, X        
      if(outcome=="binary"){
        z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = zsd, a = -Inf, b = 0)
        z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = zsd, a = 0, b = Inf)
      } else{
        z <- Y
      }
      
      # Compute posterior variance of V
      Q_0 <- diag(tau2_v, d*R)
      prec_0 <- solve(Q_0)
      yv <- diag(1/y_var, d*R)
      Var <- chol2inv(chol(prec_0 + (yv)%*%crossprod(Xb, Xb))) ### crossprod == X'X
      
      # Alt method to help with convergence issues
      s <- 10^-10
      Var_alt <- chol2inv(chol(s*(prec_0 + (yv)%*%crossprod(Xb, Xb))))
      # Compute posterior mean of V
      Mean <- Var %*% (prec_0 %*% Vs_0 + (yv)%*%crossprod(Xb, z)) 
      alt <- function(){
        c(mvrnorm(1, Mean, Var_alt * s))
      }
      Vs <- tryCatch( {c(mvrnorm(1, Mean, Var))}, #<---- Var * 10^-10
                      error = function(e){
                        alt()
                      }
      )
      
      # v Variance fixed at 1
      tau2_v <- 1
      Vs_r <- matrix(Vs,ncol = R)
      
      # Store the draws
      iter <- paste("t:",t,sep='')
      w_list[[iter]] <- Ws_r
      tauv_chain[t, ] <- tau2_v
      tau1_chain[t, ] <- tau2_1
      v_chain[[iter]] <- Vs_r
      prod_chain[t, ] <- as.vector(Ws_r %*% t(Vs_r))
      for(n in 1:N){
        for(r in 1:R){
          Xr1[,r] <-  X[n,,] %*% Vs_r[,r]
        }
        probs_train[n] <- if(outcome=="binary"){pnorm(as.vector(Ws_r) %*% as.vector(Xr1))}else{as.vector(Ws_r) %*% as.vector(Xr1)}
      }
      
      probtrain_chain[t, ] <- probs_train
      
      if(zsd==1){
        y_var <- zsd
      }else{
        b_param <- sum((Y - probs_train)^2)
        y_var <- 1/rgamma(1,0.001 + N/2, 0.001 + b_param)
      }
      y_var_chain[t, ] <- y_var
      setTxtProgressBar(pb,t)
    }
    
  }
  
  # Get posterior mean of v,w, product of v and w, test class probabilities

  post_tau2_msmw <- tau2_chain[-(1:burn_in), ]
  post_prod_msmw <- prod_chain[-(1:burn_in), ]
  post_probtrain <- probtrain_chain[-(1:burn_in), ]
  post_y_var     <- y_var_chain[-(1:burn_in), ]
  y_var_est <- median(post_y_var)
  
  Post_B <- colMeans(post_prod_msmw)
  
  
  if(outcome=="binary"){
    class_tru <- Y
    misclas <- class_tru - round(colMeans(post_probtrain))
    train_mse <- sum(abs(misclas))/N
  } else{
    class_tru <- Y
    misclas <- class_tru - colMeans(post_probtrain)
    train_mse <- sum((misclas)^2)/sum((class_tru^2))
  }
  
  taus2 <- colMeans(post_tau2_msmw)   
  
  
  output <- list(Post_B,train_mse,taus2,y_var_est,post_y_var,w_list,v_list)
  names(output) <- c("CoefMatrix","TrainErr","Tau2","VarEst","VarChain","Ws","Vs")
  return(output)
}



predict.MSMW <- function(Xtest,outcome="binary",R, Ws, Vs,burn_in=1000){
  burn <- burn_in - 1 
  Nsim <- (length(Ws)-burn)
  pb = txtProgressBar(min = 0, max = Nsim, initial = 0,style = 3) 
  Xx <- Xtest
  Xxr1 <- matrix(ncol = R)
  Nn <- dim(Xtest)[1]
  probtest_chain <- matrix(0, nrow = Nsim, ncol = Nn)   
  if(R==1){
    for(t in 1:Nsim){
      W_s <- Ws[[t]]
      V_s <- Vs[[t]]
      for(n in 1:Nn){
        probs_test[n] <- if(outcome=="binary"){pnorm(t(W_s) %*% Xx[n,,] %*% t(t(V_s)))}else{t(W_s) %*% Xx[n,,] %*% t(t(V_s))}
      }
      probtest_chain[t, ] <- probs_test
      
      setTxtProgressBar(pb,t)
    }
  }else{
    for(t in 1:(length(Ws)-burn)){
      Ws_r <- matrix(Ws[[t]],ncol = R)
      Vs_r <- matrix(Vs[[t]],ncol = R)
      for(n in 1:Nn){
        for(r in 1:R){
          Xxr1[,r] <-  Xx[n,,] %*% Vs_r[,r]
        }
        probs_test[n] <- if(outcome=="binary"){pnorm(as.vector(Ws_r) %*% as.vector(Xxr1))}else{as.vector(Ws_r) %*% as.vector(Xxr1)}
      }
      probtest_chain[t, ] <- probs_test
      
      setTxtProgressBar(pb,t)
    }
  }
  if(outcome=="binary"){
    y_prediction <- round(colMeans(probtest_chain))
  } else{
    y_prediction <- colMeans(probtest_chain)
  }
  output <- list(y_prediction)
  names(output) <- c("prediction")
  return(output)
}



