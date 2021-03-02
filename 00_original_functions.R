

# Original functions ----


# Libraries ----
library(tidyverse)
library(tidymodels)
library(tidyquant)

library(mgcv)     # function1
library(np)       # function1

library(survey)   # function2_2019
library(sampling) # gedata_test


# <><> function1 (the 6 estimators) ----
# fM_A ... A mean
# fM_B ... B mean
# fP_A ... PMIE (from sample A; uses lm())
# fKS  ... NPMIEK
# fKS1 ... (not used)
# fGAM ... NPMIEG
# fPW  ... PWE
# 1.0 - The sample mean from sample A ----
fM_A <- function(dat){
  
  y       <- dat[,1]
  x1      <- dat[,2]
  sIA     <- dat[,4]
  syA     <- y[sIA == 1]
  sx1A    <- x1[sIA == 1]
  etheta1 <- mean(syA)
  etheta2 <- mean(syA[sx1A > 2])
  etheta  <- c(etheta1,etheta2)
  
  return(etheta)
}

# 2.0 The naive estimator from sample B ----
# (sample mean)
fM_B <- function(dat){
  
  y       <- dat[,1]
  x1      <- dat[,2]
  sIB     <- dat[,5]
  syB     <- y[sIB==1]
  sx1B    <- x1[sIB==1]
  etheta1 <- mean(syB)
  etheta2 <- mean(syB[sx1B>2])
  etheta  <- c(etheta1,etheta2)
  
  return(etheta)
}

# 3.0 - The PMIE from sample A using lm ----
# (PMIE = parametric mass imputation estimator)
fP_A <- function(dat){
  
  y    <- dat[,1]
  x1   <- dat[,2]
  x2   <- dat[,3]
  sIA  <- dat[,4]
  sIB  <- dat[,5]
  sx1A <- x1[sIA == 1]
  sx1B <- x1[sIB == 1]
  sx2A <- x2[sIA == 1]
  sx2B <- x2[sIB == 1]
  syB  <- y[sIB == 1]
  
  # Fit model by using sample B
  m     <- lm(syB ~ sx1B + sx2B)
  ebeta <- m$coefficients
  
  # Mass imputed values by using sample A
  iyA     <- ebeta[1] + ebeta[2]*sx1A + ebeta[3]*sx2A
  etheta1 <- mean(iyA)
  etheta2 <- mean(iyA[sx1A > 2])
  etheta  <- c(etheta1,etheta2)
  
  return(etheta)
}

# 4.0 - NPMIEG ---- 
# modeling_method can be "GAM", "RANDOM_FOREST", "XGBOOST", "NNET"
fGAM <- function(dat){
  
  y      <- dat[,1]
  x1     <- dat[,2]
  x2     <- dat[,3]
  # Add x3
  
  # Sample Indicators for A and B
  sIA    <- dat[,4]
  sIB    <- dat[,5]
  
  # 
  sx1A   <- x1[sIA == 1]
  sx1B   <- x1[sIB == 1]
  
  # 
  sx2A   <- x2[sIA == 1]
  sx2B   <- x2[sIB == 1]
  
  # ADD sx3A and B
  
  # 
  syB    <- y[sIB == 1]
  
  # 
  datB   <- cbind(syB,sx1B,sx2B)
  datB2  <- as.data.frame(datB)
  
  datXA  <- cbind(sx1A,sx2A)
  datXA2 <- as.data.frame(datXA)
  
  colnames(datXA2) <- c('sx1B','sx2B')
  
  # Fit model by using sample B ----
  fit <- gam(syB ~ s(sx1B) + s(sx2B),
             data = datB2)
  
  iyA     <- predict(fit,datXA2)
  etheta1 <- mean(iyA)
  etheta2 <- mean(iyA[sx1A > 2])
  etheta  <- c(etheta1,etheta2)
  
  return(etheta)
}



# 5.0 - NPMIEK ----
# What's the diff between fKS and fKS1?? ----
fKS1 <- function(dat){
  
  y       <- dat[,1]
  x1      <- dat[,2]
  x2      <- dat[,3]
  sIA     <- dat[,4]
  sIB     <- dat[,5]
  w       <- dat[,6]
  N       <- length(y)
  # sx1A  <- x1[sIA == 1]
  # sx2A  <- x2[sIA == 1]
  # syB   <- y[sIB == 1]
  # sxB   <- cbind(x1[sIB == 1],x2[sIB == 1])
  # sdatA <- cbind(y[sIA == 1],x1[sIA == 1],x2[sIA == 1],w[sIA == 1])
  sdatB   <- cbind(y[sIB == 1],x1[sIB == 1],x2[sIB == 1])
  sdatB   <- as.data.frame(sdatB)
  
  colnames(sdatB) <- c('syB','sx1B','sx2B')
  
  sdatXA <- cbind(x1[sIA == 1],x2[sIA == 1])
  sdatXA <- as.data.frame(sdatXA)
  
  colnames(sdatXA) <- c('sx1B','sx2B')
  
  # nA <- dim(sdatA)[1]
  # nB <- dim(sdatB)[1]
  
  # Model estimator from B ----
  
  # eh1 <- 1.06*sd(sdatB[,2])*nB^(-1/5)
  # eh2 <- 1.06*sd(sdatB[,3])*nB^(-1/5)
  # eh  <- 1.5*N^(-1/3)
  eh1   <- 1.5*N^(-1/3)
  eh2   <- 1.5*N^(-1/3)
  
  
  bw     <- npregbw(syB ~ sx1B + sx2B,
                    data = sdatB,
                    bws = c(eh1,eh2),
                    bandwidth.compute=F)
  modelB <- npreg(bws = bw, exdat = sdatXA)
  iyA    <- modelB$mean
  
  # a1    <- kronecker(sdatA[,2],sdatB[,2],FUN="-")
  # a2    <- matrix(a1,nA,nB,byrow=T)
  # a3    <- dnorm(a2 / eh)
  # denom <- apply(a3,1,sum)
  # b0    <- rep(0,nA)
  # b1    <- kronecker(b0,sdatB[,1],FUN="+")
  # b2    <- matrix(b1,nA,nB,byrow=T)
  # b3    <- b2*a3
  # numo  <- apply(b3,1,sum)
  # iyA   <- numo / denom
  
  etheta1 <- mean(iyA)
  etheta2 <- mean(iyA[sx1A>2])
  etheta  <- c(etheta1,etheta2)
  
  return(etheta)
}

# * NPMIEK again? ----
fKS <- function(dat){
  
  y     <- dat[,1]
  x1    <- dat[,2]
  x2    <- dat[,3]
  sIA   <- dat[,4]
  sIB   <- dat[,5]
  w     <- dat[,6]
  N     <- length(y)
  sx1A  <- x1[sIA == 1]
  sx2A  <- x2[sIA == 1]
  
  sdatA <- cbind(y[sIA == 1],x1[sIA == 1],x2[sIA == 1],w[sIA == 1])
  sdatB <- cbind(y[sIB == 1],x1[sIB == 1],x2[sIB == 1])
  nA    <- dim(sdatA)[1]
  nB    <- dim(sdatB)[1]
  
  
  # Model estimator from B ----
  
  # Calculate bandwidth ----
  
  # bandwidth1 ----
  # eh1 <- 1.06*sd(sdatB[,2])*nB^{-1/5}
  # eh2 <- 1.06*sd(sdatB[,3])*nB^{-1/5}
  
  # bandwidth2 ----
  # eh1 <- 1.5*nB^-1/3
  # eh2 <- 1.5*nB^-1/3
  
  # bandwidth3
  eh1   <- 1.5*N^(-1/3)
  eh2   <- 1.5*N^(-1/3)
  
  a1_1  <- kronecker(sdatA[,2],sdatB[,2],FUN="-")
  a2_1  <- matrix(a1_1,nA,nB,byrow=T)
  a3_1  <- dnorm(a2_1 / eh1)
  a1_2  <- kronecker(sdatA[,3],sdatB[,3],FUN="-")
  a2_2  <- matrix(a1_2,nA,nB,byrow=T)
  a3_2  <- dnorm(a2_2 / eh2)
  a3    <- a3_1*a3_2
  denom <- apply(a3,1,sum)
  b0    <- rep(0,nA)
  b1    <- kronecker(b0,sdatB[,1],FUN="+")
  b2    <- matrix(b1,nA,nB,byrow=T)
  b3    <- b2*a3
  numo  <- apply(b3,1,sum)
  iyA   <- numo / denom
  
  etheta1 <- mean(iyA)
  etheta2 <- mean(iyA[sx1A > 2])
  etheta  <- c(etheta1,etheta2)
  
  return(etheta)
}

# 6.0 - PWE ----

fPW <- function(dat){
  
  y      <- dat[,1]
  x1     <- dat[,2]
  x2     <- dat[,3]
  sIA    <- dat[,4]
  sIB    <- dat[,5]
  syB    <- y[sIB == 1]
  sx1B   <- x1[sIB == 1]
  w      <- dat[,6]
  nA     <- table(sIA)[2]
  nB     <- table(sIB)[2]
  sdat_A <- cbind(y[sIA == 1],x1[sIA == 1],x2[sIA == 1],w[sIA == 1])
  sdat_B <- cbind(y[sIB == 1],x1[sIB == 1],x2[sIB == 1],rep(1,nB))
  
  # Combine sample A with sample B ----
  z        <- c(rep(1,nB),rep(0,nA))
  sdat_com <- cbind(rbind(sdat_B,sdat_A),z)
  sDat_com <- as.data.frame(sdat_com)
  
  colnames(sDat_com) <- c('y','x1','x2','w','z')
  
  l_M  <- glm(z ~ x1 + x2,
              data=sDat_com,
              family=binomial(logit))
  ep   <- l_M$fitted.values
  epB1 <- ep[z == 1]
  epB0 <- 1 - epB1
  
  log.dA   <- log(sdat_A[,4] - 1)
  fit.dA   <- lm(log.dA ~ sdat_A[,2] + sdat_A[,3])
  sigma.dA <- sqrt(sum(fit.dA$residuals^2) / fit.dA$df.residual)
  ebetaA   <- fit.dA$coefficients
  mu.dB    <- ebetaA[1] + ebetaA[2]*sdat_B[,2] + ebetaA[3]*sdat_B[,3]
  # mu.dA  <- fit.dA$fitted.values
  # dA.hat <- exp(mu.dA + sigma.dA^2 / 2) + 1
  dB.hat   <- exp(mu.dB + sigma.dA^2 / 2) + 1
  
  pw0 <- dB.hat*epB0 / epB1
  pws <- pw0*sum(w[sIA == 1]) / sum(pw0)
  
  etheta1 <- sum(syB*pws) / sum(pws)
  etheta2 <- sum(syB[sx1B > 2]*pws[sx1B > 2]) / sum(pws[sx1B > 2])
  etheta  <- c(etheta1,etheta2)
  
  return(etheta)
}



# <><> function_2_2019 (variance estimation) ----
# fboot ... bootstrap
fboot <- function(dat,theta0,B_B){
  
  y    <- dat[,1]
  x1   <- dat[,2]
  x2   <- dat[,3]
  w    <- dat[,6]
  sIA  <- dat[,4]
  sIB  <- dat[,5]
  sx1A <- x1[sIA==1]
  sx1B <- x1[sIB==1]
  sx2A <- x2[sIA==1]
  sx2B <- x2[sIB==1]
  syA  <- y[sIA==1]
  syB  <- y[sIB==1]
  swA  <- w[sIA==1]
  swB  <- w[sIB==1]
  nA   <- length(sx1A)
  nB   <- length(sx1B)
  N    <- length(y)
  etheta0_P_A <- fP_A(dat)
  etheta0_GAM <- fGAM(dat)
  etheta0_KS  <- fKS(dat)
  etheta0 <- c(etheta0_P_A,etheta0_GAM,etheta0_KS)
  datB    <- cbind(syB,sx1B,sx2B,swB)
  datB2   <- as.data.frame(datB)
  datA    <- cbind(syA,sx1A,sx2A,swA)
  datA2   <- as.data.frame(datA)
  
  # Bootstrap procedure ----
  S_etheta <- NULL
  iter <- 0
  repeat{
    iter <- iter + 1
    # * Bootstrap sample from B ----
    s_id_B  <- sample(1:nB,nB,replace=T)
    s_datB  <- datB2[s_id_B,]
    s_datB2 <- cbind(s_datB[,1:3],rep(0,nA),rep(1,nB),s_datB[,4])
    
    colnames(s_datB2) <- c('y','x1','x2','sIA','sIB','w')
    
    # * Bootstrap sample from A ----
    id_P    <- rep(1:nA,20)
    s_id_A  <- sample(id_P,nA)
    s_datA  <- datA2[s_id_A,]
    s_datA2 <- cbind(s_datA[,1:3],rep(1,nA),rep(0,nB),s_datA[,4])
    
    colnames(s_datA2) <- c('y','x1','x2','sIA','sIB','w')
    
    s_dat <- rbind(s_datB2,s_datA2)
    
    s_etheta_P_A <- fP_A(s_dat)
    s_etheta_GAM <- fGAM(s_dat)
    s_etheta_KS  <- fKS(s_dat)
    s_etheta     <- c(s_etheta_P_A,s_etheta_GAM,s_etheta_KS)
    
    S_etheta <- cbind(S_etheta,s_etheta)
    
    if (iter == B_B) break
  }
  
  eV <- apply((S_etheta-etheta0)^2,1,mean)
  
  qz <- qnorm(0.975)
  
  LB <- etheta0 - qz*sqrt(eV)
  UB <- etheta0 + qz*sqrt(eV)
  AL <- UB - LB
  CI <- as.numeric(LB < theta0 & theta0 < UB)
  
  return(c(eV,AL,CI))
}



# <><> gedata (generation of simulated data) ----
# gedata
## fsA
## fsrs
gedata <- function(B,N,nA,nB,id_m){
  
  # Super population model ----
  x1 <- rnorm(B*N,mean=2,sd=1)
  x2 <- rnorm(B*N,mean=2,sd=1)
  
  # For linear model (index = 1) ----
  if (id_m == 1){
    epsilon <- rnorm(B*N)
    y       <- 0.3 + 2*x1 + 2*x2 + epsilon
  }
  
  # For nonlinear model (index = 2) ----
  if (id_m == 2){
    epsilon   <- rnorm(B*N)
    # epsilon <- rexp(B*N) - 1
    # epsilon <- rchisq(B*N,df=2) - 2
    y         <- 0.3 + 0.5*x1^2 + 0.5*x2^2 + epsilon
  }
  
  # Selecting sample A ----
  M  <- rep(1:N,B)
  M2 <- matrix(M,B,N,byrow=T)
  
  # Function for SRS ----
  # !! Where is tt defined or passed?? ----
  fsA <- function(tt){
    sI0      <- rep(0,length(tt))
    id0      <- sample(tt,nA)
    sI0[id0] <- 1
    return(sI0)
  }
  
  MsIA0 <- apply(M2,1,fsA)
  sIA   <- as.numeric(MsIA0)
  
  # Selecting sample B ----
  n1         <- 0.7*nB
  n2         <- 0.3*nB
  st         <- rep(1,B*N)
  st[x1 > 2] <- 2
  ST         <- matrix(st,B,N,byrow=T)
  
  a1 <- cbind(ST,M2)
  a2 <- as.numeric(t(a1))
  a3 <- array(a2,c(N,2,B))
  
  # Function for stratified SRS ----
  # !! Where is dat defined or passed?? ----
  fsrs <- function(dat){
    
    st0   <- dat[,1]
    id0   <- dat[,2]
    samp1 <- sample(id0[st0==1],n1)
    samp2 <- sample(id0[st0==2],n2)
    sI0   <- rep(0,N)
    sI0[c(samp1,samp2)] = 1
    
    return(sI0)
  }
  
  MsIB0 <- apply(a3,3,fsrs)
  sIB   <- as.numeric(MsIB0)
  w     <- rep((N/nA),(B*N))
  
  My   <- matrix(y,B,N,byrow=T)
  Mx1  <- matrix(x1,B,N,byrow=T)
  Mx2  <- matrix(x2,B,N,byrow=T)
  MsIA <- matrix(sIA,B,N,byrow=T)
  MsIB <- matrix(sIB,B,N,byrow=T)
  Mw   <- matrix(w,B,N,byrow=T)
  
  # <> 6 columns of gedata are here ----
  M   <- cbind(My,Mx1,Mx2,MsIA,MsIB,Mw)
  M   <- t(M)
  aM  <- as.numeric(M)
  Res <- array(aM,c(N,6,B))
  
  if (id_m == 1){
    theta0_1 <- 8.3
    theta0_2 <- mean(y[x1 > 2])
  }
  else{
    theta0_1 <- 5.3
    theta0_2 <- mean(y[x1 > 2])
  }
  
  # !! What are theta0_1 & theta0_2?? ----
  res <- list(Res,theta0_1,theta0_2)
  
  return(res)
}
# dat1_500 <- gedata(1000,10000,500,500,1)
# dat2_500 <- gedata(1000,10000,500,500,2)
# dat1_1000 <- gedata(1000,10000,500,1000,1)
# dat2_1000 <- gedata(1000,10000,500,1000,2)



# <><> respoint (performance metrics) ----
# FRES
# FVAR
## Fboot
FRES <- function(indat,modeling_method = "GAM"){
  
  dat    <- indat[[1]]
  THETA0 <- c(indat[[2]],indat[[3]])
  
  # 1. The sample mean from sample A (Mean A) ----
  res_MA    <- apply(dat,3,fM_A)
  bias_MA   <- res_MA - THETA0
  m_bias_MA <- apply(bias_MA,1,mean)
  rb_MA     <- m_bias_MA / THETA0
  var_MA    <- apply(res_MA,1,var)
  se_MA     <- sqrt(var_MA)
  rse_MA    <- se_MA / THETA0
  mse_MA    <- m_bias_MA^2 + var_MA
  rrmse_MA  <- sqrt(mse_MA) / THETA0
  Res_MA    <- cbind(rb_MA,rse_MA,rrmse_MA)
  
  # 2. The naive estimator (sample mean) from sample B (Mean B) ----
  res_MB    <- apply(dat,3,fM_B)
  bias_MB   <- res_MB - THETA0
  m_bias_MB <- apply(bias_MB,1,mean)
  rb_MB     <- m_bias_MB / THETA0
  var_MB    <- apply(res_MB,1,var)
  se_MB     <- sqrt(var_MB)
  rse_MB    <- se_MB / THETA0
  mse_MB    <- m_bias_MB^2 + var_MB
  rrmse_MB  <- sqrt(mse_MB) / THETA0
  Res_MB    <- cbind(rb_MB,rse_MB,rrmse_MB)
  
  # 3. The PMIE from sample A by using linear model ----
  res_P     <- apply(dat,3,fP_A)
  bias_P    <- res_P - THETA0
  m_bias_P  <- apply(bias_P,1,mean)
  rb_P      <- m_bias_P / THETA0
  var_P     <- apply(res_P,1,var)
  se_P      <- sqrt(var_P)
  rse_P     <- se_P / THETA0
  mse_P     <- m_bias_P^2 + var_P
  rrmse_P   <- sqrt(mse_P) / THETA0
  Res_P     <- cbind(rb_P,rse_P,rrmse_P)
  
  # 4. The PWE from sample B ----
  res_PW    <- apply(dat,3,fPW)
  bias_PW   <- res_PW - THETA0
  m_bias_PW <- apply(bias_PW,1,mean)
  rb_PW     <- m_bias_PW / THETA0
  var_PW    <- apply(res_PW,1,var)
  se_PW     <- sqrt(var_PW)
  rse_PW    <- se_PW / THETA0
  mse_PW    <- m_bias_PW^2 + var_PW
  rrmse_PW  <- sqrt(mse_PW) / THETA0
  Res_PW    <- cbind(rb_PW,rse_PW,rrmse_PW)
  
  # 5. The PMIEK ----
  res_KS    <- apply(dat,3,fKS)
  bias_KS   <- res_KS - THETA0
  m_bias_KS <- apply(bias_KS,1,mean)
  rb_KS     <- m_bias_KS / THETA0
  var_KS    <- apply(res_KS,1,var)
  se_KS     <- sqrt(var_KS)
  rse_KS    <- se_KS / THETA0
  mse_KS    <- m_bias_KS^2 + var_KS
  rrmse_KS  <- sqrt(mse_KS) / THETA0
  Res_KS    <- cbind(rb_KS,rse_KS,rrmse_KS)
  
  # 6. The PMIEG ----
  # Modified fGAM to include modeling_method argument
  res_GAM    <- apply(dat,3,fGAM)
  bias_GAM   <- res_GAM - THETA0
  m_bias_GAM <- apply(bias_GAM,1,mean)
  rb_GAM     <- m_bias_GAM / THETA0
  var_GAM    <- apply(res_GAM,1,var)
  se_GAM     <- sqrt(var_GAM)
  rse_GAM    <- se_GAM / THETA0
  mse_GAM    <- m_bias_GAM^2 + var_GAM
  rrmse_GAM  <- sqrt(mse_GAM) / THETA0
  Res_GAM    <- cbind(rb_GAM,rse_GAM,rrmse_GAM)
  
  RES <- rbind(Res_MA,Res_MB,Res_P,Res_PW,Res_KS,Res_GAM)
  RES <- round(RES,4)
  
  rownames(RES) <- c(c('Mean_A','Domain Mean_A'),
                     c('Mean_B','Domain Mean_B'),
                     c('Mean_P','Domain Mean_P'),
                     c('Mean_PW','Domain Mean_PW'),
                     c('Mean_KS','Domain Mean_KS'),
                     c('Mean_GAM','Domain Mean_GAM'))
  
  colnames(RES) <- c('RB','RSE','RRMSE')
  
  return(RES)
}

# Plot model metrics ----
plot_metrics <- function(res) {
  res %>% 
    filter(!grepl("Domain",estimator)) %>% 
    pivot_longer(cols = RB:RRMSE) %>% 
    mutate(estimator = as_factor(estimator) %>% fct_rev) %>% 
    
    ggplot(aes(estimator,value,fill=estimator)) +
    facet_wrap(~name) +
    geom_col() + theme_tq() + scale_fill_tq() +
    coord_flip()
}
# RES1 <- FRES(list(dat1_500[[1]][,,1:1000],dat1_500[[2]],dat1_500[[3]]))
# 
# RES2 <- FRES(list(dat2_500[[1]][,,1:1000],dat2_500[[2]],dat2_500[[3]]))
# # ADD dat2_500_mod
# RES2_mod <- FRES(list(dat2_500_mod[[1]][,,1:1000],
#                       dat2_500_mod[[2]],
#                       dat2_500_mod[[3]]))
# 
# write.csv(RES1,file='R:/research/Kim/non-probability/mass imputation/submission/JSSAM/revision/simulation/results/RES1_500_bw3.csv')
# write.csv(RES2,file='R:/research/Kim/non-probability/mass imputation/submission/JSSAM/revision/simulation/results/RES2_500_bw3.csv')
# 
# RES1 <- FRES(list(dat1_1000[[1]][,,1:1000],dat1_1000[[2]],dat1_1000[[3]]))
# RES2 <- FRES(list(dat2_1000[[1]][,,1:1000],dat2_1000[[2]],dat2_1000[[3]]))
# 
# write.csv(RES1,file='R:/research/Kim/non-probability/mass imputation/submission/JSSAM/revision/simulation/results/RES1_1000_bw3.csv')
# write.csv(RES2,file='R:/research/Kim/non-probability/mass imputation/submission/JSSAM/revision/simulation/results/RES2_1000_bw3.csv')

# Variance estimation ----
# B_B <- 500
FVAR <- function(indat,B_B){
  
  Theta0 <- c(indat[[2]],
              indat[[3]],
              indat[[2]],
              indat[[3]],
              indat[[2]],
              indat[[3]])
  
  Fboot <- function(dat){
    return(fboot(dat,Theta0,B_B))
  }
  
  RES_var <- apply(indat[[1]],3,Fboot)
  
  return(RES_var)
}
