###--------------------------------------------------
# Latent Class Linear Mixed Models
# Continuous response
# K classes
###--------------------------------------------------

### Change this location to read the JAGS code 
setwd("HERE")

savepath <- "OUT/PATH"

###--------------------------------------------------
### Packages 
library(ggplot2)
library(rjags)
library(lcmm)
library(label.switching)  #(version 1.3)       
library(RColorBrewer)
library(mnormt)
library(dclone) # To run MCMC in
library(snow)   # several cores
library(R2jags)
library(abind)
library(MCMCvis)
library(parallelly)
library(loo) ### package to compute WAIC and LOO
library(foreach)
library(doParallel)

detectCores()

###--------------------------------------------------
### Truncated Normal I[Z < trb]
rnormright <- function(trb,mu,sig){
  rp <- pnorm(trb, mean=mu, sd=sig)
  u <- rp*runif(1)
  q <- qnorm(u, mean=mu, sd=sig)
  if(!is.finite(q)){ q = trb }
  return(q)
}
###--------------------------------------------------
###--------------------------------------------------

### JAGS Specifications 

n.adapt = 30000
n.update = 30000
n.iter = 30000
n.thin = 15


### JOB Parameters
if( !exists("TT") ) TT <- 8  # number of time for repetitions 
if( !exists("n") ) n <- 200   # number of subjects 
if( !exists("K2") ) K2 <- 2   # number of latent classes being estimated
if( !exists("TID") ) TID = 0 # chunk ID

K1 <- 3   # number of latent classes simulating  (always 3 in this setting, otherwise, additional specs are needed below)

cat("n=",n,"; K2=", K2, "; TT=", TT, "; TID", TID, sep = "")


#### ---- INDEPENDENT VARS SIMULATION 
set.seed(1234)
{
### Subjects and time
id <- rep(1:n, each=TT)   # identify subjects
time <- rep(1:TT, times=n)  # time
time2 <- time^2 

### Indicator for the beginning of observations for each subject for long tables format 
offset <- c((0:n)*TT+1)  

### Fixed effects
beta0 <- 16
beta <- c(0.2,0.5)#,0.5)   # coefficients for fixed effects X1, X2, etc. 
# Age: gama(59.7,8.8)
X1 <- rep(rgamma(n,shape=46.024,rate=0.771),each=TT)   # covariates fixed effects 
# Sex: female 53%
X2 <- rep(rbinom(n,prob=0.53,size=1),each=TT)   # covariates for fixed effects  
#X3 <- rnorm(n*TT,mean=1,sd=1)   # covariates fixed effects 
X = cbind(X1,X2)#,X3)   # fixed effects design matrix 

### Random effects
Gama0 <- cbind(c(0.5,0),c(0,0.5))   # variance-covariance matrix for random effects 
gama <- rmnorm(n,0,Gama0)  # random effects 
Z = cbind(rep(1,n*TT),time)   # random effects design matrix 

### Latent classes 
### { clases K1=3 simula, K=K2=2 estima
alpha <- matrix(0,K1,2)   # KxQ vector of coefficients for the latent class membership 
alpha[,1] <- c(0,-0.2,0.1)   # coefficients for covariate V1  
alpha[,2] <- c(0,0.2,-0.4)   # coefficients for covariate V2  
#alpha[,3] <- c(0.5,0)   # coefficients for covariate V3  
# BMI: gama(29.7,4.7)
### }

### { Codigo WAIC LOO
V1 = rgamma(n,shape=39.931,rate=1.344)   # covariate for the latent class membership 
# WOMAC: gama(23.5,16.9)
V2 = rgamma(n,shape=1.933,rate=0.082)    # covariate for the latent class membership 
#V3 = rnorm(n,mean=2,sd=1)   # covariate for the latent class membership 
V = cbind(V1,V2)   # latent class membership design matrix 
### } 

### Class-specific fixed effects
# lambda <- matrix(0,K,3)   # KxP vector of coefficients of class-specific fixed effects  
# U = cbind(rep(1,n*TT),time,time2)   ### class-specific fixed effects design matrix 
### { clases K1=3 simula, K=K2=2 estima
lambda <- matrix(0,K1,3)   # K1xP vector of coefficients of class-specific fixed effects  
lambda[,1] <- c(    0, 4.8,  5.2)  # covariate for the class-specific fixed effects U1
lambda[,2] <- c(-0.05,-0.5, -1.2)   # covariate for the class-specific fixed effects U2  
lambda[,3] <- c(-0.01,-0.01,-0.01)   # covariate for the class-specific fixed effects U3  
U = cbind(rep(1,n*TT),time,time2)   ### class-specific fixed effects design matrix
### }

### Variance 
tau2 <- 0.70   # variance of the latent classes 

### Measurement errors
### { clases K1=3 simula, K=K2=2 estima
sigma2 <- c(1.20, 0.1, 0.3)   # variance of the measurement errors
### }

#### ----
}

###--------------------------------------------------
### Functions for the Results
{
  etiquetar <- function(claseVer,claseEst){
    claseVer2 = claseVer
    claseEst2 = claseEst
    claseEstNew <- rep(NA,length(claseEst))
    tab <- table(claseVer2,claseEst2)
    labVer <- as.numeric(rownames(tab))
    labEst <- as.numeric(colnames(tab))
    if(length(labEst)==1){
      claseEstNew = claseEst
    }
    if(length(labEst)>1){
      for(j in 1:length(labVer)){
        tab <- table(claseVer2,claseEst2)
        if(length(tab)>0){
          labVer <- as.numeric(rownames(tab))
          labEst <- as.numeric(colnames(tab))
          id <- which(max(tab)==tab)[1]
          idV = id%%length(labVer)
          idE = floor((id+length(labEst)-1)/length(labEst))
          if(idV==0){
            idV=length(labVer)
          }
          claseEstNew[claseEst==labEst[idE]] = labVer[idV]
          claseEst2[claseEst2==labEst[idE]] = NA
          claseVer2[claseVer2==labVer[idV]] = NA
        }
      }
    }
    return(claseEstNew)
  }
  
  ###--------------------------------------------------
  
  ###--------------------------------------------------
  ### Functions for LABEL SWITCHING proposal
  
  # Compute necessary information that will be additional input to the label.switching package.  
  # Define the complete log-likelihood function. 
  
  complete.loglikelihood <- function(iter, zclass, data2,
                                     mcmc.pars, mcmc.beta0, mcmc.beta, mcmc.alpha, 
                                     mcmc.tau2, mcmc.sigma2) { 
    
    n <- data2$n
    K <- data2$K
    Yerror <- data2$Yerror
    XX <- data2$X
    ZZ <- data2$Z
    UU <- data2$U
    VV <- data2$V
    offset <- data2$offset
    
    pars = mcmc.pars[iter,,]
    beta0 = mcmc.beta0[iter]
    beta = mcmc.beta[iter,]
    alpha = mcmc.alpha[iter,,]
    tau2 = mcmc.tau2[iter]
    sigma2 = mcmc.sigma2[iter]
    
    lambda <- pars[,1:ncol(UU)]
    m <- length(Yerror)
    ppi <- matrix(NA,n,K)
    w <- matrix(NA,n,K)
    eta <- rep(NA,n) 
    eta2 <- rep(NA,n) 
    logw <- rep(NA,n)
    logl <- rep(NA,n)
    for(j in 1:n){
      for(k in 1:K){ 
        ppi[j,k] <- exp(VV[j,]%*%alpha[k,])
      } 
      w[j,] <- ppi[j,]/sum(ppi[j,])
      logw[j] <- log(w[j,zclass[j]])  
      
      idx <- (offset[j]):(offset[j+1]-1)
      eta2[idx] <- beta0 + XX[idx,]%*%beta + UU[idx,]%*%lambda[zclass[j],] 
      for(i in offset[j]){
        media <- eta2[i]
        eta[i] <- media
        varianza <- tau2 + sigma2[zclass[j]]
        logl[i] <- dnorm(Yerror[i],mean=media,sd=sqrt(varianza), log=T)
      }
      for(i in (offset[j]+1):(offset[j+1]-1)){
        beta2 <- (eta[i-1]-eta2[i])/sqrt(tau2)
        z2 <- pnorm(beta2)
        medtrunc <- sqrt(tau2)*(-dnorm(beta2))/z2 ### media de una normal truncada 
        medaux <- ifelse(is.finite(medtrunc),medtrunc,0)  
        media <- eta2[i]+medaux 
        eta[i] <- media
        vartrunc <- (1+ (-beta2*dnorm(beta2)/z2) -((-dnorm(beta2)/z2)^2)) ### varianza de una normal truncada
        varaux <- tau2*ifelse(is.finite(vartrunc) & vartrunc>0,vartrunc,1)
        varianza <- varaux + sigma2[zclass[j]]
        logl[i] <- dnorm(Yerror[i], mean=media, sd=sqrt(varianza), log=T)
      }
    } 
    return(sum(logw)+sum(logl,na.rm=TRUE))
  }
  
  
  ### computing allocation probabilities for stephens method and ECR-ITERATIVE-2
  fn.probabilities <- function(iter, data2, TT, mcmc.pars, 
                               mcmc.tau2, mcmc.beta0, mcmc.beta, mcmc.alpha, mcmc.sigma2){   
    n <- data2$n
    K <- data2$K
    m <- length(data2$Yerror)  
    Yerror <- data2$Yerror
    XX <- data2$X
    ZZ <- data2$Z
    UU <- data2$U
    VV <- data2$V
    offset <- data2$offset
    
    proba <- array(data = NA, dim = c(1, n, K))
    tau2 <- mcmc.tau2[iter]
    lambda <- mcmc.pars[iter,,]
    beta0 <- mcmc.beta0[iter]
    beta <- mcmc.beta[iter,]
    alpha <- mcmc.alpha[iter,,]
    sigma2 <- mcmc.sigma2[iter,]
    ppi <- matrix(NA,n,K)
    w <- matrix(NA,n,K)
    like <- array(NA,dim=c(m,K))
    eta <- array(NA,dim=c(m,K))
    eta2 <- array(NA,dim=c(m,K))
    num = array(NA,dim=c(n,K))
    for(j in 1:n){
      idx <- (offset[j]):(offset[j+1]-1) 
      for(k in 1:K){
        eta2[idx,k] <- beta0 + XX[idx,]%*%beta + UU[idx,]%*%lambda[k,] 
        ppi[j,k] <- exp(VV[j,]%*%alpha[k,])
        for(i in offset[j]){
          media <- eta2[i,k]
          eta[i,k] <- media
          varianza <- tau2 + sigma2[k]
          like[i,k] <- dnorm(Yerror[i],mean=media,sd=sqrt(varianza), log=FALSE)
        } ### END i in offset[j]
        for(i in (offset[j]+1):(offset[j+1]-1)){ 
          beta2 <- (eta[i-1,k]-eta2[i,k])/sqrt(tau2)
          z2 <- pnorm(beta2)
          medtrunc <- sqrt(tau2)*(-dnorm(beta2))/z2 ### media de una normal truncada 
          medaux <- ifelse(is.finite(medtrunc),medtrunc,0)  
          media <- eta2[i,k]+medaux 
          eta[i,k] <- media
          vartrunc <- (1+ (-beta2*dnorm(beta2)/z2) -((-dnorm(beta2)/z2)^2)) ### varianza de una normal truncada
          varaux <- tau2*ifelse(is.finite(vartrunc) & vartrunc>0,vartrunc,1)
          varianza <- varaux + sigma2[k]
          like[i,k] <- dnorm(Yerror[i], mean=media, sd=sqrt(varianza), log=FALSE)
        } ### END i in (offset[j]+1):(offset[j+1]-1)
        num[j,k] <- prod(like[(offset[j]):(offset[j+1]-1),k])
      }  ### END k in 1:K
      w[j,] <- ppi[j,]/sum(ppi[j,])
      proba[1, j, ] <- (num[j,]*w[j,])/sum(num[j,]*w[j,])     
    } ### END j in 1:n
    #  } ### END iter in 1:m2
    proba[is.na(proba)]=0
    return(proba)
  }
  
  ###--------------------------------------------------
  
  ###--------------------------------------------------
  ### Functions for LABEL SWITCHING lcmm
  
  # Compute necessary information that will be additional input to the label.switching package.  
  # Define the complete log-likelihood function. 
  
  complete.loglikelihood.lcmm <- function(iter, zclass, data2lcmm,
                                          mcmc.pars, mcmc.beta0, mcmc.beta, mcmc.alpha, 
                                          mcmc.tau2) { 
    n <- data2lcmm$n
    K <- data2lcmm$K
    Wtrue <- data2lcmm$Wtrue
    XX <- data2lcmm$X
    ZZ <- data2lcmm$Z
    UU <- data2lcmm$U
    VV <- data2lcmm$V
    offset <- data2lcmm$offset
    
    pars = mcmc.pars[iter,,]
    beta0 = mcmc.beta0[iter]
    beta = mcmc.beta[iter,]
    alpha = mcmc.alpha[iter,,]
    tau2 = mcmc.tau2[iter]
    
    lambda <- pars[,1:ncol(UU)]
    m <- length(Wtrue)
    ppi <- matrix(NA,n,K)
    w <- matrix(NA,n,K)
    logw <- rep(NA,n)
    logl <- rep(NA,n)
    for(j in 1:n){
      for(k in 1:K){ 
        ppi[j,k] <- exp(VV[j,]%*%alpha[k,])
      } 
      w[j,] <- ppi[j,]/sum(ppi[j,])
      logw[j] <- log(w[j,zclass[j]])  
      idx <- (offset[j]):(offset[j+1]-1)
      media <- beta0 + XX[idx,]%*%beta + UU[idx,]%*%lambda[zclass[j],] 
      varianza <- diag(tau2,length(idx)) 
      varianza <- (varianza+t(varianza))/2
      logl[j] <- dmnorm(Wtrue[idx], mean=t(media),varcov=varianza, log=T)
    } 
    return(sum(logw)+sum(logl,na.rm=TRUE))
  }
  
  # computing allocation probabilities for stephens method and ECR-ITERATIVE-2
  fn.probabilities.lcmm <- function(iter, data2lcmm, TT, 
                                    mcmc.pars, mcmc.tau2, mcmc.beta0, mcmc.beta, mcmc.alpha){   
    ###   like <- rep(0, data2lcmm$offset[n+1]-1)
    n <- data2lcmm$n
    K <- data2lcmm$K
    m <- length(data2lcmm$Wtrue)
    #  like <- rep(0, n)
    Wtrue <- data2lcmm$Wtrue
    XX <- data2lcmm$X
    ZZ <- data2lcmm$Z
    UU <- data2lcmm$U
    VV <- data2lcmm$V
    offset <- data2lcmm$offset
    proba <- array(NA, dim = c(1, n, K))
    tau2 <- mcmc.tau2[iter]
    lambda <- mcmc.pars[iter,,]
    beta0 <- mcmc.beta0[iter]
    beta <- mcmc.beta[iter,]
    alpha <- mcmc.alpha[iter,,]
    ppi <- matrix(NA,n,K)
    w <- matrix(NA,n,K)
    like <- array(NA,dim=c(n,K))
    for(j in 1:n){
      idx <- (offset[j]):(offset[j+1]-1)
      for(k in 1:K){
        ppi[j,k] <- exp(VV[j,]%*%alpha[k,])
        media <- beta0 + XX[idx,]%*%beta + UU[idx,]%*%lambda[k,] 
        varianza <- diag(tau2,length(idx)) 
        varianza <- (varianza+t(varianza))/2
        like[j,k] <- dmnorm(Wtrue[idx], mean=t(media),varcov=varianza)
      }
      w[j,] <- ppi[j,]/sum(ppi[j,])
      proba[1, j, ] <- (like[j,]*w[j,])/sum(like[j,]*w[j,])     
    } 
    #  }
    proba[is.na(proba)]=0
    return(proba)
  }
  ###--------------------------------------------------
  
  ###--------------------------------------------------
  ### computing complete loglikelihoods in order to find pivot
  fn.zmapindex <- function(zclass, logLvec, m2){
    iter <- 1
    mapindex <- 1
    maxL <- logLvec[iter]
    for(iter in 2:m2) {
      logL <- logLvec[iter]
      if(logL > maxL) {
        maxL <- logL
        mapindex <- iter
      } 
    }
    print(paste("complete likelihood pivot = ", mapindex))
    mapindex <- mapindex
    zmap <- zclass[mapindex, ]
    return(list("mapindex"=mapindex,"zmap"=zmap,"maxL"=maxL))
  }
  ###--------------------------------------------------
  
}
###--------------------------------------------------

###--------------------------------------------------
### Initialization values for JAGS
###--------------------------------------------------
{
  inits1 <- function(){	list(   
    "beta0" = rnorm(1,0,0.1) ,
    "beta" = rnorm(data2$L,0,0.1) ,
    "sigma2_prior" = rep(1,data2$K) , ### { clases K1=3 simula, K=K2=2 estima } 
    "Wtrue"=rnorm(length(data2$Z[,2]),mean=0,sd=0.1)-data2$Z[,2] 
    ### (creciente) usar el tiempo para proponer una prior
  )	}
  
  param1 = c("beta0","beta", "lambda", "alpha",
             "tau2", "sigma2"
  )
  
  param2 = c("beta0","beta", "lambda", "alpha", 
             "tau2", "sigma2", 
             "g" 
  )
  
  inits3 <- function(){	list(   
    "beta0" = rnorm(1,0,0.1) ,
    "beta" = rnorm(data3$L,0,0.1) ,
    "sigma2_prior" = rep(1,data3$K) , ### { clases K1=3 simula, K=K2=2 estima }
    "gama" = matrix(rnorm(data3$n*data3$M,0,0.01),data3$n,data3$M) ,
    "Wtrue"= rnorm(length(data3$Z[,2]),mean=-data3$Z[,2],sd=0.1) ### (decreciente) usar el tiempo para proponer una prior 
  )	}
  
  param3 = c("beta0","beta","lambda",
             "tau2", "sigma2",
             "Gama0", "invGama0" 
  )
  
  ###--------------------------------------------------
  
  inits1lcmm <- function(){	list(   
    "beta0" = rnorm(1,0,0.1) ,
    "beta" = rnorm(data2lcmm$L,0,0.1) 
  )	}
  
  param1lcmm = c("beta0","beta", "lambda", "alpha",
                 "tau2"
  )
  
  param2lcmm = c("beta0","beta", "lambda", "alpha", 
                 "tau2", 
                 "g"
  )
  
  inits3lcmm <- function(){	list(   
    "beta0" = rnorm(1,0,0.1) ,
    "beta" = rnorm(data3lcmm$L,0,0.1) ,
    "gama" = matrix(rnorm(data3lcmm$n*data3lcmm$M,0,0.01),data3lcmm$n,data3lcmm$M) 
  )	}
  
  param3lcmm = c("beta0","beta","lambda", 
                 "tau2", 
                 "Gama0", "invGama0"
  ) 
}

###--------------------------------------------------
### LABEL SWITCHING Specifications

set <- c("STEPHENS",   ### p is required for STEPHENS
         "PRA", 
         "ECR", 
         "ECR-ITERATIVE-1", 
         "ECR-ITERATIVE-2",   ### p is required for ECR-ITERATIVE-2
         "AIC"#,
         #         "SJW",  ### data is required for SJW
         #         "DATA-BASED"  ### data is required for DATA-BASED
)

###--------------------------------------------------

###--------------------------------------------------
### Save values for simulations
{
  SIM = 20
  
  ### { clases K1=3 simula, K=K2=2 estima
  CLASES = array(NA,dim=c(SIM,K1*K2))
  CLASESlcmm = array(NA,dim=c(SIM,K1*K2))
  ### }
  
  if( K2 == 2 ){
    ### { clases K1=3 simula, K=K2=2 estima   
    PARAM = array(NA,dim=c(SIM,20,5))
    PARAMlcmm = array(NA,dim=c(SIM,18,5))
    ### OJO estas dimensiones de 21 y 19 para estimar K=2 clases
    ### }
  }else if( K2 == 3 ){
    ### { clases K1=3 simula, K=K2=3 estima   
    PARAM = array(NA,dim=c(SIM,24,5))
    PARAMlcmm = array(NA,dim=c(SIM,21,5))
    ### OJO estas dimensiones de 25 y 22 para K=3 deben cambiarse
    ### en caso de cambiar el numero de parametros y clases K
    ### }
  }else if( K2 == 4 ){
    ### { clases K1=3 simula, K=K2=4 estima   
    PARAM = array(NA,dim=c(SIM,28,5))
    PARAMlcmm = array(NA,dim=c(SIM,24,5))
    ### OJO estas dimensiones de 29 y 25 para estimar K=4 clases
    ### }
  }
  
  DIC_all = array(NA,dim=c(SIM,1))
  DIClcmm_all = array(NA,dim=c(SIM,1))
  DIC1_all = array(NA,dim=c(SIM,1))
  DIC1lcmm_all = array(NA,dim=c(SIM,1))
  
  ### { Codigo WAIC LOO
  WAIC_all = array(NA,dim=c(SIM,2))
  WAIClcmm_all = array(NA,dim=c(SIM,2))
  LOO_all = array(NA,dim=c(SIM,2))
  LOOlcmm_all = array(NA,dim=c(SIM,2))
  ELPD_all = array(NA,dim=c(SIM,2))
  ELPDlcmm_all = array(NA,dim=c(SIM,2))
  ### OJO: guarda criterio de informacion = -2 deviance 
  ### }
  
  if( K2 == 2 ){
    ### { clases K1=3 simula, K=K2=2 estima   
    colnames(PARAM) = c("Gama0[1,1]","Gama0[2,1]","Gama0[1,2]","Gama0[2,2]",
                        "beta[1]","beta[2]","beta0",
                        "invGama0[1,1]","invGama0[2,1]","invGama0[1,2]","invGama0[2,2]", 
                        "lambda[1,1]","lambda[2,1]",
                        "lambda[1,2]","lambda[2,2]",
                        "lambda[1,3]","lambda[2,3]",
                        "sigma2[1]","sigma2[2]","tau2") 
    ### } 
    ### { clases K1=3 simula, K=K2=2 estima   
    colnames(PARAMlcmm) = c("Gama0[1,1]","Gama0[2,1]","Gama0[1,2]","Gama0[2,2]",
                            "beta[1]","beta[2]","beta0",
                            "invGama0[1,1]","invGama0[2,1]","invGama0[1,2]","invGama0[2,2]", 
                            "lambda[1,1]","lambda[2,1]",
                            "lambda[1,2]","lambda[2,2]",
                            "lambda[1,3]","lambda[2,3]",
                            "tau2") 
    ### } 
  }else if ( K2==3 ){
    ### OJO: para estimar K=3 clases simulando de K=3 clases
    # > unlist(jagsfit3$BUGSoutput[["mean"]])
    # Gama01       Gama02       Gama03       Gama04        beta0 beta1 
    # beta2     deviance    invGama01    invGama02    invGama03 
    # invGama04      lambda1      lambda2      lambda3      lambda4 
    # lambda5      lambda6      lambda7      lambda8      lambda9 
    # sigma21      sigma22      sigma23         tau2 
    
    ### { clases K1=3 simula, K=K2=3 estima   
    colnames(PARAM) = c("Gama0[1,1]","Gama0[2,1]","Gama0[1,2]","Gama0[2,2]",
                        "beta[1]","beta[2]","beta0",
                        "invGama0[1,1]","invGama0[2,1]","invGama0[1,2]","invGama0[2,2]",
                        "lambda[1,1]","lambda[2,1]","lambda[3,1]",
                        "lambda[1,2]","lambda[2,2]","lambda[3,2]",
                        "lambda[1,3]","lambda[2,3]","lambda[3,3]",
                        "sigma2[1]","sigma2[2]","sigma2[3]","tau2")
    ## }
    ## { clases K1=3 simula, K=K2=3 estima
    colnames(PARAMlcmm) = c("Gama0[1,1]","Gama0[2,1]","Gama0[1,2]","Gama0[2,2]",
                            "beta[1]","beta[2]","beta0",
                            "invGama0[1,1]","invGama0[2,1]","invGama0[1,2]","invGama0[2,2]",
                            "lambda[1,1]","lambda[2,1]","lambda[3,1]",
                            "lambda[1,2]","lambda[2,2]","lambda[3,2]",
                            "lambda[1,3]","lambda[2,3]","lambda[3,3]",
                            "tau2")
    ## }
  }else if ( K2==4 ){
    ### { clases K1=3 simula, K=K2=4 estima   
    colnames(PARAM) = c("Gama0[1,1]","Gama0[2,1]","Gama0[1,2]","Gama0[2,2]",
                        "beta[1]","beta[2]","beta0",
                        "invGama0[1,1]","invGama0[2,1]","invGama0[1,2]","invGama0[2,2]",
                        "lambda[1,1]","lambda[2,1]","lambda[3,1]","lambda[4,1]",
                        "lambda[1,2]","lambda[2,2]","lambda[3,2]","lambda[4,2]",
                        "lambda[1,3]","lambda[2,3]","lambda[3,3]","lambda[4,3]",
                        "sigma2[1]","sigma2[2]","sigma2[3]","sigma2[4]","tau2")
    ## }
    ## { clases K1=3 simula, K=K2=4 estima
    colnames(PARAMlcmm) = c("Gama0[1,1]","Gama0[2,1]","Gama0[1,2]","Gama0[2,2]",
                            "beta[1]","beta[2]","beta0",
                            "invGama0[1,1]","invGama0[2,1]","invGama0[1,2]","invGama0[2,2]",
                            "lambda[1,1]","lambda[2,1]","lambda[3,1]","lambda[4,1]",
                            "lambda[1,2]","lambda[2,2]","lambda[3,2]","lambda[4,2]",
                            "lambda[1,3]","lambda[2,3]","lambda[3,3]","lambda[4,3]",
                            "tau2")
    ## }
  }
  
  ### time points 1,2,3,4,5
  # Yerror[-(offset[-(n+1)]+TT-1)]
  ### times points 2,3,4,5,6
  # Yerror[-(offset[-(n+1)])]
  CLASEtrue = array(NA,dim=c(SIM,K1))
  
  # sum(Wtrue[-(offset[-(n+1)]+TT-1)]>=Wtrue[-(offset[-(n+1)])]) ### todas No-crecientes
  # sum(Wtrue[-(offset[-(n+1)]+TT-1)]<Wtrue[-(offset[-(n+1)])]) ### ninguna creciente 
  # sum(Yerror[-(offset[-(n+1)]+TT-1)]>=Yerror[-(offset[-(n+1)])])
  # sum(Yerror[-(offset[-(n+1)]+TT-1)]<Yerror[-(offset[-(n+1)])])
  
  
  ### { clases K1=3 simula, K=K2=2 estima }
  MONOTON = array(NA,dim=c(SIM,K1*2)) 
  colnames(MONOTON) = c("dec1","dec2","dec3","inc1","inc2","inc3")
}

### Begin FOR simulation --------------------------------------------------
#### Makeforeach
num_cores.iters <- 8 # for 32 cores in cl1 OR change to 10 for 40 cores in cl1
cl0 <- parallel::makeCluster(num_cores.iters)
registerDoParallel(cl0)

n.chains = 4 

### Performs initial model fitting and label switching
ptm <- Sys.time()

res.SIM1 <- foreach( sim = 1:SIM, .combine = 'c', .inorder = FALSE, .multicombine = TRUE, .verbose =  TRUE, .packages = c("label.switching", "rjags", "dclone", "R2jags", "abind", "MCMCvis"), .errorhandling = "pass") %dopar% {
  ###--------------------------------------------------
  # sim = 1
  
  cl <- parallelly::makeClusterPSOCK(n.chains, tries = 20L, delay = 10.0, setup_strategy = "sequential")
  tmp <- parallel::clusterEvalQ(cl, library(dclone))
  
  ## To compute DIC
  tmp <- parLoadModule(cl, "dic")
  tmp <- parLoadModule(cl, "glm")
  
  set.seed(1457+TID+sim*50) # make the simulation reproducible and to vary within each node (TID) and simulation
  
  ### simulate responses
  ### Classify subjects in latent classes 
  # ppi <- matrix(0, nrow=n, ncol=K1)
  # w <- matrix(0, nrow=n, ncol=K1)
  # g <- rep(0,n)   # classes 
  # 
  # for(j in 1:n){
  #   for(k in 1:K1){
  #     ppi[j,k] <- exp(V[j,]%*%alpha[k,])
  #   }
  #   for(k in 1:K1){
  #     w[j,k] <- ppi[j,k]/sum(ppi[j,])   # probability for classes
  #   }
  #   g[j] <- sample((1:K1),1,prob=w[j,])   # class
  # }
  
  ### A bit faster than the above
  ppi <- exp(V %*% t(alpha))
  w <- ppi/rowSums(ppi)
  g0 <- sapply(1:n, function(j){sample((1:K1),1,prob=w[j,])})
  
  clase <- g0
  
  ### Linear predictor and response 
  eta <- matrix(0,n*TT)   # linear predictor
  Wtrue <- rep(NA,n*TT)   # response variable 
  Yerror <- rep(NA,n*TT)
  for(j in 1:n){		
    for(i in offset[j]){	
      eta[i] <- beta0 + Z[i,]%*%gama[j,] + U[i,]%*%lambda[g0[j],] + X[i,]%*%beta 
      Wtrue[i] <- rnorm(1, eta[i], sqrt(tau2))
      Yerror[i] <- rnorm(1, Wtrue[i], sqrt(sigma2[g0[j]]))
    }
    for(i in (offset[j]+1):(offset[j+1]-1)){
      eta[i] <- beta0 + Z[i,]%*%gama[j,] + U[i,]%*%lambda[g0[j],] + X[i,]%*%beta 
      Wtrue[i] <- rnormright(Wtrue[i-1],eta[i],sqrt(tau2)) 
      Yerror[i] <- rnorm(1, Wtrue[i], sqrt(sigma2[g0[j]]))
    }
  }
  
  clase_cat = factor(clase)
  ### { clases K1=3 simula, K=K2=2 estima }
  levels(clase_cat) = as.character(c(1:K1))
  # CLASEtrue[sim,] 
  CLASEtrue = table(clase)
  # table(clase)/n
  
  MONOTON = c(sum(Yerror[-(offset[-(n+1)]+TT-1)][rep(clase,each=5)==1] >= Yerror[-(offset[-(n+1)])][rep(clase,each=5)==1]), 
              sum(Yerror[-(offset[-(n+1)]+TT-1)][rep(clase,each=5)==2] >= Yerror[-(offset[-(n+1)])][rep(clase,each=5)==2]), 
              sum(Yerror[-(offset[-(n+1)]+TT-1)][rep(clase,each=5)==3] >= Yerror[-(offset[-(n+1)])][rep(clase,each=5)==3]),
              sum(Yerror[-(offset[-(n+1)]+TT-1)][rep(clase,each=5)==1] < Yerror[-(offset[-(n+1)])][rep(clase,each=5)==1]), 
              sum(Yerror[-(offset[-(n+1)]+TT-1)][rep(`clase`,each=5)==2] < Yerror[-(offset[-(n+1)])][rep(clase,each=5)==2]), 
              sum(Yerror[-(offset[-(n+1)]+TT-1)][rep(clase,each=5)==3] < Yerror[-(offset[-(n+1)])][rep(clase,each=5)==3]))
  
  ###--------------------------------------------------
  ### JAGS
  
  ### datos de entrada para el modelo JAGS
  data2 <- list(Yerror=Yerror, R2w=(max(Yerror[offset[-(n+1)]])-min(Yerror[offset[-(n+1)]])),   ### response variable 
                X=X, L=ncol(X), zeros.beta=rep(0,ncol(X)), prior.betaB=diag(ncol(X))*0.01,   # fixed effects
                Z=Z, M=ncol(Z), zeros.gama=t(rep(0,ncol(Z))), R.u=diag(1,ncol(Z)),   # random effects
                U=U, P=ncol(U), zeros.lambda=rep(0,ncol(U)), prior.lambdaD=diag(ncol(U))*0.01,   #  class-specific fixed effects 
                V=V, Q=ncol(V),   # latent class membership
                n=n,   # number of subjects 
                offset=offset,   # Indicator for the beginning of observations for each subject for long tables format 
                ### { clases K1=3 simula, K=K2=2 estima               
                K=K2)   ### number of classes ### 
  ### } ### number of classes 
  
  parJagsModel(cl, name="fit1", file="lclmm_Kclasses_error_NOincrease.jag", 
               data=data2, inits=inits1, n.chains=n.chains, n.adapt=n.adapt)
  
  parUpdate(cl, "fit1", n.iter=n.update, thin=n.thin)
  
  sample2 <- parCodaSamples(cl, model="fit1", variable.names=param2, 
                            n.iter=n.iter, thin=n.thin)
  
  ###--------------------------------------------------
  
  ###--------------------------------------------------
  
  ### datos de entrada que necesitamos para el modelo en JAGS 
  ### datos de entrada para el modelo JAGS
  data2lcmm <- list(Wtrue=Yerror, R2w=(max(Yerror[offset[-(n+1)]])-min(Yerror[offset[-(n+1)]])),  ### response variable 
                    X=X, L=ncol(X), zeros.beta=rep(0,ncol(X)), prior.betaB=diag(ncol(X))*0.01,   # fixed effects
                    Z=Z, M=ncol(Z), zeros.gama=t(rep(0,ncol(Z))), R.u=diag(1,ncol(Z)),   # random effects
                    U=U, P=ncol(U), zeros.lambda=rep(0,ncol(U)), prior.lambdaD=diag(ncol(U))*0.01,   #  class-specific fixed effects 
                    V=V, Q=ncol(V),   # latent class membership
                    n=n,   # number of subjects 
                    offset=offset,   # Indicator for the beginning of observations for each subject for long tables format 
                    ### { clases K1=3 simula, K=K2=2 estima 
                    K=K2)     ### number of classes 
  ### }  ### number of classes 
  
  ### Necesitaremos estas cadenas para el label switching
  #sample2lcmm <- coda.samples(fit1lcmm, param2lcmm, n.iter=n.iter, thin=n.thin)
  
  parJagsModel(cl, name="fit1lcmm", file="lclmm_Kclasses_lclmm.jag", 
               data=data2lcmm, inits=inits1lcmm, n.chains=n.chains, n.adapt=n.adapt)
  
  parUpdate(cl, "fit1lcmm", n.iter=n.update, thin=n.thin)
  
  sample2lcmm <- parCodaSamples(cl, model="fit1lcmm", variable.names=param2lcmm, 
                                n.iter=n.iter, thin=n.thin)
  
  ###--------------------------------------------------
  ###--------------------------------------------------
  ### LABEL SWITCHING 
  ###------------------------------------------------
  
  ### DATA and LISTS for Label Switching 
  n <- data2$n   # # sample size 
  K <- data2$K   # # number of components, classes, to estimate
  # produce an MCMC sample 
  m1 <- n.iter/n.thin #  # define number of iterations for each chain
  m2 <- n.chains*n.iter/n.thin #  # define number of iterations for total chains 
  # write output to a convenient format for label.switching package: 
  # MCMC parameters should be saved as m\times K\times J array
  #  # different types of parameters:
  mcmc.pars <- array(data = NA, dim = c(m2, data2$K, data2$P))
  mcmc.alpha <- array(data = NA, dim = c(m2, data2$K, data2$Q))
  mcmc.beta0 <- array(data = NA, dim = c(m2))
  mcmc.beta <- array(data = NA, dim = c(m2, data2$L))
  mcmc.tau2 <- array(data = NA, dim = c(m2))
  mcmc.sigma2 <- array(data = NA, dim = c(m2,data2$K))
  zclass <- array(0,dim=c(m2,n))   ### identifica la clase latente a la que pertenece el sujeto
  
  for(chain in 1:n.chains){
    for(h1 in 1:data2$K){
      for(h2 in 1:data2$P){
        mcmc.pars[(chain-1)*m1+1:m1, h1,h2] <- sample2[[chain]][,paste0("lambda[",h1,",",h2,"]")] ### lambda
      }
      for(h2 in 1:data2$Q){
        mcmc.alpha[(chain-1)*m1+1:m1, h1,h2] <- sample2[[chain]][,paste0("alpha[",h1,",",h2,"]")] ### alpha
      }
      mcmc.sigma2[(chain-1)*m1+1:m1, h1] <- sample2[[chain]][,paste0("sigma2[",h1,"]")] ### alpha
    }
    for(h1 in 1:data2$L){
      mcmc.beta[(chain-1)*m1+1:m1, h1] <- sample2[[chain]][,paste0("beta[",h1,"]")] ### beta
    }
    mcmc.tau2[(chain-1)*m1+1:m1 ] <- sample2[[chain]][,"tau2"] ### tau2
    mcmc.beta0[(chain-1)*m1+1:m1 ] <- sample2[[chain]][,"beta0"] ### beta0
    for(i in 1:n){
      zclass[(chain-1)*m1+1:m1,i] <- sample2[[chain]][,paste0("g[",i,"]")]
    }
  }
  
  # Compute necessary information that will be additional input to the
  # label.switching package.  Define the complete log-likelihood function 
  # esta verosimilitud es aproximada
  
  ### computing complete loglikelihoods (in parallel) in order to find pivot
  nil <- parallel::clusterEvalQ(cl, library(mnormt))
  logLvec <- parallel::parSapply(cl, 1:m2, complete.loglikelihood, zclass=zclass, data2 = data2,
                                 mcmc.pars = mcmc.pars, mcmc.beta0 = mcmc.beta0, mcmc.beta = mcmc.beta, mcmc.alpha = mcmc.alpha, 
                                 mcmc.tau2 = mcmc.tau2, mcmc.sigma2=mcmc.sigma2)
  
  # computing complete loglikelihoods in order to find pivot
  zmapindex <- fn.zmapindex(zclass, logLvec, m2)
  mapindex <- zmapindex$mapindex
  zmap <- zmapindex$zmap 
  maxL <- zmapindex$maxL 
  
  # computing allocation probabilities for stephens method and ECR-ITERATIVE-2
  piter <- parallel::parLapply(cl, 1:m2, fn.probabilities, data2=data2, TT = TT, 
                               mcmc.pars = mcmc.pars, mcmc.tau2 = mcmc.tau2, 
                               mcmc.beta0 = mcmc.beta0, mcmc.beta = mcmc.beta, mcmc.alpha = mcmc.alpha, 
                               mcmc.sigma2 = mcmc.sigma2) 
  p <- do.call(abind, list(...=piter, along=1))
  
  # Run label.switching command using all methods.
  
  # The pivot for default ECR algorithm will be the allocation `mapindex` that
  # corresponds to the maximum of complete logL estimate: zpivot=z[mapindex,].
  
  # The pivot for PRA algorithm will be the parameters that correspond to the same
  # iteration: prapivot = mcmc.pars[mapindex,,].
  
  # The SJW method will be initialized using this iteration as well: sjwinit =
  # mapindex.  The complete log-likelihood is defined as: complete =
  # complete.normal.loglikelihood
  
  
  ls <- label.switching(method = set, 
                        zpivot = zmap, 
                        z = zclass, 
                        K = data2$K, ### { clases K1=3 simula, K=K2=2 estima }
                        prapivot = mcmc.pars[mapindex, , ], 
                        p = p, 
                        # complete = complete.loglikelihood, 
                        mcmc = mcmc.pars, 
                        # data = data2$Yerror, 
                        constraint=2, ### 3=pendiente del tiempo
                        sjwinit = mapindex
  )
  
  ###--------------------------------------------------
  
  ### DATA and LISTS for Label Switching 
  n <- data2lcmm$n   # # sample size 
  K <- data2lcmm$K   # # number of components, classes
  # produce an MCMC sample 
  m1 <- n.iter/n.thin #  # define number of iterations for each chain
  m2 <- n.chains*n.iter/n.thin #  # define number of iterations for total chains 
  # write output to a convenient format for label.switching package: 
  # MCMC parameters should be saved as m\times K\times P array
  # # different types of parameters:
  mcmc.pars <- array(data = NA, dim = c(m2, data2lcmm$K, data2lcmm$P)) ### { clases K1=3 simula, K=K2=2 estima } 
  mcmc.alpha <- array(data = NA, dim = c(m2, data2lcmm$K,data2lcmm$Q)) ### { clases K1=3 simula, K=K2=2 estima } 
  mcmc.beta0 <- array(data = NA, dim = c(m2))
  mcmc.beta <- array(data = NA, dim = c(m2, data2lcmm$L))
  mcmc.tau2 <- array(data = NA, dim = c(m2))
  zclass <- array(0,dim=c(m2,n))   ### identifica la clase latente a la que pertenece el sujeto
  for(chain in 1:n.chains){
    for(h1 in 1:data2lcmm$K){
      for(h2 in 1:data2lcmm$P){
        mcmc.pars[(chain-1)*m1+1:m1, h1,h2] <- sample2lcmm[[chain]][,paste0("lambda[",h1,",",h2,"]")] ### lambda
      }
      for(h2 in 1:data2lcmm$Q){
        mcmc.alpha[(chain-1)*m1+1:m1, h1,h2] <- sample2lcmm[[chain]][,paste0("alpha[",h1,",",h2,"]")] ### alpha
      }
    }
    for(h1 in 1:data2lcmm$L){
      mcmc.beta[(chain-1)*m1+1:m1, h1] <- sample2lcmm[[chain]][,paste0("beta[",h1,"]")] ### beta
    }
    mcmc.beta0[(chain-1)*m1+1:m1 ] <- sample2lcmm[[chain]][,"beta0"] ### beta0
    mcmc.tau2[(chain-1)*m1+1:m1 ] <- sample2lcmm[[chain]][,"tau2"] ### tau2
    for(i in 1:n){
      zclass[(chain-1)*m1+1:m1,i] <- sample2lcmm[[chain]][,paste0("g[",i,"]")]
    }
  }
  
  ### computing complete loglikelihoods in order to find pivot
  
  ### computing complete loglikelihoods (in parallel) in order to find pivot
  logLvec.lcmm <- parallel::parSapply(cl, 1:m2, complete.loglikelihood.lcmm, zclass=zclass, data2lcmm = data2lcmm,
                                      mcmc.pars = mcmc.pars, mcmc.beta0 = mcmc.beta0, mcmc.beta = mcmc.beta, mcmc.alpha = mcmc.alpha, 
                                      mcmc.tau2 = mcmc.tau2)
  
  # Function below finds the pivot
  zmapindex.lcmm <- fn.zmapindex(zclass, logLvec.lcmm, m2)  
  mapindex <- zmapindex.lcmm$mapindex
  zmap <- zmapindex.lcmm$zmap 
  maxL <- zmapindex.lcmm$maxL 
  
  # computing allocation probabilities for stephens method and ECR-ITERATIVE-2
  # no_cores <- detectCores()-1 # Calculate the number of cores
  # cl <- makeCluster(no_cores) # Initiate cluster
  # nil <- clusterEvalQ(cl, library(mnormt))
  # clusterExport(cl,c("data2lcmm", "TT", 
  #                    "mcmc.pars", "mcmc.tau2", 
  #                    "mcmc.beta", "mcmc.alpha", "mcmc.Gama0")) 
  piter <- parallel::parLapply(cl, 1:m2, fn.probabilities.lcmm, data2lcmm=data2lcmm, TT = TT, 
                               mcmc.pars = mcmc.pars, mcmc.tau2 = mcmc.tau2, 
                               mcmc.beta0 = mcmc.beta0, mcmc.beta = mcmc.beta, mcmc.alpha = mcmc.alpha) 
  plcmm <- do.call(abind, list(...=piter, along=1))
  
  
  
  lslcmm <- label.switching(method = set, 
                            zpivot = zmap, 
                            z = zclass, 
                            K = data2lcmm$K, ### { clases K1=3 simula, K=K2=2 estima 
                            prapivot = mcmc.pars[mapindex, , ], 
                            p = plcmm, 
                            # complete = complete.loglikelihood.lcmm, 
                            mcmc = mcmc.pars, 
                            # data = data2lcmm$Wtrue, 
                            constraint=2, ### 2=pendiente del tiempo
                            sjwinit = mapindex
  )
  
  parallel::stopCluster(cl)
  
  ###--------------------------------------------------
  ### USAR JAGS y la CLASE ESTIMADA por label.switching 
  
  ### Elegir las clases de alguno de los metodos
  ### STEPHENS, PRA, ECR, ECR-ITERATIVE-1, ECR-ITERATIVE-2, AIC
  
  # table(clase,lcmm_error_class2$pprob$class)
  # table(clase,lcmm_error_class3$pprob$class)
  # table(clase,lcmm_error_class4$pprob$class)
  
  claseSTE = etiquetar(clase,ls$clusters["STEPHENS",])
  clasePRA = etiquetar(clase,ls$clusters["PRA",])
  claseECR = etiquetar(clase,ls$clusters["ECR",])
  claseEC1 = etiquetar(clase,ls$clusters["ECR-ITERATIVE-1",])
  claseEC2 = etiquetar(clase,ls$clusters["ECR-ITERATIVE-2",])
  claseAIC = etiquetar(clase,ls$clusters["AIC",])
  
  claseEC1_cat = factor(claseEC1)
  levels(claseEC1_cat) = as.character(c(1:K))
  
  # CLASES[sim,] 
  CLASES = as.vector(table(clase_cat,claseEC1_cat))
  
  
  ###--------------------------------------------------
  
  claselcmmSTE = etiquetar(clase,lslcmm$clusters["STEPHENS",])
  claselcmmPRA = etiquetar(clase,lslcmm$clusters["PRA",])
  claselcmmECR = etiquetar(clase,lslcmm$clusters["ECR",])
  claselcmmEC1 = etiquetar(clase,lslcmm$clusters["ECR-ITERATIVE-1",])
  claselcmmEC2 = etiquetar(clase,lslcmm$clusters["ECR-ITERATIVE-2",])
  claselcmmAIC = etiquetar(clase,lslcmm$clusters["AIC",])
  
  claselcmmEC1_cat = factor(claselcmmEC1)
  levels(claselcmmEC1_cat) = as.character(c(1:K))
  
  # CLASESlcmm[sim,] 
  CLASESlcmm = as.vector(table(clase_cat,claselcmmEC1_cat))
  
  data1 <- data.frame("id"=id, "Wtrue"=Wtrue, "Yerror"=Yerror,
                      "time"=time, "time2"=time2, 
                      "timejit"=time+rep(jitter(rep(0,n),factor=10),each=TT), 
                      "x1"=X1, "x2"=X2, #"x3"=X3,
                      "v1"=V1[id], "v2"=V2[id], 
                      "g"=clase[id])
  
  return( list(list(CLASEtrue = CLASEtrue, MONOTON = MONOTON,
                    CLASES = CLASES, CLASESlcmm =CLASESlcmm,
                    data1 = data1, ls = ls, lslcmm = lslcmm)) )
  ###--------------------------------------------------
} 
# str(res.SIM1, 2)

### End FOR simulation
cat("Done with Step 1; Elapsed time: ", difftime(Sys.time(), ptm,units = "hours")," hours\n", sep = "")

stopCluster(cl0)

cl1 <- parallel::makeCluster(32)
registerDoParallel(cl1)

### Performs MCMC considering label switching with the chains disaggregated across all cores
res.SIM2a <- foreach( sim2 = 1:(SIM*n.chains), .combine = 'c', .inorder = FALSE, .multicombine = TRUE, .verbose =  TRUE, .packages = c("rjags", "R2jags"), .errorhandling = "pass" ) %dopar% {
  
  idx <- sim2%%SIM + 1
  
  idxdt <- res.SIM1[[idx]]
  dt <- idxdt$data1
  
  set.seed(1457+TID+sim2*50)
  
  data3 <- list(Yerror=dt$Yerror, R2w=(max(dt$Yerror)-min(dt$Yerror)), 
                X=X, L=ncol(X), zeros.beta=rep(0,ncol(X)), prior.betaB=diag(ncol(X))*0.01,   ### fixed effects
                Z=Z, M=ncol(Z), zeros.gama=t(rep(0,ncol(Z))), R.u=diag(1,ncol(Z)),   ### random effects
                U=U, P=ncol(U), zeros.lambda=rep(0,ncol(U)), prior.lambdaD=diag(ncol(U))*0.01,   ### class-specific fixed effects 
                n=n,
                offset=offset,   # Indicator for the beginning of observations for each subject for long tables format 
                K=K2,   ### number of classes 
                g=idxdt$ls$clusters["ECR-ITERATIVE-1",])   ### indicator for classes
  attach(data3)
  jagsfit3 <- tryCatch(jags(
    data = c(names(data3), "data3", "n.iter", "n.thin"), 
    inits = inits3, 
    parameters.to.save = c(param3, "LogLik"), 
    n.iter = n.iter, 
    n.chains = 1,
    n.thin = n.thin,
    jags.module = c("glm","dic"), DIC=TRUE, ### default
    progress.bar = "none",
    model.file="lclmm_Kclasses_error_NOincrease_labelswitching.jag"), error = function(e){message(e); return(NULL)})
  
  
  data3lcmm <- list(Wtrue=dt$Yerror, R2w=(max(dt$Yerror)-min(dt$Yerror)),   ### response variable
                    X=X, L=ncol(X), zeros.beta=rep(0,ncol(X)), prior.betaB=diag(ncol(X))*0.01,   ### fixed effects
                    Z=Z, M=ncol(Z), zeros.gama=t(rep(0,ncol(Z))), R.u=diag(1,ncol(Z)),   ### random effects
                    U=U, P=ncol(U), zeros.lambda=rep(0,ncol(U)), prior.lambdaD=diag(ncol(U))*0.01,   ### class-specific fixed effects 
                    n=n,
                    offset=offset,   # Indicator for the beginning of observations for each subject for long tables format 
                    K=K2,   ### number of classes 
                    g=idxdt$lslcmm$clusters["ECR-ITERATIVE-1",])   ### indicator for classes
  attach(data3lcmm)
  jagsfit3lcmm <- tryCatch(jags(
    data = c(names(data3lcmm), "data3lcmm", "n.iter", "n.thin"), 
    inits = inits3lcmm, 
    parameters.to.save = c(param3lcmm, "LogLik"), 
    n.iter = n.iter, 
    n.chains = 1,
    n.thin = n.thin,
    jags.module = c("glm","dic"), DIC=TRUE, ### default
    progress.bar = "none",
    model.file="lclmm_Kclasses_lclmm_labelswitching.jag"), error = function(e){message(e); return(NULL)})
  
  return( list(list(sim2=sim2, idx=idx, 
                    jagsfit3 = jagsfit3, jagsfit3lcmm = jagsfit3lcmm)) )
}

cat("Done with Step 2; Elapsed time: ", difftime(Sys.time(), ptm,units = "hours")," hours\n", sep = "")
# str(res.SIM2a, 2)
classres.SIM2a <- sapply(res.SIM2a, class)

idxgroups <- tryCatch(sapply(res.SIM2a[classres.SIM2a=="list"],`[[`, 2), error = function(e){message(e); return(NULL)})

if( is.null(idxgroups) ){
  save(res.SIM2a, file=paste0(savepath,"/res.SIMa_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".Rdata"))
  stop("Issue with extracting indices from res.SIM2a.")
}

### Integrates the chain information across samples
res.SIM2 <- foreach( idxsim = 1:SIM, .combine = 'c', .inorder = FALSE, .multicombine = TRUE, .verbose =  TRUE, .packages = c("rjags", "R2jags", "abind", "loo"), .errorhandling = "pass") %dopar% {
  
  res <- res.SIM2a[classres.SIM2a=="list"][which(idxgroups==idxsim)]
  n.chains0 <- sum(idxgroups==idxsim)
  
  #####  jagsfit3
  ### combine chains following jags.parallel
  result <- NULL
  model <- NULL
  for (ch in 1:n.chains0) {
    result <- abind(result, res[[ch]]$jagsfit3$BUGSoutput$sims.array, 
                    along = 2)
    model[[ch]] <- res[[ch]]$jagsfit3$model
  }
  result <- R2jags:::as.bugs.array2(result, model.file="lclmm_Kclasses_error_NOincrease_labelswitching.jag", program = "jags", DIC = TRUE, n.iter = n.iter, n.burnin = floor(n.iter/2), n.thin = n.thin)
  
  out <- list(model = model, BUGSoutput = result, parameters.to.save = param3, model.file = "lclmm_Kclasses_error_NOincrease_labelswitching.jag", n.iter = n.iter, DIC = TRUE)
  class(out) <- c("rjags.parallel", "rjags")
  
  ### same as coda output
  # print(MCMCsummary(sample3lcmm))
  # print(MCMCsummary(sample3lcmm, params = param3lcmm, digits=3))
  # plot(sample3lcmm[,!grepl("deviance|LogLik", varnames(sample3lcmm))])
  # MCMCplot(sample3lcmm, params = param3lcmm[-length(param3)])
  
  sample3 <- as.mcmc.list(as.mcmc(out))
  all <- as.matrix(sample3, iters = TRUE, chains=TRUE)
  pD <- var(all[,"deviance"])/2
  DICall <- mean(all[,"deviance"]) + pD
  # cat("Sample3:")
  # print(c(DIC=DICall,pD=pD))
  
  loglik_prop <- out$BUGSoutput$sims.list$LogLik
  waic_prop <- waic(loglik_prop)
  loo_prop <- loo(loglik_prop)
  
  DIC1 = out$BUGSoutput$DIC
  WAIC_ic = waic_prop$estimates[rownames(waic_prop$estimates)=="waic",]
  LOO_ic = loo_prop$est[rownames(loo_prop$est)=="looic",]
  ELPDLOO_ic = loo_prop$est[rownames(loo_prop$est)=="elpd_loo",]
  
  PARAM = cbind(summary(sample3[,!grepl("deviance|LogLik", varnames(sample3))])[1][[1]][,1:2],
                summary(sample3[,!grepl("deviance|LogLik", varnames(sample3))])[2][[1]][,c(1,3,5)])
  
  
  #####  jagsfit3lcmm
  ### combine chains following jags.parallel
  resultlcmm <- NULL
  modellcmm <- NULL
  for (ch in 1:n.chains0) {
    resultlcmm <- abind(resultlcmm, res[[ch]]$jagsfit3lcmm$BUGSoutput$sims.array, 
                        along = 2)
    modellcmm[[ch]] <- res[[ch]]$jagsfit3lcmm$model
  }
  resultlcmm <- R2jags:::as.bugs.array2(resultlcmm, model.file="lclmm_Kclasses_lclmm_labelswitching.jag", program = "jags", DIC = TRUE, n.iter = n.iter, n.burnin = floor(n.iter/2), n.thin = n.thin)
  
  outlcmm <- list(model = modellcmm, BUGSoutput = resultlcmm, parameters.to.save = param3lcmm, model.file = "lclmm_Kclasses_lclmm_labelswitching.jag", n.iter = n.iter, DIC = TRUE)
  class(outlcmm) <- c("rjags.parallel", "rjags")
  
  ### same as coda output
  sample3lcmm <- as.mcmc.list(as.mcmc(outlcmm))
  
  # print(MCMCsummary(sample3lcmm))
  # print(MCMCsummary(sample3lcmm, params = param3lcmm, digits=3))
  # plot(sample3lcmm[,!grepl("deviance|LogLik", varnames(sample3lcmm))])
  # MCMCplot(sample3lcmm, params = param3lcmm[-length(param3)])
  
  all <- as.matrix(sample3lcmm, iters = TRUE, chains=TRUE)
  pD <- var(all[,"deviance"])/2
  DIClcmmall <- mean(all[,"deviance"]) + pD
  # cat("Sample3:")
  # print(c(DIC=DIClcmmall,pD=pD))
  
  ### { Codigo WAIC LOO
  loglik_lcmm <- outlcmm$BUGSoutput$sims.list$LogLik
  waic_lcmm <- waic(loglik_lcmm)
  loo_lcmm <- loo(loglik_lcmm)
  
  DIC1lcmm = outlcmm$BUGSoutput$DIC
  WAIC_lcmm_ic = waic_lcmm$estimates[rownames(waic_lcmm$estimates)=="waic",]
  LOO_lcmm_ic = loo_lcmm$est[rownames(loo_lcmm$est)=="looic",]
  ELPDLOO_lcmm_ic = loo_lcmm$est[rownames(loo_lcmm$est)=="elpd_loo",]
  ### } 
  
  # PARAMlcmm[sim,,1:2] = summary(sample3lcmm)[1][[1]][,1:2]
  # PARAMlcmm[sim,,3:5] = summary(sample3lcmm)[2][[1]][,c(1,3,5)]
  PARAMlcmm = cbind(summary(sample3lcmm[,!grepl("deviance|LogLik", varnames(sample3lcmm))])[1][[1]][,1:2],
                    summary(sample3lcmm[,!grepl("deviance|LogLik", varnames(sample3lcmm))])[2][[1]][,c(1,3,5)])
  
  return( list(list(idx=idxsim, PARAM = PARAM, PARAMlcmm = PARAMlcmm,
                    jagsfit3 = out, jagsfit3lcmm = outlcmm,
                    sample3 = sample3, sample3lcmm = sample3lcmm,
                    DIC = DICall, DIClcmm = DIClcmmall,
                    DIC1 = DIC1, DIC1lcmm = DIC1lcmm,
                    ELPD = ELPDLOO_ic, ELPDlcmm = ELPDLOO_lcmm_ic,
                    LOO = LOO_ic, LOOlcmm = LOO_lcmm_ic,
                    WAIC = WAIC_ic, WAIClcmm = WAIC_lcmm_ic)) )
}

cat("Done with Step 3; Elapsed time: ", difftime(Sys.time(), ptm,units = "hours")," hours\n", sep = "")
# str(res.SIM2, 2)
stopCluster(cl1)

save(res.SIM1, res.SIM2, idxgroups, file=paste0(savepath,"/res.SIM_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".Rdata"))

for(sim in 1:SIM){
  cat("****************************\n",
      "****** Iteration: ", sim," *******\n",
      "****************************\n", sep="")
  MONOTON[sim,] = res.SIM1[[sim]]$MONOTON
  CLASEtrue[sim,] = res.SIM1[[sim]]$CLASEtrue
  CLASES[sim,] = res.SIM1[[sim]]$CLASES
  CLASESlcmm[sim,] = res.SIM1[[sim]]$CLASESlcmm
  PARAM[sim,,] = res.SIM2[[sim]]$PARAM
  PARAMlcmm[sim,,] = res.SIM2[[sim]]$PARAMlcmm
  DIC_all[sim,1] = res.SIM2[[sim]]$DIC
  DIClcmm_all[sim,1] = res.SIM2[[sim]]$DIClcmm
  DIC1_all[sim,1] = res.SIM2[[sim]]$DIC1
  DIC1lcmm_all[sim,1] = res.SIM2[[sim]]$DIC1lcmm
  WAIC_all[sim,] = res.SIM2[[sim]]$WAIC
  WAIClcmm_all[sim,] = res.SIM2[[sim]]$WAIClcmm
  LOO_all[sim,] = res.SIM2[[sim]]$LOO
  LOOlcmm_all[sim,] = res.SIM2[[sim]]$LOOlcmm
  ELPD_all[sim,] = res.SIM2[[sim]]$ELPD
  ELPDlcmm_all[sim,] = res.SIM2[[sim]]$ELPDlcmm
  
  cat("**** Proposal ****\n")
  print(MCMCsummary(res.SIM2[[sim]]$sample3, params = param3))
  cat("**** LCLMM ****\n")
  print(MCMCsummary(res.SIM2[[sim]]$sample3lcmm, params = param3lcmm))
}
###--------------------------------------------------

write.csv(CLASEtrue, paste0(savepath,"/CLASEtrue_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".csv"))
write.csv(MONOTON, paste0(savepath,"/MONOTON_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".csv"))
write.csv(CLASES, paste0(savepath,"/CLASES_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".csv"))
write.csv(CLASESlcmm, paste0(savepath,"/CLASESlcmm_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".csv"))
write.csv(apply(PARAM,c(2,3),"mean"), paste0(savepath,"/PARAM_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".csv"))
write.csv(apply(PARAMlcmm,c(2,3),"mean"), paste0(savepath,"/PARAMlcmm_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".csv"))

write.csv(DIC_all, paste0(savepath,"/DIC_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".csv"))
write.csv(DIClcmm_all, paste0(savepath,"/DIClcmm_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".csv"))

write.csv(DIC1_all, paste0(savepath,"/DIC1_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".csv"))
write.csv(DIC1lcmm_all, paste0(savepath,"/DIC1lcmm_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".csv"))

### { Code WAIC LOO
write.csv(WAIC_all, paste0(savepath,"/WAIC_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".csv"))
write.csv(WAIClcmm_all, paste0(savepath,"/WAIClcmm_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".csv"))
write.csv(LOO_all, paste0(savepath,"/LOO_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".csv"))
write.csv(LOOlcmm_all, paste0(savepath,"/LOOlcmm_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".csv"))
write.csv(ELPD_all, paste0(savepath,"/ELPD_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".csv"))
write.csv(ELPDlcmm_all, paste0(savepath,"/ELPDlcmm_n=", n,"_K2=",K2,"_TT=",TT,"_TID",TID,".csv"))
### }


cat("Job completed; Elapsed time: ", difftime(Sys.time(), ptm,units = "hours")," hours\n", sep = "")

