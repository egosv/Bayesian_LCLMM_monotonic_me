###--------------------------------------------------
# Latent Class Linear Mixed Models
# Growth Mixture Models
# Continuous response
# K classes
###--------------------------------------------------

### Change this path to read the JAGS code 
setwd("HERE")

savepath <- "SAVEPATH" # where the results are stored

load(file="OAI_data_long.RData")
### loads data.frame OAI_KL_JSW_data obtained from file 1_Compile_KLData.R

### Set K (defaults to 2 if not specified)
if( !exists("K") ){
  K <- 2
}
print(K)

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

library(foreach)
library(doParallel)
### Creating clusters to run a chain in each core
detectCores()
n.cores <- 40 ##detectCores()-2
#options(mc.cores = parallel::detectCores()-2)
options(mc.cores = n.cores)

cl0 <- makeCluster(n.cores-1)

n.chains = 5 # detectCores()-1 # 10
cl <- makeCluster(n.chains)
tmp <- clusterEvalQ(cl, library(dclone))

## To compute DIC
tmp <- parLoadModule(cl, "dic")

###--------------------------------------------------
### JAGS Specifications 
# 
n.adapt = 30000
n.update = 30000
n.iter = 30000
n.thin = 15

# n.adapt = 500
# n.update = 500
# n.iter = 500
# n.thin = 2


###--------------------------------------------------
### Functions for LABEL SWITCHING proposal
{
# Compute necessary information that will be additional input to the label.switching package.  
# Define the complete log-likelihood function. 

complete.loglikelihood <- function(iter, zclass, data2,
                                   mcmc.pars, mcmc.beta0, mcmc.beta, mcmc.alpha, 
                                   mcmc.tau2, mcmc.sigma2, mcmc.Gama0) { 
  
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
  Gama0 = mcmc.Gama0[iter,,]
  sigma2 = mcmc.sigma2[iter]
  
  lambda <- pars[,1:ncol(UU)]
  m <- length(Yerror)
  ###  logl <- rep(0, m)
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
      ###      eta[i] <- beta0 + XX[i,]%*%beta + UU[i,]%*%lambda[zclass[j],]
      media <- eta2[i]
      eta[i] <- media
###      varianza <- ZZ[i,]%*%Gama0%*%ZZ[i,] + tau2 + sigma2[zclass[j]]
      varianza <- tau2 + sigma2[zclass[j]]
      logl[i] <- dnorm(Yerror[i],mean=media,sd=sqrt(varianza), log=T)
    }
    for(i in (offset[j]+1):(offset[j+1]-1)){
      ###      eta[i] <- beta0 + XX[i,]%*%beta + UU[i,]%*%lambda[zclass[j],]
      beta2 <- (eta[i-1]-eta2[i])/sqrt(tau2)
      z2 <- pnorm(beta2)
      medtrunc <- sqrt(tau2)*(-dnorm(beta2))/z2 ### media de una normal truncada 
      medaux <- ifelse(is.finite(medtrunc),medtrunc,0)  
      media <- eta2[i]+medaux 
      eta[i] <- media
      vartrunc <- (1+ (-beta2*dnorm(beta2)/z2) -((-dnorm(beta2)/z2)^2)) ### varianza de una normal truncada
      varaux <- tau2*ifelse(is.finite(vartrunc) & vartrunc>0,vartrunc,1)
###      varianza <- ZZ[i,]%*%Gama0%*%ZZ[i,] + varaux + sigma2[zclass[j]]
      varianza <- varaux + sigma2[zclass[j]]
      logl[i] <- dnorm(Yerror[i], mean=media, sd=sqrt(varianza), log=T)
    }
  } 
  return(sum(logw)+sum(logl,na.rm=TRUE))
}


### computing allocation probabilities for stephens method and ECR-ITERATIVE-2
fn.probabilities <- function(iter, data2, TT, mcmc.pars, 
                             mcmc.tau2, mcmc.beta0, mcmc.beta, mcmc.alpha, mcmc.sigma2, mcmc.Gama0){   
  n <- data2$n
  K <- data2$K
  m <- length(data2$Yerror)  
  #  like <- rep(0, m)
  Yerror <- data2$Yerror
  XX <- data2$X
  ZZ <- data2$Z
  UU <- data2$U
  VV <- data2$V
  offset <- data2$offset
  
  proba <- array(data = NA, dim = c(1, n, K))
#  for(iter in 1:m2){
    tau2 <- mcmc.tau2[iter]
    lambda <- mcmc.pars[iter,,]
    beta0 <- mcmc.beta0[iter]
    beta <- mcmc.beta[iter,]
    alpha <- mcmc.alpha[iter,,]
    Gama0 <- mcmc.Gama0[iter,,]
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
          ###          eta[i,k] <- beta0 + XX[i,]%*%beta + UU[i,]%*%lambda[k,]
          media <- eta2[i,k]
          eta[i,k] <- media
###          varianza <- ZZ[i,]%*%Gama0%*%ZZ[i,] + tau2 + sigma2[k]
          varianza <- tau2 + sigma2[k]
          like[i,k] <- dnorm(Yerror[i],mean=media,sd=sqrt(varianza), log=FALSE)
        } ### END i in offset[j]
        for(i in (offset[j]+1):(offset[j+1]-1)){
          ###          eta[i,k] <- beta0 + XX[i,]%*%beta + UU[i,]%*%lambda[k,]
          
          beta2 <- (eta[i-1,k]-eta2[i,k])/sqrt(tau2)
          z2 <- pnorm(beta2)
          medtrunc <- sqrt(tau2)*(-dnorm(beta2))/z2 ### media de una normal truncada 
          medaux <- ifelse(is.finite(medtrunc),medtrunc,0)  
          media <- eta2[i,k]+medaux 
          eta[i,k] <- media
          vartrunc <- (1+ (-beta2*dnorm(beta2)/z2) -((-dnorm(beta2)/z2)^2)) ### varianza de una normal truncada
          varaux <- tau2*ifelse(is.finite(vartrunc) & vartrunc>0,vartrunc,1)
###          varianza <- ZZ[i,]%*%Gama0%*%ZZ[i,] + varaux + sigma2[k]
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
                                        mcmc.tau2, mcmc.Gama0) { 
  
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
  Gama0 = mcmc.Gama0[iter,,]
  
  lambda <- pars[,1:ncol(UU)]
  m <- length(Wtrue)
  ### logl <- rep(0, m)
  ppi <- matrix(NA,n,K)
  w <- matrix(NA,n,K)
  ###   eta <- rep(NA,m) 
  #   eta <- rep(NA,n) 
  logw <- rep(NA,n)
  logl <- rep(NA,n)
  for(j in 1:n){
    for(k in 1:K){ 
      ppi[j,k] <- exp(VV[j,]%*%alpha[k,])
    } 
    w[j,] <- ppi[j,]/sum(ppi[j,])
    logw[j] <- log(w[j,zclass[j]])  
    #      for(i in (offset[j]):(offset[j+1]-1)){
    #         eta[i] <- beta0 + XX[i,]%*%beta + UU[i,]%*%lambda[zclass[j],] 
    ###         logl[i] <- dnorm(Wtrue[i],mean=eta[i],sd=sqrt(tau2), log=T)
    #      }
    idx <- (offset[j]):(offset[j+1]-1)
    media <- beta0 + XX[idx,]%*%beta + UU[idx,]%*%lambda[zclass[j],] 
    ###      media <- beta0 + XX[idx,]%*%beta + ZZ[idx,]%*%gama.est[j,] + UU[idx,]%*%lambda[zclass[j],] 
    varianza <- diag(tau2,length(idx)) ####+ ZZ[idx,]%*%Gama0%*%t(ZZ[idx,]) 
    ###      varianza <- diag(tau2,length(idx)) 
    varianza <- (varianza+t(varianza))/2
    logl[j] <- dmnorm(Wtrue[idx], mean=t(media),varcov=varianza, log=T)
  } 
  return(sum(logw)+sum(logl,na.rm=TRUE))
}

# computing allocation probabilities for stephens method and ECR-ITERATIVE-2
fn.probabilities.lcmm <- function(iter, data2lcmm, TT, 
                                  mcmc.pars, mcmc.tau2, mcmc.beta0, mcmc.beta, mcmc.alpha, mcmc.Gama0){   
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
#  for(iter in 1:m2){
    tau2 <- mcmc.tau2[iter]
    lambda <- mcmc.pars[iter,,]
    beta0 <- mcmc.beta0[iter]
    beta <- mcmc.beta[iter,]
    alpha <- mcmc.alpha[iter,,]
    Gama0 <- mcmc.Gama0[iter,,]
    ppi <- matrix(NA,n,K)
    w <- matrix(NA,n,K)
    #      like <- array(NA,dim=c(m,K))
    like <- array(NA,dim=c(n,K))
    #     eta <- array(NA,dim=c(m,K))
    #      num = array(NA,dim=c(n,K))
    for(j in 1:n){
      idx <- (offset[j]):(offset[j+1]-1)
      for(k in 1:K){
        ppi[j,k] <- exp(VV[j,]%*%alpha[k,])
        media <- beta0 + XX[idx,]%*%beta + UU[idx,]%*%lambda[k,] 
        ###            media <- beta0 + XX[idx,]%*%beta + ZZ[idx,]%*%gama.est[j,] + UU[idx,]%*%lambda[k,] 
        varianza <- diag(tau2,length(idx)) ###+ ZZ[idx,]%*%Gama0%*%t(ZZ[idx,])  
        ###            varianza <- diag(tau2,length(idx)) 
        varianza <- (varianza+t(varianza))/2
        like[j,k] <- dmnorm(Wtrue[idx], mean=t(media),varcov=varianza)
        
        #            for(i in (offset[j]):(offset[j+1]-1)){
        #               eta[i,k] <- beta0 + XX[i,]%*%beta + UU[i,]%*%lambda[k,]
        #               like[i,k] <- dnorm(Wtrue[i],mean=eta[i,k],sd=sqrt(tau2))
        #            }
        #            num[j,k] <- prod(like[(offset[j]):(offset[j+1]-1),k])
      }
      w[j,] <- ppi[j,]/sum(ppi[j,])
      #         proba[iter, j, ] <- (num[j,]*w[j,])/sum(num[j,]*w[j,])     
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


### Data pre-processing
outcome <- "MCMJSWsum"

colskeep <- c("SEX","V00COHORT","ID","VISIT","MCMJSWsum","MCMJSWmin","TIME","BMI_BL","AGE_BL", "WOMTSmax_BL","BMI","WOMTSmax") # 

data0 <- OAI_KL_JSW_data[complete.cases(OAI_KL_JSW_data[,colskeep]),colskeep]

data0 <- droplevels(subset(data0, V00COHORT=="1: Progression"))
# table(data0$VISIT)

### Keep records with complete visits data
IDs_complete <- do.call(rbind, by(data0, INDICES = list(data0$ID), function(dt){
  data.frame(ID=dt$ID[1], NRECORDS=nrow(dt), RECORDEDVISITS=paste0(dt$VISIT, collapse="-"))
}))
# table(IDs_complete$NRECORDS)
# table(IDs_complete$RECORDEDVISITS)
# table(IDs_complete$RECORDEDVISITS,IDs_complete$NRECORDS)

data0a <- subset( data0, ID %in% IDs_complete$ID[IDs_complete$NRECORDS==7] )

### Convert to numeric
data0b <- data.frame(data0a[,-c(1:2)],  model.matrix(~SEX, data0a)[,-1, drop=FALSE])
data0b$TIMEyr <- data0b$TIME/12
data0b$TIMEyr2 <- data0b$TIMEyr^2

### calculate offset
offset <- c(1)
for(i in 2:nrow(data0b)){
  if( data0b$ID[i]==data0b$ID[i-1] ){
    next
  } else{
    offset <- c(offset,i)
  }
}
offset <- c(offset,nrow(data0b)+1)


### data needed for the JAGS models
Yerror = data0b[,outcome]
X = cbind(data0b$AGE_BL, data0b$SEXM) # overall fixed effects
Z = cbind(rep(1,nrow(data0b)),data0b$TIMEyr) ### random effects
U = cbind(rep(1,nrow(data0b)),data0b$TIMEyr,data0b$TIMEyr2)   ### class-specific fixed effects
V = cbind(data0b$BMI, data0b$WOMTSmax) # class-membership effects

n = length(unique(data0b$ID))
TT = 7

ID0 <- data.frame(ID=1:n)
rownames(ID0) <- unique(data0b$ID)
ID1 = ID0[as.character(data0b$ID),"ID"]

###--------------------------------------------------

data1 <- with(data0b, data.frame("id"=ID1, "Yerror"=Yerror, 
                                 "time"=TIMEyr, "time2"=TIMEyr2, 
                                 "x1"=AGE_BL, "x2"=SEXM, 
                      "v1"=BMI, "v2"=WOMTSmax) )

str(data1)


###--------------------------------------------------

inits1 <- function(){	list(   
  "beta0" = rnorm(1,0,0.1) ,
  "beta" = rnorm(data2$L,0,0.1) ,
  "gama" = matrix(rnorm(data2$n*data2$M,0,0.1),data2$n,data2$M) ,
  "sigma2_prior" = rep(1,data2$K) ,
  "Wtrue"=rnorm(length(data2$Z[,2]),mean=0,sd=0.1)-data2$Z[,2] 
# change line above to below in case parent-node error appears
#"Wtrue"=rnorm(length(data2$Z[,2]),mean=0,sd=0.099)-data2$Z[,2] 
)	}

param1 = c("beta0","beta", "lambda", "alpha",
           "tau2", "sigma2", 
           "Gama0"
)

param2 = c("beta0","beta", "lambda", "alpha", 
           "tau2", "sigma2", 
           "Gama0","gam0",
           "g" 
)


inits3 <- function(){	list(   
  "beta0" = rnorm(1,0,0.1) ,
  "beta" = rnorm(data3$L,0,0.1) ,
  "gama0" = matrix(rnorm(data3$n*data3$M,0,0.01),data3$n,data3$M) ,
  "tau" = 1 ,
"sigma2_prior" = rep(1,data3$K) ,
  "Wtrue"= rnorm(length(data3$Z[,2]),mean=-data3$Z[,2],sd=0.1) ### # change line above to below in case parent-node error appears
#"Wtrue"= rnorm(length(data3$Z[,2]),mean=-data3$Z[,2],sd=0.099) 
  )	}

param3 = c("beta0","beta","lambda",
           "tau2", "sigma2",
           "Gama0", "invGama0", "deviance"
)

###--------------------------------------------------

inits1lcmm <- function(){	list(   
  "beta0" = rnorm(1,0,0.1) ,
  "beta" = rnorm(data2lcmm$L,0,0.1) ,
  "gama" = matrix(rnorm(data2lcmm$n*data2lcmm$M,0,0.1),data2lcmm$n,data2lcmm$M)
)	}

param1lcmm = c("beta0","beta", "lambda", "alpha",
               "tau2", "Gama0",
               "gam0"
)

param2lcmm = c("beta0","beta", "lambda", "alpha", 
               "tau2", "Gama0", 
               "gam0",
               "g"
)

inits3lcmm <- function(){	list(   
  "beta0" = rnorm(1,0,0.1) ,
  "beta" = rnorm(data3lcmm$L,0,0.1) ,
  "gama" = matrix(rnorm(data3lcmm$n*data3lcmm$M,0,0.01),data3lcmm$n,data3lcmm$M) ,
  "tau" = 1
)	}

param3lcmm = c("beta0","beta","lambda", "tau2", 
               "Gama0", "invGama0", "deviance"  
) 

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

### JAGS

### datos de entrada para el modelo JAGS
data2 <- list(Yerror=Yerror, R2w=(max(Yerror[offset[-(n+1)]])-min(Yerror[offset[-(n+1)]])),   ### response variable 
              X=X, L=ncol(X), zeros.beta=rep(0,ncol(X)), prior.betaB=diag(ncol(X))*0.1,   # fixed effects
              Z=Z, M=ncol(Z), zeros.gama=t(rep(0,ncol(Z))), R.u=diag(1,ncol(Z)),   # random effects
              U=U, P=ncol(U), zeros.lambda=rep(0,ncol(U)), prior.lambdaD=diag(ncol(U))*0.1,   #  class-specific fixed effects 
              V=V, Q=ncol(V),   # latent class membership
              n=n,   # number of subjects 
              offset=offset,   # Indicator for the beginning of observations for each subject for long tables format 
              K=K)   ### number of classes 

#fit1 <- jags.model("lclmm_Kclasses_error_NOincrease.bug", 
#                   data=data2, inits=inits1, n.chains=n.chains,n.adapt=n.adapt) 
#update(fit1,n.iter=n.update)

tmp <- parJagsModel(cl, name="fit1", file="lclmm_Kclasses_error_NOincrease.jag", 
             data=data2, inits=inits1, n.chains=n.chains, n.adapt=n.adapt)

tmp <- parUpdate(cl, "fit1", n.iter=n.update, thin=n.thin)

### Revisar si las cadenas convergen, aunque se tenga label switching
#sample1 <- coda.samples(fit1, param1, n.iter=n.iter, thin=n.thin)

#sample1 <- parCodaSamples(cl, model="fit1", variable.names=param1, 
#                            n.iter=n.iter, thin=n.thin)

#plot(sample1)
#summary(sample1)

### Necesitaremos estas cadenas para el label switching
#sample2 <- coda.samples(fit1, param2, n.iter=n.iter, thin=n.thin)
# attributes(sample2)
#plot(sample2)
#summary(sample2)

sample2 <- parCodaSamples(cl, model="fit1", variable.names=param2, 
                          n.iter=n.iter, thin=n.thin)

###--------------------------------------------------

###--------------------------------------------------

### data for JAGS model
data2lcmm <- list(Wtrue=Yerror, R2w=(max(Yerror[offset[-(n+1)]])-min(Yerror[offset[-(n+1)]])),  ### response variable 
                  X=X, L=ncol(X), zeros.beta=rep(0,ncol(X)), prior.betaB=diag(ncol(X))*0.1,   # fixed effects
                  Z=Z, M=ncol(Z), zeros.gama=t(rep(0,ncol(Z))), R.u=diag(1,ncol(Z)),   # random effects
                  U=U, P=ncol(U), zeros.lambda=rep(0,ncol(U)), prior.lambdaD=diag(ncol(U))*0.1,   #  class-specific fixed effects 
                  V=V, Q=ncol(V),   # latent class membership
                  n=n,   # number of subjects 
                  offset=offset,   # Indicator for the beginning of observations for each subject for long tables format 
                  K=K)     ### number of classes 

tmp <- parJagsModel(cl, name="fit1lcmm", file="lclmm_Kclasses_lclmm.jag", 
             data=data2lcmm, inits=inits1lcmm, n.chains=n.chains, n.adapt=n.adapt)

tmp <- parUpdate(cl, "fit1lcmm", n.iter=n.update, thin=n.thin)

sample2lcmm <- parCodaSamples(cl, model="fit1lcmm", variable.names=param2lcmm, 
                          n.iter=n.iter, thin=n.thin)

###--------------------------------------------------
###--------------------------------------------------
### LABEL SWITCHING 
###--------------------------------------------------
### LABEL SWITCHING 
### es necesario ordenar las salidas de R que se obtienen en "sample2" para aplicar el "label.switching"
### algunas cosas se deben modificar, dependiendo de los datos: matrices Z (parametro lambda)
### se usa una verosimilitud aproximada (porque la exacta depende de las variables latentes)
### se usa una probabilidad de pertenencia a la clase aproximada 


### DATA and LISTS for Label Switching 
n <- data2$n   # # sample size 
K <- data2$K   # # number of components, classes
# produce an MCMC sample 
m1 <- n.iter/n.thin #  # define number of iterations for each chain
m2 <- n.chains*n.iter/n.thin #  # define number of iterations for total chains 
# write output to a convenient format for label.switching package: 
# MCMC parameters should be saved as m\times K\times J array
#  # different types of parameters:
mcmc.pars <- array(data = NA, dim = c(m2, K, data2$P))
mcmc.alpha <- array(data = NA, dim = c(m2, K,data2$Q))
mcmc.beta0 <- array(data = NA, dim = c(m2))
mcmc.beta <- array(data = NA, dim = c(m2, data2$L))
mcmc.tau2 <- array(data = NA, dim = c(m2))
mcmc.Gama0 <- array(data = NA, dim = c(m2, data2$M,data2$M))
mcmc.sigma2 <- array(data = NA, dim = c(m2,K))
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
  for(h1 in 1:data2$M){  
    for(h2 in 1:data2$M){
      mcmc.Gama0[(chain-1)*m1+1:m1, h1,h2] <- sample2[[chain]][,paste0("Gama0[",h1,",",h2,"]")] ### Gama0
    }  
  }
  for(i in 1:n){
    zclass[(chain-1)*m1+1:m1,i] <- sample2[[chain]][,paste0("g[",i,"]")]
  }
}

# Compute necessary information that will be additional input to the
# label.switching package.  Define the complete log-likelihood function 
# esta verosimilitud es aproximada

### computing complete loglikelihoods (in parallel) in order to find pivot
nil <- clusterEvalQ(cl0, library(mnormt))
logLvec <- parSapply(cl0, 1:m2, complete.loglikelihood, zclass=zclass, data2 = data2,
                     mcmc.pars = mcmc.pars, mcmc.beta0 = mcmc.beta0, mcmc.beta = mcmc.beta, mcmc.alpha = mcmc.alpha, 
                     mcmc.tau2 = mcmc.tau2, mcmc.sigma2=mcmc.sigma2, mcmc.Gama0 = mcmc.Gama0)

# computing complete loglikelihoods in order to find pivot
zmapindex <- fn.zmapindex(zclass, logLvec, m2)
mapindex <- zmapindex$mapindex
zmap <- zmapindex$zmap 
maxL <- zmapindex$maxL 

# computing allocation probabilities for stephens method and ECR-ITERATIVE-2
# nil <- clusterEvalQ(cl, library(mnormt))
piter <- parLapply(cl0, 1:m2, fn.probabilities, data2=data2, TT = TT, 
                   mcmc.pars = mcmc.pars, mcmc.tau2 = mcmc.tau2, 
                   mcmc.beta0 = mcmc.beta0, mcmc.beta = mcmc.beta, mcmc.alpha = mcmc.alpha, 
                   mcmc.sigma2 = mcmc.sigma2, 
                   mcmc.Gama0 = mcmc.Gama0) 
p <- do.call(abind, list(...=piter, along=1))

# Run label.switching command using all methods.

# The pivot for default ECR algorithm will be the allocation `mapindex` that
# corresponds to the maximum of complete logL estimate: zpivot=z[mapindex,].

# The pivot for PRA algorithm will be the parameters that correspond to the same
# iteration: prapivot = mcmc.pars[mapindex,,].

# The SJW method will be initialized using this iteration as well: sjwinit =
# mapindex.  The complete log-likelihood is defined as: complete =
# complete.normal.loglikelihood

mcmc.parsNOincr <- mcmc.pars

ls <- label.switching(method = set, 
                      zpivot = zmap, 
                      z = zclass, 
                      K = K, 
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
mcmc.pars <- array(data = NA, dim = c(m2, K, data2lcmm$P))
mcmc.alpha <- array(data = NA, dim = c(m2, K,data2lcmm$Q))
mcmc.beta0 <- array(data = NA, dim = c(m2))
mcmc.beta <- array(data = NA, dim = c(m2, data2lcmm$L))
mcmc.tau2 <- array(data = NA, dim = c(m2))
mcmc.Gama0 <- array(data = NA, dim = c(m2, data2lcmm$M,data2lcmm$M))
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
  for(h1 in 1:data2lcmm$M){  
    for(h2 in 1:data2lcmm$M){
      mcmc.Gama0[(chain-1)*m1+1:m1, h1,h2] <- sample2lcmm[[chain]][,paste0("Gama0[",h1,",",h2,"]")] ### Gama0
    }  
  }
  for(i in 1:n){
    zclass[(chain-1)*m1+1:m1,i] <- sample2lcmm[[chain]][,paste0("g[",i,"]")]
  }
}

### computing complete loglikelihoods (in parallel) in order to find pivot
# nil <- clusterEvalQ(cl, library(mnormt))
logLvec.lcmm <- parSapply(cl0, 1:m2, complete.loglikelihood.lcmm, zclass=zclass, data2lcmm = data2lcmm,
                     mcmc.pars = mcmc.pars, mcmc.beta0 = mcmc.beta0, mcmc.beta = mcmc.beta, mcmc.alpha = mcmc.alpha, 
                     mcmc.tau2 = mcmc.tau2, mcmc.Gama0 = mcmc.Gama0)

# Function below finds the pivot
zmapindex.lcmm <- fn.zmapindex(zclass, logLvec.lcmm, m2)  
mapindex <- zmapindex.lcmm$mapindex
zmap <- zmapindex.lcmm$zmap 
maxL <- zmapindex.lcmm$maxL 

# computing allocation probabilities for stephens method and ECR-ITERATIVE-2
piter <- parLapply(cl0, 1:m2, fn.probabilities.lcmm, data2lcmm=data2lcmm, TT = TT, 
                   mcmc.pars = mcmc.pars, mcmc.tau2 = mcmc.tau2, 
                   mcmc.beta0 = mcmc.beta0, mcmc.beta = mcmc.beta, mcmc.alpha = mcmc.alpha, 
                   mcmc.Gama0 = mcmc.Gama0) 
plcmm <- do.call(abind, list(...=piter, along=1))

mcmc.parslcmm <- mcmc.pars

lslcmm <- label.switching(method = set, 
                          zpivot = zmap, 
                          z = zclass, 
                          K = K, 
                          prapivot = mcmc.pars[mapindex, , ], 
                          p = plcmm, 
                          # complete = complete.loglikelihood.lcmm, 
                          mcmc = mcmc.pars, 
                          # data = data2lcmm$Wtrue, 
                          constraint=2, ### 2=pendiente del tiempo
                          sjwinit = mapindex
)

###--------------------------------------------------


###--------------------------------------------------
### USAR JAGS y la CLASE ESTIMADA por label.switching 

### Elegir las clases de alguno de los metodos
### STEPHENS, PRA, ECR, ECR-ITERATIVE-1, ECR-ITERATIVE-2, AIC

data3 <- list(Yerror=Yerror, R2w=(max(Yerror)-min(Yerror)), 
              X=X, L=ncol(X), zeros.beta=rep(0,ncol(X)), prior.betaB=diag(ncol(X))*0.1,   ### fixed effects
              Z=Z, M=ncol(Z), zeros.gama=t(rep(0,ncol(Z))), R.u=diag(1,ncol(Z)),   ### random effects
              U=U, P=ncol(U), zeros.lambda=rep(0,ncol(U)), prior.lambdaD=diag(ncol(U))*0.1,   ### class-specific fixed effects 
              n=n,
              offset=offset,   # Indicator for the beginning of observations for each subject for long tables format 
              K=K,   ### number of classes 
              g=ls$clusters["STEPHENS",] )   ### indicator for classes

tmp <- parJagsModel(cl, name="fit3a", file="lclmm_Kclasses_error_NOincrease_labelswitching.jag",
             data=data3, inits=inits3, n.chains=n.chains)

tmp <- parUpdate(cl, "fit3a", n.iter=n.adapt, thin=n.thin)

sample3aNOincr <- parCodaSamples(cl, model="fit3a", variable.names=param3,
                              n.iter=n.iter, thin=n.thin)


data3$g <- ls$clusters["ECR-ITERATIVE-1",]

tmp <- parJagsModel(cl, name="fit3b", file="lclmm_Kclasses_error_NOincrease_labelswitching.jag",
             data=data3, inits=inits3, n.chains=n.chains)

tmp <- parUpdate(cl, "fit3b", n.iter=n.adapt, thin=n.thin)

sample3bNOincr <- parCodaSamples(cl, model="fit3b", variable.names=param3,
                           n.iter=n.iter, thin=n.thin)


all <- as.matrix(sample3aNOincr, iters = TRUE, chains=TRUE)
pD <- var(all[,"deviance"])/2
DICall <- mean(all[,"deviance"]) + pD
cat("Proposal Sample3a:")
c(DIC=DICall,pD=pD)

all <- as.matrix(sample3bNOincr, iters = TRUE, chains=TRUE)
pD <- var(all[,"deviance"])/2
DIC <- mean(all[,"deviance"]) + pD
cat("Proposal Sample3b:")
c(DIC=DIC,pD=pD)

MCMCsummary(sample3aNOincr, digits=3)
MCMCsummary(sample3bNOincr, digits=3)

summary(sample3aNOincr)
summary(sample3bNOincr)

labswi <- ls
sample3a <- sample3aNOincr
sample3b <- sample3bNOincr
save(labswi, sample3a, sample3b, file = paste0(savepath,'/Sample3_K=',K,"_error_NOincreasing_OAI.RData"))

### If using R2jags
# list2env(data3, envir=.GlobalEnv)
# 
# jagsfit3 <- jags.parallel(
#   data=c(names(data3), "data3", "n.iter", "n.chains", "n.thin"), 
#   inits=inits3, 
#   parameters.to.save=param3, 
#   n.iter=n.iter, 
#   n.chains=n.chains,
#   n.cluster=n.chains,
#   n.thin = n.thin,
#   model.file="lclmm_Kclasses_error_NOincrease_labelswitching.jag")
# 
# 
# ### nice summary
# print(jagsfit3)
# 
# ### Same as coda.samples output
# sample3 <- as.mcmc.list(as.mcmc(jagsfit3))
# 
# ### Extract DIC
# DIC[sim,] = jagsfit3$BUGSoutput$DIC
# 
# 
# PARAM[sim,,1] = unlist(jagsfit3$BUGSoutput$mean)
# PARAM[sim,,2] = unlist(jagsfit3$BUGSoutput$sd)
# PARAM[sim,,3] = unlist(jagsfit3$BUGSoutput$median)

###--------------------------------------------------

data3lcmm <- list(Wtrue=Yerror, R2w=(max(Yerror)-min(Yerror)),   ### response variable
                  X=X, L=ncol(X), zeros.beta=rep(0,ncol(X)), prior.betaB=diag(ncol(X))*0.1,   ### fixed effects
                  Z=Z, M=ncol(Z), zeros.gama=t(rep(0,ncol(Z))), R.u=diag(1,ncol(Z)),   ### random effects
                  U=U, P=ncol(U), zeros.lambda=rep(0,ncol(U)), prior.lambdaD=diag(ncol(U))*0.1,   ### class-specific fixed effects 
                  n=n,
                  offset=offset,   # Indicator for the beginning of observations for each subject for long tables format 
                  K=K,   ### number of classes 
                  g=lslcmm$clusters["STEPHENS",])   ### indicator for classes

tmp <- parJagsModel(cl, name="fit3alcmm", file="lclmm_Kclasses_lclmm_labelswitching.jag",
             data=data3lcmm, inits=inits3lcmm, n.chains=n.chains)

tmp <- parUpdate(cl, "fit3alcmm", n.iter=n.adapt, thin=n.thin)

sample3alcmm <- parCodaSamples(cl, model="fit3alcmm", variable.names=param3lcmm,
                          n.iter=n.iter, thin=n.thin)


data3lcmm$g <- ls$clusters["ECR-ITERATIVE-1",]

tmp <- parJagsModel(cl, name="fit3blcmm", file="lclmm_Kclasses_lclmm_labelswitching.jag",
             data=data3lcmm, inits=inits3lcmm, n.chains=n.chains)

tmp <- parUpdate(cl, "fit3blcmm", n.iter=n.adapt, thin=n.thin)

sample3blcmm <- parCodaSamples(cl, model="fit3blcmm", variable.names=param3lcmm,
                           n.iter=n.iter, thin=n.thin)



all <- as.matrix(sample3alcmm, iters = TRUE, chains=TRUE)
pD <- var(all[,"deviance"])/2
DICall <- mean(all[,"deviance"]) + pD
cat("LCMM Sample3a:")
c(DIC=DICall,pD=pD)

all <- as.matrix(sample3blcmm, iters = TRUE, chains=TRUE)
pD <- var(all[,"deviance"])/2
DIC <- mean(all[,"deviance"]) + pD
cat("LCMM Sample3b:")
c(DIC=DIC,pD=pD)

MCMCsummary(sample3alcmm, digits=3)
MCMCsummary(sample3blcmm, digits=3)

summary(sample3alcmm)
summary(sample3blcmm)

labswi <- lslcmm
sample3a <- sample3alcmm
sample3b <- sample3blcmm
save(labswi, sample3a, sample3b, file = paste0(savepath,'/Sample3_K=',K,"_lclmm_OAI.RData"))

### If using R2jags
# list2env(data3lcmm, envir=.GlobalEnv)
# 
# jagsfit3.lcmm <- jags.parallel(
#   data = c(names(data3lcmm), "data3lcmm", "n.iter", "n.chains", "n.thin"), 
#   inits = inits3lcmm, 
#   parameters.to.save = param3lcmm, 
#   n.iter = n.iter, 
#   n.chains = n.chains,
#   n.cluster = n.chains,
#   n.thin = n.thin,
#   model.file="lclmm_Kclasses_lclmm_labelswitching.jag")
# 
# 
# ### nice summary 
# print(jagsfit3.lcmm) 
# 
# ### Same as coda.samples output
# sample3.lcmm <- as.mcmc.list(as.mcmc(jagsfit3.lcmm))
# 

stopCluster(cl)


