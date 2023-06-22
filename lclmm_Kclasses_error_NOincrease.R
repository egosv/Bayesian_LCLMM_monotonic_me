###--------------------------------------------------
# Latent Class Linear Mixed Models
# Growth Mixture Models
# Continuous response
# K classes
###--------------------------------------------------

### Change this address to read the JAGS code 
setwd("HERE")

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

library(foreach)
library(doParallel)
### Creating clusters to run a chain in each core
detectCores()
n.cores <- detectCores()-2
#options(mc.cores = parallel::detectCores()-2)
options(mc.cores = n.cores)

n.chains = n.cores# detectCores()-1 # 10
cl <- makeCluster(n.chains, type="SOCK")
tmp <- clusterEvalQ(cl, library(dclone))


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
###--------------------------------------------------
### Functions for the Results

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
      medtrunc <- sqrt(tau2)*(-dnorm(beta2))/z2 ### mean of truncated normal 
      medaux <- ifelse(is.finite(medtrunc),medtrunc,0)  
      media <- eta2[i]+medaux 
      eta[i] <- media
      vartrunc <- (1+ (-beta2*dnorm(beta2)/z2) -((-dnorm(beta2)/z2)^2)) ### variance of the truncated normal 
      varaux <- tau2*ifelse(is.finite(vartrunc) & vartrunc>0,vartrunc,1)
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
    tau2 <- mcmc.tau2[iter]
    lambda <- mcmc.pars[iter,,]
    beta0 <- mcmc.beta0[iter]
    beta <- mcmc.beta[iter,]
    alpha <- mcmc.alpha[iter,,]
    Gama0 <- mcmc.Gama0[iter,,]
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

###--------------------------------------------------
### SIMULATE DATA
# gama <- function(media,desvest){
#   b = media/(desvest^2)
#   a = b*media
#   return(c(a,b))
# }

set.seed(12345)

TT <- 8   # number of time for repetitions 
n <- 200   # number of subjects 
K <- 3   # number of latent classes  

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
alpha <- matrix(0,K,2)   # KxQ vector of coefficients for the latent class membership 
alpha[,1] <- c(0,-0.2,0.1)   # coefficients for covariate V1  
alpha[,2] <- c(0,0.2,-0.4)   # coefficients for covariate V2  
#alpha[,3] <- c(0.5,0)   # coefficients for covariate V3  
# BMI: gama(29.7,4.7)
V1 = rep(rgamma(n,shape=39.931,rate=1.344),each=TT)   # covariate for the latent class membership 
# WOMAC: gama(23.5,16.9)
V2 = rep(rgamma(n,shape=1.933,rate=0.082),each=TT)   # covariate for the latent class membership 
#V3 = rnorm(n,mean=2,sd=1)   # covariate for the latent class membership 
V = cbind(V1,V2)   # latent class membership design matrix 

### Class-specific fixed effects
# lambda <- matrix(0,K,3)   # KxP vector of coefficients of class-specific fixed effects  
# U = cbind(rep(1,n*TT),time,time2)   ### class-specific fixed effects design matrix 
lambda <- matrix(0,K,3)   # KxP vector of coefficients of class-specific fixed effects  
lambda[,1] <- c(    0, 4.8,  5.2)  # covariate for the class-specific fixed effects U1
lambda[,2] <- c(-0.05,-0.5, -1.2)   # covariate for the class-specific fixed effects U2  
lambda[,3] <- c(-0.01,-0.01,-0.01)   # covariate for the class-specific fixed effects U2  
U = cbind(rep(1,n*TT),time,time2)   ### class-specific fixed effects design matrix 

### Variance 
tau2 <- 0.70   # variance of the latent classes 

### Measurement errors
sigma2 <- c(1.20, 0.1, 0.3)   # variance of the measurement errors

###--------------------------------------------------
### JAGS Specifications 

n.chains = K
n.adapt = 10000 
n.update = 10000 
n.iter = 10000 
n.thin = 5 

###--------------------------------------------------

inits1 <- function(){	list(   
  "beta0" = rnorm(1,0,0.1) ,
  "beta" = rnorm(data2$L,0,0.1) ,
  "gama" = matrix(rnorm(data2$n*data2$M,0,0.1),data2$n,data2$M) ,
#  "lambda" = matrix(rnorm(data2$K*data2$P,0,0.1),data2$K,data2$P) #,
  #  "tau" = rgamma(1,1,1) ,
  "sigma2_prior" = rep(1,data2$K) , #rgamma(data2$K,0.1,0.1) , 
  #  "invGama0" = diag(ncol(data2$Z)) 
  "Wtrue"=rnorm(length(data2$Z[,2]),mean=0,sd=0.1)-data2$Z[,2] 
  ### (creciente) usar el tiempo para proponer una prior
)	}

param1 = c("beta0","beta", "lambda", "alpha",
           "tau2", "sigma2", 
           "Gama0",
           "gam0"
)

param2 = c("beta0","beta", "lambda", "alpha", 
           "tau2", "sigma2", 
           "Gama0",
           "gam0",
           "g" 
)

inits3 <- function(){	list(   
  "beta0" = rnorm(1,0,0.1) ,
  "beta" = rnorm(data3$L,0,0.1) ,
  "gama0" = matrix(rnorm(data3$n*data3$M,0,0.01),data3$n,data3$M) ,
#  "lambda" = matrix(rnorm(data3$K*data3$P,0,0.01),data3$K,data3$P) ,
  "tau" = 1 ,
"sigma2_prior" = rep(1,data3$K) , #rgamma(data3$K,1,1), 
  "Wtrue"= rnorm(length(data3$Z[,2]),mean=-data3$Z[,2],sd=0.1) ### (decreciente) usar el tiempo para proponer una prior 
  )	}

param3 = c("beta0","beta","lambda",
           "tau2", "sigma2",
           "Gama0", "invGama0" 
)

###--------------------------------------------------

inits1lcmm <- function(){	list(   
  "beta0" = rnorm(1,0,0.1) ,
  "beta" = rnorm(data2lcmm$L,0,0.1) ,
  "gama" = matrix(rnorm(data2lcmm$n*data2lcmm$M,0,0.1),data2lcmm$n,data2lcmm$M) #,
  #  ,"lambda" = t(matrix(c(quantile(Wtrue[offset[-(n+1)]],c(.25,.5,.75)),quantile(Wtrue[offset[-(n+1)]+1]-Wtrue[offset[-(n+1)]],c(.25,.5,.75)),rep(0,data2lcmm$P)),data2lcmm$K,data2lcmm$P))
  #,"lambda" = matrix(rnorm(data2lcmm$K*data2lcmm$P,0,0.1),data2lcmm$K,data2lcmm$P) ###,
  ###  "invGama0" = diag(ncol(data2lcmm$Z)) 
)	}

param1lcmm = c("beta0","beta", "lambda", "alpha",
               "tau2", 
               "Gama0",
               "gam0"
)

param2lcmm = c("beta0","beta", "lambda", "alpha", 
               "tau2", 
               "Gama0", 
               "gam0",
               "g"
)

inits3lcmm <- function(){	list(   
  "beta0" = rnorm(1,0,0.1) ,
  "beta" = rnorm(data3lcmm$L,0,0.1) ,
  "gama" = matrix(rnorm(data3lcmm$n*data3lcmm$M,0,0.01),data3lcmm$n,data3lcmm$M) ,
#  "lambda" = matrix(rnorm(data3lcmm$K*data3lcmm$P,0,0.01),data3lcmm$K,data3lcmm$P) ,
  "tau" = 1 #,
  #"rho.u"=0.5, "tau.u1"=1, "tau.u2"=1
)	}

param3lcmm = c("beta0","beta","lambda", "tau2", 
               "Gama0", "invGama0"  
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

###--------------------------------------------------
###--------------------------------------------------

### Classify subjects in latent classes 
ppi <- matrix(0,nrow=n,ncol=K)
w <- matrix(0,nrow=n,ncol=K)
g <- rep(0,n)   # classes 
for(j in 1:n){
  for(k in 1:K){
    ppi[j,k] <- exp(V[j,]%*%alpha[k,])
  }
  for(k in 1:K){
    w[j,k] <- ppi[j,k]/sum(ppi[j,])   # probability for classes
  }
  g[j] <- sample((1:K),1,prob=w[j,])   # class
}
table(g)
clase <- g


### Linear predictor and response 
eta <- matrix(0,n*TT)   # linear predictor
Wtrue <- rep(NA,n*TT)   # response variable 
Yerror <- rep(NA,n*TT)
for(j in 1:n){		
  for(i in offset[j]){	
    eta[i] <- beta0 + Z[i,]%*%gama[j,] + U[i,]%*%lambda[g[j],] + X[i,]%*%beta 
    Wtrue[i] <- rnorm(1, eta[i], sqrt(tau2))
    Yerror[i] <- rnorm(1, Wtrue[i], sqrt(sigma2[g[j]]))
  }
  for(i in (offset[j]+1):(offset[j+1]-1)){
    eta[i] <- beta0 + Z[i,]%*%gama[j,] + U[i,]%*%lambda[g[j],] + X[i,]%*%beta 
    Wtrue[i] <- rnormright(Wtrue[i-1],eta[i],sqrt(tau2)) 
    Yerror[i] <- rnorm(1, Wtrue[i], sqrt(sigma2[g[j]]))
  }
}

clase_cat = factor(clase)
levels(clase_cat) = as.character(c(1:K))
table(clase)
# table(clase)/n

# Monotone patterns decreasing
c(sum(Yerror[-(offset[-(n+1)]+TT-1)][rep(clase,each=5)==1] >= Yerror[-(offset[-(n+1)])][rep(clase,each=5)==1]), 
                     sum(Yerror[-(offset[-(n+1)]+TT-1)][rep(clase,each=5)==2] >= Yerror[-(offset[-(n+1)])][rep(clase,each=5)==2]), 
                     sum(Yerror[-(offset[-(n+1)]+TT-1)][rep(clase,each=5)==3] >= Yerror[-(offset[-(n+1)])][rep(clase,each=5)==3])) 
# Monotone patterns increasing
c(sum(Yerror[-(offset[-(n+1)]+TT-1)][rep(clase,each=5)==1] < Yerror[-(offset[-(n+1)])][rep(clase,each=5)==1]), 
                     sum(Yerror[-(offset[-(n+1)]+TT-1)][rep(clase,each=5)==2] < Yerror[-(offset[-(n+1)])][rep(clase,each=5)==2]), 
                     sum(Yerror[-(offset[-(n+1)]+TT-1)][rep(clase,each=5)==3] < Yerror[-(offset[-(n+1)])][rep(clase,each=5)==3])) 


###--------------------------------------------------

options(mc.cores = n.cores)
cl <- makeCluster(n.chains, type="SOCK")
tmp <- clusterEvalQ(cl, library(dclone))

###--------------------------------------------------
### JAGS

### data for the JAGS model 
data2 <- list(Yerror=Yerror, R2w=(max(Yerror[offset[-(n+1)]])-min(Yerror[offset[-(n+1)]])),   ### response variable 
              X=X, L=ncol(X), zeros.beta=rep(0,ncol(X)), prior.betaB=diag(ncol(X))*0.01,   # fixed effects
              Z=Z, M=ncol(Z), zeros.gama=t(rep(0,ncol(Z))), R.u=diag(1,ncol(Z)),   # random effects
              U=U, P=ncol(U), zeros.lambda=rep(0,ncol(U)), prior.lambdaD=diag(ncol(U))*0.01,   #  class-specific fixed effects 
              V=V, Q=ncol(V),   # latent class membership
              n=n,   # number of subjects 
              offset=offset,   # Indicator for the beginning of observations for each subject for long tables format 
              K=K)   ### number of classes 

#fit1 <- jags.model("lclmm_Kclasses_error_NOincrease.bug", 
#                   data=data2, inits=inits1, n.chains=n.chains,n.adapt=n.adapt) 
#update(fit1,n.iter=n.update)

parJagsModel(cl, name="fit1", file="lclmm_Kclasses_error_NOincrease.jag", 
             data=data2, inits=inits1, n.chains=n.chains, n.adapt=n.adapt)

parUpdate(cl, "fit1", n.iter=n.update, thin=n.thin)

### Review convergence 
#sample1 <- coda.samples(fit1, param1, n.iter=n.iter, thin=n.thin)

#sample1 <- parCodaSamples(cl, model="fit1", variable.names=param1, 
#                            n.iter=n.iter, thin=n.thin)

#plot(sample1)
#summary(sample1)

### We need chains for the label switching
#sample2 <- coda.samples(fit1, param2, n.iter=n.iter, thin=n.thin)
# attributes(sample2)
#plot(sample2)
#summary(sample2)

sample2 <- parCodaSamples(cl, model="fit1", variable.names=param2, 
                          n.iter=n.iter, thin=n.thin)

###--------------------------------------------------

###--------------------------------------------------

### data as input for the JAGS model 
data2lcmm <- list(Wtrue=Yerror, R2w=(max(Yerror[offset[-(n+1)]])-min(Yerror[offset[-(n+1)]])),  ### response variable 
                  X=X, L=ncol(X), zeros.beta=rep(0,ncol(X)), prior.betaB=diag(ncol(X))*0.01,   # fixed effects
                  Z=Z, M=ncol(Z), zeros.gama=t(rep(0,ncol(Z))), R.u=diag(1,ncol(Z)),   # random effects
                  U=U, P=ncol(U), zeros.lambda=rep(0,ncol(U)), prior.lambdaD=diag(ncol(U))*0.01,   #  class-specific fixed effects 
                  V=V, Q=ncol(V),   # latent class membership
                  n=n,   # number of subjects 
                  offset=offset,   # Indicator for the beginning of observations for each subject for long tables format 
                  K=K)     ### number of classes 



#fit1lcmm <- jags.model("lclmm_Kclasses_lclmm.bug", 
#                       data=data2lcmm, inits=inits1lcmm, n.chains=n.chains,n.adapt=n.adapt) 

#update(fit1lcmm,n.iter=n.update)


### Need chans for the label switching
#sample2lcmm <- coda.samples(fit1lcmm, param2lcmm, n.iter=n.iter, thin=n.thin)

parJagsModel(cl, name="fit1lcmm", file="lclmm_Kclasses_lclmm.jag", 
             data=data2lcmm, inits=inits1lcmm, n.chains=n.chains, n.adapt=n.adapt)

parUpdate(cl, "fit1lcmm", n.iter=n.update, thin=n.thin)

### Review convergence

#sample1lcmm <- parCodaSamples(cl, model="fit1lcmm", variable.names=param1lcmm, 
#                          n.iter=n.iter, thin=n.thin)

#plot(sample1lcmm)

sample2lcmm <- parCodaSamples(cl, model="fit1lcmm", variable.names=param2lcmm, 
                          n.iter=n.iter, thin=n.thin)

###--------------------------------------------------
###--------------------------------------------------
### LABEL SWITCHING 
###--------------------------------------------------
### LABEL SWITCHING 
### It is needed to order the output of codes in "sample2" to apply the "label.switching"
### according to the data, what changes is: matrices Z (parameter lambda)
### likelihood approximated (because it depends on latent variables)
### it is used a probability of pertenence to the latent classes 


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
# likelihood approximated

### computing complete loglikelihoods (in parallel) in order to find pivot
nil <- clusterEvalQ(cl, library(mnormt))
logLvec <- parSapply(cl, 1:m2, complete.loglikelihood, zclass=zclass, data2 = data2,
                     mcmc.pars = mcmc.pars, mcmc.beta0 = mcmc.beta0, mcmc.beta = mcmc.beta, mcmc.alpha = mcmc.alpha, 
                     mcmc.tau2 = mcmc.tau2, mcmc.sigma2=mcmc.sigma2, mcmc.Gama0 = mcmc.Gama0)

# computing complete loglikelihoods in order to find pivot
zmapindex <- fn.zmapindex(zclass, logLvec, m2)
mapindex <- zmapindex$mapindex
zmap <- zmapindex$zmap 
maxL <- zmapindex$maxL 

# computing allocation probabilities for stephens method and ECR-ITERATIVE-2
nil <- clusterEvalQ(cl, library(mnormt))
piter <- parLapply(cl, 1:m2, fn.probabilities, data2=data2, TT = TT, 
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
                      
#ls$similarity  #similarity of single best clusterings:      
#ls$timings  #timings for each method
#table(ls$permutations)

# # retrieve the number of observations assigned to each cluster:
# frequency <- apply(ls$clusters, 1, function(y) {
#   freq <- numeric(K)
#   for (j in 1:K) {
#     freq[j] <- length(which(y == j))
#   }
#   return(freq)
# })
# rownames(frequency) <- 1:K
# frequency

# t <- 5
# thinning <- t + seq(1, m2, by = t) - 1  # #select every t-th iteration for plotting the mcmc 
# colors <- brewer.pal(K, name = "Set1")  # #define color pallete     
# pal <- colorRampPalette(colors)


# plot the simulated means (label switching phenomenon)
# image.width <- 6
# image.height <- 6
# pointsize <- 20
# mymargin <- c(4, 4, 2, 0.5)
# par(mar = mymargin)
# par(mfrow=c(4,2))
# matplot(mcmc.pars[thinning, , 1], type = "l", lty = 1, col = pal(K), ylab = "means", 
#         xlab = "iteration", main = "Raw MCMC output")#dev.off()
# 
# # Reordered outputs
# matplot(permute.mcmc(mcmc.pars, ls$permutations$ECR)$output[thinning, , 1], type = "l", 
#         xlab = "iteration", main = "ECR", ylab = "means", lty = 1, col = pal(K))
# matplot(permute.mcmc(mcmc.pars, ls$permutations$"ECR-ITERATIVE-2")$output[thinning, , 1], type = "l", 
#         xlab = "iteration", main = "ECR-iterative-2", ylab = "means", lty = 1, col = pal(K))
# matplot(permute.mcmc(mcmc.pars, ls$permutations$"ECR-ITERATIVE-1")$output[thinning, , 1], type = "l", 
#         xlab = "iteration", main = "ECR-iterative-1", ylab = "means", lty = 1, col = pal(K))
# matplot(permute.mcmc(mcmc.pars, ls$permutations$PRA)$output[thinning, , 1], type = "l", 
#         xlab = "iteration", main = "PRA", ylab = "means", lty = 1, col = pal(K))
# matplot(permute.mcmc(mcmc.pars, ls$permutations$STEPHENS)$output[thinning, , 1], 
#         type = "l", xlab = "iteration", main = "STEPHENS", ylab = "means", lty = 1, col = pal(K))
# #matplot(permute.mcmc(mcmc.pars, ls$permutations$SJW)$output[thinning, , 1], type = "l", 
# #        xlab = "iteration", main = "SJW", ylab = "means", lty = 1, col = pal(K))
# matplot(permute.mcmc(mcmc.pars, ls$permutations$AIC)$output[thinning, , 1], type = "l", 
#         xlab = "iteration", main = "AIC", ylab = "means", lty = 1, col = pal(K))
# #matplot(permute.mcmc(mcmc.pars, ls$permutations$"DATA-BASED")$output[thinning, , 1], type = "l", 
# #        xlab = "iteration", main = "DATA-BASED", ylab = "means", lty = 1, col = pal(K))
# 
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
zclass <- array(0,dim=c(m2,n))   ### identify the latent class 
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

### computing complete loglikelihoods in order to find pivot

### computing complete loglikelihoods (in parallel) in order to find pivot
nil <- clusterEvalQ(cl, library(mnormt))
logLvec.lcmm <- parSapply(cl, 1:m2, complete.loglikelihood.lcmm, zclass=zclass, data2lcmm = data2lcmm,
                     mcmc.pars = mcmc.pars, mcmc.beta0 = mcmc.beta0, mcmc.beta = mcmc.beta, mcmc.alpha = mcmc.alpha, 
                     mcmc.tau2 = mcmc.tau2, mcmc.Gama0 = mcmc.Gama0)

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
nil <- clusterEvalQ(cl, library(mnormt))
piter <- parLapply(cl, 1:m2, fn.probabilities.lcmm, data2lcmm=data2lcmm, TT = TT, 
                   mcmc.pars = mcmc.pars, mcmc.tau2 = mcmc.tau2, 
                   mcmc.beta0 = mcmc.beta0, mcmc.beta = mcmc.beta, mcmc.alpha = mcmc.alpha, 
                   mcmc.Gama0 = mcmc.Gama0) 
plcmm <- do.call(abind, list(...=piter, along=1))


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

stopCluster(cl)

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

table(clase_cat,claseEC1_cat)


data3 <- list(Yerror=Yerror, R2w=(max(Yerror)-min(Yerror)), 
              X=X, L=ncol(X), zeros.beta=rep(0,ncol(X)), prior.betaB=diag(ncol(X))*0.01,   ### fixed effects
              Z=Z, M=ncol(Z), zeros.gama=t(rep(0,ncol(Z))), R.u=diag(1,ncol(Z)),   ### random effects
              U=U, P=ncol(U), zeros.lambda=rep(0,ncol(U)), prior.lambdaD=diag(ncol(U))*0.01,   ### class-specific fixed effects 
              n=n,
              offset=offset,   # Indicator for the beginning of observations for each subject for long tables format 
              K=K,   ### number of classes 
              g=claseEC1 )   ### indicator for classes

list2env(data3, envir=.GlobalEnv)

#fit3 <- jags.model("lclmm_Kclasses_error_NOincrease_labelswitching.bug", 
#                   data3,inits=inits3, n.chains=n.chains, n.adapt=n.adapt)
#update(fit3,n.iter=n.update)
#sample3 <- coda.samples(fit3, param3, n.iter=n.iter*3, thin=n.thin*3)

# parJagsModel(cl, name="fit3", file="lclmm_Kclasses_error_NOincrease_labelswitching.jag", 
#              data=data3, inits=inits3, n.chains=n.chains)
# parUpdate(cl, "fit3", n.iter=n.adapt, thin=n.thin)
# 
# sample3 <- parCodaSamples(cl, model="fit3", variable.names=param3, 
#                               n.iter=n.iter, thin=n.thin)
#   

jagsfit3 <- jags.parallel(
  data=c(names(data3), "data3", "n.iter", "n.chains", "n.thin"), 
  inits=inits3, 
  parameters.to.save=param3, 
  n.iter=n.iter, 
  n.chains=n.chains,
  n.cluster=n.chains,
  n.thin = n.thin,
  model.file="lclmm_Kclasses_error_NOincrease_labelswitching.jag")


### nice summary
print(jagsfit3)

### Same as coda.samples output
sample3 <- as.mcmc.list(as.mcmc(jagsfit3))

### Extract DIC
jagsfit3$BUGSoutput$DIC


#summary(sample3)
#plot(sample3)

jagsfit3$BUGSoutput

###--------------------------------------------------


claselcmmSTE = etiquetar(clase,lslcmm$clusters["STEPHENS",])
claselcmmPRA = etiquetar(clase,lslcmm$clusters["PRA",])
claselcmmECR = etiquetar(clase,lslcmm$clusters["ECR",])
claselcmmEC1 = etiquetar(clase,lslcmm$clusters["ECR-ITERATIVE-1",])
claselcmmEC2 = etiquetar(clase,lslcmm$clusters["ECR-ITERATIVE-2",])
claselcmmAIC = etiquetar(clase,lslcmm$clusters["AIC",])

claselcmmEC1_cat = factor(claselcmmEC1)
levels(claselcmmEC1_cat) = as.character(c(1:K))



data3lcmm <- list(Wtrue=Yerror, R2w=(max(Yerror)-min(Yerror)),   ### response variable
                  X=X, L=ncol(X), zeros.beta=rep(0,ncol(X)), prior.betaB=diag(ncol(X))*0.01,   ### fixed effects
                  Z=Z, M=ncol(Z), zeros.gama=t(rep(0,ncol(Z))), R.u=diag(1,ncol(Z)),   ### random effects
                  U=U, P=ncol(U), zeros.lambda=rep(0,ncol(U)), prior.lambdaD=diag(ncol(U))*0.01,   ### class-specific fixed effects 
                  n=n,
                  offset=offset,   # Indicator for the beginning of observations for each subject for long tables format 
                  K=K,   ### number of classes 
                  g=claselcmmEC1)   ### indicator for classes

# parJagsModel(cl, name="fit3lcmm", file="lclmm_Kclasses_lclmm_labelswitching.jag", 
#              data=data3lcmm, inits=inits3lcmm, n.chains=n.chains)
# parUpdate(cl, "fit3lcmm", n.iter=n.adapt, thin=n.thin)
# 
# sample3lcmm <- parCodaSamples(cl, model="fit3lcmm", variable.names=param3lcmm, 
#                           n.iter=n.iter, thin=n.thin)

list2env(data3lcmm, envir=.GlobalEnv)

jagsfit3.lcmm <- jags.parallel(
  data = c(names(data3lcmm), "data3lcmm", "n.iter", "n.chains", "n.thin"), 
  inits = inits3lcmm, 
  parameters.to.save = param3lcmm, 
  n.iter = n.iter, 
  n.chains = n.chains,
  n.cluster = n.chains,
  n.thin = n.thin,
  model.file="lclmm_Kclasses_lclmm_labelswitching.jag")


### nice summary 
print(jagsfit3.lcmm) 


### Same as coda.samples output
sample3.lcmm <- as.mcmc.list(as.mcmc(jagsfit3.lcmm))

#summary(sample3.lcmm)
#plot(sample3.lcmm)


###--------------------------------------------------

###-------------------------------------------------- 

###--------------------------------------------------

###--------------------------------------------------
data1 <- data.frame("id"=id, "Wtrue"=Wtrue, "Yerror"=Yerror,
                    "time"=time, "time2"=time2, 
                    "timejit"=time+rep(jitter(rep(0,n),factor=10),each=TT), 
                    "x1"=X1, "x2"=X2, #"x3"=X3,
                    "v1"=V1, "v2"=V2, 
                    "g"=g)
str(data1)

p1 <- ggplot(data1, aes(x=timejit, y=Wtrue, group=id, colour=factor(g[id]))) +
  theme_bw()+
  geom_point(cex=0.8) +
  geom_line(lwd=0.1) +
  stat_smooth(data=data1, aes(x=timejit, y=Wtrue,group=g[id]),
              method="loess",se=TRUE,level=0.95,lwd=1,alpha=0.5)#,color="black")
p1
p2 <- ggplot(data1, aes(x=timejit, y=Yerror, group=id, colour=factor(g[id]))) +
  theme_bw()+
  geom_point(cex=0.8) +
  geom_line(lwd=0.1) +
  stat_smooth(data=data1, aes(x=timejit, y=Yerror, group=g[id]),
              method="loess",se=TRUE,level=0.95,lwd=1,alpha=0.5)#,color="black")
p2
###--------------------------------------------------

###--------------------------------------------------

###--------------------------------------------------

g.diff = rep(NA,length(clase))
for(i1 in 1:length(table(clase))){
  for(i2 in 1:length(table(ls$clusters["ECR",]))){
    g.diff[clase==i1 & ls$clusters["ECR",]==i2] = paste("obs",i1,"-","est",i2)
  }
}

data4 <- data.frame("id"=id, "Wtrue"=Wtrue, "Yerror"=Yerror,
                    "time"=time,
                    "id"=id,"g.diff"=g.diff)

p3 <- ggplot(data4, aes(x=time, y=Wtrue, group=id, colour=factor(g.diff[id]))) +
  theme_bw()+
  geom_point() +
  geom_line() +
  stat_smooth(data=data4, aes(x=time, y=Wtrue,group=g[id]),
              method="loess",se=TRUE,level=0.95,lwd=1,alpha=0.5,color="black")
p3
p4 <- ggplot(data4, aes(x=time, y=Yerror, group=id, colour=factor(g.diff[id]))) +
  theme_bw()+
  geom_point() +
  geom_line() +
  stat_smooth(data=data4, aes(x=time, y=Yerror, group=g[id]),
              method="loess",se=TRUE,level=0.95,lwd=1,alpha=0.5,color="black")
p4
###-------------------------------------------------- 
###-------------------------------------------------- 

###-------------------------------------------------- 


