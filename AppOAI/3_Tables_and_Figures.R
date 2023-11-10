###--------------------------------------------------
# Latent Class Linear Mixed Models
# Continuous response
# K classes
# Generate tables and figures
###--------------------------------------------------
### Packages needed
library(reportRmd)
library(xtable)
library(ggplot2)
library(rjags)
library(R2jags)
library(label.switching)  #(version 1.3)       
library(RColorBrewer)
library(patchwork)
# remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(dplyr)
library(loo) # WAIC, and LOO


###--------------------------------------------------
setwd("HERE")

load(file="OAI_data_long.RData")
### loads data.frame OAI_KL_JSW_data obtained from file 1_Compile_KLData.R

savepath <- "SAVEPATH" # where the results from "2_Kclasses_TwoModels.R" are stored

###--------------------------------------------------

###--------------------------------------------------
### PROCESAR DATOS
outcome <- "MCMJSWsum"

colskeep <- c("SEX","V00COHORT","ID","VISIT","MCMJSWsum","MCMJSWmin","TIME","BMI_BL","AGE_BL","BMI","WOMTSmax","WOMTSmax_BL") # 

data0 <- OAI_KL_JSW_data[complete.cases(OAI_KL_JSW_data[,colskeep]),colskeep]

# data0 <- droplevels(subset(data0, V00COHORT=="2: Incidence"))
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

Yerror = data0b[,outcome]
X = cbind(data0b$AGE_BL, data0b$SEXM, data0b$BMI, data0b$WOMTSmax)
Z = cbind(rep(1,nrow(data0b)),data0b$TIMEyr)
U = cbind(rep(1,nrow(data0b)),data0b$TIMEyr,data0b$TIMEyr2)   
V = cbind(data0b$AGE_BL, data0b$SEXM, data0b$BMI)[data0b$VISIT=="V00", ]

n = length(unique(data0b$ID))
TT = 7

ID0 <- data.frame(ID=1:n)
rownames(ID0) <- unique(data0b$ID)
ID1 = ID0[as.character(data0b$ID),"ID"]

###--------------------------------------------------

data1 <- with(data0b, data.frame("id"=ID1, "Yerror"=Yerror, 
                                 "time"=TIMEyr, "time2"=TIMEyr2, 
                                 "x1"=AGE_BL, "x2"=SEXM, 
                                 "x3"=BMI, "x4"=WOMTSmax) )
str(data1)

### Baseline summary (Table 1):
tab <- reportRmd::rm_covsum(droplevels(subset(data0a, TIME==0)), c(outcome,"AGE_BL","SEX", "BMI", "WOMTSmax", "V00COHORT"),tableOnly = TRUE)
print(xtable::xtable(tab), include.rownames=F)



###--------------------------------------------------

### Figure 1
ggplot(data1, aes(x=time, y=Yerror, group=id)) +
  geom_point(size=3, alpha=0.3) + geom_line(alpha=0.3) +
  xlab("Time (years)") + ylab("Total MCMJSW (mm)") +
  scale_x_continuous(breaks=c(0:9)) +
  theme_classic(base_size = 17)

if( !file.exists(paste0(savepath,"/FigMCMJSWsum_trajectories.png")) ){
  ggsave(filename = "FigMCMJSWsum_trajectories.png", width = 9, height = 5,units = "in", dpi = 500)
}


# Figure 3
plotbyclass <- function(labswiobj, labswimethod) {
  data1$ga <- labswiobj$clusters[labswimethod,data1$id]
  data1$Class <- factor(data1$ga)
  proptbl <- prop.table(table(data1$Class))
  labls <- paste0(names(proptbl), " (",sprintf("%.1f",100*proptbl),")")
  data1$Class <- factor(data1$ga, labels = labls)
  p1 <- ggplot(data1, aes(x=time, y=Yerror, group=id, colour=Class)) +
    geom_point(size=1.5, alpha=0.3) + 
    geom_line(alpha=0.1, linetype=2) +
    xlab("Time (years)") + ylab("Total MCMJSW (mm)") +
    scale_x_continuous(breaks=c(0:9)) +
    theme_classic(base_size = 17) + 
    scale_color_manual(values = mycols) + 
    stat_smooth(data=data1, aes(x=time, y=Yerror, 
                                group=Class, fill=Class),
                method="loess", se=FALSE, level=0.95, lwd=1)
  return(p1)
}

### Plots by classes
for( K in 2:5 ){
  # K=2
  mycols <- brewer.pal(5, "Set1")[1:K]
  load(file = paste0(savepath,'/Sample3_K=',K,"_error_NOincreasing_OAI.RData"))
  # loads labswi, sample3a, sample3b, 
  assign(paste0("labswi_NOincr_k",K),labswi)
  assign(paste0("sample3a_NOincr_k",K),sample3a)
  assign(paste0("sample3b_NOincr_k",K),sample3b)
  assign(paste0("jagsfit3a_k",K),jagsfit3a)
  assign(paste0("jagsfit3b_k",K),jagsfit3b)
  
  load(file = paste0(savepath,'/Sample3_K=',K,"_lclmm_OAI.RData"))
  # loads labswi, sample3a, sample3b
  assign(paste0("labswi_lclmm_k",K),labswi)
  assign(paste0("sample3a_lclmm_k",K),sample3a)
  assign(paste0("sample3b_lclmm_k",K),sample3b)
  assign(paste0("jagsfit3alcmm_k",K),jagsfit3alcmm)
  assign(paste0("jagsfit3blcmm_k",K),jagsfit3blcmm)
  
  ### Make plots for each method
  for( methodid in c("STEPHENS","ECR-ITERATIVE-1") )
  {
    # methodid <- "STEPHENS"
    methodlabel <- ifelse(methodid=="STEPHENS","ste","ecr1")
    pA <- plotbyclass(get(paste0("labswi_lclmm_k",K)), methodid)
    pB <- plotbyclass(get(paste0("labswi_NOincr_k",K)), methodid)
    if( K==2 ) {
      pA <- pA + ggtitle("LCLMM")
      pB <- pB + ggtitle("Proposal")
    }
    assign(paste0("plot_k",K,"_",methodlabel), pA+pB)
  }
}

pecr1 <- plot_k2_ecr1 / plot_k3_ecr1 / plot_k4_ecr1 / plot_k5_ecr1

if( !file.exists(paste0(savepath,"/FigMCMJSWsum_classes_ECR1.png")) ){
  ggsave(pecr1, filename = paste0(savepath,"/FigMCMJSWsum_classes_ECR1.png"), width = 12, height = 18, units = "in", dpi = 300)
}

## Figure 4 - Sankey plots
Classes_df <- data.frame()
for( K in 2:5 ){
  # K=2
  Classes_dfK <- data.frame(
    get(paste0("labswi_NOincr_k",K))$clusters["STEPHENS",],
    get(paste0("labswi_NOincr_k",K))$clusters["ECR-ITERATIVE-1",],
    get(paste0("labswi_lclmm_k",K))$clusters["STEPHENS",],
    get(paste0("labswi_lclmm_k",K))$clusters["ECR-ITERATIVE-1",])
  names(Classes_dfK) <- paste(rep(c("Noincr","LCLMM"),each=2),rep(c("STEP","ECR1"),2),K, sep="_")
  if( K==2 ){
    Classes_df <- Classes_dfK
  }else Classes_df <- cbind(Classes_df, Classes_dfK)
}

mycols <- brewer.pal(5, "Set1")

for( mod in c("Noincr","LCLMM") ){
  for( lbsw in c("STEP","ECR1") ){
    df <- Classes_df |>
      make_long(paste0(mod,"_",lbsw,"_2"), 
                paste0(mod,"_",lbsw,"_3"), 
                paste0(mod,"_",lbsw,"_4"),
                paste0(mod,"_",lbsw,"_5"))
    
    psank <- ggplot(df, aes(x = x,                        
                            next_x = next_x,
                            node = node,
                            next_node = next_node,        
                            fill = factor(node),
                            label = node)) +
      geom_sankey(flow.alpha = 0.5, node.color = 1) +
      geom_sankey_label(size = 3.5, color = 1, fill = "white") +
      xlab("") +
      # scale_fill_viridis_d() +
      scale_fill_manual(values = mycols) + 
      scale_x_discrete(labels = paste0("K=",2:5)) +
      theme_sankey(base_size = 17) +
      theme(legend.position = "none")
    
    if( mod=="LCLMM" ) {
      psank <- psank + ggtitle("LCLMM")
    }else{
      psank <- psank + ggtitle("Proposal")
    }
    assign(paste0("plot_",mod,"_",lbsw), psank)
  }
}

psankeyscr1 <- plot_LCLMM_ECR1 + plot_Noincr_ECR1

if( !file.exists(paste0(savepath,"/FigSankey_ECR1.png")) ){
  ggsave(psankeyscr1, filename = paste0(savepath,"/FigSankey_ECR1.png"), width = 12, height = 7, units = "in", dpi = 300)
}

# Table 5 -  WAIC and LOO
ICres <- data.frame()
for( K in 2:5 ){
  # K=2
  BUGSout <- get(paste0("jagsfit3alcmm_k",K))$BUGSoutput
  loglik_ <- BUGSout$sims.list$LogLik
  waic_ <- waic(loglik_)
  loo_ <- loo(loglik_)
  WAIClclmma = waic_$estimates[rownames(waic_$estimates)=="waic",1]
  LOOlclmma = loo_$est[rownames(loo_$est)=="looic",1]
  
  BUGSout <- get(paste0("jagsfit3blcmm_k",K))$BUGSoutput
  loglik_ <- BUGSout$sims.list$LogLik
  waic_ <- waic(loglik_)
  loo_ <- loo(loglik_)
  WAIClclmmb = waic_$estimates[rownames(waic_$estimates)=="waic",1]
  LOOlclmmb = loo_$est[rownames(loo_$est)=="looic",1]
  
  BUGSout <- get(paste0("jagsfit3a_k",K))$BUGSoutput
  loglik_ <- BUGSout$sims.list$LogLik
  waic_ <- waic(loglik_)
  loo_ <- loo(loglik_)
  WAICnoincra = waic_$estimates[rownames(waic_$estimates)=="waic",1]
  LOOnoincra = loo_$est[rownames(loo_$est)=="looic",1]
  
  BUGSout <- get(paste0("jagsfit3b_k",K))$BUGSoutput
  loglik_ <- BUGSout$sims.list$LogLik
  waic_ <- waic(loglik_)
  loo_ <- loo(loglik_)
  WAICnoincrb = waic_$estimates[rownames(waic_$estimates)=="waic",1]
  LOOnoincrb = loo_$est[rownames(loo_$est)=="looic",1]
  
  dtK <- data.frame(K=K, 
                    LS=c("STEPHENS", "ECR-1"), 
                    WAIC.LCLMM = c(WAIClclmma, WAIClclmmb),
                    WAIC.Prop =c(WAICnoincra, WAICnoincrb),
                    LOO.LCLMM = c(LOOlclmma, LOOlclmmb),
                    LOO.Prop =c(LOOnoincra, LOOnoincrb)
  )
  ICres <- rbind(ICres, dtK)
}

print(xtable::xtable(ICres[ICres$LS=="ECR-1",-2], digits = 2), include.rownames=F)


### Table 5 - summaries for K=3 (ECR-1)
K=3
tabA1 <- summary(get(paste0("sample3b_lclmm_k",K)))$statistics
tabA2 <- summary(get(paste0("sample3b_lclmm_k",K)))$quantiles
tabB1 <- summary(get(paste0("sample3b_NOincr_k",K)))$statistics
tabB2 <- summary(get(paste0("sample3b_NOincr_k",K)))$quantiles
tabA <- data.frame(tabA1, CR=paste0("(",apply(formatC(tabA2[,c(1,5)], digits=3, format = "fg"), 1, paste0, collapse=","),")"))
tabB <- data.frame(tabB1,  CR=paste0("(",apply(formatC(tabB2[,c(1,5)], digits=3, format = "fg"), 1, paste0, collapse=","),")"))
# rownames(tabA) %in% rownames(tabB)
# rownames(tabB) %in% rownames(tabA)
Parameter <- rownames(tabB)
## reorder to show betas first
betasidx <- which(grepl("beta", Parameter))
llidx <- which(grepl("LogLik", Parameter))
Parameter <- Parameter[c(betasidx,(1:length(Parameter))[-c(betasidx,llidx)])]
tabsummK3 <- data.frame(Parameter, 
                        tabA[match(Parameter,rownames(tabA)),c(1:2,5)], 
                        tabB[match(Parameter,rownames(tabB)),c(1:2,5)])

print(xtable::xtable(tabsummK3, digits = 3), include.rownames=F)





