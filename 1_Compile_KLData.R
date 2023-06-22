setwd("OAICompleteData_ASCII") # Location of OAICompleteData_ASCII, which can be obtained (with permission) from https://nda.nih.gov/oai/  


#### Obtain X-ray dates ----
for( file in grep("XRay[[:digit:]]+.txt", dir(),value = T) ){
  filei <- read.table(paste0("",file),sep="|",header=TRUE,na.strings = c("",".: Missing Form/Incomplete Workbook"), stringsAsFactors = F)
  names(filei) <- toupper(names(filei))
  XR_col <- grep("V[[:digit:]]{2}XRDATE",names(filei),value=T)
  filt_col <- grep("V[[:digit:]]{2}EXAMTP",names(filei),value=T)
  assign(paste0("XRay",gsub("XRay(\\d+).txt","\\1",file)),filei[filei[,filt_col]=="Bilateral PA Fixed Flexion Knee",c("ID",XR_col)])
}
### ** Merge all visits ----
Xray <- XRay00
Xray$V00TIME <- 0
for( i in c('01','03','05','06','08','10') ){ #i <- '01'
  filei <- get(paste0("XRay",i))
  XR_col <- grep("V[[:digit:]]{2}XRDATE",names(filei),value=T)
  Xray <- merge(Xray, filei, by=c("ID"),all=TRUE) 
  Xray$Time <- as.numeric(as.Date(Xray[,XR_col], format = "%m/%d/%Y") - as.Date(Xray[,"V00XRDATE"], format = "%m/%d/%Y"))/(365.25/12)
  names(Xray)[which(names(Xray)=="Time")] <- paste0("V",i,"TIME")
}
rm(list=grep("XRay[[:digit:]]{2}",ls(),value=T))

###** Put data in long format ----
Dates <- grep("V[[:digit:]]{2}XRDATE",names(Xray),value=T)
Times <- grep("V[[:digit:]]{2}TIME",names(Xray),value=T)
Xray_long <- reshape(Xray, varying = list(Dates,Times), idvar = "ID", timevar = "VISIT", times = paste0("V",c('00','01','03','05','06','08','10')), direction = "long")
Xray_long <- Xray_long[order(Xray_long$ID),]
names(Xray_long)[3:4] <- c("DATE","TIME")
table(Xray_long$VISIT)

### Clinical (Time-varying covariates) ----
for( file in grep("AllClinical(00|01|03|05|06|08|10).txt", dir(),value = T) ){ # file <- grep("AllClinical[[:digit:]]+.txt", dir(),value = T)[2]
  filei <- read.table(paste0("",file),sep="|",header=TRUE,na.strings = c("",".: Missing Form/Incomplete Workbook"),fill = FALSE, quote = "")
  ### Follow-up visit date: V00FVDATE
  ### columns to obtain form the dataset
  names(filei) <- toupper(names(filei))
  nams <- names(filei)
  ### WOMAC Total score, Pain score, Stiffness score, Disability score
  womac_cols <- grep("WOMTS|WOMKP|WOMSTF|WOMADL", nams, value=T)
  bmi_cols <- grep("BMI", nams, value=T)
  date_cols <- grep("DATE", nams, value=T)
  age_cols <- grep("V[[:digit:]]{2}AGE", nams, value=T)
  kneerep_cols <- grep("KRS(L|R)|ART(L|R)", nams, value=T)
  assign(paste0("allclinical",gsub("AllClinical(\\d+).txt","\\1",file)),filei[,c("ID","VERSION",date_cols,bmi_cols,age_cols,womac_cols,kneerep_cols)])
}
### ** Merge all visits ----
### put them all together (same time-points as the KL measures)
allclinical <- allclinical00
for( i in c('01','03','05','06','08','10') ){
  filei <- get(paste0("allclinical",i))
  allclinical <- merge(allclinical,filei[,-2],by=c("ID"),all.x=TRUE) 
}
rm(list=grep("allclinical[[:digit:]]{2}",ls(),value=T))

# Rename some variables to maintain name format with other time points
names(allclinical)[which(names(allclinical)=="P01BMI")] <- "V00BMI"
names(allclinical)[which(names(allclinical)=="V00EVDATE")] <- "V00FVDATE"
allclinical$BMI_BL <- allclinical$V00BMI
allclinical$AGE_BL <- allclinical$V00AGE
allclinical$WOMTSmax_BL <- pmax(allclinical$V00WOMTSL,allclinical$V00WOMTSR)

### Knee replacement or arthroscopy data
names(allclinical)[which(names(allclinical)=="P01KRSR")] <- "V00KRSR12"
names(allclinical)[which(names(allclinical)=="P01KRSL")] <- "V00KRSL12"
allclinical$V00KRSR12 <- ifelse(allclinical$V00KRSR12=="1: Yes", 1, 0)
allclinical$V00KRSL12 <- ifelse(allclinical$V00KRSL12=="1: Yes", 1, 0)
allclinical$V01KRSR12 <- ifelse(allclinical$V01KRSR12=="1: Yes", 1, 0)
allclinical$V01KRSL12 <- ifelse(allclinical$V01KRSL12=="1: Yes", 1, 0)
### No info on KRS for visit s 6 or 8
allclinical$V06KRSR12 <- allclinical$V08KRSR12 <- NA
allclinical$V06KRSL12 <- allclinical$V08KRSL12 <- NA

### Arthroscopy data
names(allclinical)[which(names(allclinical)=="P01ARTR")] <- "V00ARTR12"
names(allclinical)[which(names(allclinical)=="P01ARTL")] <- "V00ARTL12"
allclinical$V00ARTR12 <- ifelse(allclinical$V00ARTR12=="1: Yes", 1, 0)
allclinical$V00ARTL12 <- ifelse(allclinical$V00ARTL12=="1: Yes", 1, 0)
allclinical$V01ARTR12 <- ifelse(allclinical$V01ARTR12=="1: Yes", 1, 0)
allclinical$V01ARTL12 <- ifelse(allclinical$V01ARTL12=="1: Yes", 1, 0)

### ** Put data in long format ----
Dates <- grep("V[[:digit:]]{2}FVDATE",names(allclinical),value=T)
BMIs <- grep("V[[:digit:]]{2}BMI$",names(allclinical),value=T)
Ages <- grep("V[[:digit:]]{2}AGE$",names(allclinical),value=T)
WOMACs_TSL <- grep("V[[:digit:]]{2}WOMTSL$",names(allclinical),value=T)
WOMACs_TSR <- grep("V[[:digit:]]{2}WOMTSR$",names(allclinical),value=T)

WOMACs_FNL <- grep("V[[:digit:]]{2}WOMADLL$",names(allclinical),value=T)
WOMACs_FNR <- grep("V[[:digit:]]{2}WOMADLR$",names(allclinical),value=T)

WOMACs_PNL <- grep("V[[:digit:]]{2}WOMKPL$",names(allclinical),value=T)
WOMACs_PNR <- grep("V[[:digit:]]{2}WOMKPR$",names(allclinical),value=T)

KRS_L <- grep("V[[:digit:]]{2}KRSL12",names(allclinical),value=T)
KRS_L <- KRS_L[order(KRS_L)]
KRS_R <- grep("V[[:digit:]]{2}KRSR12",names(allclinical),value=T)
KRS_R <- KRS_R[order(KRS_R)]
ART_L <- grep("V[[:digit:]]{2}ARTL12",names(allclinical),value=T)
ART_R <- grep("V[[:digit:]]{2}ARTR12",names(allclinical),value=T)

allclinical_long <- reshape(allclinical[,c("ID",grep("_BL",names(allclinical),value = T),Dates,BMIs,Ages,WOMACs_TSL,WOMACs_TSR,WOMACs_FNL,WOMACs_FNR,WOMACs_PNL,WOMACs_PNR,KRS_L,KRS_R,ART_L,ART_R)], varying = list(Dates,BMIs,Ages,WOMACs_TSL,WOMACs_TSR,WOMACs_FNL,WOMACs_FNR,WOMACs_PNL,WOMACs_PNR,KRS_L,KRS_R,ART_L,ART_R), idvar = "ID", timevar = "VISIT", times = paste0("V",c('00','01','03','05','06','08','10')), drop=, direction = "long")
allclinical_long <- allclinical_long[order(allclinical_long$ID),]
names(allclinical_long)[-c(1:5)] <- c("FVDATE","BMI","AGE","WOMTSL","WOMTSR","WOMFNL","WOMFNR","WOMPNL","WOMPNR","KRSL","KRSR","ARTL","ARTR")

KRAT_cols <- c("KRSL","KRSR","ARTL","ARTR")
## Assume missing data on KRK and ART is 0 (No)
# summary(allclinical_long[,KRAT_cols])
allclinical_long[,KRAT_cols][is.na(allclinical_long[,KRAT_cols])] <- 0

### Merge TVcovars and X-ray dates ----
covars_long <- merge(Xray_long, allclinical_long, by=c("ID","VISIT"))
covars_long$WOMTSmax <- pmax(covars_long$WOMTSL,covars_long$WOMTSR)
rm(Xray_long, allclinical_long)
# table(covars_long$VISIT)
# summary(covars_long)

### X-ray KL ----
for( file in grep("kxr_sq_bu[[:digit:]]+.txt", dir(),value = T) ){
  filei <- read.table(paste0("",file),sep="|",header=TRUE,na.strings = c("",".: Missing Form/Incomplete Workbook"))
  names(filei) <- toupper(names(filei))
  KL_col <- grep("V[[:digit:]]{2}XRKL",names(filei),value=T)
  # assign(paste0("kxr_sq_bu",gsub("kxr_sq_bu(\\d+).txt","\\1",file)),filei[,c("ID","SIDE","READPRJ",KL_col)])
  names(filei)[names(filei)==KL_col] <- "XRKL"
  filei$VISIT <- paste0("V",gsub("kxr_sq_bu(\\d+).txt","\\1",file))
  assign(paste0("kxr_sq_bu",gsub("kxr_sq_bu(\\d+).txt","\\1",file)),filei[,c("ID","VISIT","SIDE","READPRJ","XRKL")])
}
rm(filei)

### recode project for the 42 in months 72 and 96 per instructions
kxr_sq_bu08$READPRJ[kxr_sq_bu08$READPRJ==42] <- 37
kxr_sq_bu10$READPRJ[kxr_sq_bu10$READPRJ==42] <- 37

### ** Stack all visits ----
Xry_sq <- subset(kxr_sq_bu00, READPRJ==15)
for( i in c('01','03','05','06','08','10') ){
  filei <- get(paste0("kxr_sq_bu",i))
  projid <- ifelse(i %in% c('08','10'), 37, 15)
  Xry_sq <- rbind(Xry_sq, subset(filei,READPRJ==projid))
}
rm(list=grep("kxr_sq_bu[[:digit:]]{2}",ls(),value=T))
# length(unique(Xry_sq$ID))
# table(Xry_sq$VISIT,Xry_sq$SIDE)

### ** Take maximum score per side ----
Xry_sq_max0_long <- aggregate(XRKL~ID+VISIT,Xry_sq,FUN = max, na.rm=TRUE)
# dim(Xry_sq_max_long)
### Separate L and R sides
Xry_sq_LR_long <- reshape(Xry_sq, timevar = "SIDE", idvar = c("ID","VISIT"), drop = "READPRJ", direction = "wide")
names(Xry_sq_LR_long)[3:4] <- c("XRKLR","XRKLL")
# dim(Xry_sq_max_long)
# 17 records with neither left nor right knee info

Xry_sq_max_long <- merge(Xry_sq_max0_long, Xry_sq_LR_long, by=c("ID","VISIT"), all.x=T)
# Xry_sq_max_long[is.na(Xry_sq_max_long$XRKL),]

### select BL measure
Xry_sq_maxBL <- subset(Xry_sq_max_long, VISIT=="V00")
names(Xry_sq_maxBL)[3:5] <- c("XRKL_BL","XRKLR_BL","XRKLL_BL")
Xry_sq_max_long <- merge(Xry_sq_max_long, Xry_sq_maxBL[,c(1,3:5)], by="ID")
### 38 records lost due to missing baseline measure
## length(unique(Xry_sq_max_long$ID))

Xry_sq_max_long <- Xry_sq_max_long[order(Xry_sq_max_long$ID,Xry_sq_max_long$VISIT),]
### table(Xry_sq_max_long$VISIT)
### table(Xry_sq_max_long$VISIT[Xry_sq_max_long$XRKL_BL%in% 0:1])

### X-ray JSW ----
for( fileval in grep("kxr_qjsw_duryea[[:digit:]]+.txt", dir(),value = T) ){
  # fileval <- grep("kxr_qjsw_duryea[[:digit:]]+.txt", dir(),value = T)[1]
    filei <- read.table(paste0("",fileval),sep="|",header=TRUE,na.strings = c("",".: Missing Form/Incomplete Workbook"))
    names(filei) <- toupper(names(filei))
  KL_col <- grep("V[[:digit:]]{2}MCMJSW",names(filei),value=T)
  names(filei)[names(filei)==KL_col] <- "MCMJSW"
  filei$VISIT <- paste0("V",gsub("kxr_qjsw_duryea(\\d+).txt","\\1",fileval))
  assign(paste0("kxr_qjsw_duryea",gsub("kxr_qjsw_duryea(\\d+).txt","\\1",fileval)),filei[,c("ID","VISIT","SIDE","READPRJ","MCMJSW")])
}
rm(filei)


### ** Stack all visits ----
Xry_qjsw <- kxr_qjsw_duryea00
for( i in c('01','03','05','06','08','10') ){
  filei <- get(paste0("kxr_qjsw_duryea",i))
  Xry_qjsw <- rbind(Xry_qjsw, filei)
}
rm(list=grep("kxr_qjsw_duryea[[:digit:]]{2}",ls(),value=T))
# length(unique(Xry_qjsw$ID))
# table(Xry_qjsw$VISIT, Xry_qjsw$SIDE)

### ** Take summaries score per side ----
Xry_qjsw_min_long <- aggregate(MCMJSW~ID+VISIT, Xry_qjsw, FUN = min, na.rm=TRUE)
names(Xry_qjsw_min_long)[3] <- "MCMJSWmin"
Xry_qjsw_sum_long <- aggregate(MCMJSW~ID+VISIT, Xry_qjsw, FUN = sum, na.rm=TRUE)
names(Xry_qjsw_sum_long)[3] <- "MCMJSWsum"
Xry_qjsw_avg_long <- aggregate(MCMJSW~ID+VISIT, Xry_qjsw, FUN = mean, na.rm=TRUE)
names(Xry_qjsw_avg_long)[3] <- "MCMJSWavg"
Xry_qjsw_long0 <- merge(Xry_qjsw_min_long, merge(Xry_qjsw_sum_long, Xry_qjsw_avg_long))
# dim(Xry_qjsw_max_long)
### Separate L and R sides
Xry_qjsw_LR_long <- reshape(Xry_qjsw, timevar = "SIDE", idvar = c("ID","VISIT"), drop = "READPRJ", direction = "wide")
names(Xry_qjsw_LR_long)[3:4] <- c("MCMJSWR","MCMJSWL")
# dim(Xry_qjsw_LR_long)

Xry_qjsw_long <- merge(Xry_qjsw_long0, Xry_qjsw_LR_long, by=c("ID","VISIT"), all.x=T, all.y=T)

### select BL measure
Xry_qjsw_BL <- subset(Xry_qjsw_long, VISIT=="V00")
names(Xry_qjsw_BL)[-c(1:2)] <- paste0(names(Xry_qjsw_BL)[-c(1:2)],"_BL")
Xry_qjsw_long1 <- merge(Xry_qjsw_long, Xry_qjsw_BL[,c(1,3:5)], by="ID")
## length(unique(Xry_qjsw_long1$ID))

Xry_qjsw_long1 <- Xry_qjsw_long1[order(Xry_qjsw_long1$ID,Xry_qjsw_long1$VISIT),]
### table(Xry_qjsw_long1$VISIT)

### demographics ----
enrollees <- read.table("Enrollees.txt",sep="|",header=TRUE,na.strings = c("",".: Missing Form/Incomplete Workbook"))
# table(enrollees$P02SEX , useNA="ifany")
enrollees$SEX <- factor(ifelse(enrollees$P02SEX=="1: Male", "M", "F"))
# summary(enrollees$SEX )

# table(enrollees$P02RACE, enrollees$P02HISP, useNA="ifany")
enrollees$P02RACENH <- factor(mapply(function(race, hisp){
  if( !is.na(hisp) & hisp == "1: Yes" ) return(3)
  if( !is.na(race) & race == "1: White or Caucasian" ){
    return(1) 
  }else if( !is.na(race) & race == "2: Black or African American" ){
    return(2)
  }else return(4)
  }, race=enrollees$P02RACE, hisp=enrollees$P02HISP))
# table(enrollees$P02RACENH , useNA="ifany")
enrollees$RACE <- factor(ifelse(enrollees$P02RACENH==1, "White", "Non-white"), levels=c("White", "Non-white"))
# table(enrollees$RACE , useNA="ifany")
# table(enrollees$V00COHORT , useNA="ifany")


### Put all data together ----
### merge with dates, time, and time-varying covariates
covars_long1 <- merge(covars_long, Xry_sq_max_long, by=c('ID',"VISIT"), all.x=T)
covars_long1 <- merge(covars_long1, Xry_qjsw_long1, by=c('ID',"VISIT"), all.x=T)
OAI_KL_JSW_data <- merge(enrollees[, c("ID", "SEX", "RACE", "V00COHORT")], covars_long1, by="ID")
# summary(OAI_KL_JSW_data)
# length(unique(OAI_KL_JSW_data$ID))
# table(OAI_KL_JSW_data$VISIT,useNA="ifany")
OAI_KL_JSW_data$BMI_BL_cat <- droplevels(factor(sapply(OAI_KL_JSW_data$BMI_BL,function(x){
  if(is.na(x)) return(NA)
  if(x>=30) return("Obese")
  if(x>=25 & x<30) return("Overweight")
  if(x>=18.5 & x<25) return("Healthy weight")
  if(x<18.5) return("Underweight")
}),levels=c("Healthy weight","Underweight","Overweight","Obese")))
### table(OAI_KL_JSW_data$BMI_BL_cat,useNA = "ifany")
OAI_KL_JSW_data$AGE_BL_cat <- factor(ifelse(OAI_KL_JSW_data$AGE_BL<=60, "<=60", ">60"))
### table(OAI_KL_JSW_data$AGE_BL_cat,useNA = "ifany")

### calculate an approximate time for one participant with missing XRDATE (use FVDATE instead of XRDATE as proxy)
# Xray[Xray$ID==9496386, ]
# OAI_KL_JSW_data[OAI_KL_JSW_data$ID==9496386, ]
OAI_KL_JSW_data[is.na(OAI_KL_JSW_data$TIME),"DATE"] <- "08/14/2006"
OAI_KL_JSW_data[is.na(OAI_KL_JSW_data$TIME),"TIME"] <- as.numeric(as.Date("08/14/2006", format = "%m/%d/%Y") - as.Date("06/14/2005", format = "%m/%d/%Y"))/(365.25/12)
# length(unique(OAI_KL_JSW_data$ID))


### save data.frame above
save(OAI_KL_JSW_data, file="OAI_data_long.RData")




