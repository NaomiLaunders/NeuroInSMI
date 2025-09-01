#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Minimally adjusted results schizophrenia, bipolar and other psychoses vs control adjusting for AGE, SEX, ETHNICITY, REGION, CALENDAR YEAR
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all.names = TRUE))

library(dplyr)
library(tableone)
library(tidyselect)
library(tidyr)
library(ggplot2)
library(lubridate)
library(binom)
library(tidyverse)
library(data.table)
library(sandwich)
library(lmtest)
library(gridExtra)
library(RColorBrewer) 
library(cowplot)
library(forestplot)
library(miceadds)

#Bring caseID back in#
load("/Data/ReRunCleanBaseFile.Rdata")
CaseID<-select(MatchedLK, patid, case_patid)
rm(MatchedLK)
load("/Neuro/Data/ImpAnalysis.Rdata")


load("/Neuro/Data/FinalNeuroData.Rdata")

MatchedLK<-merge(x=MatchedLK, y=CaseID, by="patid", all.x=TRUE, all.y=FALSE)

#Change to list of datasets
datlist<-mids2datlist(ImputedData)

####Create three data sets####
SchizOrig<-subset(MatchedLK, diagn=="Schizophrenia")
BipolarOrig<-subset(MatchedLK, diagn=="Bipolar")
OtherOrig<-subset(MatchedLK, diagn=="Other")

ControlSchiz<-subset(MatchedLK, SMI==0 & case_patid %in% SchizOrig$patid)
ControlBip<-subset(MatchedLK, SMI==0 & case_patid %in% BipolarOrig$patid)
ControlOther<-subset(MatchedLK, SMI==0 & case_patid %in% OtherOrig$patid)

SchizAll<-rbind(SchizOrig, ControlSchiz)
BipolarAll<-rbind(BipolarOrig, ControlBip)
OtherAll<-rbind(OtherOrig, ControlOther)

rm(SchizOrig, ControlSchiz, BipolarOrig, ControlBip, OtherOrig, ControlOther)

#Create an empty dataset
#Do the first
Neuro1 <- lapply(datlist, FUN=function(data){
  glm.cluster( data=data, formula=SMI ~ Neuro_MS_prior+gender+region+AgeAtDiag+ethnicity+calendaryear,
               cluster=data$pracid, family=binomial(link="logit") )
} )

betas <- lapply(Neuro1, FUN=function(rr){ coef(rr) } )
vars <- lapply(Neuro1, FUN=function(rr){ vcov(rr) } )

# conduct statistical inference
Neuro1<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
Neuro1$OR<-exp(Neuro1$results)
Neuro1$Lower<-exp(Neuro1$`(lower`)
Neuro1$Upper<-exp(Neuro1$`upper)`)

####Min adjusted in B5####
Neuro<-vars_select(names(ImputedData$data), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('Prior')))

####Run for b5 schizophrenia####

#Select only those in the imputed data who have schizophrenia
ImpSchiz<-subset_datlist(datlist, index=1:5,
                         subset=datlist[[1]]$patid %in% SchizAll$patid)
#Create an empty data set
b5Schiz<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpSchiz, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  b5Schiz<-rbind(b5Schiz, Neuro2)
}

b5Schiz<-rownames_to_column(b5Schiz, var="NeuroType")
b5Schiz<-subset(b5Schiz, startsWith(b5Schiz$NeuroType, "Neuro"))
b5Schiz$Time<-"-5"
b5Schiz$Diagn<-"Schizophrenia"
b5Schiz<-select(b5Schiz, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars)

####Run for b5 bipolar####

#Select only those in the imputed data who have bipolar
ImpBP<-subset_datlist(datlist, index=1:5,
                      subset=datlist[[1]]$patid %in% BipolarAll$patid)

#Create an empty data set
b5BP<-Neuro1[NULL,]
#Wrap it in loop


for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpBP, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  b5BP<-rbind(b5BP, Neuro2)
}

b5BP<-rownames_to_column(b5BP, var="NeuroType")
b5BP<-subset(b5BP, startsWith(b5BP$NeuroType, "Neuro"))
b5BP$Time<-"-5"
b5BP$Diagn<-"Bipolar disorder"
b5BP<-select(b5BP, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars)

#### Run for b5 other psychoses####

#Select only those in the imputed data who have other psychoses
ImpOther<-subset_datlist(datlist, index=1:5,
                         subset=datlist[[1]]$patid %in% OtherAll$patid)

#Create an empty data set
b5Other<-Neuro1[NULL,]
#Wrap it in loop
for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpOther, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  b5Other<-rbind(b5Other, Neuro2)
}

b5Other<-rownames_to_column(b5Other, var="NeuroType")
b5Other<-subset(b5Other, startsWith(b5Other$NeuroType, "Neuro"))
b5Other$Time<-"-5"
b5Other$Diagn<-"Other Psychoses"
b5Other<-select(b5Other, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars, Neuro)

####Create one dataset for b5####
b5<-rbind(b5Schiz, b5BP, b5Other)

write.csv(b5, file = "/Neuro/Outputs/MinAdjDiagB5.csv")


####Min Adjusted, diag, in -3####
Neuro<-vars_select(names(ImputedData$data), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('B4')))


####Run for b3 schizophrenia####

#Select only those in the imputed data who have schizophrenia
ImpSchiz<-subset_datlist(datlist, index=1:5,
                         subset=datlist[[1]]$patid %in% SchizAll$patid)
#Create an empty data set
b3Schiz<-Neuro1[NULL,]
#Wrap it in loop


for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpSchiz, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  b3Schiz<-rbind(b3Schiz, Neuro2)
}

b3Schiz<-rownames_to_column(b3Schiz, var="NeuroType")
b3Schiz<-subset(b3Schiz, startsWith(b3Schiz$NeuroType, "Neuro"))
b3Schiz$Time<-"-3"
b3Schiz$Diagn<-"Schizophrenia"
b3Schiz<-select(b3Schiz, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars)

####Run for b3 bipolar####

#Select only those in the imputed data who have bipolar
ImpBP<-subset_datlist(datlist, index=1:5,
                      subset=datlist[[1]]$patid %in% BipolarAll$patid)

#Create an empty data set
b3BP<-Neuro1[NULL,]

#Wrap it in loop

b3<-Neuro1[NULL,]

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpBP, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  b3BP<-rbind(b3BP, Neuro2)
}

b3BP<-rownames_to_column(b3BP, var="NeuroType")
b3BP<-subset(b3BP, startsWith(b3BP$NeuroType, "Neuro"))
b3BP$Time<-"-3"
b3BP$Diagn<-"Bipolar disorder"
b3BP<-select(b3BP, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars)

####Run for b3 other psychoses####

#Select only those in the imputed data who have other psychoses
ImpOther<-subset_datlist(datlist, index=1:5,
                         subset=datlist[[1]]$patid %in% OtherAll$patid)

#Create an empty data set
b3Other<-Neuro1[NULL,]

#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpOther, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  b3Other<-rbind(b3Other, Neuro2)
}

b3Other<-rownames_to_column(b3Other, var="NeuroType")
b3Other<-subset(b3Other, startsWith(b3Other$NeuroType, "Neuro"))
b3Other$Time<-"-3"
b3Other$Diagn<-"Other Psychoses"
b3Other<-select(b3Other, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars, Neuro)

####Create one dataset for b3####
b3<-rbind(b3Schiz, b3BP, b3Other)

write.csv(b3, file = "/Neuro/Outputs/MinAdjDiagB3.csv")


####Minimally Adjusted, diag, in -1####
Neuro<-vars_select(names(ImputedData$data), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('B2')))

####Run for b1 schizophrenia####

#Select only those in the imputed data who have schizophrenia
ImpSchiz<-subset_datlist(datlist, index=1:5,
                         subset=datlist[[1]]$patid %in% SchizAll$patid)
#Create an empty data set
b1Schiz<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpSchiz, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  b1Schiz<-rbind(b1Schiz, Neuro2)
}

b1Schiz<-rownames_to_column(b1Schiz, var="NeuroType")
b1Schiz<-subset(b1Schiz, startsWith(b1Schiz$NeuroType, "Neuro"))
b1Schiz$Time<-"-1"
b1Schiz$Diagn<-"Schizophrenia"
b1Schiz<-select(b1Schiz, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars)

####Run for b1 bipolar####

#Select only those in the imputed data who have bipolar
ImpBP<-subset_datlist(datlist, index=1:5,
                      subset=datlist[[1]]$patid %in% BipolarAll$patid)

#Create an empty data set
b1BP<-Neuro1[NULL,]

#Wrap it in loop

b1<-Neuro1[NULL,]

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpBP, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  b1BP<-rbind(b1BP, Neuro2)
}

b1BP<-rownames_to_column(b1BP, var="NeuroType")
b1BP<-subset(b1BP, startsWith(b1BP$NeuroType, "Neuro"))
b1BP$Time<-"-1"
b1BP$Diagn<-"Bipolar disorder"
b1BP<-select(b1BP, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars)

####Run for b1 other psychoses####

#Select only those in the imputed data who have other psychoses
ImpOther<-subset_datlist(datlist, index=1:5,
                         subset=datlist[[1]]$patid %in% OtherAll$patid)

#Create an empty data set
b1Other<-Neuro1[NULL,]

#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpOther, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  b1Other<-rbind(b1Other, Neuro2)
}

b1Other<-rownames_to_column(b1Other, var="NeuroType")
b1Other<-subset(b1Other, startsWith(b1Other$NeuroType, "Neuro"))
b1Other$Time<-"-1"
b1Other$Diagn<-"Other Psychoses"
b1Other<-select(b1Other, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars, Neuro)

####Create one dataset for b1####
b1<-rbind(b1Schiz, b1BP, b1Other)

write.csv(b1, file = "/Neuro/Outputs/MinAdjDiagB1.csv")


####Min Adjusted, Diag, at index - NOTE renamed index b0 as b1 used above for -1 year####
Neuro<-vars_select(names(ImputedData$data), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('B1')))

####Run for b0 index schizophrenia####

#Select only those in the imputed data who have schizophrenia
ImpSchiz<-subset_datlist(datlist, index=1:5,
                         subset=datlist[[1]]$patid %in% SchizAll$patid)
#Create an empty data set
b0Schiz<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpSchiz, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  b0Schiz<-rbind(b0Schiz, Neuro2)
}

b0Schiz<-rownames_to_column(b0Schiz, var="NeuroType")
b0Schiz<-subset(b0Schiz, startsWith(b0Schiz$NeuroType, "Neuro"))
b0Schiz$Time<-"0"
b0Schiz$Diagn<-"Schizophrenia"
b0Schiz<-select(b0Schiz, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars)

####Run for b0 index bipolar####

#Select only those in the imputed data who have bipolar
ImpBP<-subset_datlist(datlist, index=1:5,
                      subset=datlist[[1]]$patid %in% BipolarAll$patid)

#Create an empty data set
b0BP<-Neuro1[NULL,]

#Wrap it in loop

b0<-Neuro1[NULL,]


for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpBP, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  b0BP<-rbind(b0BP, Neuro2)
}

b0BP<-rownames_to_column(b0BP, var="NeuroType")
b0BP<-subset(b0BP, startsWith(b0BP$NeuroType, "Neuro"))
b0BP$Time<-"0"
b0BP$Diagn<-"Bipolar disorder"
b0BP<-select(b0BP, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars)

####Run for b0 index other psychoses####

#Select only those in the imputed data who have other psychoses
ImpOther<-subset_datlist(datlist, index=1:5,
                         subset=datlist[[1]]$patid %in% OtherAll$patid)

#Create an empty data set
b0Other<-Neuro1[NULL,]

#Wrap it in loop


for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpOther, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  b0Other<-rbind(b0Other, Neuro2)
}

b0Other<-rownames_to_column(b0Other, var="NeuroType")
b0Other<-subset(b0Other, startsWith(b0Other$NeuroType, "Neuro"))
b0Other$Time<-"0"
b0Other$Diagn<-"Other Psychoses"
b0Other<-select(b0Other, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars, Neuro)

####Create one dataset for b0 index####
b0<-rbind(b0Schiz, b0BP, b0Other)

write.csv(b0, file = "/Neuro/Outputs/MinAdjDiagB0.csv")


####Minimally Adj Diag: A1####
ActA1<-subset(MatchedLK, Act_a1==1)

Neuro<-vars_select(names(ImputedData$data), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('A1')))

####Schizophrenia A1####
#Select only those in the imputed data who appear in the descriptive data as active

ImpA1Schiz<-subset_datlist(datlist, index=1:5,
                           subset=datlist[[1]]$patid %in% ActA1$patid & datlist[[1]]$patid %in% SchizAll$patid)

#Create an empty data set
A1Schiz<-Neuro1[NULL,]

#Wrap it in loop


for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpA1Schiz, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  A1Schiz<-rbind(A1Schiz, Neuro2)
}

A1Schiz<-rownames_to_column(A1Schiz, var="NeuroType")
A1Schiz<-subset(A1Schiz, startsWith(A1Schiz$NeuroType, "Neuro"))
A1Schiz$Time<-"+1"
A1Schiz$Diagn<-"Schizophrenia"
A1Schiz<-select(A1Schiz, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars)


####Run for a1 bipolar####

#Select only those in the imputed data who have bipolar
ImpA1BP<-subset_datlist(datlist, index=1:5,
                        subset=datlist[[1]]$patid %in% BipolarAll$patid)

#Create an empty data set
A1BP<-Neuro1[NULL,]

#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpA1BP, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  A1BP<-rbind(A1BP, Neuro2)
}

A1BP<-rownames_to_column(A1BP, var="NeuroType")
A1BP<-subset(A1BP, startsWith(A1BP$NeuroType, "Neuro"))
A1BP$Time<-"+1"
A1BP$Diagn<-"Bipolar disorder"
A1BP<-select(A1BP, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars)

####Run for A1  other psychoses####

#Select only those in the imputed data who have other psychoses
ImpA1Other<-subset_datlist(datlist, index=1:5,
                           subset=datlist[[1]]$patid %in% OtherAll$patid)

#Create an empty data set
A1Other<-Neuro1[NULL,]

#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpA1Other, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  A1Other<-rbind(A1Other, Neuro2)
}

A1Other<-rownames_to_column(A1Other, var="NeuroType")
A1Other<-subset(A1Other, startsWith(A1Other$NeuroType, "Neuro"))
A1Other$Time<-"+1"
A1Other$Diagn<-"Other Psychoses"
A1Other<-select(A1Other, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars, Neuro)

####Create one dataset for a1####
a1<-rbind(A1Schiz, A1BP, A1Other)

write.csv(a1, file = "/Neuro/Outputs/MinAdjDiagA1.csv")


####Minimally Adj diag: A3####
ActA3<-subset(MatchedLK, Act_a3==1)

Neuro<-vars_select(names(ImputedData$data), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('A3')))

####Schizophrenia A3####
#Select only those in the imputed data who appear in the descriptive data as active 

ImpA3Schiz<-subset_datlist(datlist, index=1:5,
                           subset=datlist[[1]]$patid %in% ActA3$patid & datlist[[1]]$patid %in% SchizAll$patid)

#Create an empty data set
A3Schiz<-Neuro1[NULL,]

#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpA3Schiz, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  A3Schiz<-rbind(A3Schiz, Neuro2)
}

A3Schiz<-rownames_to_column(A3Schiz, var="NeuroType")
A3Schiz<-subset(A3Schiz, startsWith(A3Schiz$NeuroType, "Neuro"))
A3Schiz$Time<-"+3"
A3Schiz$Diagn<-"Schizophrenia"
A3Schiz<-select(A3Schiz, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars)

####Run for a3 bipolar####

#Select only those in the imputed data who have bipolar
ImpA3BP<-subset_datlist(datlist, index=1:5,
                        subset=datlist[[1]]$patid %in% BipolarAll$patid)

#Create an empty data set
A3BP<-Neuro1[NULL,]

#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpA3BP, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  A3BP<-rbind(A3BP, Neuro2)
}

A3BP<-rownames_to_column(A3BP, var="NeuroType")
A3BP<-subset(A3BP, startsWith(A3BP$NeuroType, "Neuro"))
A3BP$Time<-"+3"
A3BP$Diagn<-"Bipolar disorder"
A3BP<-select(A3BP, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars)

####Run for A3  other psychoses####

#Select only those in the imputed data who have other psychoses
ImpA3Other<-subset_datlist(datlist, index=1:5,
                           subset=datlist[[1]]$patid %in% OtherAll$patid)

#Create an empty data set
A3Other<-Neuro1[NULL,]

#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpA3Other, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  A3Other<-rbind(A3Other, Neuro2)
}

A3Other<-rownames_to_column(A3Other, var="NeuroType")
A3Other<-subset(A3Other, startsWith(A3Other$NeuroType, "Neuro"))
A3Other$Time<-"+3"
A3Other$Diagn<-"Other Psychoses"
A3Other<-select(A3Other, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars, Neuro)

####Create one dataset for a3####
a3<-rbind(A3Schiz, A3BP, A3Other)

write.csv(a3, file = "/Neuro/Outputs/MinAdjDiagA3.csv")

####Minimally Adj diag: A5####
ActA5<-subset(MatchedLK, Act_a5==1)

Neuro<-vars_select(names(ImputedData$data), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('A5')))

####Schizophrenia A5####
#Select only those in the imputed data who appear in the descriptive data as active
ImpA5Schiz<-subset_datlist(datlist, index=1:5,
                           subset=datlist[[1]]$patid %in% ActA5$patid & datlist[[1]]$patid %in% SchizAll$patid)

#Create an empty data set
A5Schiz<-Neuro1[NULL,]

#Wrap it in loop


for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpA5Schiz, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  A5Schiz<-rbind(A5Schiz, Neuro2)
}

A5Schiz<-rownames_to_column(A5Schiz, var="NeuroType")
A5Schiz<-subset(A5Schiz, startsWith(A5Schiz$NeuroType, "Neuro"))
A5Schiz$Time<-"+5"
A5Schiz$Diagn<-"Schizophrenia"
A5Schiz<-select(A5Schiz, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars)

####Run for a5 bipolar####

#Select only those in the imputed data who have bipolar
ImpA5BP<-subset_datlist(datlist, index=1:5,
                        subset=datlist[[1]]$patid %in% BipolarAll$patid)

#Create an empty data set
A5BP<-Neuro1[NULL,]

#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpA5BP, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  A5BP<-rbind(A5BP, Neuro2)
}

A5BP<-rownames_to_column(A5BP, var="NeuroType")
A5BP<-subset(A5BP, startsWith(A5BP$NeuroType, "Neuro"))
A5BP$Time<-"+5"
A5BP$Diagn<-"Bipolar disorder"
A5BP<-select(A5BP, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars)

####Run for A5 other psychoses####

#Select only those in the imputed data who have other psychoses
ImpA5Other<-subset_datlist(datlist, index=1:5,
                           subset=datlist[[1]]$patid %in% OtherAll$patid)

#Create an empty data set
A5Other<-Neuro1[NULL,]

#Wrap it in loop


for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpA5Other, FUN=function(data){
    glm.cluster(data=data, as.formula( paste("SMI ~", Neuro[i], "+gender+region+AgeAtDiag+ethnicity+calendaryear")),
                cluster=data$pracid, family=binomial(link="logit") )
  } )
  
  betas <- lapply(Neuro2, FUN=function(rr){ coef(rr) } )
  vars <- lapply(Neuro2, FUN=function(rr){ vcov(rr) } )
  
  # conduct statistical inference
  Neuro2<-as.data.frame(summary(pool_mi(qhat=betas, u=vars)))
  Neuro2$OR<-exp(Neuro2$results)
  Neuro2$Lower<-exp(Neuro2$`(lower`)
  Neuro2$Upper<-exp(Neuro2$`upper)`)
  A5Other<-rbind(A5Other, Neuro2)
}

A5Other<-rownames_to_column(A5Other, var="NeuroType")
A5Other<-subset(A5Other, startsWith(A5Other$NeuroType, "Neuro"))
A5Other$Time<-"+5"
A5Other$Diagn<-"Other Psychoses"
A5Other<-select(A5Other, NeuroType, OR, Lower, Upper, Time, Diagn)

rm(betas, Neuro2, vars, Neuro)

####Create one dataset for a5####
a5<-rbind(A5Schiz, A5BP, A5Other)

write.csv(a5, file = "/Neuro/Outputs/MinAdjDiagMinA5.csv")


####Graphs for each neurological condition####

#Create a dataframe with all our odds ratios in
Counts<-rbind(b5, b3, b1, b0, a1, a3, a5)

save(Counts, file= "/Neuro/Data/MINAdjDiagResults.Rdata")

#Have a numeric variable for time as well as sometimes easier on graphs
Counts$ContTime[Counts$Time=="-5"]<- -5
Counts$ContTime[Counts$Time=="-3"]<- -3
Counts$ContTime[Counts$Time=="-1"]<- -1
Counts$ContTime[Counts$Time=="0"]<- 0
Counts$ContTime[Counts$Time=="+1"]<- 1
Counts$ContTime[Counts$Time=="+3"]<- 3
Counts$ContTime[Counts$Time=="+5"]<- 5

Counts$Diagn[Counts$Diagn=="Other Psychoses"]<-"Other psychoses"
Counts$Diagn<-factor(Counts$Diagn, levels=c("Schizophrenia", "Bipolar disorder", "Other psychoses"))

#Set CI as 7 if greater than 7
Counts$Upper[Counts$Upper>7]<-7

CountsMS<-subset(Counts, startsWith(NeuroType, "Neuro_MS"))

MS<-ggplot()+
  geom_line(data = CountsMS, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=2)+
  geom_point(data = CountsMS, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=5)+
  geom_ribbon(data = CountsMS, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Diagn, fill=Diagn), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  scale_fill_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  scale_color_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs(y = "Partially adjusted odds ratio", title="Multiple sclerosis", x="Time before and after index date (years)")+
  scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

MS

#Save for poster
ggsave("/Neuro/MinAdjMS.tiff", plot=MS, device="tiff", dpi=300)

####Graph for Cerebrovascular disease####
CountsCereb<-subset(Counts, startsWith(NeuroType, "Neuro_Cereb"))

Cereb<-ggplot()+
  geom_line(data = CountsCereb, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=2)+
  geom_point(data = CountsCereb, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=5)+
  geom_ribbon(data = CountsCereb, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Diagn, fill=Diagn), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  scale_fill_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  scale_color_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs(y = "Partially adjusted odds ratio", title="Cerebrovascular disease", x="Time before and after index date (years)")+
  scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

Cereb

#Save for poster
ggsave("/Neuro/MinAdjCereb.tiff", plot=Cereb, device="tiff", dpi=300)

####Graph for Dementia####
CountsDem<-subset(Counts, startsWith(NeuroType, "Neuro_Dementia"))

Dem<-ggplot()+
  geom_line(data = CountsDem, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=2)+
  geom_point(data = CountsDem, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=5)+
  geom_ribbon(data = CountsDem, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Diagn, fill=Diagn), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  scale_fill_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  scale_color_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs(y = "Partially adjusted odds ratio", title="Dementia", x="Time before and after index date (years)")+
  scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+  
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
Dem
#Save for poster
ggsave("/Neuro/MinAdjDem.tiff", plot=Dem, device="tiff", dpi=300)

####Graph for Ataxia####
CountsAtax<-subset(Counts, startsWith(NeuroType, "Neuro_Ataxia"))

Ataxia<-ggplot()+
  geom_line(data = CountsAtax, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=2)+
  geom_point(data = CountsAtax, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=5)+
  geom_ribbon(data = CountsAtax, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Diagn, fill=Diagn), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  scale_fill_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  scale_color_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs(y = "Partially adjusted odds ratio", title="Ataxia", x="Time before and after index date (years)")+
  scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+  
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
Ataxia
#Save for poster
ggsave("/Neuro/MinAdjAtaxia.tiff", plot=Ataxia, device="tiff", dpi=300)
####Graph for Epilepsy####
CountsEpilepsy <-subset(Counts, startsWith(NeuroType, "Neuro_Epilepsy"))

Epi<-ggplot()+
  geom_line(data = CountsEpilepsy, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=2)+
  geom_point(data = CountsEpilepsy, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=5)+
  geom_ribbon(data = CountsEpilepsy, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Diagn, fill=Diagn), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  scale_fill_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  scale_color_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs(y = "Partially adjusted odds ratio", title="Epilepsy", x="Time before and after index date (years)")+
  scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+ 
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

Epi

#Save for poster
ggsave("/Neuro/MinAdjEpi.tiff", plot=Epi, device="tiff", dpi=300)


####Graph for Parkinsons####
CountsPD <-subset(Counts, startsWith(NeuroType, "Neuro_Parkinsons"))

Park<-ggplot()+
  geom_line(data = CountsPD, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=2)+
  geom_point(data = CountsPD, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=5)+
  geom_ribbon(data = CountsPD, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Diagn, fill=Diagn), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  scale_fill_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  scale_color_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs(y = "Partially adjusted odds ratio", title="Parkinson's disease", x="Time before and after index date (years)")+
  scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

Park

#Save for poster
ggsave("/Neuro/MinAdjPark.tiff", plot=Park, device="tiff", dpi=300)


####Graph for Paralysis####
CountsParalysis <-subset(Counts, startsWith(NeuroType, "Neuro_Paralysis"))

Para<-ggplot()+
  geom_line(data = CountsParalysis, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=2)+
  geom_point(data = CountsParalysis, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=5)+
  geom_ribbon(data = CountsParalysis, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Diagn, fill=Diagn), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  scale_fill_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  scale_color_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs(y = "Partially adjusted odds ratio", title="Paralysis", x="Time before and after index date (years)")+
  scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+ 
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

Para

#Save for poster
ggsave("/Neuro/MinAdjPara.tiff", plot=Para, device="tiff", dpi=300)

####Graph for Cerebral Palsy####
CountsCP <-subset(Counts, startsWith(NeuroType, "Neuro_CP"))

CP<-ggplot()+
  geom_line(data = CountsCP, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=2)+
  geom_point(data = CountsCP, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=5)+
  geom_ribbon(data = CountsCP, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Diagn, fill=Diagn), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  scale_fill_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  scale_color_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs(y = "Partially adjusted odds ratio", title="Cerebral Palsy", x="Time before and after index date (years)")+
  scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

CP

#Save for poster
ggsave("/Neuro/MinAdjCP.tiff", plot=CP, device="tiff", dpi=300)

####Graph for CSF disorders####
CountsCSF <-subset(Counts, startsWith(NeuroType, "Neuro_CSF"))

CSF<-ggplot()+
  geom_line(data = CountsCSF, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=2)+
  geom_point(data = CountsCSF, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=5)+
  geom_ribbon(data = CountsCSF, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Diagn, fill=Diagn), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  scale_fill_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  scale_color_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs(y = "Partially adjusted odds ratio", title="CSF disorders", x="Time before and after index date (years)")+
  scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+ 
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
CSF
#Save for poster
ggsave("/Neuro/MinAdjCSF.tiff", plot=CSF, device="tiff", dpi=300)

####Graph for Spinal Cord disorders####
CountsSpinal <-subset(Counts, startsWith(NeuroType, "Neuro_Spinal"))

Spinal<-ggplot()+
  geom_line(data = CountsSpinal, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=2)+
  geom_point(data = CountsSpinal, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=5)+
  geom_ribbon(data = CountsSpinal, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Diagn, fill=Diagn), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  scale_fill_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  scale_color_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs(y = "Partially adjusted odds ratio", title="Spinal cord disorders", x="Time before and after index date (years)")+
  scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+ 
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

Spinal

ggsave("/Neuro/MinAdjSpinal.tiff", plot=Spinal, device="tiff", dpi=300)
####Graph for Parkinsonism other disorders####
###note this one may take some playing around to fit - range of confidence intervals is to 176!
# CountsParkinsonOther <-subset(Counts,startsWith(NeuroType, "Neuro_ParkinsonOther"))
# 
# ParkO<-ggplot()+
#   geom_line(data = CountsParkinsonOther, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=2)+
#   geom_point(data = CountsParkinsonOther, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=5)+
#   geom_ribbon(data = CountsParkinsonOther, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Diagn, fill=Diagn), alpha=0.3)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 14))+
#   theme(axis.text.y = element_text(size = 14))+
#   theme(axis.title.x=element_text(size=16))+
#   theme(axis.title.y=element_text(size=16))+
#   theme(title=element_text(size=16))+
#   scale_fill_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
#   scale_color_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
#   theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
#   labs(y = "Partially adjusted odds ratio", title="Parkinsonism other disorders", x="Time before and after index date (years)")+
#   scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+   
#   scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
# ParkO
# ggsave("/Neuro/MinAdjParkO.tiff", plot=ParkO, device="tiff", dpi=300)
####Graph for Disorders of nerve root, plexus or peripheral nerves####
CountsPerip <-subset(Counts, startsWith(NeuroType, "Neuro_Peripheral"))

Periph<-ggplot()+
  geom_line(data = CountsPerip, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=2)+
  geom_point(data = CountsPerip, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=5)+
  geom_ribbon(data = CountsPerip, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Diagn, fill=Diagn), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  scale_fill_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  scale_color_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs(y = "Partially adjusted odds ratio", title="Peripheral nerve disorders", x="Time before and after index date (years)")+
  scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+    
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
Periph
ggsave("/Neuro/MinAdjPeriph.tiff", plot=Periph, device="tiff", dpi=300)
####Graph for -	Encephalopathy####
# CountsEnceph <-subset(Counts, startsWith(NeuroType, "Neuro_Enceph"))
# 
# ggplot()+
#   geom_line(data = CountsEnceph, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=2)+
#   geom_point(data = CountsEnceph, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=5)+
#   geom_ribbon(data = CountsEnceph, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Diagn, fill=Diagn), alpha=0.3)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 14))+
#   theme(axis.text.y = element_text(size = 14))+
#   theme(axis.title.x=element_text(size=16))+
#   theme(axis.title.y=element_text(size=16))+
#   theme(title=element_text(size=16))+
#   scale_fill_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
#   scale_color_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
#   theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
#   labs( , y = "Odds ratio", title="Encephalopathy", x="Time before and after index date (years)")+
#   scale_y_continuous(limits = c(0.5, 21), breaks=c(0.5, 1, seq(2, 21, by=2)))+   
#   scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

#may need to adjust break here - highest CI is 20.89 

####Graph for Movement disorders (other) ####
CountsMovement <-subset(Counts, startsWith(NeuroType, "Neuro_Movement"))

Move<-ggplot()+
  geom_line(data = CountsMovement, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=2)+
  geom_point(data = CountsMovement, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=5)+
  geom_ribbon(data = CountsMovement, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Diagn, fill=Diagn), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  scale_fill_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  scale_color_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
  theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs(y = "Partially adjusted odds ratio", title="Movement disorders, other", x="Time before and after index date (years)")+
  scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
Move

ggsave("/Neuro/MinAdjMove.tiff", plot=Move, device="tiff", dpi=300)
####Graph for Motor neurone disease ####
# CountsMND <-subset(Counts, startsWith(NeuroType, "Neuro_MND"))
# 
# MND<-ggplot()+
#   geom_line(data = CountsMND, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=2)+
#   geom_point(data = CountsMND, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=5)+
#   geom_ribbon(data = CountsMND, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Diagn, fill=Diagn), alpha=0.3)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 14))+
#   theme(axis.text.y = element_text(size = 14))+
#   theme(axis.title.x=element_text(size=16))+
#   theme(axis.title.y=element_text(size=16))+
#   theme(title=element_text(size=16))+
#   scale_fill_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
#   scale_color_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
#   theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
#   labs(y = "Partially adjusted odds ratio", title="Motor Neurone Disease", x="Time before and after index date (years)")+
#   scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+
#   scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
# MND

ggsave("/Neuro/MinAdjMND.tiff", plot=MND, device="tiff", dpi=300)
#may need to adjust break here - highest CI is 42.03

####Graph for Autonomic Nervous System Disorders ####
# CountsAuto <-subset(Counts, startsWith(NeuroType, "Neuro_Autonomic"))
# 
# Auto<-ggplot()+
#   geom_line(data = CountsAuto, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=2)+
#   geom_point(data = CountsAuto, aes(x=ContTime, y=OR, group=Diagn, color=Diagn), size=5)+
#   geom_ribbon(data = CountsAuto, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Diagn, fill=Diagn), alpha=0.3)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 14))+
#   theme(axis.text.y = element_text(size = 14))+
#   theme(axis.title.x=element_text(size=16))+
#   theme(axis.title.y=element_text(size=16))+
#   theme(title=element_text(size=16))+
#   scale_fill_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
#   scale_color_manual(values = c("dodgerblue2", "firebrick2", "gold1" ))+
#   theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
#   labs(y = "Partially adjusted odds ratio", title="Autonomic Nervous System Disorders", x="Time before and after index date (years)")+
#   scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+ 
#   scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
# 
# Auto
# 
# ggsave("/Neuro/MinAdjAuto.tiff", plot=Auto, device="tiff", dpi=300)
#may need to adjust break here - highest CI is 147