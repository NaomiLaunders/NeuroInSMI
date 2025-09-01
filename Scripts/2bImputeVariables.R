#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Impute missing BMI and ethnicity variables
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(mice)
library(miceadds)
library(mitools)
library(lubridate)
library(tidyverse)
library(tidyselect)

rm(list = ls(all.names = TRUE))

load("/Neuro/Data/FinalNeuroData.Rdata")

MatchedLK$ethnicity<-as.factor(MatchedLK$ethnicity)
MatchedLK$BMI5<-as.factor(MatchedLK$BMI5)
MatchedLK$Smoke5<-as.factor(MatchedLK$Smoke5)
MatchedLK$Drugs5<-as.factor(MatchedLK$Drugs5)
MatchedLK$Alc5<-as.factor(MatchedLK$Alc5)
MatchedLK$region<-as.factor(MatchedLK$region)

Vars<-vars_select(names(MatchedLK), starts_with('Neuro_', ignore.case = TRUE))

MatchedLK[,Vars]<-lapply(MatchedLK[,Vars], as.factor)


ToImpute<-dplyr::select(MatchedLK, patid, pracid, SMI, gender, diagn, region, AgeAtDiag, ethnicity, Smoke5, BMI5, Alc5, Drugs5, calendaryear, all_of(Vars))
sapply(ToImpute[, c(1:13)], class)

ImputedData<-mice(data=ToImpute, m=5, seed=500)

summary(ImputedData$imp$ethnicity)
summary(ImputedData$imp$BMI5)
ImputedData$method
Pred<-ImputedData$predictorMatrix
save(ImputedData, file="/Neuro/Data/ImpAnalysis.Rdata")

comp_imp<-complete(ImputedData,"long")
save(comp_imp, file="/Neuro/Data/Full_Imp.Rdata")
