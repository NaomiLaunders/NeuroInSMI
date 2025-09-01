#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Unadjusted results - schizophrenia, bipolar and other psychoses vs control
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
load("/Neuro/Data/FinalNeuroData.Rdata")

MatchedLK<-merge(x=MatchedLK, y=CaseID, by="patid", all.x=TRUE, all.y=FALSE)

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
Neuro1 <- glm( data=SchizAll, formula=SMI ~ Neuro_MS_prior,
               family=binomial(link="logit") )


Sand<-vcovCL(Neuro1, cluster =SchizAll$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(Neuro1, vcov=Sand)))

Neuro1<-tidy(Neuro1)
Neuro1$OR<-exp(Neuro1$estimate)
Neuro1<-cbind(Neuro1, CIFin)

####Un adjusted in B5####
Neuro<-vars_select(names(MatchedLK), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('Prior')))

####Run for b5 schizophrenia####

#Create an empty data set
b5Schiz<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=SchizAll, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =SchizAll$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))

    
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  b5Schiz<-rbind(b5Schiz, Neuro2)
}

b5Schiz<-rownames_to_column(b5Schiz, var="NeuroType")
b5Schiz<-subset(b5Schiz, startsWith(b5Schiz$NeuroType, "Neuro"))
b5Schiz$Time<-"-5"
b5Schiz<-select(b5Schiz, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2)

####Run for b5 bipolar####

#Create an empty data set
b5BP<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=BipolarAll, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =BipolarAll$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  b5BP<-rbind(b5BP, Neuro2)
}

b5BP<-rownames_to_column(b5BP, var="NeuroType")
b5BP<-subset(b5BP, startsWith(b5BP$NeuroType, "Neuro"))
b5BP$Time<-"-5"
b5BP<-select(b5BP, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2)

#### Run for b5 other psychoses####

#Create an empty data set
b5Other<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=OtherAll, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =OtherAll$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  b5Other<-rbind(b5Other, Neuro2)
}

b5Other<-rownames_to_column(b5Other, var="NeuroType")
b5Other<-subset(b5Other, startsWith(b5Other$NeuroType, "Neuro"))
b5Other$Time<-"-5"
b5Other<-select(b5Other, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2, Neuro)

####Create one dataset for b5####
b5<-rbind(b5Schiz, b5BP, b5Other)

write.csv(b5, file = "/Neuro/Outputs/UnAdjDiagB5.csv")


####Un Adjusted, diag, in -3####
Neuro<-vars_select(names(MatchedLK), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('B4')))

####Run for b3 schizophrenia####

#Create an empty data set
b3Schiz<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=SchizAll, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =SchizAll$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  b3Schiz<-rbind(b3Schiz, Neuro2)
}

b3Schiz<-rownames_to_column(b3Schiz, var="NeuroType")
b3Schiz<-subset(b3Schiz, startsWith(b3Schiz$NeuroType, "Neuro"))
b3Schiz$Time<-"-3"
b3Schiz<-select(b3Schiz, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2)

####Run for b3 bipolar####

#Create an empty data set
b3BP<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=BipolarAll, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =BipolarAll$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  b3BP<-rbind(b3BP, Neuro2)
}

b3BP<-rownames_to_column(b3BP, var="NeuroType")
b3BP<-subset(b3BP, startsWith(b3BP$NeuroType, "Neuro"))
b3BP$Time<-"-3"
b3BP<-select(b3BP, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2)

#### Run for b3 other psychoses####

#Create an empty data set
b3Other<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=OtherAll, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =OtherAll$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  b3Other<-rbind(b3Other, Neuro2)
}

b3Other<-rownames_to_column(b3Other, var="NeuroType")
b3Other<-subset(b3Other, startsWith(b3Other$NeuroType, "Neuro"))
b3Other$Time<-"-3"
b3Other<-select(b3Other, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2, Neuro)

####Create one dataset for b3####
b3<-rbind(b3Schiz, b3BP, b3Other)

write.csv(b3, file = "/Neuro/Outputs/UnAdjDiagB3.csv")


####Un Adjusted, diag, in -1####
Neuro<-vars_select(names(MatchedLK), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('B2')))

####Run for b1 schizophrenia####

#Create an empty data set
b1Schiz<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=SchizAll, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =SchizAll$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  b1Schiz<-rbind(b1Schiz, Neuro2)
}

b1Schiz<-rownames_to_column(b1Schiz, var="NeuroType")
b1Schiz<-subset(b1Schiz, startsWith(b1Schiz$NeuroType, "Neuro"))
b1Schiz$Time<-"-1"
b1Schiz<-select(b1Schiz, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2)

####Run for b1 bipolar####

#Create an empty data set
b1BP<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=BipolarAll, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =BipolarAll$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  b1BP<-rbind(b1BP, Neuro2)
}

b1BP<-rownames_to_column(b1BP, var="NeuroType")
b1BP<-subset(b1BP, startsWith(b1BP$NeuroType, "Neuro"))
b1BP$Time<-"-1"
b1BP<-select(b1BP, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2)

#### Run for b1 other psychoses####

#Create an empty data set
b1Other<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=OtherAll, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =OtherAll$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  b1Other<-rbind(b1Other, Neuro2)
}

b1Other<-rownames_to_column(b1Other, var="NeuroType")
b1Other<-subset(b1Other, startsWith(b1Other$NeuroType, "Neuro"))
b1Other$Time<-"-1"
b1Other<-select(b1Other, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2, Neuro)

####Create one dataset for b1####
b1<-rbind(b1Schiz, b1BP, b1Other)

write.csv(b1, file = "/Neuro/Outputs/UnAdjDiagB1.csv")

####UnAdjusted, diag, at index - NOTE renamed index b0 as b1 used above for -1 year####
Neuro<-vars_select(names(MatchedLK), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('B1')))

####Run for b0 index schizophrenia####

#Create an empty data set
b0Schiz<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=SchizAll, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =SchizAll$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  b0Schiz<-rbind(b0Schiz, Neuro2)
}

b0Schiz<-rownames_to_column(b0Schiz, var="NeuroType")
b0Schiz<-subset(b0Schiz, startsWith(b0Schiz$NeuroType, "Neuro"))
b0Schiz$Time<-"0"
b0Schiz<-select(b0Schiz, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2)

####Run for b0 index bipolar####

#Create an empty data set
b0BP<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=BipolarAll, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =BipolarAll$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  b0BP<-rbind(b0BP, Neuro2)
}

b0BP<-rownames_to_column(b0BP, var="NeuroType")
b0BP<-subset(b0BP, startsWith(b0BP$NeuroType, "Neuro"))
b0BP$Time<-"0"
b0BP<-select(b0BP, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2)

#### Run for b0 index other psychoses####

#Create an empty data set
b0Other<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=OtherAll, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =OtherAll$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  b0Other<-rbind(b0Other, Neuro2)
}

b0Other<-rownames_to_column(b0Other, var="NeuroType")
b0Other<-subset(b0Other, startsWith(b0Other$NeuroType, "Neuro"))
b0Other$Time<-"0"
b0Other<-select(b0Other, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2, Neuro)

####Create one dataset for b0 index####
b0<-rbind(b0Schiz, b0BP, b0Other)

write.csv(b0, file = "/Neuro/Outputs/UnAdjDiagB0.csv")

####UnAdj diag: A1####
ActA1<-subset(MatchedLK, Act_a1==1)

Neuro<-vars_select(names(MatchedLK), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('A1')))

####Schizophrenia A1####
#Select only those in the imputed data who appear in the descriptive data as active 

Acta1Schiz<-subset(SchizAll, patid %in% ActA1$patid)

#Create an empty data set
a1Schiz<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=Acta1Schiz, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =Acta1Schiz$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  a1Schiz<-rbind(a1Schiz, Neuro2)
}

a1Schiz<-rownames_to_column(a1Schiz, var="NeuroType")
a1Schiz<-subset(a1Schiz, startsWith(a1Schiz$NeuroType, "Neuro"))
a1Schiz$Time<-"+1"
a1Schiz<-select(a1Schiz, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2)

####Run for a1 bipolar####

#Select only those in the imputed data who have bipolar
Acta1BP<-subset(BipolarAll, patid %in% ActA1$patid)

#Create an empty data set
a1BP<-Neuro1[NULL,]
#Wrap it in loop
for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=Acta1BP, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =Acta1BP$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  a1BP<-rbind(a1BP, Neuro2)
}

a1BP<-rownames_to_column(a1BP, var="NeuroType")
a1BP<-subset(a1BP, startsWith(a1BP$NeuroType, "Neuro"))
a1BP$Time<-"+1"
a1BP<-select(a1BP, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2)

####Run for A1  other psychoses####

#Select only those in the imputed data who have other psychoses
Acta1Other<-subset(OtherAll, patid %in% ActA1$patid)

#Create an empty data set
a1Other<-Neuro1[NULL,]
#Wrap it in loop
for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=Acta1Other, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =Acta1Other$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  a1Other<-rbind(a1Other, Neuro2)
}

a1Other<-rownames_to_column(a1Other, var="NeuroType")
a1Other<-subset(a1Other, startsWith(a1Other$NeuroType, "Neuro"))
a1Other$Time<-"+1"
a1Other<-select(a1Other, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2, Neuro)

####Create one dataset for a1####
a1<-rbind(a1Schiz, a1BP, a1Other)

write.csv(a1, file = "/Neuro/Outputs/UnAdjDiagA1.csv")

####UnAdj diag: A3####
ActA3<-subset(MatchedLK, Act_a3==1)

Neuro<-vars_select(names(MatchedLK), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('A3')))

####Schizophrenia A3####

#Select only those in the imputed data who appear in the descriptive data as active
Acta3Schiz<-subset(SchizAll, patid %in% ActA3$patid)

#Create an empty data set
a3Schiz<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=Acta3Schiz, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =Acta3Schiz$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  a3Schiz<-rbind(a3Schiz, Neuro2)
}

a3Schiz<-rownames_to_column(a3Schiz, var="NeuroType")
a3Schiz<-subset(a3Schiz, startsWith(a3Schiz$NeuroType, "Neuro"))
a3Schiz$Time<-"+3"
a3Schiz<-select(a3Schiz, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2)

####Run for a3 bipolar####

#Select only those in the imputed data who have bipolar
Acta3BP<-subset(BipolarAll, patid %in% ActA3$patid)

#Create an empty data set
a3BP<-Neuro1[NULL,]
#Wrap it in loop
for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=Acta3BP, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =Acta3BP$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  a3BP<-rbind(a3BP, Neuro2)
}

a3BP<-rownames_to_column(a3BP, var="NeuroType")
a3BP<-subset(a3BP, startsWith(a3BP$NeuroType, "Neuro"))
a3BP$Time<-"+3"
a3BP<-select(a3BP, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2)

####Run for A3  other psychoses####

#Select only those in the imputed data who have other psychoses
Acta3Other<-subset(OtherAll, patid %in% ActA3$patid)

#Create an empty data set
a3Other<-Neuro1[NULL,]
#Wrap it in loop
for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=Acta3Other, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =Acta3Other$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  a3Other<-rbind(a3Other, Neuro2)
}

a3Other<-rownames_to_column(a3Other, var="NeuroType")
a3Other<-subset(a3Other, startsWith(a3Other$NeuroType, "Neuro"))
a3Other$Time<-"+3"
a3Other<-select(a3Other, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2, Neuro)

####Create one dataset for a3####
a3<-rbind(a3Schiz, a3BP, a3Other)

write.csv(a3, file = "/Neuro/Outputs/UnAdjDiagA3.csv")

####Un Adj diag: A5####
ActA5<-subset(MatchedLK, Act_a5==1)

Neuro<-vars_select(names(MatchedLK), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('A5')))

####Schizophrenia A5####

#Select only those in the imputed data who appear in the descriptive data as active
Acta5Schiz<-subset(SchizAll, patid %in% ActA5$patid)

#Create an empty data set
a5Schiz<-Neuro1[NULL,]
#Wrap it in loop

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=Acta5Schiz, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =Acta5Schiz$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  a5Schiz<-rbind(a5Schiz, Neuro2)
}

a5Schiz<-rownames_to_column(a5Schiz, var="NeuroType")
a5Schiz<-subset(a5Schiz, startsWith(a5Schiz$NeuroType, "Neuro"))
a5Schiz$Time<-"+5"
a5Schiz<-select(a5Schiz, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2)

####Run for a5 bipolar####

#Select only those in the imputed data who have bipolar
Acta5BP<-subset(BipolarAll, patid %in% ActA5$patid)

#Create an empty data set
a5BP<-Neuro1[NULL,]
#Wrap it in loop
for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=Acta5BP, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =Acta5BP$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  a5BP<-rbind(a5BP, Neuro2)
}

a5BP<-rownames_to_column(a5BP, var="NeuroType")
a5BP<-subset(a5BP, startsWith(a5BP$NeuroType, "Neuro"))
a5BP$Time<-"+5"
a5BP<-select(a5BP, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2)

####Run for A5  other psychoses####

#Select only those in the imputed data who have other psychoses
Acta5Other<-subset(OtherAll, patid %in% ActA5$patid)

#Create an empty data set
a5Other<-Neuro1[NULL,]
#Wrap it in loop
for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=Acta5Other, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =Acta5Other$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  a5Other<-rbind(a5Other, Neuro2)
}

a5Other<-rownames_to_column(a5Other, var="NeuroType")
a5Other<-subset(a5Other, startsWith(a5Other$NeuroType, "Neuro"))
a5Other$Time<-"+5"
a5Other<-select(a5Other, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

rm(CIFin, Sand, Neuro2, Neuro)


####Create one dataset for a5####
a5<-rbind(a5Schiz, a5BP, a5Other)

write.csv(a5, file = "/Neuro/Outputs/UnAdjDiagUnA5.csv")

#Create a dataframe with all our odds ratios in
Counts<-rbind(b5, b3, b1, b0, a1, a3, a5)

save(Counts, file= "/Neuro/Data/UnAdjDiagResults.Rdata")


