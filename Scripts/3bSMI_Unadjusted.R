#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Unadjusted results SMI vs no SMI 
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
library(broom)

load("/Neuro/Data/FinalNeuroData.Rdata")

#### Unadjusted, SMI, in -5####
#Do the first
Neuro1 <- glm(data=MatchedLK, formula=SMI ~ Neuro_MS_prior,
               family=binomial(link="logit") )

Sand<-vcovCL(Neuro1, cluster =MatchedLK$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
CIFin<-cbind(exp(coefci(Neuro1, vcov=Sand)))

Neuro1<-tidy(Neuro1)
Neuro1$OR<-exp(Neuro1$estimate)
Neuro1<-cbind(Neuro1, CIFin)

#### Unadjusted, SMI, in -5####
#Create an empty data set
b5<-Neuro1[NULL,]
#Wrap it in loop

Neuro<-vars_select(names(MatchedLK), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('Prior')))

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=MatchedLK, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )

  Sand<-vcovCL(Neuro2, cluster =MatchedLK$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  # conduct statistical inference
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  b5<-rbind(b5, Neuro2)
}

b5<-rownames_to_column(b5, var="NeuroType")
b5<-subset(b5, startsWith(b5$NeuroType, "Neuro"))
b5$Time<-"-5"
b5<-select(b5, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

write.csv(b5, file = "/Neuro/Outputs/UnAdjSMIB5.csv")
rm(CIFin, Sand, Neuro2, Neuro)

####Unadjusted, SMI, in -3####
b3<-Neuro1[NULL,]

Neuro<-vars_select(names(MatchedLK), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('B4')))

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=MatchedLK, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =MatchedLK$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  # conduct statistical inference
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  b3<-rbind(b3, Neuro2)
}

b3<-rownames_to_column(b3, var="NeuroType")
b3<-subset(b3, startsWith(b3$NeuroType, "Neuro"))
b3$Time<-"-3"
b3<-select(b3, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

write.csv(b3, file = "/Neuro/Outputs/UnAdjSMIB3.csv")
rm(CIFin, Sand, Neuro2, Neuro)

####Unadjusted, SMI, in -1####
b1<-Neuro1[NULL,]

Neuro<-vars_select(names(MatchedLK), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('B2')))

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=MatchedLK, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =MatchedLK$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  # conduct statistical inference
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  b1<-rbind(b1, Neuro2)
}

b1<-rownames_to_column(b1, var="NeuroType")
b1<-subset(b1, startsWith(b1$NeuroType, "Neuro"))
b1$Time<-"-1"
b1<-select(b1, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

write.csv(b1, file = "/Neuro/Outputs/UnAdjSMIB1.csv")
rm(CIFin, Sand, Neuro2, Neuro)

####Unadjusted, SMI, at index - NOTE renamed index b0 as b1 used above for -1 year####
b0<-Neuro1[NULL,]

Neuro<-vars_select(names(MatchedLK), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('B1')))

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=MatchedLK, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =MatchedLK$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  # conduct statistical inference
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  b0<-rbind(b0, Neuro2)
}

b0<-rownames_to_column(b0, var="NeuroType")
b0<-subset(b0, startsWith(b0$NeuroType, "Neuro"))
b0$Time<-"0"
b0<-select(b0, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

write.csv(b0, file = "/Neuro/Outputs/UnAdjSMIB0index.csv")
rm(CIFin, Sand, Neuro2, Neuro)

####For those after outcome we need to limit to those who are active#####

####Unadjusted SMI: A1####
ActA1<-subset(MatchedLK, Act_a1==1)
#Select only those in the imputed data who appear in the descriptive data as active

a1<-Neuro1[NULL,]

Neuro<-vars_select(names(MatchedLK), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('A1')))

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=ActA1, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =ActA1$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  # conduct statistical inference
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  a1<-rbind(a1, Neuro2)
}

a1<-rownames_to_column(a1, var="NeuroType")
a1<-subset(a1, startsWith(a1$NeuroType, "Neuro"))
a1$Time<-"+1"
a1<-select(a1, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

write.csv(a1, file = "/Neuro/Outputs/UnAdjSMIA1.csv")
rm(CIFin, Sand, Neuro2, Neuro)

####Unadjusted SMI: A3####

ActA3<-subset(MatchedLK, Act_a3==1)
#Select only those in the imputed data who appear in the descriptive data as active

a3<-Neuro1[NULL,]

Neuro<-vars_select(names(MatchedLK), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('A3')))

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=ActA3, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =ActA3$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  # conduct statistical inference
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  a3<-rbind(a3, Neuro2)
}

a3<-rownames_to_column(a3, var="NeuroType")
a3<-subset(a3, startsWith(a3$NeuroType, "Neuro"))
a3$Time<-"+3"
a3<-select(a3, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

write.csv(a3, file = "/Neuro/Outputs/UnAdjSMIA3.csv")
rm(CIFin, Sand, Neuro2, Neuro)

####Unadjusted SMI: A5####
ActA5<-subset(MatchedLK, Act_a5==1)
#Select only those in the imputed data who appear in the descriptive data as active

a5<-Neuro1[NULL,]

Neuro<-vars_select(names(MatchedLK), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('A5')))

for (i in (1:length(Neuro))) {
  Neuro2 <- glm(data=ActA5, as.formula( paste("SMI ~", Neuro[i])),
                family=binomial(link="logit") )
  
  Sand<-vcovCL(Neuro2, cluster =ActA5$pracid, type = NULL, sandwich = TRUE, fix = FALSE)
  CIFin<-cbind(exp(coefci(Neuro2, vcov=Sand)))
  
  # conduct statistical inference
  Neuro2<-tidy(Neuro2)
  Neuro2$OR<-exp(Neuro2$estimate)
  Neuro2<-cbind(Neuro2, CIFin)
  a5<-rbind(a5, Neuro2)
}

a5<-rownames_to_column(a5, var="NeuroType")
a5<-subset(a5, startsWith(a5$NeuroType, "Neuro"))
a5$Time<-"+5"
a5<-select(a5, NeuroType, OR, Lower=`2.5 %`, Upper=`97.5 %`, Time)

write.csv(a5, file = "/Neuro/Outputs/UnAdjSMIA5.csv")
rm(CIFin, Sand, Neuro2, Neuro)

#### can stop running code here ####

####Graphs - Unadjusted SMI vs no SMI####
#Create a dataframe with all our odds ratios in
Counts<-rbind(b5, b3, b1, b0, a1, a3, a5)

#Have a numeric variable for time as well as sometimes easier on graphs
Counts$ContTime[Counts$Time=="-5"]<- -5
Counts$ContTime[Counts$Time=="-3"]<- -3
Counts$ContTime[Counts$Time=="-1"]<- -1
Counts$ContTime[Counts$Time=="Diagnosis"]<- 0
Counts$ContTime[Counts$Time=="+1"]<- 1
Counts$ContTime[Counts$Time=="+3"]<- 3
Counts$ContTime[Counts$Time=="+5"]<- 5

####Graph for MS####
CountsMS<-subset(Counts, Neurotype=="Neuro_MS")

Neuro<-ggplot()+
  geom_line(data = CountsMS, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsMS, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsMS, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs( x="", y = "Odds ratio", title="Multiple sclerosis")+
  scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

####Graph for Cerebrovascular disease####
CountsCereb<-subset(Counts, Neurotype=="Neuro_Cereb")

Neuro<-ggplot()+
  geom_line(data = CountsCereb, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsCereb, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsCereb, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs( x="", y = "Odds ratio", title="Cerebrovascular disease")+
  scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

####Graph for Dementia####
CountsDem<-subset(Counts, Neurotype=="Neuro_Dementia")

Neuro<-ggplot()+
  geom_line(data = CountsDem, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsDem, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsDem, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs( x="", y = "Odds ratio", title="Dementia")+
  scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

####Graph for Ataxia####
CountsAtax<-subset(Counts, Neurotype=="Neuro_Ataxia")

Neuro<-ggplot()+
  geom_line(data = CountsAtax, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsAtax, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsAtax, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs( x="", y = "Odds ratio", title="Ataxia")+
  scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

####Graph for Epilepsy####
CountsEpilepsy <-subset(Counts, Neurotype=="Neuro_Epilepsy")

Neuro<-ggplot()+
  geom_line(data = CountsEpilepsy, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsEpilepsy, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsEpilepsy, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs( x="", y = "Odds ratio", title="Epilepsy")+
  scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

####Graph for Parkinsons####
CountsPD <-subset(Counts, Neurotype=="Neuro_Parkinsons")

Neuro<-ggplot()+
  geom_line(data = CountsPD, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsPD, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsPD, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs( x="", y = "Odds ratio", title="Parkinson's disease")+
  scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

####Graph for Paralysis####
CountsParalysis <-subset(Counts, Neurotype=="Neuro_Paralysis")

Neuro<-ggplot()+
  geom_line(data = CountsParalysis, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsParalysis, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsParalysis, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs( x="", y = "Odds ratio", title="Paralysis")+
  scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

####Graph for Cerebral Palsy####
CountsCP <-subset(Counts, Neurotype=="Neuro_CP")

Neuro<-ggplot()+
  geom_line(data = CountsCP, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsCP, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsCP, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs( x="", y = "Odds ratio", title="Cerebral Palsy")+
  scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

####Graph for CSF disorders####
CountsCSF <-subset(Counts, Neurotype=="Neuro_CSF")

Neuro<-ggplot()+
  geom_line(data = CountsCSF, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsCSF, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsCSF, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs( x="", y = "Odds ratio", title="CSF disorders")+
  scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

####Graph for Spinal Cord disorders####
CountsSpinal <-subset(Counts, Neurotype=="Neuro_Spinal")

Neuro<-ggplot()+
  geom_line(data = CountsSpinal, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsSpinal, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsSpinal, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs( x="", y = "Odds ratio", title="Spinal Cord disorders")+
  scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

####Graph for Parkinsonism other disorders####
CountsParkinsonOther <-subset(Counts, Neurotype=="Neuro_ParkinsonOther")

Neuro<-ggplot()+
  geom_line(data = CountsParkinsonOther, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsParkinsonOther, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsParkinsonOther, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs( x="", y = "Odds ratio", title="Parkinsonism other disorders")+
  scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

####Graph for Disorders of nerve root, plexus or peripheral nerves####
CountsPerip <-subset(Counts, Neurotype=="Neuro_Peripheral")

Neuro<-ggplot()+
  geom_line(data = CountsPerip, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsPerip, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsPerip, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs( x="", y = "Odds ratio", title="Disorders of nerve root, plexus or peripheral nerves")+
  scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

####Graph for -	Encephalopathy####
CountsEnceph <-subset(Counts, Neurotype=="Neuro_Enceph")

Neuro<-ggplot()+
  geom_line(data = CountsEnceph, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsEnceph, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsEnceph, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs( x="", y = "Odds ratio", title="Encephalopathy")+
  scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

####Graph for Movement disorders (other) ####
CountsMovement <-subset(Counts, Neurotype=="Neuro_Movement")

Neuro<-ggplot()+
  geom_line(data = CountsMovement, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsMovement, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsMovement, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs( x="", y = "Odds ratio", title="Movement disorders, other")+
  scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

####Graph for Motor neurone disease ####
CountsMND <-subset(Counts, Neurotype=="Neuro_MND")

Neuro<-ggplot()+
  geom_line(data = CountsMND, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsMND, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsMND, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs( x="", y = "Odds ratio", title="Motor Neurone Disease")+
  scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

####Graph for Autonomic Nervous System Disorders ####
CountsAuto <-subset(Counts, Neurotype=="Neuro_Autonomic")

Neuro<-ggplot()+
  geom_line(data = CountsAuto, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsAuto, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsAuto, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs( x="", y = "Odds ratio", title="Autonomic Nervous System Disorders")+
  scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

