#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Minimally adjusted results SMI vs no SMI adjusting for AGE, SEX, ETHNICITY, REGION, CALENDAR YEAR
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

####Minimally adjusted results SMI vs no SMI - AGE, SEX, ETHNICITY, REGION, CALENDAR YEAR####

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

load("/Neuro/Data/ImpAnalysis.Rdata")

#Change to list of datasets
datlist<-mids2datlist(ImputedData)

#### Minimally Adjusted, SMI, in -5####
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

#We can put this into a function to save time 

#These are all minimally adjusted now. I looked at SMI vs no SMI first, and then by diagnosis.

#### Minimally Adjusted, SMI, in -5####
#Create an empty data set
b5<-Neuro1[NULL,]
#Wrap it in loop

Neuro<-vars_select(names(ImputedData$data), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('Prior')))

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(datlist, FUN=function(data){
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
  b5<-rbind(b5, Neuro2)
}

b5<-rownames_to_column(b5, var="NeuroType")
b5<-subset(b5, startsWith(b5$NeuroType, "Neuro"))
b5$Time<-"-5"
b5<-select(b5, NeuroType, OR, Lower, Upper, Time)

write.csv(b5, file = "/Neuro/Outputs/MinAdjSMIB5.csv")
rm(betas, Neuro2, vars, Neuro)

####Minimally Adjusted, SMI, in -3####
b3<-Neuro1[NULL,]

Neuro<-vars_select(names(ImputedData$data), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('B4')))

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(datlist, FUN=function(data){
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
  b3<-rbind(b3, Neuro2)
}

b3<-rownames_to_column(b3, var="NeuroType")
b3<-subset(b3, startsWith(b3$NeuroType, "Neuro"))
b3$Time<-"-3"
b3<-select(b3, NeuroType, OR, Lower, Upper, Time)

write.csv(b3, file = "/Neuro/Outputs/MinAdjSMIB3.csv")
rm(betas, Neuro2, vars, Neuro)

####Minimally Adjusted, SMI, in -1####
b1<-Neuro1[NULL,]

Neuro<-vars_select(names(ImputedData$data), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('B2')))

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(datlist, FUN=function(data){
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
  b1<-rbind(b1, Neuro2)
}

b1<-rownames_to_column(b1, var="NeuroType")
b1<-subset(b1, startsWith(b1$NeuroType, "Neuro"))
b1$Time<-"-1"
b1<-select(b1, NeuroType, OR, Lower, Upper, Time)

write.csv(b1, file = "/Neuro/Outputs/MinAdjSMIB1.csv")
rm(betas, Neuro2, vars, Neuro)


####Minimally Adjusted, SMI, at index - NOTE renamed index b0 as b1 used above for -1 year####
b0<-Neuro1[NULL,]

Neuro<-vars_select(names(ImputedData$data), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('B1')))

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(datlist, FUN=function(data){
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
  b0<-rbind(b0, Neuro2)
}

b0<-rownames_to_column(b0, var="NeuroType")
b0<-subset(b0, startsWith(b0$NeuroType, "Neuro"))
b0$Time<-"0"
b0<-select(b0, NeuroType, OR, Lower, Upper, Time)

write.csv(b0, file = "/Neuro/Outputs/MinAdjSMIB0index.csv")
rm(betas, Neuro2, vars, Neuro)


####For those after outcome we need to limit to those who are active#####
#Work out who is active from the descriptive data
load("/Neuro/Data/FinalNeuroData.Rdata")

####Min Adj SMI: A1####
ActA1<-subset(MatchedLK, Act_a1==1)
#Select only those in the imputed data who appear in the descriptive data as active
ImpA1<-subset_datlist(datlist, index=1:5,
                      subset=datlist[[1]]$patid %in% ActA1$patid)

a1<-Neuro1[NULL,]

Neuro<-vars_select(names(ImputedData$data), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('A1')))

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpA1, FUN=function(data){
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
  a1<-rbind(a1, Neuro2)
}

a1<-rownames_to_column(a1, var="NeuroType")
a1<-subset(a1, startsWith(a1$NeuroType, "Neuro"))
a1$Time<-"+1"
a1<-select(a1, NeuroType, OR, Lower, Upper, Time)

write.csv(a1, file = "/Neuro/Outputs/MinAdjSMIA1.csv")
rm(betas, Neuro2, vars, Neuro)


####Minimally Adj SMI: A3####
ActA3<-subset(MatchedLK, Act_a3==1)
#Select only those in the imputed data who appear in the descriptive data as active
ImpA3<-subset_datlist(datlist, index=1:5,
                      subset=datlist[[1]]$patid %in% ActA3$patid)

a3<-Neuro1[NULL,]

Neuro<-vars_select(names(ImputedData$data), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('A3')))

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpA3, FUN=function(data){
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
  a3<-rbind(a3, Neuro2)
}

a3<-rownames_to_column(a3, var="NeuroType")
a3<-subset(a3, startsWith(a3$NeuroType, "Neuro"))
a3$Time<-"+3"
a3<-select(a3, NeuroType, OR, Lower, Upper, Time)

write.csv(a3, file = "/Neuro/Outputs/MinAdjSMIA3.csv")
rm(betas, Neuro2, vars, Neuro)


####Minimally Adj SMI: A5####
ActA5<-subset(MatchedLK, Act_a5==1)
#Select only those in the imputed data who appear in the descriptive data as active
ImpA5<-subset_datlist(datlist, index=1:5,
                      subset=datlist[[1]]$patid %in% ActA5$patid)

a5<-Neuro1[NULL,]

Neuro<-vars_select(names(ImputedData$data), starts_with('Neuro_', ignore.case = TRUE)&(ends_with('A5')))

for (i in (1:length(Neuro))) {
  Neuro2 <- lapply(ImpA5, FUN=function(data){
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
  a5<-rbind(a5, Neuro2)
}

a5<-rownames_to_column(a5, var="NeuroType")
a5<-subset(a5, startsWith(a5$NeuroType, "Neuro"))
a5$Time<-"+5"
a5<-select(a5, NeuroType, OR, Lower, Upper, Time)

write.csv(a5, file = "/Neuro/Outputs/MinAdjSMIA5.csv")
rm(betas, Neuro2, vars, Neuro)

####Graphs - minimally adjusted SMI vs no SMI####

#Create a dataframe with all our odds ratios in
Counts<-rbind(b5, b3, b1, b0, a1, a3, a5)

save(Counts, file="/Neuro/Data/MinSMIResults.Rdata")

#Have a numeric variable for time as well as sometimes easier on graphs
Counts$ContTime[Counts$Time=="-5"]<- -5
Counts$ContTime[Counts$Time=="-3"]<- -3
Counts$ContTime[Counts$Time=="-1"]<- -1
Counts$ContTime[Counts$Time=="Diagnosis"]<- 0
Counts$ContTime[Counts$Time=="+1"]<- 1
Counts$ContTime[Counts$Time=="+3"]<- 3
Counts$ContTime[Counts$Time=="+5"]<- 5

Counts$SMI<-"SMI"
CountsPO<-Counts
CountsPO$Upper[CountsPO$Upper>20]<-20
Counts$Upper[Counts$Upper>7]<-7

####Graph for MS####
# CountsMS<-subset(Counts, NeuroType=="Neuro_MS")
# 
# Neuro<-ggplot()+
#   geom_line(data = CountsMS, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
#   geom_point(data = CountsMS, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
#   geom_ribbon(data = CountsMS, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 14))+
#   theme(axis.text.y = element_text(size = 14))+
#   theme(axis.title.x=element_text(size=16))+
#   theme(axis.title.y=element_text(size=16))+
#   theme(title=element_text(size=16))+
#   theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
#   labs( x="", y = "Odds ratio", title="Multiple sclerosis")+
#   scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
#   scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
# 
# ####Graph for Cerebrovascular disease####
# CountsCereb<-subset(Counts, NeuroType=="Neuro_Cereb")
# 
# Neuro<-ggplot()+
#   geom_line(data = CountsCereb, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
#   geom_point(data = CountsCereb, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
#   geom_ribbon(data = CountsCereb, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 14))+
#   theme(axis.text.y = element_text(size = 14))+
#   theme(axis.title.x=element_text(size=16))+
#   theme(axis.title.y=element_text(size=16))+
#   theme(title=element_text(size=16))+
#   theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
#   labs( x="", y = "Odds ratio", title="Cerebrovascular disease")+
#   scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
#   scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
# 
# ####Graph for Dementia####
# CountsDem<-subset(Counts, NeuroType=="Neuro_Dementia")
# 
# Neuro<-ggplot()+
#   geom_line(data = CountsDem, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
#   geom_point(data = CountsDem, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
#   geom_ribbon(data = CountsDem, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 14))+
#   theme(axis.text.y = element_text(size = 14))+
#   theme(axis.title.x=element_text(size=16))+
#   theme(axis.title.y=element_text(size=16))+
#   theme(title=element_text(size=16))+
#   theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
#   labs( x="", y = "Odds ratio", title="Dementia")+
#   scale_y_continuous(limits = c(0.5, 6), breaks=seq(0.5,6, by=0.5), trans = "log2")+   
#   scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
# 
# #may need to change breaks to 1 - test run on data first
# 
# ####Graph for Ataxia####
# CountsAtax<-subset(Counts, NeuroType=="Neuro_Ataxia")
# 
# Neuro<-ggplot()+
#   geom_line(data = CountsAtax, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
#   geom_point(data = CountsAtax, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
#   geom_ribbon(data = CountsAtax, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 14))+
#   theme(axis.text.y = element_text(size = 14))+
#   theme(axis.title.x=element_text(size=16))+
#   theme(axis.title.y=element_text(size=16))+
#   theme(title=element_text(size=16))+
#   theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
#   labs( x="", y = "Odds ratio", title="Ataxia")+
#   scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
#   scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
# 
# ####Graph for Epilepsy####
# CountsEpilepsy <-subset(Counts, NeuroType=="Neuro_Epilepsy")
# 
# Neuro<-ggplot()+
#   geom_line(data = CountsEpilepsy, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
#   geom_point(data = CountsEpilepsy, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
#   geom_ribbon(data = CountsEpilepsy, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 14))+
#   theme(axis.text.y = element_text(size = 14))+
#   theme(axis.title.x=element_text(size=16))+
#   theme(axis.title.y=element_text(size=16))+
#   theme(title=element_text(size=16))+
#   theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
#   labs( x="", y = "Odds ratio", title="Epilepsy")+
#   scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
#   scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
# 
# ####Graph for Parkinsons####
# CountsPD <-subset(Counts, NeuroType=="Neuro_Parkinsons")
# 
# Neuro<-ggplot()+
#   geom_line(data = CountsPD, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
#   geom_point(data = CountsPD, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
#   geom_ribbon(data = CountsPD, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 14))+
#   theme(axis.text.y = element_text(size = 14))+
#   theme(axis.title.x=element_text(size=16))+
#   theme(axis.title.y=element_text(size=16))+
#   theme(title=element_text(size=16))+
#   theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
#   labs( x="", y = "Odds ratio", title="Parkinson's disease")+
#   scale_y_continuous(limits = c(0.5, 5), breaks=seq(0.5,5, by=0.5), trans = "log2")+   
#   scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
# 
# ####Graph for Paralysis####
# CountsParalysis <-subset(Counts, NeuroType=="Neuro_Paralysis")
# 
# Neuro<-ggplot()+
#   geom_line(data = CountsParalysis, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
#   geom_point(data = CountsParalysis, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
#   geom_ribbon(data = CountsParalysis, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 14))+
#   theme(axis.text.y = element_text(size = 14))+
#   theme(axis.title.x=element_text(size=16))+
#   theme(axis.title.y=element_text(size=16))+
#   theme(title=element_text(size=16))+
#   theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
#   labs( x="", y = "Odds ratio", title="Paralysis")+
#   scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
#   scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
# 
# ####Graph for Cerebral Palsy####
# CountsCP <-subset(Counts, NeuroType=="Neuro_CP")
# 
# Neuro<-ggplot()+
#   geom_line(data = CountsCP, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
#   geom_point(data = CountsCP, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
#   geom_ribbon(data = CountsCP, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 14))+
#   theme(axis.text.y = element_text(size = 14))+
#   theme(axis.title.x=element_text(size=16))+
#   theme(axis.title.y=element_text(size=16))+
#   theme(title=element_text(size=16))+
#   theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
#   labs( x="", y = "Odds ratio", title="Cerebral Palsy")+
#   scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
#   scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
# 
# ####Graph for CSF disorders####
# CountsCSF <-subset(Counts, NeuroType=="Neuro_CSF")
# 
# Neuro<-ggplot()+
#   geom_line(data = CountsCSF, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
#   geom_point(data = CountsCSF, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
#   geom_ribbon(data = CountsCSF, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 14))+
#   theme(axis.text.y = element_text(size = 14))+
#   theme(axis.title.x=element_text(size=16))+
#   theme(axis.title.y=element_text(size=16))+
#   theme(title=element_text(size=16))+
#   theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
#   labs( x="", y = "Odds ratio", title="CSF disorders")+
#   scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
#   scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
# 
# ####Graph for Spinal Cord disorders####
# CountsSpinal <-subset(Counts, NeuroType=="Neuro_Spinal")
# 
# Neuro<-ggplot()+
#   geom_line(data = CountsSpinal, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
#   geom_point(data = CountsSpinal, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
#   geom_ribbon(data = CountsSpinal, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 14))+
#   theme(axis.text.y = element_text(size = 14))+
#   theme(axis.title.x=element_text(size=16))+
#   theme(axis.title.y=element_text(size=16))+
#   theme(title=element_text(size=16))+
#   theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
#   labs( x="", y = "Odds ratio", title="Spinal Cord disorders")+
#   scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
#   scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

####Graph for Parkinsonism other disorders####
CountsParkinsonOther <-subset(CountsPO, startsWith(NeuroType, "Neuro_ParkinsonOther"))

ParkO<-ggplot()+
  geom_line(data = CountsParkinsonOther, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsParkinsonOther, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsParkinsonOther, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  scale_fill_manual(values = c("firebrick2"))+
  scale_color_manual(values = c("firebrick2"))+
  theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs(y = "Partially adjusted odds ratio", title="Parkinsonism other disorders", x="Time before and after index date (years)")+
  scale_y_continuous(limits = c(0, 20), breaks=c(0.5, 1, seq(2, 20, by=1)))+
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
ParkO
ggsave("/Neuro/MinAdjParkO.tiff", plot=ParkO, device="tiff", dpi=300)

####Graph for Disorders of nerve root, plexus or peripheral nerves####
# CountsPerip <-subset(Counts, NeuroType=="Neuro_Peripheral")
# 
# Neuro<-ggplot()+
#   geom_line(data = CountsPerip, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
#   geom_point(data = CountsPerip, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
#   geom_ribbon(data = CountsPerip, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 14))+
#   theme(axis.text.y = element_text(size = 14))+
#   theme(axis.title.x=element_text(size=16))+
#   theme(axis.title.y=element_text(size=16))+
#   theme(title=element_text(size=16))+
#   theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
#   labs( x="", y = "Odds ratio", title="Disorders of nerve root, plexus or peripheral nerves")+
#   scale_y_continuous(limits = c(0.5, 4), breaks=seq(0.5,4, by=0.5), trans = "log2")+   
#   scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

# ####Graph for -	Encephalopathy####
# CountsEnceph <-subset(Counts, NeuroType=="Neuro_Enceph")
# 
# Neuro<-ggplot()+
#   geom_line(data = CountsEnceph, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
#   geom_point(data = CountsEnceph, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
#   geom_ribbon(data = CountsEnceph, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.5)+
#   theme_classic()+
#   theme(axis.text.x = element_text(size = 14))+
#   theme(axis.text.y = element_text(size = 14))+
#   theme(axis.title.x=element_text(size=16))+
#   theme(axis.title.y=element_text(size=16))+
#   theme(title=element_text(size=16))+
#   theme(legend.position = c("none"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
#   labs( x="", y = "Odds ratio", title="Encephalopathy")+
#   scale_y_continuous(limits = c(0.5, 9), breaks=seq(0.5,9, by=1), trans = "log2")+   
#   scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))

####Graph for Movement disorders (other) ####
CountsMovement <-subset(Counts, startsWith(NeuroType, "Neuro_Movement"))

Move<-ggplot()+
  geom_line(data = CountsMovement, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsMovement, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsMovement, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  scale_fill_manual(values = c("firebrick2"))+
  scale_color_manual(values = c("firebrick2"))+
  theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs(y = "Partially adjusted odds ratio", title="Movement disorders, other", x="Time before and after index date (years)")+
  scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
Move
ggsave("/Neuro/MinAdjMove.tiff", plot=Move, device="tiff", dpi=300)


####Graph for Motor neurone disease ####
CountsMND <-subset(Counts, startsWith(NeuroType, "Neuro_MND"))

MND<-ggplot()+
  geom_line(data = CountsMND, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsMND, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsMND, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  scale_fill_manual(values = c("firebrick2"))+
  scale_color_manual(values = c("firebrick2"))+
  theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs(y = "Partially adjusted odds ratio", title="Motor Neurone Disease", x="Time before and after index date (years)")+
  scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
MND
ggsave("/Neuro/MinAdjMND.tiff", plot=MND, device="tiff", dpi=300)

####Graph for Autonomic Nervous System Disorders ####
CountsAuto <-subset(Counts, startsWith(NeuroType, "Neuro_Autonomic"))

Auto<-ggplot()+
  geom_line(data = CountsAuto, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=2)+
  geom_point(data = CountsAuto, aes(x=ContTime, y=OR, group=SMI, color=SMI), size=5)+
  geom_ribbon(data = CountsAuto, aes(x=ContTime, ymin = Lower, ymax = Upper, group=SMI, fill=SMI), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  scale_fill_manual(values = c("firebrick2"))+
  scale_color_manual(values = c("firebrick2"))+
  theme(legend.position = c("bottom"),legend.title = element_blank(), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs(y = "Partially adjusted odds ratio", title="Autonomic Nervous System Disorders", x="Time before and after index date (years)")+
  scale_y_continuous(limits = c(0, 7), breaks=c(0.5, 1, seq(2, 7, by=1)))+
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
Auto
ggsave("/Neuro/MinAdjAuto.tiff", plot=Auto, device="tiff", dpi=300)


####Graph for all neuro ####
Counts$Neuro[startsWith(Counts$NeuroType, "Neuro_MS")]<-"Multiple Sclerosis"
Counts$Neuro[startsWith(Counts$NeuroType, "Neuro_Cereb")]<-"Cerebrovascular disease"
Counts$Neuro[startsWith(Counts$NeuroType, "Neuro_Dementia")]<-"Dementia"
Counts$Neuro[startsWith(Counts$NeuroType, "Neuro_Ataxia")]<-"Ataxia"
Counts$Neuro[startsWith(Counts$NeuroType, "Neuro_Epilepsy")]<-"Epilepsy"
Counts$Neuro[startsWith(Counts$NeuroType, "Neuro_Parkinsons")]<-"Parkinson's disease"
Counts$Neuro[startsWith(Counts$NeuroType, "Neuro_Paralysis")]<-"Paralysis"
Counts$Neuro[startsWith(Counts$NeuroType, "Neuro_CP")]<-"Cerebral Palsy"
Counts$Neuro[startsWith(Counts$NeuroType, "Neuro_CSF")]<-"CSF disorders"
Counts$Neuro[startsWith(Counts$NeuroType, "Neuro_Spinal")]<-"Spinal cord disorders"
Counts$Neuro[startsWith(Counts$NeuroType, "Neuro_ParkinsonOther")]<-"Parkinsonism other"
Counts$Neuro[startsWith(Counts$NeuroType, "Neuro_Peripheral")]<-"Disorders of nerve root, plexus or peripheral nerves"
# Counts$Neuro[startsWith(Counts$NeuroType, "Neuro_Enceph")]<-"Encephalopathy"
Counts$Neuro[startsWith(Counts$NeuroType, "Neuro_Movement")]<-"Movement Disorders (other)"
Counts$Neuro[startsWith(Counts$NeuroType, "Neuro_MND")]<-"Motor Neurone Disease"
Counts$Neuro[startsWith(Counts$NeuroType, "Neuro_Autonomic")]<-"Autonomic Nervous System Disorders"

CountsTop7<-subset(Counts, Neuro=="Multiple Sclerosis"|Neuro=="Cerebrovascular disease"| Neuro=="Dementia"|Neuro=="Epilepsy"|Neuro=="Parkinson's disease"
                   |Neuro=="Paralysis"|Neuro=="Cerebral Palsy")

All<-ggplot()+
  geom_line(data = CountsTop7, aes(x=ContTime, y=OR, group=Neuro, color=Neuro), size=2)+
  geom_point(data = CountsTop7, aes(x=ContTime, y=OR, group=Neuro, color=Neuro), size=5)+
  geom_ribbon(data = CountsTop7, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Neuro, fill=Neuro), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.title=element_blank(), legend.position=c("bottom"), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs(y = "Partially adjusted odds ratio", title="", x="Time before and after index date (years)")+
  scale_y_continuous(limits = c(0.9, 6), breaks=seq(1,6, by=1))+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))
All

#Save for poster
ggsave("/Neuro/PartAdj7.tiff", plot=All, device="tiff", dpi=300)

CountsNext8<-subset(Counts, Neuro=="CSF disorders"|Neuro=="Disorders of nerve root, plexus or peripheral nerves"| 
                      Neuro=="Spinal cord disorders"|Neuro=="Ataxia"| Neuro=="Movement Disorders (other)"|
                      Neuro=="Motor Neurone Disease"|Neuro=="Autonomic Nervous System Disorders")

ggplot()+
  geom_line(data = CountsNext8, aes(x=ContTime, y=OR, group=Neuro, color=Neuro), size=2)+
  geom_point(data = CountsNext8, aes(x=ContTime, y=OR, group=Neuro, color=Neuro), size=5)+
  geom_ribbon(data = CountsNext8, aes(x=ContTime, ymin = Lower, ymax = Upper, group=Neuro, fill=Neuro), alpha=0.3)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 14))+
  theme(axis.text.y = element_text(size = 14))+
  theme(axis.title.x=element_text(size=16))+
  theme(axis.title.y=element_text(size=16))+
  theme(title=element_text(size=16))+
  theme(legend.position=c("bottom"), legend.text = element_text(colour = "black", size = 16), axis.text.x=element_text(size=16, colour="black"), axis.title.x=element_text(size=16),axis.text.y=element_text(size=16, colour="black"), axis.title.y=element_text(size=16), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))+
  labs( x="", y = "Odds ratio", title="All Neurological Disorders")+
  scale_y_continuous(limits = c(0.1, 10), breaks=seq(1,10, by=1))+   
  scale_x_continuous(limits = c(-5, 5), breaks=c(-5, -3, -1, 0, 1, 3, 5))