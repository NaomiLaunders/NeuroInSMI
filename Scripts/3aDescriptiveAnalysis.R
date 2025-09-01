#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Descriptive analysis
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

load("/Neuro/Data/FinalNeuroData.Rdata")

MatchedLK$ethnicity<-as.factor(MatchedLK$ethnicity)
MatchedLK$BMI5<-as.factor(MatchedLK$BMI5)
MatchedLK$Smoke5<-as.factor(MatchedLK$Smoke5)
MatchedLK$Drugs5<-as.factor(MatchedLK$Drugs5)
MatchedLK$Alc5<-as.factor(MatchedLK$Alc5)
MatchedLK$region<-as.factor(MatchedLK$region)

Vars<-vars_select(names(MatchedLK), starts_with('Neuro_', ignore.case = TRUE))

MatchedLK[,Vars]<-lapply(MatchedLK[,Vars], as.factor)

MatchedLK$DiedInFU<-MatchedLK$died
MatchedLK$DiedInFU[MatchedLK$deathdate>(MatchedLK$case_index+1826.25)]<-0
table(MatchedLK$DiedInFU, MatchedLK$died)
                   
#Table 1 to include: n per group (no SMI, SMI, SMI subsets), age at SMI diagnosis, sex, ethnicity, BMI category, smoking status, drug misuse, alcohol misuse, died during follow up, age at death

MyVars<-c("AgeAtDiag", "gender", "ethnicity", "region", "BMI5", "Smoke5", "Drugs5", "Alc5", "DiedInFU", "AgeAtDeath", "FU", "Baseline")
Table1<-CreateTableOne(vars=MyVars, strata="SMI", data=MatchedLK, includeNA = TRUE)
print(Table1)

Table1Exp <- print(Table1, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE,  printToggle = FALSE, catDigits=2) 
write.csv(Table1Exp, file = "/Neuro/Outputs/TableOne.csv")

Table1b<-CreateTableOne(vars=MyVars, strata="diagn", data=MatchedLK, includeNA = TRUE)
print(Table1b)

Table1Expb <- print(Table1b, nonnormal=TRUE, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, catDigits=2)

write.csv(Table1Expb, file = "/Neuro/Outputs/TableOneb.csv")

MatchedLK<-select(MatchedLK, -DiedInFU)

####Table 2: Baseline and ever####
comorbs<-vars_select(names(MatchedLK), (starts_with('Neuro_', ignore.case = TRUE)&(!ends_with('A5'))&(!ends_with('After'))
                                         &(!ends_with('A4'))&(!ends_with('A3'))&(!ends_with('A2'))&(!ends_with('A1'))
                                         &(!ends_with('B5'))&(!ends_with('B3'))))

table2<-CreateTableOne(vars=comorbs, factorVars=comorbs, strata="SMI", data=MatchedLK, includeNA = TRUE)
Table2Exp <- print(table2, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, catDigits=2)
write.csv(Table2Exp, file = "/Neuro/Outputs/TableTwoNum.csv")

#Then for follow up we need separate tables because need to limit to those who were active in that period
####Table 2: A1####
A1<-subset(MatchedLK, Act_a1==1)
comorbsA1<-vars_select(names(MatchedLK), (starts_with('Neuro_', ignore.case = TRUE)&(ends_with('A1'))))
table2A1<-CreateTableOne(vars=comorbsA1, factorVars=comorbsA1, strata="SMI", data=A1, includeNA = TRUE)
Table2A1Exp <- print(table2A1, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, catDigits=2)
write.csv(Table2A1Exp, file = "/Neuro/Outputs/TableTwoA1.csv")

####Table 2: A3####
A3<-subset(MatchedLK, Act_a3==1)
comorbsA3<-vars_select(names(MatchedLK), (starts_with('Neuro_', ignore.case = TRUE)&(ends_with('A3'))))
table2A3<-CreateTableOne(vars=comorbsA3, factorVars=comorbsA3, strata="SMI", data=A3, includeNA = TRUE)
Table2A3Exp <- print(table2A3, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, catDigits=2)
write.csv(Table2A3Exp, file = "/Neuro/Outputs/TableTwoA3.csv")

####Table 2: A5####
A5<-subset(MatchedLK, Act_a5==1)
comorbsA5<-vars_select(names(MatchedLK), (starts_with('Neuro_', ignore.case = TRUE)&(ends_with('A5'))))
table2A5<-CreateTableOne(vars=comorbsA5, factorVars=comorbsA5, strata="SMI", data=A5, includeNA = TRUE)
Table2A5Exp <- print(table2A5, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, catDigits=2)
write.csv(Table2A5Exp, file = "/Neuro/Outputs/TableTwoA5.csv")

#Missing date
comorbsMissing<-vars_select(names(MatchedLK), (starts_with('Miss', ignore.case = TRUE)))
table2Miss<-CreateTableOne(vars=comorbsMissing, factorVars=comorbsMissing, strata="SMI", data=MatchedLK, includeNA = TRUE)
Table2MissExp <- print(table2Miss, quote = FALSE, noSpaces = TRUE, printToggle = FALSE, catDigits=2)
write.csv(Table2MissExp, file = "/Neuro/Outputs/TableTwoMiss.csv")

####Multimorbidity count####
#We could use this to look at how many groups they fall in, but I guess there wouldnt be many that have multiple neuro conditions?
MorbPrev<-select(MatchedLK, patid, SMI, count_prior, count_b5, count_b4, count_b3, count_b2, count_b1, count_a1, count_a2, count_a3, count_a4, count_a5, Act_a1, Act_a2, Act_a3, Act_a4, Act_a5)

MorbPrevbP<-MorbPrev%>%
  group_by(SMI)%>%
  mutate(pop=n())%>%
  group_by(SMI, count_prior, pop)%>%
  summarise(Sum_prior=n())%>%
  pivot_wider(names_from = count_prior, values_from=Sum_prior)%>%
  mutate(Time="prior")


MorbPrevb5<-MorbPrev%>%
  group_by(SMI)%>%
  mutate(pop=n())%>%
  group_by(SMI, count_b5, pop)%>%
  summarise(Sum_b5=n())%>%
  pivot_wider(names_from = count_b5, values_from=Sum_b5)%>%
  mutate(Time="b5")
  
MorbPrevb4<-MorbPrev%>%
  group_by(SMI)%>%
  mutate(pop=n())%>%
  group_by(SMI, count_b4 ,pop)%>%
  summarise(Sum_b4=n())%>%
  pivot_wider(names_from = count_b4, values_from=Sum_b4)%>%
  mutate(Time="b4")

MorbPrevb3<-MorbPrev%>%
  group_by(SMI)%>%
  mutate(pop=n())%>%
  group_by(SMI, count_b3 ,pop)%>%
  summarise(Sum_b3=n())%>%
  pivot_wider(names_from = count_b3, values_from=Sum_b3)%>%
  mutate(Time="b3")

MorbPrevb2<-MorbPrev%>%
  group_by(SMI)%>%
  mutate(pop=n())%>%
  group_by(SMI, count_b2 ,pop)%>%
  summarise(Sum_b2=n())%>%
  pivot_wider(names_from = count_b2, values_from=Sum_b2)%>%
  mutate(Time="b2")

MorbPrevb1<-MorbPrev%>%
  group_by(SMI)%>%
  mutate(pop=n())%>%
  group_by(SMI, count_b1 ,pop)%>%
  summarise(Sum_b1=n())%>%
  pivot_wider(names_from = count_b1, values_from=Sum_b1)%>%
  mutate(Time="b1")

MorbPreva1<-MorbPrev%>%
  group_by(SMI)%>%
  subset(Act_a1==1)%>%
  mutate(pop=n())%>%
  group_by(SMI, count_a1 ,pop)%>%
  summarise(Sum_a1=n())%>%
  pivot_wider(names_from = count_a1, values_from=Sum_a1)%>%
  mutate(Time="a1")

MorbPreva2<-MorbPrev%>%
  group_by(SMI)%>%
  subset(Act_a2==1)%>%
  mutate(pop=n())%>%
  group_by(SMI, count_a2 ,pop)%>%
  summarise(Sum_a2=n())%>%
  pivot_wider(names_from = count_a2, values_from=Sum_a2)%>%
  mutate(Time="a2")

MorbPreva3<-MorbPrev%>%
  group_by(SMI)%>%
  subset(Act_a3==1)%>%
  mutate(pop=n())%>%
  group_by(SMI, count_a3 ,pop)%>%
  summarise(Sum_a3=n())%>%
  pivot_wider(names_from = count_a3, values_from=Sum_a3)%>%
  mutate(Time="a3")

MorbPreva4<-MorbPrev%>%
  group_by(SMI)%>%
  subset(Act_a4==1)%>%
  mutate(pop=n())%>%
  group_by(SMI, count_a4 ,pop)%>%
  summarise(Sum_a4=n())%>%
  pivot_wider(names_from = count_a4, values_from=Sum_a4)%>%
  mutate(Time="a4")

MorbPreva5<-MorbPrev%>%
  group_by(SMI)%>%
  subset(Act_a5==1)%>%
  mutate(pop=n())%>%
  group_by(SMI, count_a5 ,pop)%>%
  summarise(Sum_a5=n())%>%
  pivot_wider(names_from = count_a5, values_from=Sum_a5)%>%
  mutate(Time="a5")

MorbPrevFin<-rbind(MorbPrevbP, MorbPrevb5,MorbPrevb4, MorbPrevb3, MorbPrevb2, MorbPrevb1, MorbPreva1, MorbPreva2, MorbPreva3, MorbPreva4,MorbPreva5 )
MorbPrevFin[,][is.na(MorbPrevFin[,])]<-0
MorbPrevFin<-select(MorbPrevFin, 9, 1:8, 10:12)

MorbPrevOver5<-MorbPrevFin
MorbPrevOver5$Over5<-MorbPrevOver5$'11'+MorbPrevOver5$'12'+MorbPrevOver5$'14'
MorbPrevOver5<-select(MorbPrevOver5, -c(10:12))

MorbPrevPerc<-MorbPrevOver5%>%
  group_by(Time, SMI)%>%
  mutate(Count0=`0`/pop*100,Count1=`1`/pop*100,Count2=`2`/pop*100,Count3=`3`/pop*100,Count4=`4`/pop*100,Count5=`5`/pop*100,CountOver5=`Over5`/pop*100)
MorbPrevPerc<-select(MorbPrevPerc, 1:3, 11:17)

MorbPrevPerc<-pivot_longer(MorbPrevPerc, cols=4:10, names_to = "Count", values_to = "Percent")
MorbPrevPerc$NewTime[MorbPrevPerc$Time=="b1"]<-"Diagnosis"
MorbPrevPerc$NewTime[MorbPrevPerc$Time=="b2"]<-"-1"
MorbPrevPerc$NewTime[MorbPrevPerc$Time=="b3"]<-"-2"
MorbPrevPerc$NewTime[MorbPrevPerc$Time=="b4"]<-"-3"
MorbPrevPerc$NewTime[MorbPrevPerc$Time=="b5"]<-"-4"
MorbPrevPerc$NewTime[MorbPrevPerc$Time=="prior"]<-"-5"
MorbPrevPerc$NewTime[MorbPrevPerc$Time=="a1"]<-"+1"
MorbPrevPerc$NewTime[MorbPrevPerc$Time=="a2"]<-"+2"
MorbPrevPerc$NewTime[MorbPrevPerc$Time=="a3"]<-"+3"
MorbPrevPerc$NewTime[MorbPrevPerc$Time=="a4"]<-"+4"
MorbPrevPerc$NewTime[MorbPrevPerc$Time=="a5"]<-"+5"
MorbPrevPerc$NewTime<-factor(MorbPrevPerc$NewTime, levels=c("-5", "-4", "-3", "-2", "-1", "Diagnosis", "+1", "+2", "+3", "+4", "+5"))

MorbPrevPerc$NewCount[MorbPrevPerc$Count=="Count0"]<-"0"
MorbPrevPerc$NewCount[MorbPrevPerc$Count=="Count1"]<-"1"
MorbPrevPerc$NewCount[MorbPrevPerc$Count=="Count2"]<-"2"
MorbPrevPerc$NewCount[MorbPrevPerc$Count=="Count3"]<-"3"
MorbPrevPerc$NewCount[MorbPrevPerc$Count=="Count4"]<-"4"
MorbPrevPerc$NewCount[MorbPrevPerc$Count=="Count5"]<-"5"
MorbPrevPerc$NewCount[MorbPrevPerc$Count=="CountOver5"]<-"Over 5"

MorbPrevPerc$NewCount<-factor(MorbPrevPerc$NewCount, levels=c("0", "1", "2", "3", "4", "5", "Over 5"))

MorbPrevSMI<-subset(MorbPrevPerc, SMI==1)
MorbPrevNoSMI<-subset(MorbPrevPerc, SMI==0)

SMI<-ggplot(MorbPrevSMI, aes(x=NewTime, y=Percent, fill=NewCount))+
  geom_bar(stat="identity", position=position_stack(reverse = TRUE))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 18))+
  labs( x="Time*", y = "Percentage of patients", title="SMI", fill="Number of LTCs")+
  scale_fill_brewer(palette = "Set2")+
  theme(plot.title = element_text(size = 20))+
  theme(legend.position=c(0.2, 0.2), legend.title = element_text(colour = "black", size = 18), legend.text = element_text(colour = "black", size = 18), axis.title.x=element_text(size=18),axis.text.y=element_text(size=18, colour="black"), axis.title.y=element_text(size=18), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))
  
NoSMI<-ggplot(MorbPrevNoSMI, aes(x=NewTime, y=Percent, fill=NewCount))+
  geom_bar(stat="identity", position=position_stack(reverse = TRUE))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 18))+
  labs( x="Time*", y = "Percentage of patients", title="No SMI", fill="Number of LTCs")+
  scale_fill_brewer(palette = "Set2")+
  theme(plot.title = element_text(size = 20))+
  theme(legend.position="none", legend.title = element_text(colour = "black", size = 18), legend.text = element_text(colour = "black", size = 18), axis.title.x=element_text(size=18),axis.text.y=element_text(size=18, colour="black"), axis.title.y=element_text(size=18), axis.line=element_line(color="black", size = 0.8), axis.ticks=element_line(color="black", size = 0.8))

grid.arrange(SMI, NoSMI, nrow=1)
         
write.csv(MorbPrevFin, file = "/Neuro/Outputs/Multimorb.csv")


