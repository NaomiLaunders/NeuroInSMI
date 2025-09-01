#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Find smoking status, drug and alcohol misuse status closest to five years before SMI diagnosis (or diagnosis of matched person for comparators)
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all.names = TRUE))

library(dplyr)
library(lubridate)

####Smoking####

load("/Neuro/Data/BMINeuroData.Rdata")

#Bring in file from previous data cleaning with all smoking events in 
load("/Interim Datasets/CaSmoke.Rdata")
Smoking<-subset(Smoking, patid %in% MatchedLK$patid)

length(which(is.na(Smoking$eventdate)))

#If obs date is less than 1900 or before birth or after index make it NA
length(which(Smoking$eventdate<as.Date("1900-01-01")|Smoking$eventdate>as.Date("2019-06-01")))
Smoking$eventdate[Smoking$eventdate<as.Date("1900-01-01")|Smoking$eventdate>as.Date("2019-06-01")]<-Smoking$sysdate[Smoking$eventdate<as.Date("1900-01-01")|Smoking$eventdate>as.Date("2019-06-01")]

Smoking<-subset(Smoking, eventdate>=as.Date("1900-01-01")&Smoking$eventdate<as.Date("2019-06-01"))

YoB<-select(MatchedLK, patid, yob)
Smoking<-merge(x=Smoking, y=YoB, by="patid", all.x=TRUE, all.y=FALSE)
table(Smoking$yob[Smoking$eventdate<Smoking$yob])
length(which(Smoking$eventdate<Smoking$yob))
Smoking$eventdate[Smoking$eventdate<Smoking$yob]<-Smoking$sysdate[Smoking$eventdate<Smoking$yob]
length(which(Smoking$eventdate<Smoking$yob))
Smoking<-subset(Smoking, Smoking$eventdate>=Smoking$yob)

Index<-select(MatchedLK, patid, case_index)
Smoking<-merge(x=Smoking, y=Index, by="patid", all.x=TRUE, all.y=FALSE)
length(which(Smoking$eventdate>Smoking$case_index))

#Remove all the NAs because only want those with a valid date
length(which(is.na(Smoking$eventdate)))

#Create eversmoke
EverSmoke<-subset(Smoking, currsmoke!=4)

#Take closest to 5 years prior to SMI diagnosis
#Take most recent
Smoking<-Smoking%>%
  mutate(Time5yr=case_index - 1826.25)%>%
  mutate(SmokeDistance=as.numeric(eventdate - Time5yr)/365.25)%>%
  subset(SmokeDistance<=5)%>%
  mutate(SmokeDistancePos=case_when(SmokeDistance<0 ~ SmokeDistance*-1,
                                  TRUE ~ SmokeDistance))%>%
  mutate(Priority=case_when(SmokeDistance<=0  ~1,
                            SmokeDistance>0~ 2,
                            TRUE ~ 3))%>%
  group_by(patid)%>%
  mutate(HighPriority=min(Priority))%>%
  filter(HighPriority==Priority)%>%
  mutate(ClosestSmoke=min(SmokeDistancePos))%>%
  filter(ClosestSmoke==SmokeDistancePos)%>%
  mutate(currsmoke=as.numeric(currsmoke))%>%
  mutate(Mincurrsmoke=min(currsmoke))%>%
  filter(Mincurrsmoke==currsmoke)

Smoking<-select(Smoking, patid, currsmoke, eventdate, SmokeDistancePos, SmokeDistance)
Smoking<-distinct(Smoking)


#Non smokers with a history of smoking
DateOfSmoke<-select(Smoking, patid, Closeto5=eventdate)
EverSmoke<-merge(x=EverSmoke, y=DateOfSmoke, by="patid", all.x=TRUE, all.y=FALSE)

EverSmoke<-subset(EverSmoke, eventdate<=Closeto5)
EverSmoke$EverSmoke<-1
EverSmoke<-select(EverSmoke, patid, EverSmoke)
EverSmoke<-distinct(EverSmoke)

Smoking<-merge(x=Smoking, y=EverSmoke, by="patid", all.x=TRUE, all.y=FALSE)

Smoking$currsmoke[Smoking$EverSmoke==1&Smoking$currsmoke==4]<-3

#Recode as three categories
Smoking$Smoke5[Smoking$currsmoke==1|Smoking$currsmoke==2]<-1
Smoking$Smoke5[Smoking$currsmoke==3]<-2
Smoking$Smoke5[Smoking$currsmoke==4]<-3

#Check distribution
Smoking<-Smoking %>% mutate(new_bin = cut(SmokeDistance, breaks=seq(-30, 30, 1)))

ggplot(Smoking, aes(new_bin))+
  geom_bar(color="black", fill="white")+
  theme_classic() 

#Merge to MatchedLK
Smoking<-select(Smoking, patid, Smoke5)
MatchedLK<-merge(x=MatchedLK, y=Smoking, by="patid", all.x=TRUE, all.y=FALSE)

####Drugs and alcohol####
load("/Interim Datasets/CaPatElixAllAurum.Rdata")
load("/Interim Datasets/CaPatElixAllGold.Rdata")

PatElixAll<-rename(PatElixAll, eventdate=obsdate, sysdate=enterdate, medcode=medcodeid)
PatElixGAll<-rename(PatElixGAll, medcode=V1, term=readterm)
DrugsAlc<-rbind(PatElixAll, PatElixGAll)
DrugsAlc<-subset(DrugsAlc, patid %in% MatchedLK$patid)
DrugsAlc<-subset(DrugsAlc, cat=="Alcohol abuse"|cat=="Drug abuse")

#If obs date is less than 1900 or before birth or after index make it NA
DrugsAlc$eventdate<-as.Date(DrugsAlc$eventdate, format="%d/%m/%Y")
DrugsAlc$sysdate<-as.Date(DrugsAlc$sysdate, format="%d/%m/%Y")
length(which(DrugsAlc$eventdate<as.Date("1900-01-01")|DrugsAlc$eventdate>as.Date("2019-06-01")))
DrugsAlc$eventdate[DrugsAlc$eventdate<as.Date("1900-01-01")|DrugsAlc$eventdate>as.Date("2019-06-01")]<-DrugsAlc$sysdate[DrugsAlc$eventdate<as.Date("1900-01-01")|Smoking$eventdate>as.Date("2019-06-01")]

DrugsAlc<-subset(DrugsAlc, eventdate>=as.Date("1900-01-01")&DrugsAlc$eventdate<as.Date("2019-06-01"))

YoB<-select(MatchedLK, patid, yob)
DrugsAlc<-merge(x=DrugsAlc, y=YoB, by="patid", all.x=TRUE, all.y=FALSE)
table(DrugsAlc$yob[DrugsAlc$eventdate<DrugsAlc$yob])
length(which(DrugsAlc$eventdate<DrugsAlc$yob))

DrugsAlc$eventdate[DrugsAlc$eventdate<DrugsAlc$yob]<-DrugsAlc$sysdate[DrugsAlc$eventdate<DrugsAlc$yob]
length(which(DrugsAlc$eventdate<DrugsAlc$yob))
DrugsAlc<-subset(DrugsAlc, DrugsAlc$eventdate>=DrugsAlc$yob)

Index<-select(MatchedLK, patid, case_index)
DrugsAlc<-merge(x=DrugsAlc, y=Index, by="patid", all.x=TRUE, all.y=FALSE)
length(which(DrugsAlc$eventdate>DrugsAlc$case_index))

#Remove all the NAs because only want those with a valid date
length(which(is.na(DrugsAlc$eventdate)))

####Drugs####
Drugs<-subset(DrugsAlc, cat=="Drug abuse")

#Take closest to 5 years prior to SMI diagnosis
#Take most recent
Drugs<-Drugs%>%
  mutate(Time5yr=case_index - 1826.25)%>%
  mutate(DrugsDistance=as.numeric(eventdate - Time5yr)/365.25)%>%
  subset(DrugsDistance<=5)%>%
  mutate(DrugsDistancePos=case_when(DrugsDistance<0 ~ DrugsDistance*-1,
                                    TRUE ~ DrugsDistance))%>%
  mutate(Priority=case_when(DrugsDistance<=0  ~1,
                            DrugsDistance>0 ~ 2,
                            TRUE ~ 3))%>%
  group_by(patid)%>%
  mutate(HighPriority=min(Priority))%>%
  filter(HighPriority==Priority)%>%
  mutate(ClosestDrug=min(DrugsDistancePos))%>%
  filter(ClosestDrug==DrugsDistancePos)

length(unique(Drugs$patid))

#If multiple find unique
Drugs<-select(Drugs, patid, eventdate, DrugsDistancePos, DrugsDistance)
Drugs<-distinct(Drugs)
Drugs$Drugs5<-1

#Check distribution
Drugs<-Drugs %>% mutate(new_bin = cut(DrugsDistance, breaks=seq(-30, 30, 1)))

ggplot(Drugs, aes(new_bin))+
  geom_bar(color="black", fill="white")+
  theme_classic() 

#Merge to MatchedLK
Drugs<-select(Drugs, patid, Drugs5)
MatchedLK<-merge(x=MatchedLK, y=Drugs, by="patid", all.x=TRUE, all.y=FALSE)

####Alcohol####
Alc<-subset(DrugsAlc, cat=="Alcohol abuse")

#Take closest to 5 years prior to SMI diagnosis
#Take most recent
Alc<-Alc%>%
  mutate(Time5yr=case_index - 1826.25)%>%
  mutate(AlcDistance=as.numeric(eventdate - Time5yr)/365.25)%>%
  subset(AlcDistance<=5)%>%
  mutate(AlcDistancePos=case_when(AlcDistance<0 ~ AlcDistance*-1,
                                    TRUE ~ AlcDistance))%>%
  mutate(Priority=case_when(AlcDistance<=0 ~1,
                            AlcDistance>0 ~ 2,
                            TRUE ~ 4))%>%
  group_by(patid)%>%
  mutate(HighPriority=min(Priority))%>%
  filter(HighPriority==Priority)%>%
  mutate(ClosestAlc=min(AlcDistancePos))%>%
  filter(ClosestAlc==AlcDistancePos)

length(unique(Alc$patid))

#If multiple find unique
Alc<-select(Alc, patid, eventdate, AlcDistancePos, AlcDistance)
Alc<-distinct(Alc)
Alc$Alc5<-1

#Check distribution
Alc<-Alc %>% mutate(new_bin = cut(AlcDistance, breaks=seq(-30, 30, 1)))

ggplot(Alc, aes(new_bin))+
  geom_bar(color="black", fill="white")+
  theme_classic() 

#Merge to MatchedLK
Alc<-select(Alc, patid, Alc5)
MatchedLK<-merge(x=MatchedLK, y=Alc, by="patid", all.x=TRUE, all.y=FALSE)

MatchedLK$Smoke5[is.na(MatchedLK$Smoke5)]<-3
MatchedLK$Drugs5[is.na(MatchedLK$Drugs5)]<-0
MatchedLK$Alc5[is.na(MatchedLK$Alc5)]<-0

save(MatchedLK, file="/Neuro/Data/RisksNeuroData.Rdata")
