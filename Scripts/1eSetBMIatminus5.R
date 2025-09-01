#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Find BMI closest to five years before SMI diagnosis (or diagnosis of matched person for comparators)
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all.names = TRUE))

library(dplyr)
library(lubridate)
library(ggplot2)

load("/Neuro/Data/FullNeuroData.Rdata")

####BMI listed value####

load("/Interim Datasets/LeiahPatBMIValueAurum.Rdata")

BMIAurumVal<-rename(PatEthnAll)
rm(PatEthnAll)

#Gold
load("/Interim Datasets/LeiahBMIValG.Rdata")
BMIGoldVal<-rename(BMIG)
rm(BMIG)
####BMI category####
#Aurum
load("/Interim Datasets/LeiahPatBMICatAurum.Rdata")
BMIAurumCat<-rename(PatEthnAll)
rm(PatEthnAll)

#Gold
load("/Interim Datasets/LeiahPatBMICatGold.Rdata")
BMIGoldCat<-rename(PatEthnGold)
rm(PatEthnGold)

####BMI calculated####
#Weight = aurum
load("/Interim Datasets/LeiahPatWeightAurum.Rdata")
AurumWeight<-rename(PatEthnAll)
rm(PatEthnAll)

#Weight - Gold
load("/Interim Datasets/LeiahPatWeightGold.Rdata")
GoldWeight<-rename(WeightG)
rm(WeightG)

####Height####
#Height
load("/Interim Datasets/CaPatHeightAurum.Rdata")
AurumHeight<-rename(PatEthnAll)
rm(PatEthnAll)

load("/Interim Datasets/CaPatHeightGold.Rdata")
GoldHeight<-rename(HeightG)
rm(HeightG)

####BMI####
Dates<-select(MatchedLK, patid, case_index, yob)
BMIValueA<-subset(BMIAurumVal, patid %in% Dates$patid)
BMIValueA<-merge(x=Dates, y=BMIValueA, by ="patid", all.y = FALSE, all.x=FALSE)
BMIValueA<-select(BMIValueA, -numunitid, -medcodeid, -Description)
names(BMIValueA)[names(BMIValueA)=="value"]<-"BMIValue"
names(BMIValueA)[names(BMIValueA)=="obsdate"]<-"BMIValueDate"
length(unique(BMIValueA$patid))

BMIValueG<-subset(BMIGoldVal, patid %in% Dates$patid)
BMIValueG<-merge(x=Dates, y=BMIValueG, by ="patid", all.y = FALSE, all.x=FALSE)
names(BMIValueG)[names(BMIValueG)=="BMI"]<-"BMIValue"
names(BMIValueG)[names(BMIValueG)=="BMIDate"]<-"BMIValueDate"
names(BMIValueG)[names(BMIValueG)=="sysdate"]<-"enterdate"
length(unique(BMIValueG$patid))
BMIValueG<-BMIValueG[, c(1:3,5, 6, 4) ]

BMIValue<-rbind(BMIValueA, BMIValueG)

#If an error, take sysdate
length(which(BMIValue$BMIValueDate<as.Date("1900-01-01")|BMIValue$BMIValueDate>as.Date("2019-06-01")))
BMIValue$BMIdate[BMIValue$BMIValueDate<as.Date("1900-01-01")|BMIValue$BMIValueDate>as.Date("2019-06-01")]<-BMIValue$enterdate[BMIValue$BMIdate<as.Date("1900-01-01")|BMIValue$BMIdate>as.Date("2019-06-01")]

BMIValue<-subset(BMIValue, BMIValueDate>=as.Date("1900-01-01")&BMIValueDate<as.Date("2019-06-01"))

length(which(BMIValue$BMIValueDate<BMIValue$yob))
BMIValue$BMIValueDate[BMIValue$BMIValueDate<BMIValue$yob]<-BMIValue$enterdate[BMIValue$BMIValueDate<BMIValue$yob]
length(which(BMIValue$BMIValueDate<BMIValue$yob))
BMIValue<-subset(BMIValue, BMIValueDate>=yob)

#Only take those which are on or before index and over 16.
BMIValue$BMIage<-as.numeric(format(BMIValue$BMIValueDate, "%Y"))-BMIValue$yob
summary(BMIValue$BMIage)
#Only take those as an adult
BMIValue<-subset(BMIValue, BMIage>=16)

#Take closest to 5 years prior to SMI diagnosis
#Take most recent
BMIValue<-BMIValue%>%
  mutate(Time5yr=case_index - 1826.25)%>%
  mutate(BMIDistance=as.numeric(BMIValueDate - Time5yr)/365.25)%>%
  subset(BMIDistance<=5)%>%
  mutate(BMIDistancePos=case_when(BMIDistance<0 ~ BMIDistance*-1,
                                  TRUE ~ BMIDistance),
         Priority=case_when(BMIDistance<=0 ~1,
                            BMIDistance>0 ~ 2,
                            TRUE ~ 3))%>%
  group_by(patid)%>%
  mutate(HighPriority=min(Priority))%>% #Take those before 5 years then those between 5 years and 0
  filter(HighPriority==Priority)%>%
  mutate(ClosestBMI=min(BMIDistancePos))%>% #Take the closest to 5 years
  filter(ClosestBMI==BMIDistancePos)%>%
  mutate(MaxBMIValue=max(BMIValue))%>% #Take the heaviest
  filter(MaxBMIValue==BMIValue)

length(unique(BMIValue$patid))

#If multiple on the same date take heaviest

BMIValue<-select(BMIValue, patid, BMIValue, ValPriority=HighPriority, BMIValueDate, ValBMIDistance=BMIDistance)
BMIValue<-distinct(BMIValue)
length(which(is.na(BMIValue$BMIValue)))

#BMI calculations
#Height
BMIHeightA<-subset(AurumHeight, patid %in% Dates$patid)
BMIHeightA<-merge(x=Dates, y=BMIHeightA, by ="patid", all.y = FALSE, all.x=FALSE)
BMIHeightA<-select(BMIHeightA, -numunitid, -medcodeid, -Description, -MaxHeight)
names(BMIHeightA)[names(BMIHeightA)=="value"]<-"Height"
names(BMIHeightA)[names(BMIHeightA)=="obsdate"]<-"HeightDate"

BMIHeightG<-subset(GoldHeight, patid.x %in% Dates$patid)
BMIHeightG<-merge(x=Dates, y=BMIHeightG, by.x ="patid", by.y="patid.x", all.y = FALSE, all.x=FALSE)
BMIHeightG<-select(BMIHeightG, -data1, -data2, -enttype.x, -Weight, -MaxHeight)
names(BMIHeightG)[names(BMIHeightG)=="eventdate"]<-"HeightDate"
names(BMIHeightG)[names(BMIHeightG)=="sysdate"]<-"enterdate"

BMIHeight<-rbind(BMIHeightA, BMIHeightG)
length(unique(BMIHeight$patid))

#If an error, take sysdate
length(which(BMIHeight$HeightDate<as.Date("1900-01-01")|BMIHeight$HeightDate>as.Date("2019-06-01")))
BMIHeight$HeightDate[BMIHeight$HeightDate<as.Date("1900-01-01")|BMIHeight$HeightDate>as.Date("2019-06-01")]<-BMIHeight$enterdate[BMIHeight$HeightDate<as.Date("1900-01-01")|BMIHeight$HeightDate>as.Date("2019-06-01")]

BMIHeight<-subset(BMIHeight, HeightDate>=as.Date("1900-01-01")&HeightDate<as.Date("2019-06-01"))

length(which(BMIHeight$HeightDate<BMIHeight$yob))
BMIHeight$HeightDate[BMIHeight$HeightDate<BMIHeight$yob]<-BMIHeight$enterdate[BMIHeight$HeightDate<BMIHeight$yob]
length(which(BMIHeight$HeightDate<BMIHeight$yob))
BMIHeight<-subset(BMIHeight, HeightDate>=yob)

#Check height was taken as an adult
BMIHeight$heightage<-as.numeric(format(BMIHeight$HeightDate, "%Y"))-BMIHeight$yob
table(BMIHeight$heightage)
#Only take those as an adult
BMIHeight<-subset(BMIHeight, heightage>=16)

#Check height
summary(BMIHeight$Height)
length(which(BMIHeight$Height<1))
BMIHeight<-subset(BMIHeight, Height>=1)

#Height doesnt need to be taken prior to index because doesnt change
length(which(is.na(BMIHeight$HeightDate)))
length(which(is.na(BMIHeight$Height)))
length(unique(BMIHeight$patid))

#Weight
BMIWeightA<-subset(AurumWeight, patid %in% Dates$patid)
BMIWeightA<-merge(x=Dates, y=BMIWeightA, by ="patid", all.y = FALSE, all.x=FALSE)
BMIWeightA<-select(BMIWeightA, -numunitid, -medcodeid, -Description)
names(BMIWeightA)[names(BMIWeightA)=="value"]<-"Weight"
names(BMIWeightA)[names(BMIWeightA)=="obsdate"]<-"WeightDate"

BMIWeightG<-subset(GoldWeight, patid.x %in% Dates$patid)
BMIWeightG<-merge(x=Dates, y=BMIWeightG, by.x ="patid", by.y="patid.x", all.y = FALSE, all.x=FALSE)
BMIWeightG<-select(BMIWeightG, -data1, -data2, -enttype.x)
names(BMIWeightG)[names(BMIWeightG)=="eventdate"]<-"WeightDate"
BMIWeightG$enterdate<-NA

BMIWeight<-rbind(BMIWeightA, BMIWeightG)
length(unique(BMIWeight$patid))

#If an error, take sysdate
length(which(BMIWeight$WeightDate<as.Date("1900-01-01")|BMIWeight$WeightDate>as.Date("2019-06-01")))
BMIWeight$WeightDate[BMIWeight$WeightDate<as.Date("1900-01-01")|BMIWeight$WeightDate>as.Date("2019-06-01")]<-BMIWeight$enterdate[BMIWeight$WeightDate<as.Date("1900-01-01")|BMIWeight$WeightDate>as.Date("2019-06-01")]

BMIWeight<-subset(BMIWeight, WeightDate>=as.Date("1900-01-01")&WeightDate<as.Date("2019-06-01"))

length(which(BMIWeight$WeightDate<BMIWeight$yob))
BMIWeight$WeightDate[BMIWeight$WeightDate<BMIWeight$yob]<-BMIWeight$enterdate[BMIWeight$WeightDate<BMIWeight$yob]
length(which(BMIWeight$WeightDate<BMIWeight$yob))
BMIWeight<-subset(BMIWeight, WeightDate>=yob)

#Check weight was taken as an adult
BMIWeight$weightage<-as.numeric(format(BMIWeight$WeightDate, "%Y"))-BMIWeight$yob
table(BMIWeight$weightage)
#Only take those as an adult
BMIWeight<-subset(BMIWeight, weightage>=16)

#Check weight
summary(BMIWeight$Weight)
length(which(BMIWeight$Weight<40))
BMIWeight<-subset(BMIWeight, Weight>=40)

length(which(is.na(BMIWeight$Weight)))
length(which(is.na(BMIWeight$WeightDate)))
length(unique(BMIWeight$patid))

#Take closest to 5 years
BMIWeight<-BMIWeight%>%
  mutate(Time5yr=case_index - 1826.25)%>%
  mutate(BMIDistance=as.numeric(WeightDate - Time5yr)/365.25)%>%
  subset(BMIDistance<=5)%>%
  mutate(BMIDistancePos=case_when(BMIDistance<0 ~ BMIDistance*-1,
                                  TRUE ~ BMIDistance),
         Priority=case_when(BMIDistance<=0 ~1,
                            BMIDistance>0 ~ 2,
                            TRUE ~ 3))%>%
  group_by(patid)%>%
  mutate(HighPriority=min(Priority))%>%
  filter(HighPriority==Priority)%>%
  mutate(ClosestBMI=min(BMIDistancePos))%>%
  filter(ClosestBMI==BMIDistancePos)%>%
  mutate(MaxWeight=max(Weight))%>%
  filter(MaxWeight==Weight)

length(unique(BMIWeight$patid))

BMIWeight<-select(BMIWeight, patid, Weight, WeightDate, HighPriority, BMIDistance)
BMIWeight<-distinct(BMIWeight)

#Dont keep height and weight unless both are present
BMICalc<-merge(x=BMIHeight, y=BMIWeight, by="patid", all.x=FALSE, all.y=FALSE)

#Calculate BMI
BMICalc$BMICalc<-BMICalc$Weight/(BMICalc$Height^2)
summary(BMICalc$BMICalc)
BMICalc<-subset(BMICalc, BMICalc<=80&BMICalc>=10)
summary(BMICalc$BMICalc)

length(which(is.na(BMICalc$HeightDate)))
length(which(is.na(BMICalc$WeightDate)))
BMICalc$CalcDate<-BMICalc$WeightDate
length(unique(BMICalc$patid))
length(which(is.na(BMICalc$BMICalc)))
BMICalc<-select(BMICalc, patid, BMICalc, WeightDate, CalcPriority=HighPriority, CalcBMIDistance=BMIDistance)
####BMI category####
BMICatA<-subset(BMIAurumCat, patid %in% Dates$patid)
BMICatA<-merge(x=Dates, y=BMICatA, by ="patid", all.y = FALSE, all.x=FALSE)
BMICatA<-select(BMICatA,-MedCodeId, -value)
names(BMICatA)[names(BMICatA)=="obsdate"]<-"BMICatDate"

BMICatG<-subset(BMIGoldCat, patid %in% Dates$patid)
BMICatG<-merge(x=Dates, y=BMICatG, by ="patid", all.y = FALSE, all.x=FALSE)
BMICatG<-select(BMICatG,-medcode)
names(BMICatG)[names(BMICatG)=="eventdate"]<-"BMICatDate"
names(BMICatG)[names(BMICatG)=="sysdate"]<-"enterdate"

BMICat<-rbind(BMICatA, BMICatG)
length(unique(BMICat$patid))

#If an error, take sysdate
length(which(BMICat$BMICatDate<as.Date("1900-01-01")|BMICat$BMICatDate>as.Date("2019-06-01")))
BMICat$BMICatDate[BMICat$BMICatDate<as.Date("1900-01-01")|BMICat$BMICatDate>as.Date("2019-06-01")]<-BMICat$enterdate[BMICat$BMICatDate<as.Date("1900-01-01")|BMICat$BMICatDate>as.Date("2019-06-01")]

BMICat<-subset(BMICat, BMICatDate>=as.Date("1900-01-01")&BMICatDate<as.Date("2019-06-01"))

length(which(BMICat$BMICatDate<BMICat$yob))
BMICat$BMICatDate[BMICat$BMICatDate<BMICat$yob]<-BMICat$enterdate[BMICat$BMICatDate<BMICat$yob]
length(which(BMICat$BMICatDate<BMICat$yob))
BMICat<-subset(BMICat, BMICatDate>=yob)

#Check category was taken as an adult
BMICat$BMIage<-as.numeric(format(BMICat$BMICatDate, "%Y"))-BMICat$yob
table(BMICat$BMIage)
#Only take those as an adult
BMICat<-subset(BMICat, BMIage>=16)

length(which(is.na(BMICat$BMICatDate)))
length(which(is.na(BMICat$BMICat)))


BMICat<-BMICat%>%
  mutate(Time5yr=case_index - 1826.25)%>%
  mutate(BMIDistance=as.numeric(BMICatDate - Time5yr)/365.25)%>%
  subset(BMIDistance<=5)%>%
  mutate(BMIDistancePos=case_when(BMIDistance<0 ~ BMIDistance*-1,
                                  TRUE ~ BMIDistance),
         Priority=case_when(BMIDistance<=0 ~1,
                            BMIDistance>0 ~ 2,
                            TRUE ~ 3))%>%
  group_by(patid)%>%
  mutate(HighPriority=min(Priority))%>%
  filter(HighPriority==Priority)%>%
  mutate(ClosestBMI=min(BMIDistancePos))%>%
  filter(ClosestBMI==BMIDistancePos)%>%
  mutate(MaxBMIValue=max(BMICat))%>%
  filter(MaxBMIValue==BMICat)

length(unique(BMICat$patid))

BMICat<-select(BMICat, patid, BMICat, BMICatDate, CatPriority=HighPriority, CatBMIDistance=BMIDistance)
BMICat<-distinct(BMICat)

#Combine all three to give most recent
BMIFinal<-merge(x=BMICat, y=BMIValue, by="patid", all.x=TRUE, all.y=TRUE)
BMIFinal<-merge(x=BMIFinal, y=BMICalc, by="patid", all.x=TRUE, all.y=TRUE)
length(unique(BMIFinal$patid))

BMIFinal$FinalPriority<-pmin(BMIFinal$CatPriority, BMIFinal$ValPriority, BMIFinal$CalcPriority, na.rm=TRUE)
BMIFinal$Measure<-0
BMIFinal$Measure[BMIFinal$ValPriority==BMIFinal$FinalPriority]<-"Value"
BMIFinal$Measure[BMIFinal$Measure==0&BMIFinal$CatPriority==BMIFinal$FinalPriority]<-"Category"
BMIFinal$Measure[BMIFinal$Measure==0&BMIFinal$CalcPriority==BMIFinal$FinalPriority]<-"Calculated"
table(BMIFinal$Measure, useNA="ifany")

BMIFinal<-BMIFinal%>%
  mutate(BMIDate=case_when(Measure=="Value" ~ BMIValueDate,
                                Measure=="Category" ~ BMICatDate,
                                Measure=="Calculated" ~ WeightDate,
                                TRUE ~ as.Date("9999-01-01")),
         BMIDistance=case_when(Measure=="Value" ~ ValBMIDistance,
                               Measure=="Category" ~ CatBMIDistance,
                               Measure=="Calculated" ~ CalcBMIDistance,
                               TRUE ~ (999)))%>%
  select(patid, BMICat, BMIValue, BMICalc, BMIDate, BMIDistance, Measure, FinalPriority)
       
#For those where BMI measure = value

BMIFinal$BMIGrp<-cut(BMIFinal$BMIValue, breaks=c(0, 18.5, 25, 30, Inf), right=FALSE, labels=c("<18.5", "18.5-<25", "25-<30", "30+"))

BMIFinal$BMI[BMIFinal$Measure=="Value"&BMIFinal$BMIGrp=="<18.5"]<-1
BMIFinal$BMI[BMIFinal$Measure=="Value"&BMIFinal$BMIGrp=="18.5-<25"]<-2
BMIFinal$BMI[BMIFinal$Measure=="Value"&BMIFinal$BMIGrp=="25-<30"]<-3
BMIFinal$BMI[BMIFinal$Measure=="Value"&BMIFinal$BMIGrp=="30+"]<-4

#For those where BMI measure = category

BMIFinal$BMI[BMIFinal$Measure=="Category"&BMIFinal$BMICat==1]<-1
BMIFinal$BMI[BMIFinal$Measure=="Category"&BMIFinal$BMICat==2]<-3
BMIFinal$BMI[BMIFinal$Measure=="Category"&BMIFinal$BMICat==3]<-4

#For those where BMI measure is calculated

BMIFinal$BMIGrpCalc<-cut(BMIFinal$BMICalc, breaks=c(0, 18.5, 25, 30, Inf), right=FALSE, labels=c("<18.5", "18.5-<25", "25-<30", "30+"))

BMIFinal$BMI[BMIFinal$Measure=="Calculated"&BMIFinal$BMIGrpCalc=="<18.5"]<-1
BMIFinal$BMI[BMIFinal$Measure=="Calculated"&BMIFinal$BMIGrpCalc=="18.5-<25"]<-2
BMIFinal$BMI[BMIFinal$Measure=="Calculated"&BMIFinal$BMIGrpCalc=="25-<30"]<-3
BMIFinal$BMI[BMIFinal$Measure=="Calculated"&BMIFinal$BMIGrpCalc=="30+"]<-4

table(BMIFinal$BMI, useNA="ifany")

####Check distribution####
BMIFinal<-BMIFinal %>% mutate(new_bin = cut(BMIDistance, breaks=seq(-36, 5, 1)))

ggplot(BMIFinal, aes(new_bin))+
  geom_bar(color="black", fill="white")+
  theme_classic() 

BMIMerge<-select(BMIFinal, patid, BMI5=BMIGrp, BMIDistance)
MatchedLK<-merge(x=MatchedLK, y=BMIMerge, by="patid", all.x=TRUE, all.y=FALSE)

save(MatchedLK, file="/Neuro/Data/BMINeuroData.Rdata")
