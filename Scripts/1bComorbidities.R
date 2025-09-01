#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Recode comorbidities from elixhauser and charlson
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls(all.names = TRUE))

library(dplyr) 
library(tidyr)
library(tidyselect)
library(lubridate)

load ("/Data/ReRunCleanBaseFile.Rdata")

####New comorbidity score####
#Dont include depression, dementia or psychosis
#Dont include chronic pulmonary disease, as COPD and asthma added in seperately later
#Dont include weight loss as that is a symptom
#Dont include blood loss anaemia as low prev an symptom
#Alcohol and drugs as risk factors

####Combine comorbidity indices

MatchedLK$Lau_CardiacArrythmia<-0
MatchedLK$Lau_CardiacArrythmia[MatchedLK$Elix_Cardiacarrythmia==1]<-1
MatchedLK$Lau_CongestiveHeartFailure<-0
MatchedLK$Lau_CongestiveHeartFailure[MatchedLK$Ch_Congestiveheartfailure==1|MatchedLK$Elix_Congestiveheartfailure==1]<-1
MatchedLK$Lau_MyocaridalInfarction<-0
MatchedLK$Lau_MyocaridalInfarction[MatchedLK$Ch_Myocardialinfarction==1]<-1

MatchedLK$Lau_CerebrovascularDisease<-0
MatchedLK$Lau_CerebrovascularDisease[MatchedLK$Ch_Cerebrovasculardisease==1]<-1
MatchedLK$Lau_Cancer<-0
MatchedLK$Lau_Cancer[MatchedLK$Elix_Solidtumourorleukaemia==1|MatchedLK$Elix_Lymphoma==1|MatchedLK$Ch_Anymalignancy==1|MatchedLK$Elix_Metastaticcancer==1|MatchedLK$Ch_Metastaticsolidtumour]<-1

MatchedLK$Lau_Diabetes<-0
MatchedLK$Lau_Diabetes[MatchedLK$`Ch_Diabetes(uncomplicated)`==1|MatchedLK$`Ch_Diabeteswithend-organdamage`==1|MatchedLK$`Elix_Diabetes(uncomplicated)`==1|MatchedLK$`Elix_Diabeteswithend-organdamage`==1]<-1
MatchedLK$Lau_Hypothyroid<-0
MatchedLK$Lau_Hypothyroid[MatchedLK$`Elix_Hypothyroidism`==1]<-1
MatchedLK$Lau_LiverDisease<-0
MatchedLK$Lau_LiverDisease[MatchedLK$Elix_Liverdisease==1|MatchedLK$Ch_Mildliverdisease==1|MatchedLK$Ch_Moderateorsevereliverdisease==1]<-1
MatchedLK$Lau_RenalDisease<-0
MatchedLK$Lau_RenalDisease[MatchedLK$Elix_Renaldisease==1|MatchedLK$Ch_Renaldisease==1]<-1
MatchedLK$Lau_PepticUlcer<-0
MatchedLK$Lau_PepticUlcer[MatchedLK$Elix_Pepticulcerdisease==1|MatchedLK$Ch_Pepticulcerdisease==1]<-1
MatchedLK$Lau_RheumaticDisease<-0
MatchedLK$Lau_RheumaticDisease[MatchedLK$Elix_Rheumatoidarthritisandcollagendiseases==1|MatchedLK$Ch_Rheumaticdisease==1]<-1

MatchedLK$Lau_Paralysis<-0
MatchedLK$Lau_Paralysis[MatchedLK$Ch_Hemiplegiaorparaplegia==1|MatchedLK$Elix_Paralysis==1]<-1
MatchedLK$Lau_HIV<-0
MatchedLK$Lau_HIV[MatchedLK$`Elix_HIV/AIDS`==1|MatchedLK$`Ch_HIV/AIDS`==1]<-1
MatchedLK$Lau_Hypertension<-0
MatchedLK$Lau_Hypertension[MatchedLK$`Elix_Hypertension(uncomplicated)`==1|MatchedLK$`Elix_Hypertensionwithend-organdamage`==1]<-1
MatchedLK$Lau_PeripheralVascularDisease<-0
MatchedLK$Lau_PeripheralVascularDisease[MatchedLK$Elix_Peripheralvasculardisease==1|MatchedLK$Ch_Peripheralvasculardisease==1]<-1
MatchedLK$Lau_PulmonaryCirculationDisorders<-0
MatchedLK$Lau_PulmonaryCirculationDisorders[MatchedLK$Elix_Pulmonarycirculationdisorders==1]<-1
MatchedLK$Lau_ValvularDisease<-0
MatchedLK$Lau_ValvularDisease[MatchedLK$Elix_Valvulardisease==1]<-1
MatchedLK$Lau_DeficiencyAnaemia<-0
MatchedLK$Lau_DeficiencyAnaemia[MatchedLK$Elix_Deficiencyanaemia==1]<-1
MatchedLK$Lau_BloodLossAnaemia<-0
MatchedLK$Lau_BloodLossAnaemia[MatchedLK$Elix_Bloodlossanaemia==1]<-1
MatchedLK$Lau_Coagulopathy<-0
MatchedLK$Lau_Coagulopathy[MatchedLK$Elix_Coagulopathy==1]<-1
MatchedLK$Lau_FluidsElectrolyte<-0
MatchedLK$Lau_FluidsElectrolyte[MatchedLK$Elix_Fluidandelectrolytedisorders==1]<-1

####Risk factors####
MatchedLK$Drugs<-0
MatchedLK$Drugs[MatchedLK$Elix_Drugabuse==1]<-1
MatchedLK$Alcohol<-0
MatchedLK$Alcohol[MatchedLK$Elix_Alcoholabuse==1]<-1

####Sort out date of first diagnosis####
#Recode dates
MatchedLK$Date_Lau_CardiacArrythmia<-MatchedLK$E_Min_Cardiacarrythmia
MatchedLK$Date_Lau_CongestiveHeartFailure<-pmin(MatchedLK$C_Min_Congestiveheartfailure, MatchedLK$E_Min_Congestiveheartfailure, na.rm = TRUE)
MatchedLK$Date_Lau_MyocaridalInfarction<-MatchedLK$C_Min_Myocardialinfarction

MatchedLK$Date_Lau_CerebrovascularDisease<-MatchedLK$C_Min_Cerebrovasculardisease
MatchedLK$Date_Lau_Cancer<-pmin(MatchedLK$E_Min_Solidtumourorleukaemia, MatchedLK$E_Min_Lymphoma, MatchedLK$C_Min_Anymalignancy, MatchedLK$E_Min_Metastaticcancer, MatchedLK$C_Min_Metastaticsolidtumour, na.rm=TRUE)

MatchedLK$Date_Lau_Diabetes<-pmin(MatchedLK$`C_Min_Diabetes(uncomplicated)`,MatchedLK$`C_Min_Diabeteswithend-organdamage`, MatchedLK$`E_Min_Diabetes(uncomplicated)`, MatchedLK$`E_Min_Diabeteswithend-organdamage`, na.rm=TRUE)
MatchedLK$Date_Lau_Hypothyroid<-MatchedLK$`E_Min_Hypothyroidism`
MatchedLK$Date_Lau_LiverDisease<-pmin(MatchedLK$E_Min_Liverdisease, MatchedLK$C_Min_Mildliverdisease, MatchedLK$C_Min_Moderateorsevereliverdisease, na.rm=TRUE)
MatchedLK$Date_Lau_RenalDisease<-pmin(MatchedLK$E_Min_Renaldisease, MatchedLK$C_Min_Renaldisease, na.rm=TRUE)
MatchedLK$Date_Lau_PepticUlcer<-pmin(MatchedLK$E_Min_Pepticulcerdisease, MatchedLK$C_Min_Pepticulcerdisease, na.rm=TRUE)
MatchedLK$Date_Lau_RheumaticDisease<-pmin(MatchedLK$E_Min_Rheumatoidarthritisandcollagendiseases, MatchedLK$C_Min_Rheumaticdisease, na.rm=TRUE)

MatchedLK$Date_Lau_Paralysis<-pmin(MatchedLK$C_Min_Hemiplegiaorparaplegia, MatchedLK$E_Min_Paralysis, na.rm=TRUE)
MatchedLK$Date_Lau_HIV<-pmin(MatchedLK$`E_Min_HIV/AIDS`, MatchedLK$`C_Min_HIV/AIDS`, na.rm=TRUE)
MatchedLK$Date_Lau_Hypertension<-pmin(MatchedLK$`E_Min_Hypertension(uncomplicated)`, MatchedLK$`E_Min_Hypertensionwithend-organdamage`, na.rm=TRUE)
MatchedLK$Date_Lau_PeripheralVascularDisease<-pmin(MatchedLK$E_Min_Peripheralvasculardisease, MatchedLK$C_Min_Peripheralvasculardisease, na.rm=TRUE)
MatchedLK$Date_Lau_PulmonaryCirculationDisorders<-MatchedLK$E_Min_Pulmonarycirculationdisorders
MatchedLK$Date_Lau_ValvularDisease<-MatchedLK$E_Min_Valvulardisease
MatchedLK$Date_Lau_DeficiencyAnaemia<-MatchedLK$E_Min_Deficiencyanaemia
MatchedLK$Date_Lau_BloodLossAnaemia<-MatchedLK$E_Min_Bloodlossanaemia
MatchedLK$Date_Lau_Coagulopathy<-MatchedLK$E_Min_Coagulopathy
MatchedLK$Date_Lau_FluidsElectrolyte<-MatchedLK$E_Min_Fluidandelectrolytedisorders

####Risk factors####
MatchedLK$Date_Drugs<-(MatchedLK$E_Min_Drugabuse)
MatchedLK$Date_Alcohol<-(MatchedLK$E_Min_Alcoholabuse)

#check none are missing
#Note, the algorithm already chooses event date and then system date if that is missing, or if before birth, before 1900 or after June 2019. Then rechecks 
#that it isnt after birth, before 1900 or after 2019. It then takes the minimum valid date. Missing ones will have no entry with a valid date.

####Sort out neurological disease####
load("/interim datasets/PatNeuroAll_WithDates.rdata")

PatNeuroAllUnique<-select(PatNeuroAllUnique, patid, Date_Lau_Neuro=min, Lau_NeurologicalDisease=Neuro)

MatchedLK<-merge(x=MatchedLK, y=PatNeuroAllUnique, by="patid", all.x=TRUE, all.y=FALSE)
MatchedLK$Lau_NeurologicalDisease[is.na(MatchedLK$Lau_NeurologicalDisease)]<-0

table(MatchedLK$Lau_NeurologicalDisease, MatchedLK$comparator)
prop.table(table(MatchedLK$Lau_NeurologicalDisease, MatchedLK$comparator),2)*100

####Add in COPD and Asthma ####

load("/interim datasets/PatCOPD_Hosp.rdata")

#Change to wide
PatCOPDUnique<-ungroup(PatCOPDUnique)
PatCOPD<-subset(PatCOPDUnique, cat=="COPD")
PatCOPD<-select(PatCOPD, patid, cat, obsdate, min)
PatCOPD$score<-1
PatCOPD<-spread(PatCOPD, cat, score)
PatCOPD<-rename(PatCOPD, Date_Lau_COPD=min, Lau_COPD=COPD)
PatCOPD<-select(PatCOPD, -obsdate)

PatAsthma<-subset(PatCOPDUnique, cat=="Asthma")
PatAsthma<-select(PatAsthma, patid, cat, obsdate, min)
PatAsthma$score<-1
PatAsthma<-spread(PatAsthma, cat, score)
PatAsthma<-rename(PatAsthma, Date_Lau_Asthma=min, Lau_Asthma=Asthma)
PatAsthma<-select(PatAsthma, -obsdate)

MatchedLK<-merge(x=MatchedLK, y=PatCOPD, by="patid", all.x=TRUE, all.y=FALSE)
MatchedLK<-merge(x=MatchedLK, y=PatAsthma, by="patid", all.x=TRUE, all.y=FALSE)

MatchedLK$Lau_COPD[is.na(MatchedLK$Lau_COPD)]<-0
MatchedLK$Lau_Asthma[is.na(MatchedLK$Lau_Asthma)]<-0

#Create a count of comorbidities
Lau<-vars_select(names(MatchedLK), starts_with('Lau_', ignore.case = TRUE))
Score<-select(MatchedLK, patid, all_of(Lau))

#Convert from wide to long
Score1<-Score%>%
  gather("Comorb", "Total", 2:25)

#Count score
Score2<-Score1 %>%
  group_by(patid) %>%
  summarise(Lau_Total= sum(Total))

#Merge score back to main dataset and save
MatchedLK<-merge(x=MatchedLK, y=Score2, by="patid", all.x=TRUE, all.y=TRUE)


save(MatchedLK, file="/Data/ReRunMatchedLKComorbs.Rdata")




