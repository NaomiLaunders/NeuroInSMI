#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Create clean file for analysis
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~Libraries
library(tidyverse)
library(lubridate)
library(tidyselect)


#Clear environment
rm(list = ls(all.names = TRUE))
load("/Neuro/Data/RisksNeuroData.Rdata")

#Define which fields are the conditions
Neuro<-vars_select(names(MatchedLK), starts_with('Neuro_', ignore.case = TRUE))
#define which fields are the dates
Dates<-vars_select(names(MatchedLK), starts_with('UDate_', ignore.case = TRUE))

#Just double check they are in the same order
Check<-cbind(Neuro, Dates)

#A few will have missing diagnosis so put these as "missing"
for (i in (1:length(Neuro))) {
  MatchedLK [paste0("Miss", Neuro[i])]<-if_else((MatchedLK[[paste0(Neuro[i])]]==1)&is.na(MatchedLK[[paste0(Dates[i])]]), 1, 0)
  }

####Split by years before/after diagnosis for Neuro, risk factors and death####

#Create a field that shows if they were active after diagnosis

MatchedLK$Act_a1<-if_else(MatchedLK$end-MatchedLK$case_index>0, 1, 0) # Do they enter the first year after diagnosis?
MatchedLK$Act_a2<-if_else(MatchedLK$end-MatchedLK$case_index>365.25, 1, 0) #Do they enter the second year?
MatchedLK$Act_a3<-if_else(MatchedLK$end-MatchedLK$case_index>730.5, 1, 0) # Do they enter the third year?
MatchedLK$Act_a4<-if_else(MatchedLK$end-MatchedLK$case_index>1095.75, 1, 0) #Do they enter the fourth year?
MatchedLK$Act_a5<-if_else(MatchedLK$end-MatchedLK$case_index>1461, 1, 0) #at any time after the fifth year
MatchedLK$Act_after<-if_else(MatchedLK$end-MatchedLK$case_index>1826.25, 1, 0)#at any time after the fifth year

#Then run a loop to create variables which show whether they have each disease at each time point. For those after, only include if they are "active"
for (i in (1:length(Neuro))) {       

    MatchedLK [paste0(Neuro[i], "_prior")]<-if_else(MatchedLK$case_index-MatchedLK[[paste0(Dates[i])]]>1826.25&!is.na(MatchedLK[[paste0(Dates[i])]]), 1, 0) #Are they diagnosed at 5 years before?
    MatchedLK [paste0(Neuro[i], "_b5")]<-if_else(MatchedLK$case_index-MatchedLK[[paste0(Dates[i])]]>1461&!is.na(MatchedLK[[paste0(Dates[i])]]), 1, 0) #Are they diagnosed at 4 years before
    MatchedLK [paste0(Neuro[i], "_b4")]<-if_else(MatchedLK$case_index-MatchedLK[[paste0(Dates[i])]]>1095.75&!is.na(MatchedLK[[paste0(Dates[i])]]), 1, 0) #Are they diagnosed at 3 years before
    MatchedLK [paste0(Neuro[i], "_b3")]<-if_else(MatchedLK$case_index-MatchedLK[[paste0(Dates[i])]]>730.5&!is.na(MatchedLK[[paste0(Dates[i])]]), 1, 0) #Are they diagnosed at 2 years before
    MatchedLK [paste0(Neuro[i], "_b2")]<-if_else(MatchedLK$case_index-MatchedLK[[paste0(Dates[i])]]>365.25&!is.na(MatchedLK[[paste0(Dates[i])]]), 1, 0) #Are they diagnosed at 1 year before
    MatchedLK [paste0(Neuro[i], "_b1")]<-if_else(MatchedLK$case_index-MatchedLK[[paste0(Dates[i])]]>0&!is.na(MatchedLK[[paste0(Dates[i])]]), 1, 0)#Are they diagnosed by 1 day before index?
    
    MatchedLK [paste0(Neuro[i], "_a1")]<-if_else(MatchedLK$case_index-MatchedLK[[paste0(Dates[i])]]>-365.25&!is.na(MatchedLK[[paste0(Dates[i])]])&MatchedLK$Act_a1==1, 1, 0) #Are they diagnosed at 1 year after
    MatchedLK [paste0(Neuro[i], "_a2")]<-if_else(MatchedLK$case_index-MatchedLK[[paste0(Dates[i])]]>-730.5&!is.na(MatchedLK[[paste0(Dates[i])]])&MatchedLK$Act_a2==1, 1, 0) #Are they diagnosed at 2 years after
    MatchedLK [paste0(Neuro[i], "_a3")]<-if_else(MatchedLK$case_index-MatchedLK[[paste0(Dates[i])]]>-1095.75&!is.na(MatchedLK[[paste0(Dates[i])]])&MatchedLK$Act_a3==1, 1, 0)#Are they diagnosed at 3 years after
    MatchedLK [paste0(Neuro[i], "_a4")]<-if_else(MatchedLK$case_index-MatchedLK[[paste0(Dates[i])]]>-1461&!is.na(MatchedLK[[paste0(Dates[i])]])&MatchedLK$Act_a4==1, 1, 0) #are they diagnosed at 4 years after
    MatchedLK [paste0(Neuro[i], "_a5")]<-if_else(MatchedLK$case_index-MatchedLK[[paste0(Dates[i])]]>-1826.25&!is.na(MatchedLK[[paste0(Dates[i])]])&MatchedLK$Act_a5==1, 1, 0) #Are they diagnosed at 5 years after
    MatchedLK [paste0(Neuro[i], "_after")]<-if_else(!is.na(MatchedLK[[paste0(Dates[i])]])&MatchedLK$Act_after==1, 1, 0)
    
}

####Create count of multimorbidity####
table(MatchedLK$UmbrellaTotal)
MatchedLK<-rename(MatchedLK, Count=UmbrellaTotal)

####Do the same for various time points####
Neuro<-vars_select(names(MatchedLK), ends_with('prior', ignore.case = TRUE))
Count<-select(MatchedLK, patid, all_of(Neuro))
#Convert from wide to long
Count1<-Count%>%
  gather("Neuro", "Total", 2:16)
#Count score
Count2<-Count1 %>%
  group_by(patid) %>%
  summarise(count_prior= sum(Total))
#Merge score back to main dataset and save
MatchedLK<-merge(x=MatchedLK, y=Count2, by="patid", all.x=TRUE, all.y=TRUE)

Neuro<-vars_select(names(MatchedLK), ends_with('b5', ignore.case = TRUE))
Count<-select(MatchedLK, patid, all_of(Neuro))
#Convert from wide to long
Count1<-Count%>%
  gather("Neuro", "Total", 2:16)
#Count score
Count2<-Count1 %>%
  group_by(patid) %>%
  summarise(count_b5= sum(Total))
#Merge score back to main dataset and save
MatchedLK<-merge(x=MatchedLK, y=Count2, by="patid", all.x=TRUE, all.y=TRUE)

Neuro<-vars_select(names(MatchedLK), ends_with('b4', ignore.case = TRUE))
Count<-select(MatchedLK, patid, all_of(Neuro))
#Convert from wide to long
Count1<-Count%>%
  gather("Neuro", "Total", 2:16)
#Count score
Count2<-Count1 %>%
  group_by(patid) %>%
  summarise(count_b4= sum(Total))
#Merge score back to main dataset and save
MatchedLK<-merge(x=MatchedLK, y=Count2, by="patid", all.x=TRUE, all.y=TRUE)

Neuro<-vars_select(names(MatchedLK), ends_with('b3', ignore.case = TRUE))
Count<-select(MatchedLK, patid, all_of(Neuro))
#Convert from wide to long
Count1<-Count%>%
  gather("Neuro", "Total", 2:16)
#Count score
Count2<-Count1 %>%
  group_by(patid) %>%
  summarise(count_b3= sum(Total))
#Merge score back to main dataset and save
MatchedLK<-merge(x=MatchedLK, y=Count2, by="patid", all.x=TRUE, all.y=TRUE)

Neuro<-vars_select(names(MatchedLK), ends_with('b2', ignore.case = TRUE))
Count<-select(MatchedLK, patid, all_of(Neuro))
#Convert from wide to long
Count1<-Count%>%
  gather("Neuro", "Total", 2:16)
#Count score
Count2<-Count1 %>%
  group_by(patid) %>%
  summarise(count_b2= sum(Total))
#Merge score back to main dataset and save
MatchedLK<-merge(x=MatchedLK, y=Count2, by="patid", all.x=TRUE, all.y=TRUE)

Neuro<-vars_select(names(MatchedLK), ends_with('b1', ignore.case = TRUE))
Count<-select(MatchedLK, patid, all_of(Neuro))
#Convert from wide to long
Count1<-Count%>%
  gather("Neuro", "Total", 2:16)
#Count score
Count2<-Count1 %>%
  group_by(patid) %>%
  summarise(count_b1= sum(Total))
#Merge score back to main dataset and save
MatchedLK<-merge(x=MatchedLK, y=Count2, by="patid", all.x=TRUE, all.y=TRUE)

Neuro<-vars_select(names(MatchedLK), ends_with('a1', ignore.case = TRUE))
Neuro<-Neuro[2:16] #Need to remove "the act_A1 field"
Count<-select(MatchedLK, patid, all_of(Neuro))
#Convert from wide to long
Count1<-Count%>%
  gather("Neuro", "Total", 2:16)
#Count score
Count2<-Count1 %>%
  group_by(patid) %>%
  summarise(count_a1= sum(Total))
#Merge score back to main dataset and save
MatchedLK<-merge(x=MatchedLK, y=Count2, by="patid", all.x=TRUE, all.y=TRUE)

Neuro<-vars_select(names(MatchedLK), ends_with('a2', ignore.case = TRUE))
Neuro<-Neuro[2:16]
Count<-select(MatchedLK, patid, all_of(Neuro))
#Convert from wide to long
Count1<-Count%>%
  gather("Neuro", "Total", 2:16)
#Count score
Count2<-Count1 %>%
  group_by(patid) %>%
  summarise(count_a2= sum(Total))
#Merge score back to main dataset and save
MatchedLK<-merge(x=MatchedLK, y=Count2, by="patid", all.x=TRUE, all.y=TRUE)

Neuro<-vars_select(names(MatchedLK), ends_with('a3', ignore.case = TRUE))
Neuro<-Neuro[2:16]
Count<-select(MatchedLK, patid, all_of(Neuro))
#Convert from wide to long
Count1<-Count%>%
  gather("Neuro", "Total", 2:16)
#Count score
Count2<-Count1 %>%
  group_by(patid) %>%
  summarise(count_a3= sum(Total))
#Merge score back to main dataset and save
MatchedLK<-merge(x=MatchedLK, y=Count2, by="patid", all.x=TRUE, all.y=TRUE)

Neuro<-vars_select(names(MatchedLK), ends_with('a4', ignore.case = TRUE))
Neuro<-Neuro[2:16]
Count<-select(MatchedLK, patid, all_of(Neuro))
#Convert from wide to long
Count1<-Count%>%
  gather("Neuro", "Total", 2:16)
#Count score
Count2<-Count1 %>%
  group_by(patid) %>%
  summarise(count_a4= sum(Total))
#Merge score back to main dataset and save
MatchedLK<-merge(x=MatchedLK, y=Count2, by="patid", all.x=TRUE, all.y=TRUE)

Neuro<-vars_select(names(MatchedLK), ends_with('a5', ignore.case = TRUE))
Neuro<-Neuro[2:16]
Count<-select(MatchedLK, patid, all_of(Neuro))
#Convert from wide to long
Count1<-Count%>%
  gather("Neuro", "Total", 2:16)
#Count score
Count2<-Count1 %>%
  group_by(patid) %>%
  summarise(count_a5= sum(Total))
#Merge score back to main dataset and save
MatchedLK<-merge(x=MatchedLK, y=Count2, by="patid", all.x=TRUE, all.y=TRUE)

Neuro<-vars_select(names(MatchedLK), ends_with('after', ignore.case = TRUE))
Neuro<-Neuro[2:16]
Count<-select(MatchedLK, patid, all_of(Neuro))
#Convert from wide to long
Count1<-Count%>%
  gather("Neuro", "Total", 2:16)
#Count score
Count2<-Count1 %>%
  group_by(patid) %>%
  summarise(count_after= sum(Total))
#Merge score back to main dataset and save
MatchedLK<-merge(x=MatchedLK, y=Count2, by="patid", all.x=TRUE, all.y=TRUE)

####Re-do start date so can calculate baseline time####

#Create follow up time after index
MatchedLK$FU<-MatchedLK$end-MatchedLK$case_index
MatchedLK$FU<-as.numeric(MatchedLK$FU)/365.25

#Follow up before index
MatchedLK$Baseline<-MatchedLK$case_index-MatchedLK$crd
MatchedLK$Baseline<-as.numeric(MatchedLK$Baseline)/365.25

#Age at death
MatchedLK$AgeAtDeath<-MatchedLK$deathdate-MatchedLK$dob
MatchedLK$AgeAtDeath<-as.numeric(MatchedLK$AgeAtDeath)
MatchedLK$AgeAtDeath<-MatchedLK$AgeAtDeath/365.25

#Calendar years
MatchedLK$calendaryear<-year(MatchedLK$case_index)

#Sort SMI
MatchedLK$SMI<-0
MatchedLK$SMI[MatchedLK$comparator==0]<-1
MatchedLK$SMI<-as.factor(MatchedLK$SMI)
MatchedLK<-select(MatchedLK, -comparator)

####Save analysis file####
save(MatchedLK, file="/Neuro/Data/FinalNeuroData.Rdata") 

Variables<-as.data.frame(names(MatchedLK))

write.table(Variables, "/Neuro/Outputs/FieldNames.txt")
