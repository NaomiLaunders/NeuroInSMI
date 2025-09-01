
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Append neurological disease categories and find the first instance of each
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#clear everything
rm(list=ls(all.names = TRUE))

#Load tidyverse
library(tidyverse)
library(tidyselect)
library(tableone)

####Load look ups####

CodeLists<-read.table("/Neuro/CodeLists/2606NeuroPara.txt", sep="\t", header=TRUE, quote="")
table(CodeLists$Umbrella, useNA="ifany")
CodeLists$Umbrella[CodeLists$Umbrella=='"Disorders of nerve root, plexus or peripheral nerves"']<-"Disorders of nerve root, plexus or peripheral nerves"
table(CodeLists$Group, useNA="ifany")
CodeLists$Group[CodeLists$Group=="Vascular syndromes of brain in cerebrovascular diseases\xa0\xa0"|CodeLists$Group=="Vascular syndromes of brain in cerebrovascular diseases\xe1\xe1"]<-"Vascular syndromes of brain in cerebrovascular diseases"

CodeLists$NeuroPara<-"Neuro"
CodeLists$NeuroPara[c(687:788)]<-"Para"

#Remove quotes to make matching easier
CodeLists$term<-gsub("'","",CodeLists$term)
CodeLists$term<-gsub('"',"",CodeLists$term)
CodeLists$term<-str_trim(CodeLists$term, side = c("right"))
CodeLists$Group<-str_trim(CodeLists$Group, side = c("right"))
CodeLists$Umbrella<-str_trim(CodeLists$Umbrella, side = c("right"))

CodeLists<-CodeLists%>%
  distinct()

#Merge back read codes
DictionaryA<-read.table("/Dictionary Browsers/CPRD_CodeBrowser_201904_Aurum/CPRDAurumMedical.txt", header=TRUE, fill=TRUE, sep="\t", quote = "", dec = ".")
DictionaryA$term<-tolower(DictionaryA$Term)
DictionaryA$term<-gsub("'","",DictionaryA$term)
DictionaryA$term<-gsub('"',"",DictionaryA$term)
DictionaryA$term<-str_trim(DictionaryA$term, side = c("right"))

CodeListsFin<-merge(x=CodeLists, y=DictionaryA, by="term", all.x=TRUE, all.y=FALSE)
CodeListsA<-subset(CodeListsFin, !is.na(CodeListsFin$MedCodeId))
CodeListsA<-select(CodeListsA, term, Group, Umbrella, OriginalReadCode, CleansedReadCode)

DictionaryG<-read.table("/Dictionary Browsers/CPRD_CodeBrowser_201904_GOLD/medical_CopyWithHeadings_NL.txt", header=TRUE, fill=TRUE, sep="\t", quote = "", dec = ".")
DictionaryG$term<-tolower(DictionaryG$Read.Term)
DictionaryG$term<-gsub("'","",DictionaryG$term)
DictionaryG$term<-gsub('"',"",DictionaryG$term)
DictionaryG$term<-str_trim(DictionaryG$term, side = c("right"))

CodeListsG<-subset(CodeLists, !(term %in% CodeListsA$term ))
CodeListsG<-merge(x=CodeListsG, y=DictionaryG, by="term", all.x=TRUE, all.y=FALSE)
CodeListsG$OriginalReadCode<-NA
CodeListsG<-select(CodeListsG, term, Group, Umbrella, OriginalReadCode, CleansedReadCode=Read.Code)

CodeListFin<-rbind(CodeListsA, CodeListsG)

write.csv(CodeListFin, file = "/Neuro/CodeLists/NeuroCodesForPub.csv")

####Sort SMI codes####
SMIOrig<-read.table("/Neuro/CodeLists/SMI.txt", header=TRUE, fill=TRUE, sep="\t", quote = "", dec = ".")
SMIOrig$term<-tolower(SMIOrig$Term)
SMIOrig$term<-gsub("'","",SMIOrig$term)
SMIOrig$term<-gsub('"',"",SMIOrig$term)
SMIOrig$term<-str_trim(SMIOrig$term, side = c("right"))
SMIOrig<-rename(SMIOrig, OriginalReadCode=Read.code)

SMI<-merge(x=SMIOrig, y=DictionaryA, by="OriginalReadCode", all.x=TRUE, all.y=FALSE)
SMI<-select(SMI, term=term.x, Group, OriginalReadCode, CleansedReadCode)

write.csv(CodeListFin, file = "/Neuro/CodeLists/SMICodesForPub.csv")

####Load patient observation####

load("/Interim Datasets/NeuroAll.Rdata")
load("/Interim Datasets/ElixParaAll.Rdata")
load("/Interim Datasets/ChParaAll.Rdata")
ChParaAll<-select(ChParaAll, -ChScore)

NeuroParaAll<-rbind(NeuroAll, ChParaAll, ElixParaAll)

####Load our cohort####
load("/Data/ReRunMatchedLKComorbs.Rdata")

#Are all our obsevations for our patients?#
length(which(NeuroParaAll$patid %in% MatchedLK$patid)) #Yes they are

#Remove duplicates (so same patient, day, observation)
NeuroParaAll<-NeuroParaAll%>%
  distinct()

#Remove quotations for all to make it easier to match, and remove duplicates
NeuroParaAll$term<-gsub("'","",NeuroParaAll$term)
NeuroParaAll$term<-gsub('"',"",NeuroParaAll$term)
NeuroParaAll$term<-str_trim(NeuroParaAll$term, side = c("right"))

NeuroParaAll<-NeuroParaAll%>%
  distinct()

#Code them up by term
CodedObs<-merge(x=NeuroParaAll, y=CodeLists, by="term", all.x=TRUE, all.y=FALSE)

#Check codes
table(CodedObs$NeuroPara, useNA="ifany")
table(CodedObs$Group, useNA="ifany")
table(CodedObs$Umbrella, useNA="ifany")

#remove those excluded codes
CodedObs<-subset(CodedObs, Group!="Exclude")

####FindFirstEverGroup####
length(which(is.na(CodedObs$eventdate)))

#Create a list of group
NeuroPara<-CodedObs%>%
  group_by(patid, NeuroPara)%>%
  summarise(DateNP=min(eventdate, na.rm=TRUE))%>%
  mutate(NP=1, DateNP=case_when(is.infinite(DateNP) ~ as.Date(NA),
                            TRUE ~ DateNP))%>%
  group_by(patid)%>%
  pivot_wider(names_from = NeuroPara, values_from = c(NP, DateNP))%>%
  mutate(NP_Neuro=case_when(is.na(NP_Neuro)~0,
                               TRUE ~ NP_Neuro),
         NP_Para=case_when(is.na(NP_Para)~0,
                               TRUE ~ NP_Para))

length(which(NeuroPara$NP_Neuro>0))
length(which(NeuroPara$NP_Para>0))
length(which(NeuroPara$NP_Para>0&NeuroPara$NP_Neuro>0))


Umbrella<-CodedObs%>%
  group_by(patid, Umbrella)%>%
  summarise(UDate=min(eventdate, na.rm=TRUE))%>%
  mutate(umbrella=1, UDate=case_when(is.infinite(UDate) ~ as.Date(NA),
                        TRUE ~ UDate))%>%
  group_by(patid)%>%
  pivot_wider(names_from = Umbrella, values_from = c(umbrella, UDate))


UmbrellaVars<-vars_select(names(Umbrella), starts_with('umbrella_', ignore.case = TRUE))
Umbrella<-Umbrella%>%
  mutate_at(all_of(UmbrellaVars), ~replace_na(.,0))
  
Group<-CodedObs%>%
  group_by(patid, Group)%>%
  summarise(GDate=min(eventdate, na.rm=TRUE))%>%
  mutate(group=1, GDate=case_when(is.infinite(GDate) ~ as.Date(NA),
                                 TRUE ~ GDate))%>%
  group_by(patid)%>%
  pivot_wider(names_from = Group, values_from = c(group, GDate))


GroupVars<-vars_select(names(Group), starts_with('group_', ignore.case = TRUE))
Group<-Group%>%
  mutate_at(all_of(GroupVars), ~replace_na(.,0))

#Create counts of umbrella and fine categories
CountUmbrella<-CodedObs%>%
  select(patid, Umbrella)%>%
  group_by(patid)%>%
  distinct()%>%
  summarise(UmbrellaTotal=n())

CountGroup<-CodedObs%>%
  select(patid, Group)%>%
  group_by(patid)%>%
  distinct()%>%
  summarise(GroupTotal=n())

length(which(CountGroup$GroupTotal>1))
length(which(CountUmbrella$UmbrellaTotal>1))

Umbrella<-rename(Umbrella, Neuro_MS=`umbrella_Multiple sclerosis or other white matter disorders`,
                 Neuro_Cereb=`umbrella_Cerebrovascular disease`,
                 Neuro_Dementia=`umbrella_Dementia`,
                 Neuro_Ataxia=`umbrella_Ataxic disorders`,
                 Neuro_Epilepsy=`umbrella_Epilepsy or seizures`,
                 Neuro_Parkinsons=`umbrella_Parkinson's Disease`,
                 Neuro_Paralysis =`umbrella_Paralysis`,
                 Neuro_CP=`umbrella_Cerebral palsy`,
                 Neuro_CSF = `umbrella_Structural developmental anomalies or disorders of CSF fluid pressure`,
                 Neuro_Spinal = `umbrella_Spinal cord disorders`,
                 Neuro_ParkinsonOther=`umbrella_Parkinsonism other`,
                 Neuro_Peripheral = `umbrella_Disorders of nerve root, plexus or peripheral nerves`,
                 Neuro_Movement = `umbrella_Movement disorders other`,
                 Neuro_MND = `umbrella_Motor neuron diseases or related disorders`,
                 Neuro_Autonomic = `umbrella_Disorders of autonomic nervous system`)

UmbrellaVars<-vars_select(names(Umbrella), starts_with('Neuro_', ignore.case = TRUE))
####Merge onto our data####

MatchedLK<-merge(x=MatchedLK, y=NeuroPara, by="patid", all.x=TRUE, all.y=TRUE)
MatchedLK<-merge(x=MatchedLK, y=Umbrella, by="patid", all.x=TRUE, all.y=TRUE)
MatchedLK<-merge(x=MatchedLK, y=Group, by="patid", all.x=TRUE, all.y=TRUE)
MatchedLK<-merge(x=MatchedLK, y=CountUmbrella, by="patid", all.x=TRUE, all.y=TRUE)
MatchedLK<-merge(x=MatchedLK, y=CountGroup, by="patid", all.x=TRUE, all.y=TRUE)

MatchedLK[ , GroupVars][is.na(MatchedLK[ , GroupVars])] <- 0
MatchedLK[ , UmbrellaVars][is.na(MatchedLK[ , UmbrellaVars])] <- 0
MatchedLK$NP_Neuro[is.na(MatchedLK$NP_Neuro)]<-0
MatchedLK$NP_Para[is.na(MatchedLK$NP_Para)]<-0
MatchedLK$GroupTotal[is.na(MatchedLK$GroupTotal)]<-0
MatchedLK$UmbrellaTotal[is.na(MatchedLK$UmbrellaTotal)]<-0


####Remove variables we dont need####
####remove unneccessary variables####
names(MatchedLK)

MatchedLK<-MatchedLK[, c(1,3,4,5,6,7,8,10,22,23,28, 49, 197, 198, 199, 202, 205, 260:365)]

names(MatchedLK)

MatchedLK<-select(MatchedLK, -diagnosis_date)

table(MatchedLK$NP_Neuro, MatchedLK$NP_Para)

####Rename cryptic variables####
MatchedLK<-rename(MatchedLK, ethnicity = Ethn_missing2022)

save(MatchedLK, file="/Neuro/Data/FullNeuroData.Rdata")

####Create a table of umbrella and groupings####
table1<-CreateTableOne(vars=GroupVars, factorVars=GroupVars, strata="comparator", data=MatchedLK)
Table1Exp <- print(table1, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(Table1Exp, file = "/Neuro/Outputs/GroupCat.csv")

table2<-CreateTableOne(vars=UmbrellaVars, factorVars=UmbrellaVars, strata="comparator", data=MatchedLK)
Table2Exp <- print(table2, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
write.csv(Table2Exp, file = "/Neuro/Outputs/UmbrellaCat.csv")


