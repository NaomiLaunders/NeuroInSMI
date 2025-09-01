#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Find the used neurological codes in the cohort and export them for clinical review
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#clear everything
rm(list=ls(all.names = TRUE))

#Load tidyverse
library(tidyverse)

####Find all neurological conditions for each patient####
#GOLD DATABASE#
#Load the first elixhasuer and limit to neuro
load("/Raw data/observations/NewGoldElix/Ca/Elix001.Rdata")
PatNeuroAllG<-subset(PatElix, cat=="Other neurological disorders")
rm(PatElix)

#Count how many observation files we have
Observation_files <- list.files('/Raw data/observations/NewGoldElix/Ca/')

#Create a variable that just has the file numbers at the end
num<-c(1:length(Observation_files))
num[num<100&num>9]<-paste0("0", num[num<100&num>9])
for (i in 1:9) {
  num[i]<-paste0("00", num[i])
}

#Loop through all the files, load them, and "bind" them onto the "PatNeuroAll"
for (i in (2:length(num))) {
  NewFile<-paste0("/Raw data/observations/NewGoldElix/Ca/Elix", num[i], ".Rdata")
  load(NewFile)
  PatElix<-subset(PatElix, cat=="Other neurological disorders")
  PatNeuroAllG<-rbind(x=PatNeuroAllG, y= PatElix)
  print(num[i])
}

rm(PatElix)

#tidy them up
#Format dates
PatNeuroAllG$eventdate<-as.Date(PatNeuroAllG$eventdate, format="%d/%m/%Y")
PatNeuroAllG$sysdate<-as.Date(PatNeuroAllG$sysdate, format="%d/%m/%Y")

#If obs date is less than 1900 or before birth make it NA and hope theres another obs
load("/Data/ReRunMatchedLKComorbs.Rdata")
YoB<-select(MatchedLK, patid, yob)
PatNeuroAllG$eventdate[PatNeuroAllG$eventdate<as.Date("1900-01-01")]<-NA
PatNeuroAllG<-merge(x=PatNeuroAllG, y=YoB, by="patid", all.x=TRUE, all.y=FALSE)
PatNeuroAllG$eventdate[year(PatNeuroAllG$eventdate)<PatNeuroAllG$yob]<-NA
#If after 2019 then NA
PatNeuroAllG$eventdate[PatNeuroAllG$eventdate>as.Date("2019-06-01")]<-NA

#If it's NA take system date 
PatNeuroAllG$eventdate[is.na(PatNeuroAllG$eventdate)]<-PatNeuroAllG$sysdate[is.na(PatNeuroAllG$eventdate)]

#If system date now makes it out of range set back to NA
PatNeuroAllG$eventdate[PatNeuroAllG$eventdate<as.Date("1900-01-01")]<-NA
PatNeuroAllG$eventdate[year(PatNeuroAllG$eventdate)<PatNeuroAllG$yob]<-NA
PatNeuroAllG$eventdate[PatNeuroAllG$eventdate>as.Date("2019-06-01")]<-NA

#change cat to character
PatNeuroAllG$cat<-as.character(PatNeuroAllG$cat)
save(PatNeuroAllG, file="/Interim Datasets/NeuroGAll.Rdata")

#AURUM database

#Load the first elixhasuer and limit to neuro
load("/Raw data/observations/NewAurumElix/Ca/Elix001.Rdata")
PatNeuroAllA<-subset(PatElix, cat=="Other neurological disorders")
rm(PatElix)

#Count how many observation files we have
Observation_files <- list.files('/Raw data/observations/NewAurumElix/Ca/')

#Create a variable that just has the file numbers at the end
num<-c(1:length(Observation_files))
num[num<100&num>9]<-paste0("0", num[num<100&num>9])
for (i in 1:9) {
  num[i]<-paste0("00", num[i])
}

#Loop through all the files, load them, and "bind" them onto the "PatNeuroAll"
for (i in (2:length(num))) {
  NewFile<-paste0("/Raw data/observations/NewAurumElix/Ca/Elix", num[i], ".Rdata")
  load(NewFile)
  PatElix<-subset(PatElix, cat=="Other neurological disorders")
  PatNeuroAllA<-rbind(x=PatNeuroAllA, y= PatElix)
  print(num[i])
}

rm(PatElix)

#tidy them up
#Format dates
PatNeuroAllA$eventdate<-as.Date(PatNeuroAllA$obsdate, format="%d/%m/%Y")
PatNeuroAllA$sysdate<-as.Date(PatNeuroAllA$enterdate, format="%d/%m/%Y")
PatNeuroAllA<-select(PatNeuroAllA, -obsdate, -enterdate)

#If obs date is less than 1900 or before birth make it NA and hope theres another obs
PatNeuroAllA$eventdate[PatNeuroAllA$eventdate<as.Date("1900-01-01")]<-NA
PatNeuroAllA<-merge(x=PatNeuroAllA, y=YoB, by="patid", all.x=TRUE, all.y=FALSE)
PatNeuroAllA$eventdate[year(PatNeuroAllA$eventdate)<PatNeuroAllA$yob]<-NA
#If after 2019 then NA
PatNeuroAllA$eventdate[PatNeuroAllA$eventdate>as.Date("2019-06-01")]<-NA

#If it's NA take system date 
PatNeuroAllA$eventdate[is.na(PatNeuroAllA$eventdate)]<-PatNeuroAllA$sysdate[is.na(PatNeuroAllA$eventdate)]

#If system date now makes it out of range set back to NA
PatNeuroAllA$eventdate[PatNeuroAllA$eventdate<as.Date("1900-01-01")]<-NA
PatNeuroAllA$eventdate[year(PatNeuroAllA$eventdate)<PatNeuroAllA$yob]<-NA
PatNeuroAllA$eventdate[PatNeuroAllA$eventdate>as.Date("2019-06-01")]<-NA

#Change cat to character
PatNeuroAllA$cat<-as.character(PatNeuroAllA$cat)

save(PatNeuroAllA, file="/Interim Datasets/NeuroAAll.Rdata")


####Merge the two datasets together####

#How many patients are there and how many are in our cleaned file?
length(unique(PatNeuroAllA$patid))
length(unique(PatNeuroAllG$patid))

length(which(PatNeuroAllA$patid %in% MatchedLK$patid))
length(which(PatNeuroAllG$patid %in% MatchedLK$patid))

#Only include those in our study
NeuroA<-subset(PatNeuroAllA, patid %in% MatchedLK$patid)
NeuroG<-subset(PatNeuroAllG, patid %in% MatchedLK$patid)

#rename messy field
NeuroG<-rename(NeuroG, medcodeid=V1, term=readterm)
#reorder fields
NeuroG<-NeuroG[, c(1, 4, 5, 6, 7, 2, 3, 8)]
NeuroAll<-rbind(NeuroA, NeuroG)
NeuroAll<-select(NeuroAll, -yob, -cat, -score)

save(NeuroAll, file="/Interim Datasets/NeuroAll.Rdata")

####Create a list of codes used####
#Find which ones are in other lists from my codelist
LauComorbG<-read.table("/CodeLists/Charlson/Final lookups/LaundersComorbidityGold.txt", sep="\t", header=TRUE, colClasses=(c(medcode="character")))
LauComorbA<-read.table("/CodeLists/Charlson/Final lookups/LaundersComorbidityAurum.txt", sep="\t", header=TRUE, colClasses=(c(medcodeid="character")))
LauComorbA<-select(LauComorbA, medcodeid, newcat)
LauComorbG<-select(LauComorbG, medcode, newcat)

LauComorbA<-subset(LauComorbA, newcat=="Neurological disorders")
LauComorbG<-subset(LauComorbG, newcat=="Neurological disorders")

NeuroA<-merge(x=NeuroA, y=LauComorbA, by="medcodeid", all.x=TRUE, all.y=FALSE)
NeuroG<-merge(x=NeuroG, y=LauComorbG, by.x="medcodeid", by.y="medcode", all.x=TRUE, all.y=FALSE)

table(NeuroA$cat, NeuroA$newcat, useNA = "ifany")
table(NeuroG$cat, NeuroG$newcat, useNA = "ifany")

NeuroAllCat<-rbind(NeuroA, NeuroG)


#Create a table of neurological terms with a count of how many patients
NeuroTable<-NeuroAllCat%>%
  group_by(term, newcat)%>%
  summarise(n=n())

write.csv(NeuroTable, file = "/Neuro/Outputs/NeuroTermsWithStrokeAndDementiaComplete.csv")

####Now look at paralysis and paresis####
####Gold Elixhauser####
load("/Raw data/observations/NewGoldElix/Ca/Elix001.Rdata")
PatParalysisElixG<-subset(PatElix, cat=="Paralysis")
rm(PatElix)

#Count how many observation files we have
Observation_files <- list.files('/Raw data/observations/NewGoldElix/Ca/')

#Create a variable that just has the file numbers at the end
num<-c(1:length(Observation_files))
num[num<100&num>9]<-paste0("0", num[num<100&num>9])
for (i in 1:9) {
  num[i]<-paste0("00", num[i])
}

#Loop through all the files, load them, and "bind" them onto the "PatParalysisAll"
for (i in (2:length(num))) {
  NewFile<-paste0("/Raw data/observations/NewGoldElix/Ca/Elix", num[i], ".Rdata")
  load(NewFile)
  PatElix<-subset(PatElix, cat=="Paralysis")
  PatParalysisElixG<-rbind(x=PatParalysisElixG, y= PatElix)
  print(num[i])
}

rm(PatElix)

PatParalysisElixG$eventdate<-as.Date(PatParalysisElixG$eventdate, format="%d/%m/%Y")
PatParalysisElixG$sysdate<-as.Date(PatParalysisElixG$sysdate, format="%d/%m/%Y")

#If obs date is less than 1900 or before birth make it NA and hope theres another obs
PatParalysisElixG$eventdate[PatParalysisElixG$eventdate<as.Date("1900-01-01")]<-NA
PatParalysisElixG<-merge(x=PatParalysisElixG, y=YoB, by="patid", all.x=TRUE, all.y=FALSE)
PatParalysisElixG$eventdate[year(PatParalysisElixG$eventdate)<PatParalysisElixG$yob]<-NA
#If after 2019 then NA
PatParalysisElixG$eventdate[PatParalysisElixG$eventdate>as.Date("2019-06-01")]<-NA

#If it's NA take system date 
PatParalysisElixG$eventdate[is.na(PatParalysisElixG$eventdate)]<-PatParalysisElixG$sysdate[is.na(PatParalysisElixG$eventdate)]

#If system date now makes it out of range set back to NA
PatParalysisElixG$eventdate[PatParalysisElixG$eventdate<as.Date("1900-01-01")]<-NA
PatParalysisElixG$eventdate[year(PatParalysisElixG$eventdate)<PatParalysisElixG$yob]<-NA
PatParalysisElixG$eventdate[PatParalysisElixG$eventdate>as.Date("2019-06-01")]<-NA

#change cat to character
PatParalysisElixG$cat<-as.character(PatParalysisElixG$cat)
save(PatParalysisElixG, file="/Interim Datasets/ParalysisElixAAll.Rdata")

####Gold Charlson####
load("/Raw data/observations/NewGoldCh/Ca/Ch001.Rdata")
PatParalysisChG<-subset(PatCh, cat=="Hemiplegia or paraplegia")
rm(PatCh)

#Count how many observation files we have
Observation_files <- list.files('/Raw data/observations/NewGoldCh/Ca/')

#Create a variable that just has the file numbers at the end
num<-c(1:length(Observation_files))
num[num<100&num>9]<-paste0("0", num[num<100&num>9])
for (i in 1:9) {
  num[i]<-paste0("00", num[i])
}

#Loop through all the files, load them, and "bind" them onto the "PatParalysisAll"
for (i in (2:length(num))) {
  NewFile<-paste0("/Raw data/observations/NewGoldCh/Ca/Ch", num[i], ".Rdata")
  load(NewFile)
  PatCh<-subset(PatCh, cat=="Hemiplegia or paraplegia")
  PatParalysisChG<-rbind(x=PatParalysisChG, y= PatCh)
  print(num[i])
}

rm(PatCh)

PatParalysisChG$eventdate<-as.Date(PatParalysisChG$eventdate, format="%d/%m/%Y")
PatParalysisChG$sysdate<-as.Date(PatParalysisChG$sysdate, format="%d/%m/%Y")

#If obs date is less than 1900 or before birth make it NA and hope theres another obs
PatParalysisChG$eventdate[PatParalysisChG$eventdate<as.Date("1900-01-01")]<-NA
PatParalysisChG<-merge(x=PatParalysisChG, y=YoB, by="patid", all.x=TRUE, all.y=FALSE)
PatParalysisChG$eventdate[year(PatParalysisChG$eventdate)<PatParalysisChG$yob]<-NA
#If after 2019 then NA
PatParalysisChG$eventdate[PatParalysisChG$eventdate>as.Date("2019-06-01")]<-NA

#If it's NA take system date 
PatParalysisChG$eventdate[is.na(PatParalysisChG$eventdate)]<-PatParalysisChG$sysdate[is.na(PatParalysisChG$eventdate)]

#If system date now makes it out of range set back to NA
PatParalysisChG$eventdate[PatParalysisChG$eventdate<as.Date("1900-01-01")]<-NA
PatParalysisChG$eventdate[year(PatParalysisChG$eventdate)<PatParalysisChG$yob]<-NA
PatParalysisChG$eventdate[PatParalysisChG$eventdate>as.Date("2019-06-01")]<-NA

#change cat to character
PatParalysisChG$cat<-as.character(PatParalysisChG$cat)
save(PatParalysisChG, file="/Interim Datasets/ParalysisChGAll.Rdata")

####Aurum Elixhauser####
load("/Raw data/observations/NewAurumElix/Ca/Elix001.Rdata")
PatParalysisElixA<-subset(PatElix, cat=="Paralysis")
rm(PatElix)

#Count how many observation files we have
Observation_files <- list.files('/Raw data/observations/NewAurumElix/Ca/')

#Create a variable that just has the file numbers at the end
num<-c(1:length(Observation_files))
num[num<100&num>9]<-paste0("0", num[num<100&num>9])
for (i in 1:9) {
  num[i]<-paste0("00", num[i])
}

#Loop through all the files, load them, and "bind" them onto the "PatParalysisAll"
for (i in (2:length(num))) {
  NewFile<-paste0("/Raw data/observations/NewAurumElix/Ca/Elix", num[i], ".Rdata")
  load(NewFile)
  PatElix<-subset(PatElix, cat=="Paralysis")
  PatParalysisElixA<-rbind(x=PatParalysisElixA, y= PatElix)
  print(num[i])
}

rm(PatElix)

PatParalysisElixA$eventdate<-as.Date(PatParalysisElixA$obsdate, format="%d/%m/%Y")
PatParalysisElixA$sysdate<-as.Date(PatParalysisElixA$enterdate, format="%d/%m/%Y")
PatParalysisElixA<-select(PatParalysisElixA, -obsdate, -enterdate)

#If obs date is less than 1900 or before birth make it NA and hope theres another obs
PatParalysisElixA$eventdate[PatParalysisElixA$eventdate<as.Date("1900-01-01")]<-NA
PatParalysisElixA<-merge(x=PatParalysisElixA, y=YoB, by="patid", all.x=TRUE, all.y=FALSE)
PatParalysisElixA$eventdate[year(PatParalysisElixA$eventdate)<PatParalysisElixA$yob]<-NA
#If after 2019 then NA
PatParalysisElixA$eventdate[PatParalysisElixA$eventdate>as.Date("2019-06-01")]<-NA

#If it's NA take system date 
PatParalysisElixA$eventdate[is.na(PatParalysisElixA$eventdate)]<-PatParalysisElixA$sysdate[is.na(PatParalysisElixA$eventdate)]

#If system date now makes it out of range set back to NA
PatParalysisElixA$eventdate[PatParalysisElixA$eventdate<as.Date("1900-01-01")]<-NA
PatParalysisElixA$eventdate[year(PatParalysisElixA$eventdate)<PatParalysisElixA$yob]<-NA
PatParalysisElixA$eventdate[PatParalysisElixA$eventdate>as.Date("2019-06-01")]<-NA

#change cat to character
PatParalysisElixA$cat<-as.character(PatParalysisElixA$cat)
save(PatParalysisElixA, file="/Interim Datasets/ParalysisElixAAll.Rdata")

####Aurum Charlson####
load("/Raw data/observations/NewAurumCh/Ca/Ch001.Rdata")
PatParalysisChA<-subset(PatCh, cat=="Hemiplegia or paraplegia")
rm(PatCh)

#Count how many observation files we have
Observation_files <- list.files('/Raw data/observations/NewAurumCh/Ca/')

#Create a variable that just has the file numbers at the end
num<-c(1:length(Observation_files))
num[num<100&num>9]<-paste0("0", num[num<100&num>9])
for (i in 1:9) {
  num[i]<-paste0("00", num[i])
}

#Loop through all the files, load them, and "bind" them onto the "PatParalysisAll"
for (i in (2:length(num))) {
  NewFile<-paste0("/Raw data/observations/NewAurumCh/Ca/Ch", num[i], ".Rdata")
  load(NewFile)
  PatCh<-subset(PatCh, cat=="Hemiplegia or paraplegia")
  PatParalysisChA<-rbind(x=PatParalysisChA, y= PatCh)
  print(num[i])
}

rm(PatCh)

PatParalysisChA$eventdate<-as.Date(PatParalysisChA$obsdate, format="%d/%m/%Y")
PatParalysisChA$sysdate<-as.Date(PatParalysisChA$enterdate, format="%d/%m/%Y")
PatParalysisChA<-select(PatParalysisChA, -obsdate, -enterdate)

#If obs date is less than 1900 or before birth make it NA and hope theres another obs
PatParalysisChA$eventdate[PatParalysisChA$eventdate<as.Date("1900-01-01")]<-NA
PatParalysisChA<-merge(x=PatParalysisChA, y=YoB, by="patid", all.x=TRUE, all.y=FALSE)
PatParalysisChA$eventdate[year(PatParalysisChA$eventdate)<PatParalysisChA$yob]<-NA
#If after 2019 then NA
PatParalysisChA$eventdate[PatParalysisChA$eventdate>as.Date("2019-06-01")]<-NA

#If it's NA take system date 
PatParalysisChA$eventdate[is.na(PatParalysisChA$eventdate)]<-PatParalysisChA$sysdate[is.na(PatParalysisChA$eventdate)]

#If system date now makes it out of range set back to NA
PatParalysisChA$eventdate[PatParalysisChA$eventdate<as.Date("1900-01-01")]<-NA
PatParalysisChA$eventdate[year(PatParalysisChA$eventdate)<PatParalysisChA$yob]<-NA
PatParalysisChA$eventdate[PatParalysisChA$eventdate>as.Date("2019-06-01")]<-NA

#change cat to character
PatParalysisChA$cat<-as.character(PatParalysisChA$cat)
save(PatParalysisChA, file="/Interim Datasets/ParalysisChAAll.Rdata")

NotInElixG<-subset(PatParalysisChG, !(V1 %in% PatParalysisElixG$V1))
NotInCHG<-subset(PatParalysisElixG, !(V1 %in% PatParalysisChG$V1))

Check<-subset(PatParalysisChG, !(patid %in% PatParalysisElixG$patid))
length(unique(PatParalysisChG$patid))
length(unique(PatParalysisElixG$patid))

####Merge the two Elixhauser datasets together####
#Import the analysis file from the previous version of the study
load("/Data/ReRunMatchedLKComorbs.Rdata")

#How many patients are there and how many are in our cleaned file?
length(unique(PatParalysisElixA$patid))
length(unique(PatParalysisElixG$patid))

length(which(PatParalysisElixA$patid %in% MatchedLK$patid))
length(which(PatParalysisElixG$patid %in% MatchedLK$patid))

#Only include those in our study
ElixParaA<-subset(PatParalysisElixA, patid %in% MatchedLK$patid)
ElixParaG<-subset(PatParalysisElixG, patid %in% MatchedLK$patid)

#rename messy field
ElixParaG<-rename(ElixParaG, medcodeid=V1, term=readterm)
#reorder fields
ElixParaG<-ElixParaG[, c(1, 4, 5, 6, 7, 2, 3, 8)]
ElixParaAll<-rbind(ElixParaA, ElixParaG)
ElixParaAll<-select(ElixParaAll, -yob, -cat, -score)

save(ElixParaAll, file="/Interim Datasets/ElixParaAll.Rdata")

####Create a list of codes used####
#Create a table of paralysis terms with a count of how many patients
ElixParaTable<-ElixParaAll%>%
  group_by(term)%>%
  summarise(n=n())

write.csv(ElixParaTable, file = "/Neuro/Outputs/ElixParaAll.csv")

####Merge the two charlson datasets together####
#Import the analysis file from the previous version of the study
load("/Data/ReRunMatchedLKComorbs.Rdata")

#How many patients are there and how many are in our cleaned file?
length(unique(PatParalysisChA$patid))
length(unique(PatParalysisChG$patid))

length(which(PatParalysisChA$patid %in% MatchedLK$patid))
length(which(PatParalysisChG$patid %in% MatchedLK$patid))

#Only include those in our study
ChParaA<-subset(PatParalysisChA, patid %in% MatchedLK$patid)
ChParaG<-subset(PatParalysisChG, patid %in% MatchedLK$patid)

#rename messy field
ChParaG<-rename(ChParaG, medcodeid=V1, term=readterm)
#reorder fields
ChParaG<-ChParaG[, c(1, 4, 5, 6, 7, 2, 3, 8)]
ChParaAll<-rbind(ChParaA, ChParaG)
ChParaAll<-select(ChParaAll, -yob, -cat)

save(ChParaAll, file="/Interim Datasets/ChParaAll.Rdata")

####Create a list of codes used####
#Create a table of paralysis terms with a count of how many patients
ChParaTable<-ChParaAll%>%
  group_by(term)%>%
  summarise(n=n())

write.csv(ChParaTable, file = "/Neuro/Outputs/ChParaAll.csv")

#What are the differences

NotInElix<-subset(ChParaAll, !(medcodeid %in% ElixParaAll$medcodeid))
NotInCh<-subset(ElixParaAll, !(medcodeid %in% ChParaAll$medcodeid))
