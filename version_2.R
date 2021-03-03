options( echo = T)
options( digits = 12)
options(stringsAsFactors = F)
#
library(rrtools)
# Program name
rm(list=ls())
progName <- "MedianMoMs"
# Read CSV File
##
AllData <- read.csv("./Data/StatsServicesExport.csv")
AllData$SampleDate <- as.Date(AllData$SampleDate)
AllData$DateofBirth <- as.Date(AllData$DateofBirth)
AllData$USDate <- as.Date(AllData$USDate)
dim(AllData)
AllData <- AllData[!duplicated(AllData[,"SampleNumber"]), ]
dim(AllData)
#
## DATA CLEANING
#
AllData =AllData[is.na(AllData$SampleNumber)!="",]
AllData = AllData[AllData$NumberofFetuses == 1,]
AllData = AllData[is.na(AllData$T21FinalRisk) == F,]
AllData = AllData[is.na(AllData$SampleDate) == F,]
AllData = AllData[is.na(AllData$Gestation) == F,]
AllData = AllData[is.na(AllData$IVF)==T,]

AllData = AllData[is.na(AllData$USDate) == F,]

AllData = AllData[is.na(AllData$Weight) == F,]

table(AllData$GABasedOn)

AllData = AllData[(AllData$GABasedOn == 5 | AllData$GABasedOn == 6 | AllData$GABasedOn == 7) & is.na(AllData$GABasedOn)==F,]


AllData$USDateDiff = as.numeric(AllData$SampleDate-AllData$USDate)
AllData$USGestationalAge = as.numeric(AllData$GestationalAge-AllData$USDateDiff)

AllData[,"PatientID"] = seq(nrow(AllData))
AllData[,"Lastname"] = paste("Lastname",AllData[,"PatientID"],sep="")
#

AllData[,"GA"] = AllData[,"GestationalAge"]
AllData[,"USGA"] = AllData[,"USGestationalAge"]

#DATA CLEAN: SMOKERS/NONSMOKERS
#Not known = -1
#Non-smoker = 0,
#Smoker = 1
#Ethnic Group classification
AllData$EthnicGroup = ifelse((AllData$Ethnicity=="51"), "Afro-Caribbean","Unknown")
AllData$EthnicGroup = ifelse((AllData$Ethnicity=="52"), "Asian",AllData$EthnicGroup)
AllData$EthnicGroup = ifelse((AllData$Ethnicity=="53"),"Caucasian",AllData$EthnicGroup)
AllData$EthnicGroup = ifelse((AllData$Ethnicity=="0"),"Default",AllData$EthnicGroup)
AllData$EthnicGroup = ifelse((AllData$Ethnicity=="54"),"Oriental",AllData$EthnicGroup)
AllData$EthnicGroup = ifelse((AllData$Ethnicity=="60"),"Other",AllData$EthnicGroup)
##
## Smoking Group classification
##
AllData$SmokingGroup = ifelse((AllData$Smoking=="0"), "No","Unknown")
AllData$SmokingGroup = ifelse((AllData$Smoking=="1"), "Yes",AllData$SmokingGroup)
##
##  Diabetic Group Classification
##
AllData$DiabeticGroup = ifelse((AllData$Diabetic==0), "No","Unknown")
AllData$DiabeticGroup = ifelse((AllData$Diabetic==1), "Yes", AllData$DiabeticGroup)
AllData$DiabeticGroup = ifelse((is.na(AllData$Diabetic)==T), "Unknown", AllData$DiabeticGroup)
##
## Gestational Age Group
##
AllData$GAweek = floor(AllData$GestationalAge/7)
AllData$USGAweek = floor(AllData$USGestationalAge/7)

