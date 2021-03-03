options( echo = T)
options( digits = 12)
options(stringsAsFactors = F)
#
#Install bigdata
##
# devtools::install_version("bigdata",version="0.1")
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
##
## UNLABLED CODE CHUCKS
##
defaultHCGb = function(GA){10^(-1.29035 + 0.105936 * GA - 0.00115454 * GA^2 + 0.0000036376 * GA^3)}
defaultHCGbWeight = function(Weight){10^(0.588572117 - 0.013070881*Weight + 0.0000552465*Weight^2)}
defaultPAPPA = function(GA){1.117*10^(-1.275406923 + 0.077942717*GA - 0.0002758*GA^2)}
defaultPAPPAWeight  = function(Weight){10^(1.07327971 - 0.02420418*Weight + 0.00010835*Weight^2)}

defaultAFP = function(GA){10^(-0.891564 + 0.0307963 * GA  -0.0000885561 * GA^2)}
defaultAFPWeight  = function(Weight){10^(0.47860770 - 0.010647039*Weight + 0.00004592*Weight^2)}

defaultUE3 = function(GA){1.144* 10^(-3.434651128 + 0.052359706*GA - 0.000153282*GA^2)}
defaultUE3Weight = function(Weight){10^(0.212648008 - 0.00494743*Weight + 0.000023679*Weight^2)}

defaultINHIBIN = function(GA){10^(4.02288 -0.0269949 * GA + 0.000108652 * GA^2) }
defaultINHIBINWeight  = function(Weight){0.503954 + 34.6597*(1/Weight)}

defaultINHIBA = function(GA){10^(6.64689 - 0.0726839 * GA + 0.000296817* GA^2) }
defaultINHIBAWeight  = function(Weight){0.503954 + 34.6597*(1/Weight)}

defaultNT = function(CRL){(0.857491 + 0.00553044*CRL + 0.0000512788*CRL^2)}		

defaultPlGF = function(GA){0.8*10^(1.3906 - 0.005608 * GA + 0.00008935* GA^2)}
defaultPlGFWeight  = function(Weight){10^(-0.001578 * (Weight-69))}
##
##  MORE UNLABELED CODE CHUCKS
##
AllData[,"Corr.PAPPA"] = AllData[,"PAPPA"]
AllData[AllData$SmokingGroup=="Yes","Corr.PAPPA"] = AllData[AllData$SmokingGroup=="Yes","Corr.PAPPA"]/0.82
AllData[AllData$EthnicGroup=="Afro-Caribbean","Corr.PAPPA"] = AllData[AllData$EthnicGroup=="Afro-Caribbean","Corr.PAPPA"]/1.55
AllData[AllData$EthnicGroup=="Asian","Corr.PAPPA"] = AllData[AllData$EthnicGroup=="Asian","Corr.PAPPA"]/1.08
AllData[AllData$EthnicGroup=="Oriental","Corr.PAPPA"] = AllData[AllData$EthnicGroup=="Oriental","Corr.PAPPA"]/1.09
AllData[AllData$DiabeticGroup=="Yes","Corr.PAPPA"] = AllData[AllData$DiabeticGroup=="Yes","Corr.PAPPA"]/1.02

AllData[,"Corr.HCGb"] = AllData[,"HCGb"]
AllData[AllData$SmokingGroup=="Yes","Corr.HCGb"] = AllData[AllData$SmokingGroup=="Yes","Corr.HCGb"]/0.88
AllData[AllData$EthnicGroup=="Afro-Caribbean","Corr.HCGb"] = AllData[AllData$EthnicGroup=="Afro-Caribbean","Corr.HCGb"]/1.11
AllData[AllData$EthnicGroup=="Asian","Corr.HCGb"] = AllData[AllData$EthnicGroup=="Asian","Corr.HCGb"]/0.93
AllData[AllData$EthnicGroup=="Oriental","Corr.HCGb"] = AllData[AllData$EthnicGroup=="Oriental","Corr.HCGb"]/1.06
AllData[AllData$DiabeticGroup=="Yes","Corr.HCGb"] = AllData[AllData$DiabeticGroup=="Yes","Corr.HCGb"]/0.96

AllData[,"Corr.AFP"] = AllData[,"AFP"]
AllData[AllData$EthnicGroup=="Afro-Caribbean" & AllData$GA < 98,"Corr.AFP"] = AllData[AllData$EthnicGroup=="Afro-Caribbean" & AllData$GA < 98,"Corr.AFP"]/1.23
AllData[AllData$EthnicGroup=="Afro-Caribbean" & AllData$GA >97,"Corr.AFP"] = AllData[AllData$EthnicGroup=="Afro-Caribbean" & AllData$GA > 97,"Corr.AFP"]/1.14
AllData[AllData$SmokingGroup=="Yes","Corr.AFP"] = AllData[AllData$SmokingGroup=="Yes","Corr.AFP"]/1.03
AllData[AllData$DiabeticGroup=="Yes","Corr.AFP"] = AllData[AllData$DiabeticGroup=="Yes","Corr.AFP"]/0.90

AllData[,"Corr.uE3"] = AllData[,"uE3"]
AllData[AllData$SmokingGroup=="Yes","Corr.uE3"] = AllData[AllData$SmokingGroup=="Yes","Corr.uE3"]/0.97
AllData[AllData$DiabeticGroup=="Yes","Corr.uE3"] = AllData[AllData$DiabeticGroup=="Yes","Corr.uE3"]/0.93

AllData[,"Corr.Inhibin"] = AllData[,"Inhibin"]
AllData[AllData$SmokingGroup=="Yes","Corr.Inhibin"] = AllData[AllData$SmokingGroup=="Yes","Corr.Inhibin"]/1.3
AllData[AllData$EthnicGroup=="Afro-Caribbean","Corr.Inhibin"] = AllData[AllData$EthnicGroup=="Afro-Caribbean","Corr.Inhibin"]/0.92
##
## Corr.PlGF must be made from missing *.ssc function
##
AllData[,"Corr.PlGF"] = AllData[,"PlGF"]
AllData[AllData$SmokingGroup=="Yes","Corr.PlGF"] = AllData[AllData$SmokingGroup=="Yes","Corr.PlGF"]/1.35
AllData[AllData$EthnicGroup=="Afro-Caribbean","Corr.PlGF"] = AllData[AllData$EthnicGroup=="Afro-Caribbean","Corr.PlGF"]/1.3
AllData[AllData$EthnicGroup=="Asian","Corr.PlGF"] = AllData[AllData$EthnicGroup=="Asian","Corr.PlGF"]/1.2
AllData[AllData$EthnicGroup=="Oriental","Corr.PlGF"] = AllData[AllData$EthnicGroup=="Oriental","Corr.PlGF"]/1.1
AllData[AllData$DiabeticGroup=="Yes","Corr.PlGF"] = AllData[AllData$DiabeticGroup=="Yes","Corr.PlGF"]/0.88
##

AllData[,"Re.MoM.PAPPA"]=AllData[,"Corr.PAPPA"]/(defaultPAPPA(AllData[,"GA"])*defaultPAPPAWeight(AllData[,"Weight"]))
AllData[,"Diff.PAPPA"]=AllData[,"Re.MoM.PAPPA"]-AllData[,"PAPPAMoM"]

AllData[,"Re.MoM.HCGb"]=AllData[,"Corr.HCGb"]/(defaultHCGb(AllData[,"GA"])*defaultHCGbWeight(AllData[,"Weight"]))
AllData[,"Diff.HCGb"]=AllData[,"Re.MoM.HCGb"]-AllData[,"HCGbMoM"]

AllData[,"Re.MoM.AFP"]=AllData[,"Corr.AFP"]/(defaultAFP(AllData[,"GA"])*defaultAFPWeight(AllData[,"Weight"]))
AllData[,"Diff.AFP"]=AllData[,"Re.MoM.AFP"]-AllData[,"AFPMoM"]

AllData[,"Re.MoM.uE3"]=AllData[,"Corr.uE3"]/(defaultUE3(AllData[,"GA"])*defaultUE3Weight(AllData[,"Weight"]))
AllData[,"Diff.uE3"]=AllData[,"Re.MoM.uE3"]-AllData[,"uE3MoM"]

AllData[,"Re.MoM.Inhibin"]=AllData[,"Corr.Inhibin"]/(defaultINHIBIN(AllData[,"GA"])*defaultINHIBINWeight(AllData[,"Weight"]))
AllData[,"Diff.Inhibin"]=AllData[,"Re.MoM.Inhibin"]-AllData[,"InhibinMoM"]
##
## Missing *.ssc PlGF Function
##
AllData[,"Re.MoM.PlGF"]=AllData[,"Corr.PlGF"]/(defaultPlGF(AllData[,"GA"])*defaultPlGFWeight(AllData[,"Weight"]))
AllData[,"Diff.PlGF"]=AllData[,"Re.MoM.PlGF"]-AllData[,"PlGFMoM"]
##
AllData[,"Re.MoM.NT"]=AllData[,"NT"]/defaultNT(AllData[,"CRL"])
AllData[,"Diff.NT"]=AllData[,"Re.MoM.NT"]-AllData[,"NTMoM"]

## 

AllData[,"Re.MoM.NT"]=AllData[,"NT"]/defaultNT(AllData[,"CRL"])
AllData[,"Diff.NT"]=AllData[,"Re.MoM.NT"]-AllData[,"NTMoM"]
##
## Missing timeDate Function
##
max(abs(AllData[AllData$SampleDate > timeDate("6/30/2018"),"Diff.PAPPA"]),na.rm=T)
AllData[AllData$SampleDate > timeDate("6/30/2018") & abs(AllData$Diff.PAPPA) >0.001 & is.na(AllData$Diff.PAPPA)==F ,]
max(abs(AllData[AllData$SampleDate > timeDate("6/30/2018"),"Diff.HCGb"]),na.rm=T)
AllData[AllData$SampleDate > timeDate("6/30/2018") & abs(AllData$Diff.HCGb) >0.001 & is.na(AllData$Diff.HCGb)==F ,]
max(abs(AllData[AllData$SampleDate > timeDate("6/30/2018"),"Diff.AFP"]),na.rm=T)
AllData[AllData$SampleDate > timeDate("6/30/2018") & abs(AllData$Diff.AFP) >0.001 & is.na(AllData$Diff.AFP)==F ,]
max(abs(AllData[AllData$SampleDate > timeDate("6/30/2018"),"Diff.uE3"]),na.rm=T)
AllData[AllData$SampleDate > timeDate("6/30/2018") & abs(AllData$Diff.uE3) >0.001 & is.na(AllData$Diff.uE3)==F ,]
max(abs(AllData[AllData$SampleDate > timeDate("6/30/2018"),"Diff.Inhibin"]),na.rm=T)
AllData[AllData$SampleDate > timeDate("6/30/2018") & abs(AllData$Diff.Inhibin) >0.0001 & is.na(AllData$Diff.Inhibin)==F,]
max(abs(AllData[AllData$SampleDate > timeDate("6/30/2018"),"Diff.PlGF"]),na.rm=T)
AllData[AllData$SampleDate > timeDate("6/30/2018") & abs(AllData$Diff.PlGF) >0.001 & is.na(AllData$Diff.PlGF)==F ,]
max(abs(AllData[AllData$SampleDate > timeDate("6/30/2018"),"Diff.NT"]),na.rm=T)
AllData[AllData$SampleDate > timeDate("6/30/2018") & abs(AllData$Diff.NT) >0.001 & is.na(AllData$Diff.NT)==F ,]
##
## create weight categories
##
wbreakPoints =  c( seq( from=30,to=120,by=10) ,150 )
wlabelNames = paste(wbreakPoints[-length(wbreakPoints)],wbreakPoints[-1],sep=("+ - "))
round(table(cut(AllData$Weight, breaks =  wbreakPoints, labels =wlabelNames) ) / dim(AllData)[1] *100,digits=1)
AllData[,"WeightCat"] = cut(AllData$Weight, breaks =  wbreakPoints, labels =wlabelNames)

class(AllData$SampleDate)

##
## Create collection month
##

AllData[is.na(AllData$SampleDate)==T,]
##
## Missing cut.dates Function
##
AllData[,"CollMonths"] = cut.dates(dates(substring(as.character(AllData$SampleDate),1,10)), "months")
AllData[,"USCollMonths"] = cut.dates(dates(substring(as.character(AllData$USDate),1,10)), "months")
##
## create Age categories
##
wbreakPoints =  c( seq( from=15,to=48,by=1) )
wlabelNames = paste(wbreakPoints[-length(wbreakPoints)],wbreakPoints[-1],sep=("+ - "))
#round(table(cut(AllData$Weight, breaks =  wbreakPoints, labels =wlabelNames) ) / dim(AllData)[1] *100,digits=1)
AllData[,"AgeCat"] = cut(AllData$AgeatEDD, breaks =  wbreakPoints, labels = wlabelNames)
##
##  Create basic tables contiBy and kateBy Missing Function
##
Weight = contiBy(AllData, "Weight", "GAweek",rounder=3, desc=c("Median", "Mean", "SD", "Min", "Max"), total=TRUE)
Ethnicity = kateBy(AllData, "EthnicGroup")
Smoking = kateBy(AllData, "SmokingGroup")
##
## T21 high missing MoM column
##
AllData$T21High = ifelse(AllData$T21FinalRisk > AllData$T21CutOff, 0, 1)

markers=c("HCGb","PAPPA","AFP","uE3","Inhibin","PlGF","NT")
markersFr=list()
i=1

for (marker in markers) {
  
  MoM = paste(marker,"MoM",sep="")
  markersFr[[marker]] = AllData[is.na(AllData[,eval(MoM)])==F & is.na(AllData[,eval(marker)])==F,]
  
}
##
