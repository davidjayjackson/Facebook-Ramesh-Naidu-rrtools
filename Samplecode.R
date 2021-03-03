options( echo = T)
options( digits = 12)
options(stringsAsFactors = F)
## CAN'T FIND rtftools PACKAGE
# library(rtftools)
library(rrtools)
### CAN'T FIND bigdata package
# library(bigdata)

# Program name
progName <- "MedianMoMs"   

studyPath = getwd()

# Data folder
dataPath	= paste(studyPath,"Data",sep="\\")

# Script folder
scriptPath	= paste(studyPath,"Script",sep="\\")
funcPath	= paste(scriptPath,"Functions",sep="\\")
dataBasePath = paste(scriptPath,"dataBase",sep="\\")

# Output folder
outputPath	= paste(studyPath,"Output",sep="\\")
graphPath	= paste(outputPath,"Graphs",sep="\\")

# Read study inits
#source( paste( scriptPath, "StudyInits.ssc", sep = "\\"))

# Execution started

############################################################################################################
# PART 1: Import Data and Calculate weekly medians:
############################################################################################################
# # Get the source code for creating medians for gestational age grouping
# source(paste(funcPath, "descAndMc.ssc", sep="\\") , echo = T) 
# source(paste(funcPath, "ErrorGraphs.ssc", sep="\\") , echo = T) 
# source(paste(funcPath, "makeRtf.ssc", sep="\\") , echo = T)
# source(paste(funcPath, "CalculateMoMDailyMedians.ssc", sep="\\") , echo = T) 


#############################################
### CHANGED FILE TYPE TO CSV
### 
###Data Imported:
fileName1 = paste(dataPath,"StatsServicesExport.csv", sep="\\")
AllData = importData(fileName1, stringsAsFactors = F, time.in.format =  "%d[/]%m[/]%y")
dim(AllData)
AllData=AllData[!duplicated(AllData[,"SampleNumber"]), ]

dim(AllData)
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

#Smoker:
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

#Smoking Group classification
AllData$SmokingGroup = ifelse((AllData$Smoking=="0"), "No","Unknown")
AllData$SmokingGroup = ifelse((AllData$Smoking=="1"), "Yes",AllData$SmokingGroup)

AllData$DiabeticGroup = ifelse((AllData$Diabetic==0), "No","Unknown")
AllData$DiabeticGroup = ifelse((AllData$Diabetic==1), "Yes", AllData$DiabeticGroup)
AllData$DiabeticGroup = ifelse((is.na(AllData$Diabetic)==T), "Unknown", AllData$DiabeticGroup)

AllData$GAweek = floor(AllData$GestationalAge/7)
AllData$USGAweek = floor(AllData$USGestationalAge/7)
table(AllData$BiochemGestAge)

	
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

AllData[,"Corr.PlGF"] = AllData[,"PlGF"]
AllData[AllData$SmokingGroup=="Yes","Corr.PlGF"] = AllData[AllData$SmokingGroup=="Yes","Corr.PlGF"]/1.35
AllData[AllData$EthnicGroup=="Afro-Caribbean","Corr.PlGF"] = AllData[AllData$EthnicGroup=="Afro-Caribbean","Corr.PlGF"]/1.3
AllData[AllData$EthnicGroup=="Asian","Corr.PlGF"] = AllData[AllData$EthnicGroup=="Asian","Corr.PlGF"]/1.2
AllData[AllData$EthnicGroup=="Oriental","Corr.PlGF"] = AllData[AllData$EthnicGroup=="Oriental","Corr.PlGF"]/1.1
AllData[AllData$DiabeticGroup=="Yes","Corr.PlGF"] = AllData[AllData$DiabeticGroup=="Yes","Corr.PlGF"]/0.88

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

AllData[,"Re.MoM.PlGF"]=AllData[,"Corr.PlGF"]/(defaultPlGF(AllData[,"GA"])*defaultPlGFWeight(AllData[,"Weight"]))
AllData[,"Diff.PlGF"]=AllData[,"Re.MoM.PlGF"]-AllData[,"PlGFMoM"]

AllData[,"Re.MoM.NT"]=AllData[,"NT"]/defaultNT(AllData[,"CRL"])
AllData[,"Diff.NT"]=AllData[,"Re.MoM.NT"]-AllData[,"NTMoM"]

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

#create weight categories
wbreakPoints =  c( seq( from=30,to=120,by=10) ,150 )
wlabelNames = paste(wbreakPoints[-length(wbreakPoints)],wbreakPoints[-1],sep=("+ - "))
round(table(cut(AllData$Weight, breaks =  wbreakPoints, labels =wlabelNames) ) / dim(AllData)[1] *100,digits=1)
AllData[,"WeightCat"] = cut(AllData$Weight, breaks =  wbreakPoints, labels =wlabelNames)

class(AllData$SampleDate)

#Create collection month

AllData[is.na(AllData$SampleDate)==T,]

AllData[,"CollMonths"] = cut.dates(dates(substring(as.character(AllData$SampleDate),1,10)), "months")
AllData[,"USCollMonths"] = cut.dates(dates(substring(as.character(AllData$USDate),1,10)), "months")

#create Age categories
wbreakPoints =  c( seq( from=15,to=48,by=1) )
wlabelNames = paste(wbreakPoints[-length(wbreakPoints)],wbreakPoints[-1],sep=("+ - "))
#round(table(cut(AllData$Weight, breaks =  wbreakPoints, labels =wlabelNames) ) / dim(AllData)[1] *100,digits=1)
AllData[,"AgeCat"] = cut(AllData$AgeatEDD, breaks =  wbreakPoints, labels = wlabelNames)

#Create basic tables 
Weight = contiBy(AllData, "Weight", "GAweek",rounder=3, desc=c("Median", "Mean", "SD", "Min", "Max"), total=TRUE)
Ethnicity = kateBy(AllData, "EthnicGroup")
Smoking = kateBy(AllData, "SmokingGroup")

#T21 high
AllData$T21High = ifelse(AllData$T21FinalRisk > AllData$T21CutOff, 0, 1)

markers=c("HCGb","PAPPA","AFP","uE3","Inhibin","PlGF","NT")
markersFr=list()
i=1

for (marker in markers) {
	
	MoM = paste(marker,"MoM",sep="")
	markersFr[[marker]] = AllData[is.na(AllData[,eval(MoM)])==F & is.na(AllData[,eval(marker)])==F,]

}

#Create basic tables

###Summary sample size tables
Weektable=data.frame(table(AllData$GAweek))
Weektable$GAweek=row.names(Weektable)
Weektable$N=Weektable$X1
Weektable = data.frame(t(Weektable[,c(2,3)]))
Weektable[,"Header"]=c("GA week","N")
Weektable2 = as.matrix(Weektable[,c(dim(Weektable)[2],1:dim(Weektable)[2]-1)])

Monthtable=data.frame(table(AllData$CollMonths))
Monthtable$Month=row.names(Monthtable)
Monthtable$N=Monthtable$X1
Monthtable = data.frame(t(Monthtable[,c(2,3)]))
Monthtable[,"Header"]=c("Month","N")
Monthtable2 = as.matrix(Monthtable[,c(dim(Monthtable)[2],1:dim(Monthtable)[2]-1)])

Weighttable=data.frame(table( AllData$WeightCat))
Weighttable$Weight=row.names(Weighttable)
Weighttable$N=Weighttable$X1
Weighttable = data.frame(t(Weighttable[,c(2,3)]))
Weighttable[,"Header"]=c("Weight","N")
Weighttable2 = as.matrix(Weighttable[,c(dim(Weighttable)[2],1:dim(Weighttable)[2]-1)])

###Output Information
filename = paste(customer,"Original Median Data Summary",sep=" ")
study = customer

##Raportti----------------------------------------------------------------------------
##Open Rtf-file
fun.REP.RTF.START(paste(outputPath, "\\",sep=""), filename,studyName=study, analyte="", paperHorizontal = F,  RtfFontMap = myRtfFontMap)


wmfInfo <- rtfStartGraph(file = rtfFile, width = 7, height = 4,setsize = T, font = "Cambria",pointsize = 15)
par(mfrow=c(1,1))
min(AllData[,"AgeatEDD"])
max(AllData[,"AgeatEDD"])
hist(AllData[,"AgeatEDD"],breaks=c(seq(15,50,by=1)), probability=T,xlab="Mother's age")
axis(1, at = seq(15, 50, 5)) 
title(main="Age distribution")
rtfEndGraph(rtfFile, wmfInfo)
rtfSpaceBetweenTables(rtfFile, pageBreak=F)

rtfWriteCaption("Table ",paste("Sample sizes by gestational weeks") ,rtfFile)  
printTable(rtfFile, 
	       Weektable2, 
	       columnWidth=800, 
	       isHeaderPart= F, 
	       verticalMergeColumnCount = 2
	       )
	       
rtfSpaceBetweenTables(rtfFile, pageBreak=F)

rtfWriteCaption("Table ",paste(" Sample sizes by month") ,rtfFile) 	 
printTable(rtfFile, 
	       Monthtable2, 
	       columnWidth=800, 
	       isHeaderPart= F, 
	       verticalMergeColumnCount = 2)
	       
rtfSpaceBetweenTables(rtfFile, pageBreak=F)

rtfWriteCaption("Table ",paste(" Sample sizes by weight groups") ,rtfFile) 			 
printTable(rtfFile, 
	       Weighttable2, 
	       columnWidth=800, 
	       isHeaderPart= F, 
	       verticalMergeColumnCount = 2)
	       
rtfSpaceBetweenTables(rtfFile, pageBreak=F)

rtfWriteCaption("Table ",paste(" Ethnic distribution") ,rtfFile) 
printHeaders(rtfFile, 
			 headersText = dimnames(Ethnicity)[[2]], 
			 columnWidth=1400, 
			 isHeaderPart= T) 		 
printTable(rtfFile, 
	       Ethnicity, 
	       columnWidth=1400, 
	       isHeaderPart= F, 
	       verticalMergeColumnCount = 2)
	       
rtfSpaceBetweenTables(rtfFile, pageBreak=F)

rtfWriteCaption("Table ",paste(" Smoker distribution") ,rtfFile) 
printHeaders(rtfFile, 
			 headersText = dimnames(Smoking)[[2]], 
			 columnWidth=1400, 
			 isHeaderPart= T) 	 
printTable(rtfFile, 
	       Smoking, 
	       columnWidth=1400, 
	       isHeaderPart= F, 
	       verticalMergeColumnCount = 2)
	       
rtfSpaceBetweenTables(rtfFile, pageBreak=F)

rtfWriteCaption("Table ",paste(" Maternal Weight distribution"),rtfFile) 
printHeaders(rtfFile, headersText = dimnames(Weight)[[2]], columnWidths=1000, isHeaderPart= T)      
printTable(rtfFile, Weight, columnWidths=1000, isHeaderPart= F, verticalMergeColumnCount = 1) 
		   
rtfSpaceBetweenTables(rtfFile, pageBreak=F)
AllData1Tri = AllData[AllData$GA<98,]
AllData2Tri = AllData[AllData$GA>97,]

#Plot SPR
OriginalSPR1 = tapply(AllData1Tri[,"T21High"],AllData1Tri[, "CollMonths"],function(x) round((sum(x)/length(x))*100,digits=2))
allSPR1 = round((sum(AllData1Tri[,"T21High"])/length(AllData1Tri[,"T21High"]))*100,digits=2)
OriginalSPR1[["All"]]=allSPR1

OriginalSPR2 = tapply(AllData2Tri[,"T21High"],AllData2Tri[, "CollMonths"],function(x) round((sum(x)/length(x))*100,digits=2))

allSPR2 = round((sum(AllData2Tri[,"T21High"])/length(AllData2Tri[,"T21High"]))*100,digits=2)
OriginalSPR2[["All"]]=allSPR2

wmfInfo <- rtfStartGraph(file = rtfFile, width = 7, height = 4,setsize = T, font = "Cambria",pointsize = 15)
par(mfrow=c(1,1))

nobs1 = tapply(AllData1Tri[,"T21High"],AllData1Tri[, "CollMonths"],function(x) length(x))
nobs1[["All"]]=dim(markersFr[["PAPPA"]])[1]

xnames <- names(OriginalSPR1)
yvalues = t(matrix(c(OriginalSPR1),nrow=1,byrow=T))
xvalues = t(matrix(c(1:length(xnames)),nrow=1,byrow=T))

matplot(xvalues,yvalues,type="b",pch=c(18),col=c("green"),lty=c(1), axes=F,ylab="SPR", ylim=c(0,5))
axis(1, at=1:length(xnames), labels=xnames, cex.axis=0.7, srt = 45, adj = 1)
axis(2, at = seq(0,5,by=1), labels= seq(0,5,by=1),adj=1 , cex.axis=0.7 )
title(main="First Trimester SPR")
text(seq(1, length(xnames), by=1),OriginalSPR1 + 0.5 , labels = yvalues, cex=0.5)

rtfEndGraph(rtfFile, wmfInfo)
rtfSpaceBetweenTables(rtfFile, pageBreak=F)
rtfWriteCaption("Graph ","Screening Positive Rate", rtfFile)
rtfSpaceBetweenTables(rtfFile, pageBreak=F)

wmfInfo <- rtfStartGraph(file = rtfFile, width = 7, height = 4,setsize = T, font = "Cambria",pointsize = 15)
par(mfrow=c(1,1))

nobs1 = tapply(AllData2Tri[,"T21High"],AllData2Tri[, "CollMonths"],function(x) length(x))
nobs1[["All"]]=dim(AllData2Tri)[1]

xnames <- names(OriginalSPR2)
yvalues = t(matrix(c(OriginalSPR2),nrow=1,byrow=T))
xvalues = t(matrix(c(1:length(xnames)),nrow=1,byrow=T))

matplot(xvalues,yvalues,type="b",pch=c(18),col=c("green"),lty=c(1), axes=F,ylab="SPR", ylim=c(0,10))
axis(1, at=1:length(xnames), labels=xnames, cex.axis=0.7, srt = 45, adj = 1)
axis(2, at = seq(0,20,by=1), labels= seq(0,20,by=1),adj=1 , cex.axis=0.7 )
title(main="Second Trimester SPR")
text(seq(1, length(xnames), by=1),OriginalSPR2 + 0.5 , labels = yvalues, cex=0.5)


rtfEndGraph(rtfFile, wmfInfo)
rtfSpaceBetweenTables(rtfFile, pageBreak=F)
rtfWriteCaption("Graph ","Screening Positive Rate", rtfFile)
rtfSpaceBetweenTables(rtfFile, pageBreak=F)

table(AllData$NB, AllData$CollMonths)

markers=c("HCGb","PAPPA","AFP","uE3","Inhibin","PlGF","NT")

for (marker in markers) {
	
	print(marker)
	
	yValue= paste(marker,"MoM",sep="")
	catVal = ifelse(marker=="NT", "USGAweek", "GAweek" )
	
 	graphData=markersFr[[marker]]
	
	caption=substituteString("\\*",marker,"Original * MoMs by gestational week", global = T)	
	tableMoM = contiBy(graphData, yValue, catVal,rounder=6, desc=c("Median", "Mean","SD", "Min", "Max"), total=TRUE)
	rtfWriteCaption("Table ",caption,rtfFile) ##otsikkoa voi muuttaa
	printHeaders(rtfFile, headersText = dimnames(tableMoM)[[2]], columnWidths=1000,isHeaderPart= T)     
	printTable(rtfFile, tableMoM, columnWidths=1000, isHeaderPart= F, verticalMergeColumnCount = 1) 
		   
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	
	wmfInfo <- rtfStartGraph(file = rtfFile, width = 7, height = 4,setsize = T, font = "Cambria",pointsize = 15)
	par(mfrow=c(1,1))
	errorBarStratificationMedian(valueFrame = graphData , valueFrameName = caption,  xLabel = "GA week", yValueCol=yValue , catValueCol=catVal )
	rtfEndGraph(rtfFile, wmfInfo)
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	rtfWriteCaption("Graph ",caption, rtfFile)
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	
	#Weight
	caption=substituteString("\\*",marker,"Original * MoMs by weight categories", global = T)
	tableMoM = contiBy(graphData, yValue, "WeightCat",rounder=3, desc=c("Median", "Mean","SD", "Min", "Max"), total=TRUE)
	rtfWriteCaption("Table ",caption,rtfFile) ##otsikkoa voi muuttaa
	printHeaders(rtfFile, headersText = dimnames(tableMoM)[[2]], columnWidths=1000,isHeaderPart= T)     
	printTable(rtfFile, tableMoM, columnWidths=1000, isHeaderPart= F, verticalMergeColumnCount = 1) 
		   
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	
	wmfInfo <- rtfStartGraph(file = rtfFile, width = 7, height = 4,setsize = T, font = "Cambria",pointsize = 15)
	par(mfrow=c(1,1))
	errorBarStratificationMedian(valueFrame = graphData , valueFrameName = caption,  xLabel = "Weight (kg)", yValueCol=yValue , catValueCol="WeightCat" )
	rtfEndGraph(rtfFile, wmfInfo)
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	rtfWriteCaption("Graph ",caption, rtfFile)
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	
	#Month
	caption=substituteString("\\*",marker,"Original * MoM results by month and year", global = T)
	tableMoM = contiBy(graphData, yValue, "CollMonths",rounder=3, desc=c("Median", "Mean","SD", "Min", "Max"), total=TRUE)
	rtfWriteCaption("Table ",caption,rtfFile) ##otsikkoa voi muuttaa
	printHeaders(rtfFile, headersText = dimnames(tableMoM)[[2]], columnWidths=1000,isHeaderPart= T)     
	printTable(rtfFile, tableMoM, columnWidths=1000, isHeaderPart= F, verticalMergeColumnCount = 1) 
		   
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	
	wmfInfo <- rtfStartGraph(file = rtfFile, width = 7, height = 4,setsize = T, font = "Cambria",pointsize = 15)
	par(mfrow=c(1,1))
	errorBarStratificationMedian(valueFrame = graphData , valueFrameName = caption,  xLabel = "Specimen Month", yValueCol=yValue , catValueCol="CollMonths" )
	rtfEndGraph(rtfFile, wmfInfo)
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	rtfWriteCaption("Graph ",caption, rtfFile)
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	
	#Ethnic Group
	graphData=markersFr[[marker]]
	
	
	caption=substituteString("\\*",marker,"Original * MoM results by ethnic group", global = T)
	tableMoM = contiBy(graphData, yValue, "EthnicGroup",rounder=3, desc=c("Median", "Mean","SD", "Min", "Max"), total=TRUE)
	rtfWriteCaption("Table ",caption,rtfFile) ##otsikkoa voi muuttaa
	printHeaders(rtfFile, headersText = dimnames(tableMoM)[[2]], columnWidths=1000,isHeaderPart= T)     
	printTable(rtfFile, tableMoM, columnWidths=1000, isHeaderPart= F, verticalMergeColumnCount = 1) 
		   
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	
	wmfInfo <- rtfStartGraph(file = rtfFile, width = 7, height = 4,setsize = T, font = "Cambria",pointsize = 15)
	par(mfrow=c(1,1))
	errorBarStratificationMedian(valueFrame = graphData , valueFrameName = caption,  xLabel = "Ethnic", yValueCol=yValue , catValueCol="EthnicGroup" )
	rtfEndGraph(rtfFile,wmfInfo)
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	rtfWriteCaption("Graph ",caption, rtfFile)
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	
	#Smoking Group
	
	graphData=markersFr[[marker]]
	
	caption=substituteString("\\*",marker,"Original * MoM results by smoking group", global = T)
	tableMoM = contiBy(graphData, yValue, "SmokingGroup",rounder=3, desc=c("Median", "Mean","SD", "Min", "Max"), total=TRUE)
	rtfWriteCaption("Table ",caption,rtfFile) ##otsikkoa voi muuttaa
	printHeaders(rtfFile, headersText = dimnames(tableMoM)[[2]], columnWidths=1000,isHeaderPart= T)     
	printTable(rtfFile, tableMoM, columnWidths=1000, isHeaderPart= F, verticalMergeColumnCount = 1) 
		   
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	
	wmfInfo <- rtfStartGraph(file = rtfFile, width = 7, height = 4,setsize = T, font = "Cambria",pointsize = 15)
	par(mfrow=c(1,1))
	errorBarStratificationMedian(valueFrame = graphData , valueFrameName = caption,  xLabel = "Smoking", yValueCol=yValue , catValueCol="SmokingGroup" )
	rtfEndGraph(rtfFile, wmfInfo)
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	rtfWriteCaption("Graph ",caption, rtfFile)
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	
	#Diabetic Group
	graphData=markersFr[[marker]]
	
	caption=substituteString("\\*",marker,"Original * MoM results by diabetes group", global = T)
	tableMoM = contiBy(graphData, yValue, "DiabeticGroup",rounder=3, desc=c("Median", "Mean","SD", "Min", "Max"), total=TRUE)
	rtfWriteCaption("Table ",caption,rtfFile) ##otsikkoa voi muuttaa
	printHeaders(rtfFile, headersText = dimnames(tableMoM)[[2]], columnWidths=1000,isHeaderPart= T)     
	printTable(rtfFile, tableMoM, columnWidths=1000, isHeaderPart= F, verticalMergeColumnCount = 1) 
		   
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	
	wmfInfo <- rtfStartGraph(file = rtfFile, width = 7, height = 4,setsize = T, font = "Cambria",pointsize = 15)
	par(mfrow=c(1,1))
	errorBarStratificationMedian(valueFrame = graphData , valueFrameName = caption,  xLabel = "Diabetic", yValueCol=yValue , catValueCol="DiabeticGroup" )
	rtfEndGraph(rtfFile, wmfInfo)
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	rtfWriteCaption("Graph ",caption, rtfFile)
	rtfSpaceBetweenTables(rtfFile, pageBreak=F)
	
}

rtfClose(rtfFile)
