
rm(list=ls())

library(maptools)
library(rgeos)
library(rgdal)
library(stringr)

options(stringsAsFactors=F) 


## Village shapefile
PembaVill <- readOGR(paste("output/GIS/PembaVill_NBS2012",sep=""), paste("PembaVill_NBS2012_cleaned",sep=""))

## NBS 2012 human population census
NBSpop <- read.csv("data/NBS2012_PembaVillagePop.csv")



## make sure districts match
matchDist <- match(unique(PembaVill$District_N),unique(NBSpop$District))
unique(PembaVill$District_N)[which(is.na(matchDist))]

## make sure wards match
NBSpop$matchWard <- NBSpop$Ward
NBSpop$matchWard[which(NBSpop$matchWard=="Makombeni")[3:5]] <- "Makoongwe"
NBSpop$DW <- paste(NBSpop$District,NBSpop$matchWard)
PembaVill$DW <- paste(PembaVill$District_N,PembaVill$Ward_Name)
matchDW <- match(unique(PembaVill$DW),unique(NBSpop$DW))
unique(PembaVill$DW)[which(is.na(matchDW))]
which(is.na(match(unique(NBSpop$DW),unique(PembaVill$DW)))) 

## Match villages exactly (or minus most of the weird characters anyway!)
NBSpop$matchVillage <- NBSpop$Village
NBSpop$DWV <- paste(NBSpop$District,NBSpop$matchWard,tolower(NBSpop$matchVill))
PembaVill$DWV <- paste(PembaVill$District_N,PembaVill$Ward_Name,tolower(PembaVill$Vil_Mtaa_N))
matchVill <- match(PembaVill$DWV,NBSpop$DWV) 
PembaVill$DWV[which(is.na(matchVill))] # 26 unmatched
NBSpop$DWV[which(is.na(match(tolower(NBSpop$DWV),tolower(PembaVill$DWV)))) ] # 18 unmatched

## Find approximate matches
partialMatch<-rep(NA,length(which(is.na(matchVill))))
for(i in 1:length(which(is.na(matchVill)))){
  pmatch<-agrep(PembaVill$DWV[which(is.na(matchVill))[i]],NBSpop$DWV[which(is.na(match(NBSpop$DWV,PembaVill$DWV)))],max=1)
  if(length(pmatch)==1){partialMatch[i]<-agrep(PembaVill$DWV[which(is.na(matchVill))[i]],NBSpop$DWV[which(is.na(match(NBSpop$DWV,PembaVill$DWV)))],max=1)}
}
cbind(PembaVill$DWV[which(is.na(matchVill))[which(!is.na(partialMatch))]],NBSpop$DWV[which(is.na(match(NBSpop$DWV,PembaVill$DWV)))[partialMatch[which(!is.na(partialMatch))]]])
# all matches found look ok

## Add correct village names to the NBS data
NBSpop$matchVillage[which(is.na(match(NBSpop$DWV,PembaVill$DWV)))[partialMatch[which(!is.na(partialMatch))]]]<-PembaVill$Vil_Mtaa_N[which(is.na(matchVill))[which(!is.na(partialMatch))]]

## Match villages again
NBSpop$DWV <- paste(NBSpop$District,NBSpop$matchWard,tolower(NBSpop$matchVill))
PembaVill$DWV <- paste(PembaVill$District_N,PembaVill$Ward_Name,tolower(PembaVill$Vil_Mtaa_N))
matchVill <- match(PembaVill$DWV,NBSpop$DWV) 
PembaVill$DWV[which(is.na(matchVill))] # 3 unmatched
NBSpop$DWV[which(is.na(match(tolower(NBSpop$DWV),tolower(PembaVill$DWV)))) ] # 3 unmatched
## no obvious matches though

##working at ward/shehia level anyway so not an issue


## Save matched Census data
write.csv(NBSpop,file="output/NBS2012_PembaVillagePop_matched.csv",row.names = F)

