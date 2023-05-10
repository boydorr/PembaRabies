
library(lubridate)
library(plyr)
library(dplyr)
library(rgdal)
library(stringr)

rm(list=ls())

## Pemba Ward shapefile
PembaWard <- readOGR(paste("output/GIS/PembaWard_NBS2012",sep=""), paste("PembaWard_NBS2012",sep=""))


## Transect data
transects <- read.csv("data/TransectData2013_2014.csv")[,-1]


## Find and delete duplicated records
duplicates<-transects[which(duplicated(transects) | duplicated(transects[nrow(transects):1,])[nrow(transects):1]),]
if(nrow(duplicates)>0){transects <- transects[-which(duplicated(transects)),]}


## Any villages with multiple records for different transects on the same day?
which(duplicated(transects[,c(1:4)])) #none


## Match district names between datasets
studyDist <- match(unique(transects$District),unique(PembaWard$District_N))
unique(transects$District)[which(is.na(studyDist))] 
transects$District[which(transects$District=="Chakechake")] <- "Chake Chake"


## Match ward names between datasets
studyWard <- match(unique(paste(transects$District,transects$Ward)),unique(paste(PembaWard$District_N,PembaWard$Ward_Name)))
unique(paste(transects$District,transects$Ward))[which(is.na(studyWard))] 


## Get correct amalgamated district_ward name
transects$DW <- paste(transects$District,transects$Ward,sep="_")


## Save cleaned transect data
write.csv(transects,paste("output/transects_cleaned.csv",sep=""),row.names = F)




