#
## Create village and district shapefiles for STz region
#__________________________

rm(list=ls())

library(maptools)
library(raster)
library(rgeos)
library(rgdal)
library(stringr)

options(stringsAsFactors=F) 




## Read in data 
#--------------

## Projections
crs37S <- CRS("+proj=utm +zone=37 +south +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +units=m +no_defs")
crs36S <- CRS("+proj=utm +zone=36 +south +ellps=clrk80 +towgs84=-160,-6,-302,0,0,0,0 +units=m +no_defs")
crsLL <- CRS("+proj=longlat")

## Load in Tz administrative shapefiles 
districts <- readOGR("data/GIS/TZ_District_2012","TZ_District_2012")
wards <- readOGR("data/GIS/TZ_Ward_2012","TZ_Ward_2012")
vill <- readOGR("data/GIS/NBSVill","NBS2012vill_new_v6", p4s=as.character(crs36S))



## Extract out Pemba and save
#--------------

## Pemba districts
PembaDist <- districts[which(districts$Region_Nam=="Kaskazini Pemba"|districts$Region_Nam=="Kusini Pemba"),]
plot(PembaDist)

## Pemba wards (shehia)
PembaWard <- wards[which(wards$Region_Nam=="Kaskazini Pemba"|wards$Region_Nam=="Kusini Pemba"),]
plot(PembaWard)

## Pemba villages
PembaVill <- vill[which(vill$Region_Nam=="Kaskazini Pemba"|vill$Region_Nam=="Kusini Pemba"),]
plot(PembaVill)

## Transform to UTMs
PembaVill <- spTransform(PembaVill,crs37S)
PembaWard <- spTransform(PembaWard,crs37S)
PembaDist <- spTransform(PembaDist,crs37S)

## Add areas to shapefiles
PembaVill$Area_kmsq<-gArea(PembaVill,byid = T)/1000000
PembaWard$Area_kmsq<-gArea(PembaWard,byid = T)/1000000
PembaDist$Area_kmsq<-gArea(PembaDist,byid = T)/1000000

## Save study area shapefiles
if(!dir.exists("output/GIS")){dir.create("output/GIS")}
writeOGR(PembaVill, dsn="output/GIS/PembaVill_NBS2012", "PembaVill_NBS2012", driver="ESRI Shapefile",overwrite_layer=T)
writeOGR(PembaWard, dsn="output/GIS/PembaWard_NBS2012", "PembaWard_NBS2012", driver="ESRI Shapefile",overwrite_layer=T)
writeOGR(PembaDist, dsn="output/GIS/PembaDist_NBS2012", "PembaDist_NBS2012", driver="ESRI Shapefile",overwrite_layer=T)


## Village shapefile cleaned in Q (mostly removal of stray boundary lines)
PembaVill <- readOGR(paste("output/GIS/PembaVill_NBS2012",sep=""), paste("PembaVill_NBS2012_cleaned",sep=""))
PembaVill$Area_kmsq<-gArea(PembaVill,byid = T)/1000000
writeOGR(PembaVill, dsn="output/GIS/PembaVill_NBS2012", "PembaVill_NBS2012_cleaned", driver="ESRI Shapefile",overwrite_layer=T)

