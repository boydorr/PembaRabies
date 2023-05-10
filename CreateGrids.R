#
## Create grids of study region
#__________________________

rm(list=ls())

library(maptools)
library(raster)
library(rgeos)
library(rgdal)
library(stringr)

options(stringsAsFactors=F) 

# DropboxShared <- "D:/Dropbox/SharedFolders"
# setwd("D:/Dropbox/RabiesPostdoc/Pemba")

cell_size<-1 #km


## Read in data 
#--------------

## Load in administrative shapefiles 
PembaVill <- readOGR("output/GIS/PembaVill_NBS2012","PembaVill_NBS2012_cleaned")
PembaWard <- readOGR("output/GIS/PembaWard_NBS2012","PembaWard_NBS2012")
PembaDist <- readOGR("output/GIS/PembaDist_NBS2012","PembaDist_NBS2012")



## Create Pemba grids
#--------------

## Create grid
grid <- raster(extent(PembaVill),crs=PembaVill@proj4string)
res(grid) <- cell_size*1000

##Create Pemba Grid
PembaVill$studyArea<-1
gridPoints <- SpatialPoints(rasterToPoints(grid), proj4string = PembaVill@proj4string)
values <- over(gridPoints,PembaVill)[,"studyArea"]; values[which(values==0)]<-1
PembaGrid <- grid
PembaGrid[] <- as.numeric(values)
plot(PembaGrid)

## Create 2012 village grid
PembaVill@data$VillID <- 1:nrow(PembaVill@data)
values <- over(gridPoints,PembaVill)$VillID
villGrid <- grid
villGrid[] <- values
cellsToFill <- which(is.na(villGrid[])&!is.na(PembaGrid[]))
plot(villGrid)

## Create 2012 district grid
PembaDist$DistID<-1:nrow(PembaDist)
values <- over(gridPoints,PembaDist)$DistID
distGrid <- grid
distGrid[] <- values
cellsToFill <- which(is.na(villGrid[])&!is.na(distGrid[]))
distGrid[cellsToFill]<-NA
plot(distGrid)

## Create 2012 ward grid
PembaWard$WardID<-1:nrow(PembaWard)
values <- over(gridPoints,PembaWard)$WardID
wardGrid <- grid
wardGrid[] <- values
wardGrid[cellsToFill]<-NA
plot(wardGrid)

## Create 2012 region grid
PembaVill$RegionID<-match(PembaVill$Region_Nam,sort(unique(PembaVill$Region_Nam)))
values <- over(gridPoints,PembaVill)$RegionID
regionGrid <- grid
regionGrid[] <- values
plot(regionGrid)

##Save grids
if(!dir.exists(paste("output/GIS/",cell_size^2,"kmsqGrids",sep=""))){dir.create(paste("output/GIS/",cell_size^2,"kmsqGrids",sep=""))}
writeRaster(PembaGrid,file=paste("output/GIS/",cell_size^2,"kmsqGrids/PembaGrid",cell_size^2,"kmsq.grd",sep=""),overwrite=T)
writeRaster(distGrid,file=paste("output/GIS/",cell_size^2,"kmsqGrids/distGrid",cell_size^2,"kmsq.grd",sep=""),overwrite=T)
writeRaster(villGrid,file=paste("output/GIS/",cell_size^2,"kmsqGrids/villGrid",cell_size^2,"kmsq.grd",sep=""),overwrite=T)
writeRaster(wardGrid,file=paste("output/GIS/",cell_size^2,"kmsqGrids/wardGrid",cell_size^2,"kmsq.grd",sep=""),overwrite=T)
writeRaster(regionGrid,file=paste("output/GIS/",cell_size^2,"kmsqGrids/regionGrid",cell_size^2,"kmsq.grd",sep=""),overwrite=T)


## Convert PembaGrid to polygons and crop to land 
PembaUTM <- rasterToPolygons(PembaGrid)
PembaUTM <- PembaUTM[which(PembaUTM@data$layer==1),]
plot(PembaUTM)

## Give ID to each Pemba cell
PembaUTM$cellID<-1:nrow(PembaUTM@data)
cellGrid<-PembaGrid
cellGrid[]<-over(gridPoints,PembaUTM)$cellID
cellGrid[which(is.na(cellGrid@data@values))]<-0
plot(cellGrid)

## Add additional information on cells to gridded polygon
PembaUTM$Region <- sort(unique(PembaVill$Region_Nam))[(regionGrid[which(!is.na(regionGrid@data@values))])]
PembaUTM$RegionID <- regionGrid[which(!is.na(regionGrid@data@values))]
PembaUTM$District <- PembaDist$District_N[(distGrid[which(!is.na(distGrid@data@values))])]
PembaUTM$DistrictID <- distGrid[which(!is.na(distGrid@data@values))]
PembaUTM$Ward <- PembaWard$matchVill[(wardGrid[which(!is.na(wardGrid@data@values))])]
PembaUTM$WardID <- wardGrid[which(!is.na(wardGrid@data@values))]
PembaUTM$Village <- PembaVill$Vil_Mtaa_N[(villGrid[which(!is.na(villGrid@data@values))])]
PembaUTM$VillageID <- villGrid[which(!is.na(villGrid@data@values))]

# Save cellGrid and PembaUTM for use in IBM
write.table(as.matrix(cellGrid),paste("output/Pemba_matrix_",cell_size^2,"kmsq_cellID.csv",sep=""),row.names=F,col.names=F,sep=",")
writeOGR(PembaUTM, dsn=paste("output/GIS/Pemba_gridded",cell_size^2,"kmsq",sep=""), paste("Pemba_gridded",cell_size^2,"kmsq",sep=""), driver="ESRI Shapefile", overwrite_layer=T, check_exists=T)


