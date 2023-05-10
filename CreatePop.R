
rm(list=ls())

library(maptools)
library(rgeos)
library(rgdal)
library(stringr)
library(lubridate)
library(raster)
library("RColorBrewer") 

options(stringsAsFactors=F) 

set.seed(0)



## Read in data
#__________________

## Gridded shapefile
PembaUTM <- readOGR(paste("output/GIS/Pemba_gridded1kmsq",sep=""), paste("Pemba_gridded1kmsq",sep=""))

## Ward shapefiles
PembaWard <- readOGR(paste("output/GIS/PembaWard_NBS2012",sep=""), paste("PembaWard_NBS2012",sep=""))
TzWard2002 <- readOGR(paste("data/GIS/Tz_wards_2002",sep=""), paste("Tz_wards_2002",sep="")) #2002 version for pop data

## District shapefile
PembaDist <- readOGR(paste("output/GIS/PembaDist_NBS2012",sep=""), paste("PembaDist_NBS2012",sep=""))

## Pemba Grid
cell_size <- 1
PembaGrid <- raster(paste("output/GIS/",cell_size^2,"kmsqGrids/PembaGrid",cell_size^2,"kmsq.grd",sep=""))

## NBS Human Census
NBS2012 <- read.csv("output/NBS2012_PembaVillagePop_matched.csv")

## Census data
census2012 <- read.csv("output/dogPopCensus2012_matched.csv")
census2017_2019 <- read.csv("output/dogPopCensus2017_2019_matched.csv")

## Dog Vaccination data
 vax <- read.csv("output/VaxPemba.csv")
# vax <- read.csv("output/VaxPemba_edit2014dates.csv")

## Transect Data
transects <- read.csv("output/transects_cleaned.csv")

## Pup adult ratio from Serengeti census
PAR <- 0.2562058 #pup/adult ratio

## Human population growth rate
pop2002 <- sum(TzWard2002$TOTAL[which(TzWard2002$REGION%in%c("North Pemba","South Pemba"))])
pop2012 <- 211732+195116 # from census
monthly_growth_rate <- (pop2012/pop2002)^(1/((2012-2002)*12))-1

## WorldPop data (re-projected and clipped to Pemba grid)
PembaPop15 <- raster(paste("data/GIS/WorldPop_reproject_",cell_size^2,"kmsq_2015.grd",sep=""))
PembaPop10 <- raster(paste("data/GIS/WorldPop_reproject_",cell_size^2,"kmsq_2010.grd",sep=""))

## Function for distributing individuals in cells for each village 
source("R/popMap.R")




## Time info
#__________________

## Time period of interest
years <- 2010:2020
months <- 12*length(years)
startDate <- as.Date(paste0(min(years),"-01-01"))
endDate <- as.Date(paste0(max(years),"-01-01"))

## add month & year to vaccination and transect data
vax$dateVaccination <- as.Date(vax$vacc_date,format="%d/%m/%Y")
vax$month <- month(vax$dateVaccination) + (year(vax$dateVaccination)-year(startDate))*12
vax$year <- year(vax$dateVaccination)
transects$Date <- as.Date(transects$Date,format = "%d/%m/%Y")
transects$month <- month(transects$Date) + (year(transects$Date)-year(startDate))*12
transects$year <- year(transects$Date)

## Aggregate transect data by shehia and month throw out unecessary columns
transects <- aggregate(.~DW+month,transects[,c("DW","month","year","Dogs.without.collars","Dogs.with.collars","Totalcounts")],sum)




## Add human population from NBS census to shapefile
#__________________

PembaWard$DW <- paste(PembaWard$District_N,PembaWard$Ward_Name,sep="_")
NBS2012$DW <- paste(NBS2012$District,NBS2012$matchWard,sep="_")
village_pops <- rowsum(NBS2012$Population,NBS2012$DW)
matchPop <- match(rownames(village_pops),PembaWard$DW)
matchPop2 <- match(PembaWard$DW,rownames(village_pops))
which(is.na(matchPop))
which(is.na(matchPop2))
PembaWard$humanPop2012NBS <- village_pops[matchPop2]



## Get human population in each ward in each month 
#__________________

## Set up population matrix
popMat <- matrix(nrow=nrow(PembaWard),ncol=months)
census_month <- month(as.Date("2012-08-15")) + (year(as.Date("2012-08-15"))-year(startDate))*12 # NBS census happened in August 2012
popMat[,census_month] <- PembaWard$humanPop2012NBS

## Project forward
forward <- (census_month+1):months
popMat[,forward] <- rep(popMat[,census_month],length(forward))*(1+monthly_growth_rate)^rep((1:length(forward))/12,each=nrow(PembaWard))

## Project backwards
backward <- (census_month-1):1
popMat[,backward] <- rep(popMat[,census_month],length(backward))*(1+monthly_growth_rate)^rep(-(1:length(backward))/12,each=nrow(PembaWard))

## Round to the nearest human
popMat <- round(popMat)
sum(popMat[,1])
sum(popMat[,months])



## Get human population matrix on grid
#__________________

## Exponential human pop growth rate in each grid cell based on Worldpop
#--------------

par(mar=c(4,4,1,3),mfrow=c(1,1))
PembaPopGrowth <- (PembaPop15/PembaPop10)^(1/5)-1
PembaPopGrowth[which(is.na(PembaPopGrowth[])&!is.na(PembaGrid[]))]<-0
plot(PembaPopGrowth) 


##Project a population distribution for 2012 using Worldpop
#--------------

PembaPop12 <- PembaPop10*(1+PembaPopGrowth)^2
plot(PembaPop12)


##Distribute humans into cells
#--------------

popMat_grid <- matrix(0,ncol=ncol(popMat),nrow=length(which(!is.na(PembaGrid[]))))
wardPops<-PembaWard$humanPop2012NBS
popMat_grid[,census_month]<- PopMap(PembaUTM, PembaWard, init=0, probMap=PembaPop12, wardPops=wardPops,
                                    missedVill=NA, missedVillCells=NA,villageIDs=PembaUTM$WardID)

## Project forward
forward <- (census_month+1):months
popMat_grid[,forward] <- rep(popMat_grid[,census_month],length(forward))*(1+monthly_growth_rate)^rep((1:length(forward))/12,each=nrow(popMat_grid))

## Project backwards
backward <- (census_month-1):1
popMat_grid[,backward] <- rep(popMat_grid[,census_month],length(backward))*(1+monthly_growth_rate)^rep(-(1:length(backward))/12,each=nrow(popMat_grid))


par(mfrow=c(1,1),mar=c(0,0,0,6))
popGrid<-PembaGrid
popGrid[which(!is.na(popGrid[]))]<-popMat_grid[,census_month]
colours <- colorRampPalette(c("white",brewer.pal(8,"YlOrRd")[2:8]))(100)
plot(PembaDist)
plot(popGrid,add=T,col=colours,breaks=seq(0,max(popGrid[],na.rm=T),length.out=100),legend=F)
plot(PembaDist,add=T)
plot(popGrid, breaks=seq(0,max(popGrid[],na.rm=T),length.out=100),
     legend.only=T, add=T,col=colours,
     legend.args=list(text=expression(paste("Humans/km"^2)), side=4, font=2, line=3.8, cex=1.2),
     axis.args=list(at=seq(0,max(popGrid[],na.rm=T),2000)),smallplot=c(0.75,0.76, .25,.75))




## Monthly dogs vaccinated in each ward
#__________________

## Create monthly vaccination matrix and add data
dogVaxMat <- matrix(0,nrow=nrow(PembaWard),ncol=months)
vax$matchWards <- match(vax$DW,PembaWard$DW)
for(i in 1:nrow(vax)){
  dogVaxMat[vax$matchWards[i],vax$month[i]] <- dogVaxMat[vax$matchWards[i],vax$month[i]] + vax$dogs_vaccinated[i]
}
vax[which(is.na(vax$matchWards)),] # 17 dogs to be distributed between Tumbe Mashariki & Tumbe Magharibi
dogVaxMat[which(PembaWard$Ward_Name=="Tumbe Mashariki"),vax[which(is.na(vax$matchWards)),"month"]] <- round(vax[which(is.na(vax$matchWards)),"dogs_vaccinated"]*PembaWard$humanPop2012NBS[which(PembaWard$Ward_Name=="Tumbe Mashariki")]/(PembaWard$humanPop2012NBS[which(PembaWard$Ward_Name=="Tumbe Mashariki")]+PembaWard$humanPop2012NBS[which(PembaWard$Ward_Name=="Tumbe Magharibi")]))
dogVaxMat[which(PembaWard$Ward_Name=="Tumbe Magharibi"),vax[which(is.na(vax$matchWards)),"month"]] <- round(vax[which(is.na(vax$matchWards)),"dogs_vaccinated"]*PembaWard$humanPop2012NBS[which(PembaWard$Ward_Name=="Tumbe Magharibi")]/(PembaWard$humanPop2012NBS[which(PembaWard$Ward_Name=="Tumbe Mashariki")]+PembaWard$humanPop2012NBS[which(PembaWard$Ward_Name=="Tumbe Magharibi")]))




## Obtain dog population from transects and government censuses
#__________________

## Dog Population matrix
dogPopMat_census_transects <- matrix(nrow=nrow(popMat),ncol=ncol(popMat))

## Add census information
dogPopMat_census_transects[cbind(census2012$ward_match,census2012$month)] <-  PembaWard$dogs2012[census2012$ward_match]# 2012
dogPopMat_census_transects[cbind(census2017_2019$ward_match[which(census2017_2019$year==2017)],census2017_2019$month[which(census2017_2019$year==2017)])] <-  PembaWard$dogs2017[census2017_2019$ward_match[which(census2017_2019$year==2017)]] # 2017
dogPopMat_census_transects[cbind(census2017_2019$ward_match[which(census2017_2019$year==2018)],census2017_2019$month[which(census2017_2019$year==2018)])] <-  PembaWard$dogs2018[census2017_2019$ward_match[which(census2017_2019$year==2018)]] # 2018
dogPopMat_census_transects[cbind(census2017_2019$ward_match[which(census2017_2019$year==2019)],census2017_2019$month[which(census2017_2019$year==2019)])] <-  PembaWard$dogs2019[census2017_2019$ward_match[which(census2017_2019$year==2019)]] # 2019

## Throw out transects with few observations
transects <- transects[which(transects$Totalcounts>=10),] 

## Throw out transects with 0 observed vaccinated dogs (lead to zero denominator problems)
transects <- transects[which(transects$Dogs.with.collars>0),]

## Estimate coverage 
transects$coverage <- transects$Dogs.with.collars/transects$Totalcounts

## Separate out years
transects2013 <- transects[which(transects$year==2013),]
transects2014 <- transects[which(transects$year==2014),]

## Get dogs known to be vaccinated in each village with a transect
transects2013$dogsVax <- NA
for(i in 1:nrow(transects2013)){ 
  if(length(which(vax$DW==transects2013$DW[i] & vax$month%in%c(37:49)))>0){
    transects2013$dogsVax[i] <- sum(vax$dogs_vaccinated[which(vax$DW==transects2013$DW[i] & vax$month%in%c(37:49))]) 
  }
}
transects2014$dogsVax <- NA
for(i in 1:nrow(transects2014)){ 
  if(length(which(vax$DW==transects2014$DW[i] & vax$month%in%c(59:60)))>0){
    transects2014$dogsVax[i] <- sum(vax$dogs_vaccinated[which(vax$DW==transects2014$DW[i] & vax$month%in%c(59:60))])  
  }
}


# some transects with no corresponding vax record - will remove here as can't use to get dogs estimate
transects_no_vax <- rbind(transects2013[which(is.na(transects2013$dogsVax)),],transects2014[which(is.na(transects2014$dogsVax)),])
write.table(transects_no_vax,"output/transects_with_no_corresponding_vaccination_records.csv",row.names = F,sep=",")
transects2013 <- transects2013[which(!is.na(transects2013$dogsVax)),] 
transects2014 <- transects2014[which(!is.na(transects2014$dogsVax)),]


## Estimate dog populations from each transect (accounting for pups)
## Get population size from transect data
transects2013$popEst <- round((transects2013$dogsVax/transects2013$coverage)*(1+PAR)) 
transects2014$popEst <- round((transects2014$dogsVax/transects2014$coverage)*(1+PAR))

## Add population estimates to new dog pop matrix
transects2013$matchDW <- match(transects2013$DW,PembaWard$DW)
dogPopMat_census_transects[cbind(transects2013$matchDW,transects2013$month)] <- transects2013$popEst
transects2014$matchDW <- match(transects2014$DW,PembaWard$DW)
dogPopMat_census_transects[cbind(transects2014$matchDW,transects2014$month)] <- transects2014$popEst


## Check for wards with no information 
check <- rep(NA,nrow(dogPopMat_census_transects))
for(i in 1:nrow(dogPopMat_census_transects)){
  check[i] <- length(which(!is.na(dogPopMat_census_transects[i,])))
}
which(check==0) # 0 villages where we have no information 
table(check) # some villages have 6 estimates


## Fill in gaps
maxVax<-1
for(i in 1:nrow(dogPopMat_census_transects)){
  
  done <- 0
  while(done==0){
    
    ## Project forwards and backwards from known population estimates to fill in rows
    ## For villages with just one population estimate available
    if(length(which(!is.na(dogPopMat_census_transects[i,])))==1){
      
      dogPop_i <- dogPopMat_census_transects[i,]
      
      ## Get HDR from the known estimate
      month <- which(!is.na(dogPopMat_census_transects[i,]))
      dogs_i <- dogPopMat_census_transects[i,month]
      humans_i <- popMat[i,month]
      HDR_i <- humans_i/dogs_i
      
      ## In empty months use HDR to convert humans to dogs
      dogPop_i[-month] <- round(popMat[i,-month]/rep(HDR_i,months-1))
      
      
    ## For villages with more than one population estimate available
    }else if(length(which(!is.na(dogPopMat_census_transects[i,])))>1){
      
      dogPop_i <- dogPopMat_census_transects[i,]
      
      ## months where estimates available and estimates themselves
      month <- which(!is.na(dogPopMat_census_transects[i,]))
      dogs_i <- dogPopMat_census_transects[i,month]
      
      ## calculate growth rates between estimates
      for(j in 1:(length(month))){
        
        ## Project between estimates
        if(j!=(length(month))){
          growth_i <- log(dogs_i[j+1]/dogs_i[j])/(month[j+1]-month[j])
          forward <- (month[j]+1):(month[j+1]-1)
          dogPop_i[forward] <- round(dogs_i[j]*exp(growth_i*(1:length(forward))))
        }
        
        ## Project beyond estimates
        if(j==length(month)){
          
          forward <- (month[j]+1):months
          
          ## Get HDR from the last estimate
          humans_i <- popMat[i,month[j]]
          HDR_i <- humans_i/dogs_i[j]
          
          ## In empty months use HDR to convert humans to dogs
          dogPop_i[forward] <- round(popMat[i,forward]/HDR_i)
          
        }
        
        ## Project back before first estimate
        if(j==1){
          backward <- (month[1]-1):1
          
          ## Get HDR from the first estimate
          humans_i <- popMat[i,month[1]]
          HDR_i <- humans_i/dogs_i[1]
          
          ## In empty months use HDR to convert humans to dogs
          dogPop_i[backward] <- round(popMat[i,backward]/HDR_i)
        }
      }
    }
    
    ## Check once vax data incorporated won't get coverages >maxVax
    vax_i <- dogVaxMat[i,]
    vc_i <- vax_i/dogPop_i
    if(length(which(vax_i>round(dogPop_i*maxVax)))==0){
      
      ##If we don't, break out of the while loop
      done<-1
      
      ##And record projected populations
      dogPopMat_census_transects[i,]<- dogPop_i
      
    ## If we do, add a new known population point based on the highest coverage month and stay in the loop 
    }else{
      dogPopMat_census_transects[i,which.max(vc_i)] <- round(vax_i[which.max(vc_i)]/maxVax) # new estimate of population size at month when vc is greatest
    }
  }  
}





## Get dog population by cell
#______________________

## Use village matrix to distribute dogs at the cell level for the first month
popGrid1<-PembaGrid
popGrid1[which(!is.na(popGrid1[]))]<-popMat_grid[,1]
dogPopMat_grid <- matrix(0,nrow=length(which(!is.na(PembaGrid[]))),ncol=ncol(dogPopMat_census_transects))
dogPopMat_grid[,1]<- PopMap(PembaUTM, PembaWard, probMap=popGrid1, wardPops=dogPopMat_census_transects[,1],
                       missedVill=NA, missedVillCells=NA,villageIDs = PembaUTM$WardID)

## Fill in for rest of months
diff <- t(diff(t(dogPopMat_census_transects)))
for(i in 2:ncol(dogPopMat_grid)){
  popGrid<-PembaGrid
  popGrid[which(!is.na(popGrid[]))]<-popMat_grid[,i]
  dogPopMat_grid[,i] <- PopMap(PembaUTM, PembaWard, init=dogPopMat_grid[,i-1], probMap=popGrid, wardPops=diff[,i-1],
                          missedVill=NA, missedVillCells=NA,villageIDs = PembaUTM$WardID)
}


 

##Save
#______________________

## human population by ward/month
write.table(popMat,file="output/humanPopMatPembaWard.csv",row.names = PembaWard$DW, col.names = F, sep=",")

## dog population by ward/month
write.table(dogPopMat_census_transects,file="output/dogPopMatPembaWard.csv",row.names = PembaWard$DW, col.names = F, sep=",")

## human population by cell/month
write.table(popMat_grid,file="output/humanPopMatPembaCell.csv",row.names = F, col.names = F, sep=",")

## dog population by cell/month
write.table(dogPopMat_grid,file="output/dogPopMatPembaCell.csv",row.names = F, col.names = F, sep=",")

## human population by cell/month
write.table(rowsum(popMat,PembaWard$District_N),file="output/humanPopMatPembaDist.csv",row.names = sort(PembaDist$District_N), col.names = F, sep=",")

## dog population by cell/month
write.table(rowsum(dogPopMat_census_transects,PembaWard$District_N),file="output/dogPopMatPembaDist.csv",row.names = sort(PembaDist$District_N), col.names = F, sep=",")


