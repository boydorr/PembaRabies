
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

## Ward shapefile
PembaWard <- readOGR(paste("output/GIS/PembaWard_NBS2012",sep=""), paste("PembaWard_NBS2012",sep=""))

## District shapefile
PembaDist <- readOGR(paste("output/GIS/PembaDist_NBS2012",sep=""), paste("PembaDist_NBS2012",sep=""))

## Pemba Grid
cell_size <- 1
PembaGrid <- raster(paste("output/GIS/",cell_size^2,"kmsqGrids/PembaGrid",cell_size^2,"kmsq.grd",sep=""))

## Dog Vaccination data
vax <- read.csv("output/VaxPemba.csv")


## dog population by ward/month
dogPopWardMat <- as.matrix(read.csv("output/dogPopMatPembaWard.csv",row.names = 1,header=F))

## human population by ward/month
humanPopWardMat <- as.matrix(read.csv("output/humanPopMatPembaWard.csv",row.names = 1,header=F))

## Maximum vaccination coverage
PAR <- 0.2562058 #pup/adult ratio
propPups <- 0.2039521 #proportion of dogs that are pups
maxVax <- 1-propPups # maximum vaccination coverage can reach 




## Time info
#__________________

## Time period of interest
years <- 2010:2020
months <- 12*length(years)
startDate <- as.Date(paste0(min(years),"-01-01"))
endDate <- as.Date(paste0(max(years),"-01-01"))

## add month to vax data
vax$dateVaccination <- as.Date(vax$vacc_date,format="%d/%m/%Y")
vax$month <- month(vax$dateVaccination) + (year(vax$dateVaccination)-year(startDate))*12
vax$year <- year(vax$dateVaccination)




## Monthly dogs vaccinated in each ward
#__________________

## Create monthly vaccination matrix and add data
dogVaxMat <- matrix(0,nrow=nrow(PembaWard),ncol=months)
vax$matchWards <- match(vax$DW,PembaWard$DW)
for(i in 1:nrow(vax)){
  dogVaxMat[vax$matchWards[i],vax$month[i]] <- dogVaxMat[vax$matchWards[i],vax$month[i]] + vax$dogs_vaccinated[i]
}
vax[which(is.na(vax$matchWards)),]
dogVaxMat[which(PembaWard$Ward_Name=="Tumbe Mashariki"),vax[which(is.na(vax$matchWards)),"month"]]  <- round(vax[which(is.na(vax$matchWards)),"dogs_vaccinated"]*humanPopWardMat[which(PembaWard$Ward_Name=="Tumbe Mashariki"),12+8]/(humanPopWardMat[which(PembaWard$Ward_Name=="Tumbe Mashariki"),12+8]+humanPopWardMat[which(PembaWard$Ward_Name=="Tumbe Magharibi"),12+8]))
dogVaxMat[which(PembaWard$Ward_Name=="Tumbe Magharibi"),vax[which(is.na(vax$matchWards)),"month"]]  <- round(vax[which(is.na(vax$matchWards)),"dogs_vaccinated"]*humanPopWardMat[which(PembaWard$Ward_Name=="Tumbe Magharibi"),12+8]/(humanPopWardMat[which(PembaWard$Ward_Name=="Tumbe Mashariki"),12+8]+humanPopWardMat[which(PembaWard$Ward_Name=="Tumbe Magharibi"),12+8]))




## Monthly & annual vaccination coverage in each ward
#__________________

## Get vaccination coverages in each month
vcVill <- dogVaxMat/dogPopWardMat

## Check values
length(which(vcVill>0))
length(which(vcVill>1))/length(which(vcVill>0)) 
length(which(vcVill>maxVax))/length(which(vcVill>0)) 
hist(vcVill[which(vcVill>0)],breaks="fd")
max(vcVill[which(dogPopWardMat!=0)],na.rm=T)
vcVill[which(is.nan(vcVill))] <- 0 #incase any villages have zero dogs


## Getting annual coverage by summing coverages over year isn't going to work if 2013 campaign leaked in to 2014 because the same dogs may have been vaccinated twice
## To get annual coverages for plotting, assume that Jan/Feb 2014 vaccinations actually happened in 2013
vcVillYear <- matrix(0,nrow=nrow(PembaWard),ncol=ncol(vcVill)/12)
for(i in 1:ncol(vcVillYear)){
  if(i==3){ # 2013
    vcVillYear[,i]<-rowSums(vcVill[,37:50])
  }else if(i==4){ # 2014
    vcVillYear[,i]<-rowSums(vcVill[,51:60])
  }else{
    vcVillYear[,i]<-rowSums(vcVill[,(1:12)+12*(i-1)])
  }
}
# vcVillYear[which(vcVillYear>maxVax)]<-maxVax
vcVillYear[which(vcVillYear>1)]<-1
hist(vcVillYear) 

## Vaccination Rounds
vacc_rounds <- vcVillYear[,which(colSums(vcVillYear)>0)]


## Save files
write.table(vcVill,file="output/vcMonthPembaWards.csv",row.names = PembaWard$DW,col.names=F,sep=",")
write.table(vacc_rounds,file="output/vcRoundsPembaWards.csv",row.names = PembaWard$DW,col.names=F,sep=",")
write.table(vcVillYear,file="output/vcYearPembaWards.csv",row.names = PembaWard$DW,col.names=F,sep=",")



## Annual vaccination coverage in each district
#__________________

## dogs vaccinated in each district each month
vaxDistMat <- rowsum(dogVaxMat,PembaWard$District_N)

## dog population in each district each month
dogPopDistMat <- rowsum(dogPopWardMat,PembaWard$District_N)

## vaccination coverages by district
vcDistMat <- vaxDistMat/dogPopDistMat

## vaccination coverage by district and year
vcDistYear <- matrix(0,nrow=nrow(vcDistMat),ncol=ncol(vcVill)/12,dimnames = list(rownames(vcDistMat),NULL))
for(i in 1:ncol(vcDistYear)){
  if(i==3){ # 2013
    vcDistYear[,i]<-rowSums(vcDistMat[,37:50])
  }else if(i==4){ # 2014
    vcDistYear[,i]<-rowSums(vcDistMat[,51:60])
  }else{
    vcDistYear[,i]<-rowSums(vcDistMat[,(1:12)+12*(i-1)])
  }
}

## Save file
write.table(vcDistYear,file="output/vcYearPembaDist.csv",col.names=F,sep=",")




## Monthly waning coverage 
#__________________

## Parameters
vaccine_duration <- 3 #years
annual_death_rate <- 0.448375 # from Anna's data
lambda <- exp(-(1/vaccine_duration - log(1-annual_death_rate))/12) # annual death rate is difference eq version (r) - need to convert to k (differential equation version)

## Set up matrix
vc_waning <- vcVill*dogPopWardMat

## Estimate number of dogs that remain protected each month based on pop turnover  immunity loss
## Assume no dog is re-vaccinated the same year (except in 2014, where dogs vaccinated in Jan can be revax later)
years <- rep(1:(months/12),each=12) # group columns by year
years[49:50]<-years[48] # 2013 vaccinations extended into 2014
previous_years_vax <- rep(0,nrow(vc_waning)) # initialise vector describing number of vaccinated dogs that were not vaccinated in the current year
for(i in 2:ncol(vc_waning)){ # each month
  if(years[i]>years[i-1]){previous_years_vax <- vc_waning[,i-1]} # all the currently vaxed dogs were vaxed in previous years
  popChange <- dogPopWardMat[,i]/dogPopWardMat[,i-1]
  popChange[which(is.nan(popChange))] <- 1
  popChange <- pmin(1,popChange)
  vax_this_month <- vc_waning[,i] # number of dogs to be vaccinated this month
  vc_waning[,i] <- pmin(vc_waning[,i-1]*lambda*popChange + # vaccinated dogs last month + waning
                          pmax(0,vax_this_month-previous_years_vax*lambda*popChange), # + any new vaccinations to susceptibles (dogs vax in previous years vaccinated preferentially)
                        dogPopWardMat[,i]) # can't exceed dog population
  previous_years_vax <- pmax(0, previous_years_vax*lambda*popChange - vax_this_month) # preferentially vaccinate dogs that were vaccinated before (but not this year!)
}


## Estimate vaccination coverage per month per ward
vc_waningVill<-vc_waning/dogPopWardMat
vc_waningVill[which(is.nan(vc_waningVill))] <- NA

## Waning coverage on whole island
vc_waningTotal<-colSums(vc_waning)/colSums(dogPopWardMat)

## Waning coverage in each district
districts <- (unique(PembaWard$District_N))
vaxWaningDist <-rowsum(vc_waning,PembaWard$District_N)
vcWaningDist <- vaxWaningDist/dogPopDistMat

## Save files
write.table(vc_waningTotal,file="output/dogVCWaningPemba.csv",row.names = F, col.names=F,sep=",")
write.table(vc_waning,file="output/dogVaxWaningPembaWard.csv",row.names = PembaWard$DW,col.names = F,sep=",")
write.table(vc_waningVill,file="output/dogVCWaningPembaWard.csv",row.names = PembaWard$DW,col.names = F,sep=",")
write.table(vaxWaningDist,file="output/dogVaxWaningPembaDist.csv",row.names = rownames(vcWaningDist),col.names = F,sep=",")
write.table(vcWaningDist,file="output/dogVCWaningPembaDist.csv",row.names = rownames(vcWaningDist),col.names = F,sep=",")

 


# ## Create data frame containing annual village information
# #__________________

## Sum vaccinated dogs by year
dogVaxYear <- matrix(0,nrow=nrow(PembaWard),ncol=ncol(vcVill)/12)
for(i in 1:ncol(dogVaxYear)){
  if(i==3){ # 2013
    dogVaxYear[,i]<-rowSums(dogVaxMat[,37:50]) # 2013 campaign continued into Feb 2014, so add dogs vaccinated in Jan/Feb 2014 to 2013
  }else if(i==4){ # 2014
    dogVaxYear[,i]<-rowSums(dogVaxMat[,51:60])
  }else{
    dogVaxYear[,i]<-rowSums(dogVaxMat[,(1:12)+12*(i-1)])
  }
}


## Exposures and suspect animals by year

## human exposures and animal cases (de-identified)
human <- readRDS(file = "output/humans_deid.rda")
animalCases <- readRDS(file = "output/rabid_deid.rda")

# Separate the healthy and suspected bites
humanExposures <- human[which(human$Rabid == "Yes" | human$Rabid == "Unknown"),]

## Match ward names for exposures to shapefile
humanExposures$DW <- paste(humanExposures$District,humanExposures$Ward,sep="_")
matchWard <- match(humanExposures$DW,PembaWard$DW)
which(is.na(matchWard))
humanExposures$matchShapefile <- matchWard

## Match ward names to shapefile
animalCases$DW <- paste(animalCases$District,animalCases$Ward,sep="_")
matchWard <- match(animalCases$DW,PembaWard$DW)
which(is.na(matchWard))
animalCases$matchShapefile <- matchWard



# Subset data for same time period as Vax data
humanExposures <- humanExposures[which(humanExposures$Date.bitten >= startDate),]
animalCases <- animalCases[which(animalCases$Symptoms.started >= startDate),]

 exposures_year <- matrix(0,nrow=nrow(PembaWard),ncol=ncol(dogVaxMat)/12)
 cases_year <- matrix(0,nrow=nrow(PembaWard),ncol=ncol(dogVaxMat)/12)
 for(i in 1:ncol(exposures_year)){
   for(j in 1:nrow(exposures_year)){
     exposures_year[j,i]<-length(which(year(humanExposures$Date.bitten)==years[i] & humanExposures$matchShapefile==j))
     cases_year[j,i]<-length(which(year(animalCases$Symptoms.started)==years[i] & animalCases$matchShapefile==j))
   }
 }


# ## Create and save data frame
data <- data.frame("district"=rep(PembaWard$District_N,length(years)),
                   "shehia"=rep(PembaWard$Ward_Name,length(years)),
                   "d_s"=rep(PembaWard$DW,length(years)),
                   "year"=rep(years,each=nrow(PembaWard)),
                   "dog_pop_est"=c(dogPopWardMat[,seq(7,ncol(dogPopWardMat),12)]), # use value from July of each year
                   "dogs_vacc"=c(dogVaxYear),
                   "vc_est"=c(vcVillYear),
                   "vc_waning_est"=c(vc_waningVill[,seq(7,ncol(dogPopWardMat),12)]), # use value from July of each year
                   "human_exposures"=c(exposures_year),
                   "animal_cases"=c(cases_year))
write.table(data, file="output/VillageData.csv", sep=",", row.names=F)




