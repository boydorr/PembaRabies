rm(list=ls())

library(maptools)
library(rgeos)
library(rgdal)
library(stringr)
library(lubridate)

options(stringsAsFactors = F)
set.seed(0)


## Census data
census2012 <- read.csv("data/dogPopCensus2012.csv")
census2017_2019 <- read.csv("data/Vax_Pemba_2010_to_2019.csv"); census2017_2019$vacc_date[which(census2017_2019$vacc_date=="2017|2018|2019")]

## Ward shapefile
PembaWard <- readOGR(paste("output/GIS/PembaWard_NBS2012",sep=""), paste("PembaWard_NBS2012",sep=""))
PembaWard$DW <- paste(PembaWard$District_N,PembaWard$Ward_Name,sep="_")

## Dog Vaccination data
vax <- read.csv("output/VaxPemba.csv")



## Match dog census data to shapefile
#__________________

## 2012 data
#----------

## Match district names
census2012$DISTRICT[which(census2012$DISTRICT=="Chakechake")] <- "Chake Chake"
dist_match <- match(unique(census2012$DISTRICT),unique(PembaWard$District_N))
unique(census2012$DISTRICT)[which(is.na(dist_match))] #all matching

## Match ward names
PembaWard$DW <- paste(PembaWard$District_N,PembaWard$Ward_Name,sep="_")
census2012$DW <- paste(census2012$DISTRICT,census2012$VILLAGE.SHEHIA, sep="_")
ward_match <- match(PembaWard$DW,census2012$DW)
census2012$DW[which(is.na(match(census2012$DW,PembaWard$DW)))]
PembaWard$DW[which(is.na(ward_match))]  #all matching
census2012$ward_match <- match(census2012$DW,PembaWard$DW)


## Add 2012 human population to shapefile
census2012$totalHumans <- census2012$HumanFEMALEs_vill + census2012$HumanMALES_vill
PembaWard$humanPop12 <- census2012$totalHumans[ward_match]

## Add 2012 dog population to shapefile
census2012$DogsStray_estimate_mean <- (census2012$DogsStray_estimate1 + census2012$DogsStray_estimate2)/2
census2012$totalDogs <- census2012$DogsOwned + census2012$DogsStray_estimate_mean
census2012$totalDogs[which(is.na(census2012$totalDogs))] <- census2012$DogsOwned[which(is.na(census2012$totalDogs))]
PembaWard$dogs2012 <- census2012$totalDogs[ward_match]

## add month of census (assume this was same as date of vaccination?)
vax$dateVaccination <- as.Date(vax$vacc_date,format="%d/%m/%Y")
vax$month <- month(vax$dateVaccination) + (year(vax$dateVaccination)-2010)*12
vax$year <- year(vax$dateVaccination)
census2012$month <- NA
for(i in 1:nrow(census2012)){
  vax_i <- vax[which(vax$DW==census2012$DW[i] & vax$year==2012),]
  # print(nrow(vax_i))
  # if(nrow(vax_i)>1){print(vax_i)}
  if(nrow(vax_i)>0){census2012$month[i] <- vax_i$month[1]}
}
table(census2012$month)
census2012$month[which(is.na(census2012$month))] <- 30 #assume missing campaign dates were in june

## Missing dog population data
length(which(is.na(PembaWard$dogs2012))) #19 wards without a dog population estimate
PembaWard$DW[which(is.na(PembaWard$dogs2012))]



## 2017-2019 data
#----------

## Date info
census2017_2019$vacc_date <- as.Date(census2017_2019$vacc_date,format="%d/%m/%Y")
census2017_2019$year <- year(census2017_2019$vacc_date)
census2017_2019$month <- month(census2017_2019$vacc_date) + (year(census2017_2019$vacc_date)-2010)*12

## Throw out data for years we don't need
census2017_2019 <- census2017_2019[which(census2017_2019$year %in% c(2017:2019)),]

## Match district names
census2017_2019$District[which(census2017_2019$District=="Chakechake")] <- "Chake Chake"
dist_match <- match(unique(census2017_2019$District),unique(PembaWard$District_N))
unique(census2017_2019$District)[which(is.na(dist_match))] #all matching

## Match ward names
census2017_2019$DW <- paste(census2017_2019$District,census2017_2019$Shehia, sep="_")
ward_match <- match(PembaWard$DW,census2017_2019$DW)
sort(unique(census2017_2019$DW[which(is.na(match(census2017_2019$DW,PembaWard$DW)))])) #not all matching

## Some ward fixes
census2017_2019$Shehia[which(census2017_2019$DW=="Chake Chake_Gombani")] <- "Ng'ambwa"
census2017_2019$Shehia[which(census2017_2019$DW=="Chake Chake_Kibirinzi"|census2017_2019$DW=="Chake Chake_Tibirizi")] <- "Tibirinzi"
census2017_2019$Shehia[which(census2017_2019$DW=="Chake Chake_Kilindini")] <- "Kilindi"
census2017_2019$Shehia[which(census2017_2019$DW=="Chake Chake_Kumvini")] <- "Pujini"
census2017_2019$Shehia[which(census2017_2019$DW=="Chake Chake_Mchanga mrima"|census2017_2019$DW=="Chake Chake_Mchangamrima")] <- "Mchanga Mdogo"; census2017_2019$District[which(census2017_2019$DW=="Chake Chake_Mchanga mrima"|census2017_2019$DW=="Chake Chake_Mchangamrima")] <- "Wete"
census2017_2019$Shehia[which(census2017_2019$DW=="Chake Chake_Mfikia")] <- "Mfikiwa"
census2017_2019$Shehia[which(census2017_2019$DW=="Chake Chake_Mchungwani")] <- "Michungwani"
census2017_2019$Shehia[which(census2017_2019$DW=="Chake Chake_Mjini ole"|census2017_2019$DW=="Chake Chake_Mjini Ole")]  <- "Mjini Ole"; census2017_2019$District[which(census2017_2019$DW=="Chake Chake_Mjini ole"|census2017_2019$DW=="Chake Chake_Mjini Ole")]  <- "Wete"
census2017_2019$Shehia[which(census2017_2019$DW=="Chake Chake_Ngambwa"|census2017_2019$DW=="Chake Chake_Ng'mbwa")] <- "Ng'ambwa"
census2017_2019$District[which(census2017_2019$DW=="Chake Chake_Ole")] <- "Wete"
census2017_2019$Shehia[which(census2017_2019$DW=="Micheweni_Kipange")] <- "Konde"
census2017_2019$Shehia[which(census2017_2019$DW=="Micheweni_Kiuyu")] <- "Kiuyu Mbuyuni"
census2017_2019$Shehia[which(census2017_2019$DW=="Micheweni_Msuka  Magharibi")] <- "Msuka Magharibi"
census2017_2019$Shehia[which(census2017_2019$DW=="Micheweni_Muhogoni")] <- "Mihogoni"
census2017_2019$Shehia[which(census2017_2019$DW=="Micheweni_Njuguni")] <- "Wingwi Njuguni"
census2017_2019$Shehia[which(census2017_2019$DW=="Micheweni_Shumba Vyamboni")] <- "Shumba Viamboni"
census2017_2019$Shehia[which(census2017_2019$DW=="Mkoani_Chole")] <- "Kengeja"
census2017_2019$District[which(census2017_2019$DW=="Mkoani_Dodo")] <- "Chake Chake"
census2017_2019$Shehia[which(census2017_2019$DW=="Mkoani_Kukuu")] <- "Kuukuu"
census2017_2019$Shehia[which(census2017_2019$DW=="Mkoani_Mwambe")] <- "Muambe"
census2017_2019$Shehia[which(census2017_2019$DW=="Wete_Chanjaani")] <- "Mtambwe Kaskazini"
census2017_2019$District[which(census2017_2019$DW=="Wete_Finya")] <- "Micheweni"
census2017_2019$District[which(census2017_2019$DW=="Wete_Kinyasini")] <- "Micheweni"
census2017_2019$Shehia[which(census2017_2019$DW=="Wete_Kiuyu kigongoni")] <- "Kiuyu Kigongoni"
census2017_2019$District[which(census2017_2019$DW=="Wete_Mgogoni")] <- "Micheweni"
census2017_2019$Shehia[which(census2017_2019$DW=="Wete_Mjini kiuyu")] <- "Kiuyu Minungwini"
census2017_2019$Shehia[which(census2017_2019$DW=="Wete_Mlindo")] <- "Pandani"
census2017_2019$District[which(census2017_2019$DW=="Wete_Mtemani")] <- "Micheweni"
census2017_2019$Shehia[which(census2017_2019$DW=="Wete_Mzambarau Takao")] <- "Mzambarauni Takao"
census2017_2019$Shehia[which(census2017_2019$DW=="Wete_Mzambarauni")] <- "Mzambarauni Takao"
census2017_2019$Shehia[which(census2017_2019$DW=="Wete_Selemu")] <- "Selem"
census2017_2019$District[which(census2017_2019$DW=="Chake Chake_Kipangani")] <- "Wete"
census2017_2019$Shehia[which(census2017_2019$DW=="Mkoani_Shibu")] <- "Shidi"
census2017_2019$Shehia[which(census2017_2019$DW=="Wete_Mjananza")] <- "Wingwi Mjananza"; census2017_2019$District[which(census2017_2019$DW=="Wete_Mjananza")] <- "Micheweni" # not sure
census2017_2019$Shehia[which(census2017_2019$DW=="Mkoani_Mkanyangeni")] <- "Mkanyageni"
census2017_2019$Shehia[which(census2017_2019$DW=="Mkoani_Mchakwe")] <- "Muambe"
census2017_2019$Shehia[which(census2017_2019$DW=="Wete_Kinyakani")] <- "Kinyikani"
census2017_2019$Shehia[which(census2017_2019$DW=="Micheweni_Wingwi Mtemani")] <- "Mtemani"
census2017_2019$Shehia[which(census2017_2019$DW=="Micheweni_Mapofu")] <- "Wingwi Mapofu"


## Match ward names again
census2017_2019$DW <- paste(census2017_2019$District,census2017_2019$Shehia, sep="_")
ward_match <- match(PembaWard$DW,census2017_2019$DW)
sort(unique(census2017_2019$DW[which(is.na(match(census2017_2019$DW,PembaWard$DW)))])) #all matching
PembaWard$DW[which(is.na(ward_match))] # don't have information for 5 Shehia
census2017_2019$ward_match <- match(census2017_2019$DW,PembaWard$DW)

## Multiple records for many Shehia - aggregate these.
census2017_2019 <- aggregate(.~DW+month+year+ward_match,census2017_2019[,c("DW","year","month","ward_match","dogs_vaccinated","total_dog_estimate")],sum)
length(which(duplicated(census2017_2019[,c("DW","year")])))
duplicates<-census2017_2019[which(duplicated(census2017_2019[,c("DW","year")]) | duplicated(census2017_2019[nrow(census2017_2019):1,c("DW","year")])[nrow(census2017_2019):1]),]
# no cases where duplicates were more than a month apart
for(i in 1:nrow(duplicates)){
  census2017_2019$month[as.numeric(row.names(duplicates)[i])] <-
    min(duplicates$month[which(duplicates$DW==duplicates$DW[i] & duplicates$year==duplicates$year[i])])
}
census2017_2019 <- aggregate(.~DW+month+year+ward_match,census2017_2019[,c("DW","year","month","ward_match","dogs_vaccinated","total_dog_estimate")],sum)

## Add 2017 dog population to shapefile
PembaWard$dogs2017<-NA
PembaWard$dogs2017[census2017_2019$ward_match[which(census2017_2019$year==2017)]] <- census2017_2019$total_dog_estimate[which(census2017_2019$year==2017)]

## Add 2018 dog population to shapefile
PembaWard$dogs2018<-NA
PembaWard$dogs2018[census2017_2019$ward_match[which(census2017_2019$year==2018)]] <- census2017_2019$total_dog_estimate[which(census2017_2019$year==2018)]

## Add 2019 dog population to shapefile
PembaWard$dogs2019<-NA
PembaWard$dogs2019[census2017_2019$ward_match[which(census2017_2019$year==2019)]] <- census2017_2019$total_dog_estimate[which(census2017_2019$year==2019)]



## 550 stray dogs estimated for the island in 2017-2018, assume these are
## distributed in prop with the human population
non_zero_dogs <- which(!PembaWard$DW%in%c("Mkoani_Makoongwe", "Mkoani_Kisiwa Panza", "Mkoani_Shamiani"  ))
strays <- rmultinom(1,550,prob=PembaWard$humanPop12[non_zero_dogs])
PembaWard$dogs2017[non_zero_dogs] <- as.numeric(PembaWard$dogs2017[non_zero_dogs] + strays)
PembaWard$dogs2018[non_zero_dogs] <- as.numeric(PembaWard$dogs2018[non_zero_dogs] + strays)
PembaWard$dogs2019[non_zero_dogs] <- as.numeric(PembaWard$dogs2019[non_zero_dogs] + strays)



## Save outputs
#__________________

writeOGR(PembaWard, dsn="output/GIS/PembaWard_NBS2012", "PembaWard_NBS2012", driver="ESRI Shapefile",overwrite_layer=T)
write.csv(census2012,file="output/dogPopCensus2012_matched.csv",row.names = F)
write.csv(census2017_2019,file="output/dogPopcensus2017_2019_matched.csv",row.names = F)


