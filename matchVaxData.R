#
## Fix vaccination data so that it matches to the village shapefile
#__________________________________

rm(list=ls())

library(lubridate)
library(rgeos)
library(rgdal)


## Vaccination data
vax <- read.csv("data/Vax_Pemba_2010_to_2019.csv",stringsAsFactors = F)

## Ward shapefile
PembaWard <- readOGR("output/GIS/PembaWard_NBS2012","PembaWard_NBS2012")

## Match districts
vax$District[which(vax$District=="Chakechake")] <- "Chake Chake"
match_dist <- match(unique(vax$District),unique(PembaWard$District_N))
unique(vax$District)[which(is.na(match_dist))]

## Match wards
vax$DW <- paste(vax$District,vax$Shehia,sep="_")
PembaWard$DW <- paste(PembaWard$District_N,PembaWard$Ward_Name,sep="_")
match_ward <- match(unique(vax$DW),unique(PembaWard$DW))
sort(unique(vax$DW)[which(is.na(match_ward))])
unique(PembaWard$DW)[which(is.na(match(PembaWard$DW,vax$DW)))]

## Some ward fixes
vax$Shehia[which(vax$DW=="Mkoani_Kukuu")] <- "Kuukuu"
vax$Shehia[which(vax$DW=="Mkoani_Chole")] <- "Kengeja" # Chole was originally a subvillage of Kengeja
vax$Shehia[which(vax$DW=="Mkoani_Mwambe")] <- "Muambe"
vax$Shehia[which(vax$DW=="Chake Chake_Ngambwa")] <- "Ng'ambwa"
vax$Shehia[which(vax$DW=="Chake Chake_Mfikia")] <- "Mfikiwa"
vax$Shehia[which(vax$DW=="Chake Chake_Mjini ole"|vax$DW=="Chake Chake_Mjini Ole")]  <- "Mjini Ole"; vax$District[which(vax$DW=="Chake Chake_Mjini ole"|vax$DW=="Chake Chake_Mjini Ole")]  <- "Wete"
vax$District[which(vax$DW=="Chake Chake_Ole")] <- "Wete"
vax$District[which(vax$DW=="Wete_Finya")] <- "Micheweni"
vax$District[which(vax$DW=="Wete_Kinyasini")] <- "Micheweni"
vax$District[which(vax$DW=="Wete_Mtemani")] <- "Micheweni"
vax$District[which(vax$DW=="Wete_Mgogoni")] <- "Micheweni"
vax$Shehia[which(vax$DW=="Chake Chake_Tibirizi")] <- "Tibirinzi"
vax$Shehia[which(vax$DW=="Chake Chake_Gombani")] <- "Ng'ambwa" # Chole is a village in Ng'ambwa ward
vax$Shehia[which(vax$DW=="Micheweni_Muhogoni")] <- "Mihogoni"
vax$Shehia[which(vax$DW=="Micheweni_Kipange")] <- "Konde" # Kipange was a village in Konde
vax$Shehia[which(vax$DW=="Wete_Mjini kiuyu")] <- "Kiuyu Minungwini"
vax$Shehia[which(vax$DW=="Micheweni_Kiuyu")] <- "Kiuyu Mbuyuni"
vax$Shehia[which(vax$DW=="Wete_Selemu")] <- "Selem"
vax$Shehia[which(vax$DW=="Wete_Mzambarau takao"|vax$DW=="Wete_Mzambarau Takao")] <- "Mzambarauni Takao"
vax$Shehia[which(vax$DW=="Wete_Mzambarauni")] <- "Mzambarauni Takao"
vax$Shehia[which(vax$DW=="Wete_Chanjaani")] <- "Mtambwe Kaskazini"
vax$Shehia[which(vax$DW=="Wete_Mlindo")] <- "Pandani" # not completely sure about this one
vax$Shehia[which(vax$DW=="Micheweni_Msuka  Magharibi")] <- "Msuka Magharibi"
vax$Shehia[which(vax$DW=="Micheweni_Shumba Vyamboni")] <- "Shumba Viamboni"
vax$Shehia[which(vax$DW=="Micheweni_Njuguni")] <- "Wingwi Njuguni"
vax$Shehia[which(vax$DW=="Wete_Kiuyu kigongoni")] <- "Kiuyu Kigongoni"
vax$Shehia[which(vax$DW=="Chake Chake_Kilindini")] <- "Kilindi"
vax$Shehia[which(vax$DW=="Chake Chake_Kumvini")] <- "Pujini" # not sure
vax$Shehia[which(vax$DW=="Chake Chake_Ng'mbwa")] <- "Ng'ambwa"
vax$Shehia[which(vax$DW=="Chake Chake_Mchungwani")] <- "Michungwani"
vax$District[which(vax$DW=="Mkoani_Dodo")] <- "Chake Chake"
vax$Shehia[which(vax$DW=="Micheweni_Mapofu")] <- "Wingwi Mapofu"
vax$Shehia[which(vax$DW=="Micheweni_Mjiniwingwi")] <- "Mjini Wingwi"
vax$Shehia[which(vax$DW=="Micheweni_Kipangani")] <- "Kifundi"
vax$Shehia[which(vax$DW=="Mkoani_Kengeja ")] <- "Kengeja"
vax$Shehia[which(vax$DW=="Chake Chake_Mchanga mrima"|vax$DW=="Chake Chake_Mchangamrima")] <- "Mchanga Mdogo"; vax$District[which(vax$DW=="Chake Chake_Mchanga mrima"|vax$DW=="Chake Chake_Mchangamrima")] <- "Wete" # this might be pushing it a little, not very convinced it's right
vax$Shehia[which(vax$DW=="Mkoani_Mkanyangeni")] <- "Mkanyageni"
vax$Shehia[which(vax$DW=="Mkoani_Mchakwe")] <- "Muambe"
vax$Shehia[which(vax$DW=="Wete_Kinyakani")] <- "Kinyikani"
vax$Shehia[which(vax$DW=="Micheweni_Wingwi Mtemani")] <- "Mtemani"
vax$Shehia[which(vax$DW=="Chake Chake_Kibirinzi")] <- "Tibirinzi"
vax$District[which(vax$DW=="Chake Chake_Kipangani")] <- "Wete"
vax$Shehia[which(vax$DW=="Mkoani_Shibu")] <- "Shidi"
vax$Shehia[which(vax$DW=="Wete_Mjananza")] <- "Wingwi Mjananza"; vax$District[which(vax$DW=="Wete_Mjananza")] <- "Micheweni" # not sure


## Match wards again
vax$DW <- paste(vax$District,vax$Shehia,sep="_")
PembaWard$DW <- paste(PembaWard$District_N,PembaWard$Ward_Name,sep="_")
match_ward <- match(unique(vax$DW),unique(PembaWard$DW))
unique(vax$DW)[which(is.na(match_ward))] #One record for 17 dogs in Tumbe in 2016 - not clear whether east or west - probably just have to divide between
unique(PembaWard$DW)[which(is.na(match(PembaWard$DW,vax$DW)))] # 3 Shehia in shapefile with no vaccination data: "Mkoani_Makoongwe"    "Mkoani_Kisiwa Panza" "Mkoani_Shamiani"

## Remove records where no dogs were vaccinated (or with NAs)
vax<-vax[which(vax$dogs_vaccinated>0),]

## Check for duplicates
duplicates<-vax[which(duplicated(vax[,c("vacc_date","dogs_vaccinated","Cats_vaccinated","subvillage","DW")]) | duplicated(vax[nrow(vax):1,c("vacc_date","dogs_vaccinated","Cats_vaccinated","subvillage","DW")])[nrow(vax):1]),]
duplicates <- duplicates[order(duplicates$DW),] # lots in 2019, but this is expected given the altered vaccination campaign structure

## Save cleaned vaccination data
write.csv(vax,paste("output/VaxPemba.csv",sep=""),row.names = F)

