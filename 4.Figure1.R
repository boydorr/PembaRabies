rm(list=ls())

## load all libraries
library(maptools)
library(rgeos)
library(RColorBrewer)
library(scales)
library(rgdal)
#library(GISTools)
library(gplots)
library(dplyr)
library(lubridate)
library(ggplot2)
library(stringr)
library(raster)
library(tidyr)
library(sf)
library(ggsn)
library(ggpubr)
library(grid)

# Set global options
options(stringsAsFactors=F)

# Set colour scheme
# col_pal = c("Cases"="#4A76E7", "Exposures"="#f0b928", "Deaths"="#D22254") # black", "#012583", "#E0AA1E")
col_pal = c("Cases"="#f0b928", "Exposures"="#D22254", "Deaths"="black") # black", "#012583", "#E0AA1E")

# Load functions
source("R/elapsed_months.R")
source("R/every_nth.R")
source("R/g_legend.R")

# Load raster data
cell_size <- 1
PembaGrid <- raster(paste("Output/GIS/",cell_size^2,"kmsqGrids/PembaGrid",cell_size^2,"kmsq.grd",sep=""))

# Load dog population by ward/month
dogPopWardMat <- read.csv("Output/dogPopMatPembaWard.csv", row.names=1, header=F)

# Load dog pop by cell/month (1km cell?)
dogPopWardMat_cell <- as.matrix(read.csv("Output/dogPopMatPembaCell.csv",header=F))

# Load human exposures and animal cases
human <- readRDS(file = "output/humans_deid.rda")
animalCases <- readRDS(file = "output/rabid_deid.rda")

# Load vaccination data
vcWardYear <- as.matrix(read.csv("Output/vcYearPembaWards.csv", row.names=1, header=F))

# Load waning coverage at island, district and ward level
vcWaningTotal <- as.matrix(read.csv("Output/dogVCWaningPemba.csv", header=F))
vcWaningWard <- as.matrix(read.csv("Output/dogVCWaningPembaWard.csv", row.names=1, header=F))
vcWaningDist <- as.matrix(read.csv("Output/dogVCWaningPembaDist.csv", row.names=1, header=F))

# Load clinic locations
clinic_locs <- read.csv("data/Clinic_GPS.csv", stringsAsFactors=FALSE)

# Load shapefiles
PembaWard <- read_sf("Output/GIS/PembaWard_NBS2012/PembaWard_NBS2012.shp")
PembaDist <- read_sf("Output/GIS/PembaDist_NBS2012/PembaDist_NBS2012.shp")
Tz_region <- read_sf("data/GIS/TZ_Region_2012/TZ_Region_2012.shp")

#----- Setup for script --------------------------------------------------------

#  Set start and end dates
startDate <- as.Date("2010-01-01")
endDate <- as.Date ("2021-12-31")

#  Set vaccination start date
vaxStartDate <- startDate # Vax starts in 2010
vaxEndDate <- as.Date("2019-12-31") # Vax data only until 2019, but other data runs until end of 2021

# Set vector of months of STUDY period, and calculate number of months
months_v <- seq(from=startDate, to=endDate, by="month")
n_months <- length(months_v)

# Set vector of months of VAX period, and calculate number of months
months_v_vax <- seq(from=vaxStartDate, to=vaxEndDate, by="month")
n_months_vax <- length(months_v_vax)

# Time period of interest
ys = year(startDate):year(endDate)
ys_vax = year(vaxStartDate):year(vaxEndDate)

#----- Process data ------------------------------------------------------------

# Create a new column stating if person attended hospital
human$attended_hospital <- ifelse(human$PEP.1=="true" | human$Wound.washed.at.hospital=="true" |
                                    human$Source.Hospital=="true" | human$Source.Health.facility=="true" |
                                    human$Source.Pharmacy=="true", "true", "false")
table(human$attended_hospital)

# Separate the health and suspected bites
humanAttendedHos <- human[which(human$attended_hospital=="true"),]
humanExposures <- human[which(human$Rabid == "Yes" & human$Patient.outcome != "Died"),] # Exclude those that died (counted seperately)
humanDeaths <- human[which(human$Rabid == "Yes" & human$Cause.of.death == "Rabies"),]

# Subset data for set time period
humanAttendedHos <- humanAttendedHos[which(humanAttendedHos$Date.bitten >= startDate),]; nrow(humanAttendedHos)
humanExposures <- humanExposures[which(humanExposures$Date.bitten >= startDate),]; nrow(humanExposures)
humanDeaths <- humanDeaths[which(humanDeaths$Date.bitten >= startDate),]; nrow(humanDeaths)
animalCases <- animalCases[which(animalCases$Symptoms.started >= startDate),]; nrow(animalCases)

# Get time series of bite victims and suspected animal rabies cases
ts_breaks <- seq(0, n_months, 1)
ts_suspectBites <- hist(humanExposures$month_bitten, breaks = ts_breaks, plot=F)
ts_humanDeaths <- hist(humanDeaths$month_bitten, breaks = ts_breaks, plot=F)
ts_suspectAnimals <- hist(animalCases$month_symp, breaks = ts_breaks, plot=F)
ts_allBites <- hist(humanAttendedHos$month_bitten, breaks = ts_breaks, plot=F)

# To avoid showing on map, change 0's to -5
# ts_humanDeaths$counts[which(ts_humanDeaths$counts==0)] = -5

# Subset waning totals to Vax period
vcWaningTotal_sub <- vcWaningTotal[1:n_months_vax]
vcWaningWard_sub <- vcWaningWard[,1:n_months_vax]
vcWaningDist_sub <- vcWaningDist[,1:n_months_vax]

# Calculates the waning vaccination coverage for the selected month
selected_month = 12
selected_months = seq(from=selected_month, to=n_months_vax, by=12) # Collect vector of selected months
vcWaningWard_Yr <- vcWaningWard_sub[,selected_months] # Subset matrix

#----- Panel A -----------------------------------------------------------------

# formatted much like: Pemba/figs/waningVax_cases_exp_Pemba.pdf, but human exposures
# (blue) and deaths (black) stacked, and overlaid with a line for bite patient presentations.
# The x year labels can be removed as this will be lined up with the B) animal time
# series below it, and the y-axis label will therefore change.

# Combine counts in df
exposures_df = data.frame("month_n"=1:length(ts_suspectBites$counts),
                          "Exposures"=ts_suspectBites$counts,
                          "Deaths"=ts_humanDeaths$counts)
bite_presentations = data.frame("month_n"=1:length(ts_allBites$counts),
                                "Bites"=ts_allBites$counts)

# Tranform to long dataframe
exposures_df = gather(exposures_df, exp, n, Exposures:Deaths)

# Set x-axis breaks and labels
breaks = 1:length(ts_suspectBites$counts)
labels = rep("", length(ts_suspectBites$counts))

# Create dataframe for timings bars at top
timing_bars_df = data.frame("Status" = c("PEP free", "PEP shortage", "PEP free"),
                            "xmin" = c(16, 60, 88), # April 2011 BMGF, Jan 2015 PEP shortage, May 2017 free PEP govt policy
                            "xmax" = c(60, 84, 120), # December 2014 Last date PEP available, PEP resupplied but not free, to date
                            "ymin" = c(30, 35, 39),
                            "ymax" = c(32, 37, 41))

# Plot
panel_a = ggplot() +
  annotate("rect", xmin = timing_bars_df$xmin[1], xmax = timing_bars_df$xmax[1],
           ymin = timing_bars_df$ymin[1], ymax = timing_bars_df$ymax[1], alpha = 0.8,fill = "grey") +
  annotate("rect", xmin = timing_bars_df$xmin[2], xmax = timing_bars_df$xmax[2],
           ymin = timing_bars_df$ymin[2], ymax = timing_bars_df$ymax[2], alpha = 0.5, fill = "#D22254") +
  annotate("rect", xmin = timing_bars_df$xmin[3], xmax = timing_bars_df$xmax[3],
           ymin = timing_bars_df$ymin[3], ymax = timing_bars_df$ymax[3], alpha = 0.8,fill = "grey") +
  annotate("text", label = "Free PEP (BMGF demonstration project)", x = 40, y = 34, size = 3, type = 2) + 
  annotate("text", label = "PEP shortages", x = 72, y = 39, size = 3, type = 2) + 
  annotate("text", label = "Free PEP policy", x = 105, y = 43, size = 3, type = 2) + 
  
  geom_line(data=bite_presentations, aes(x=month_n, y=Bites), size=0.7, col="grey") +
  geom_col(data=exposures_df, aes(x=month_n, y=n, fill=exp), width=1, alpha=1) +
  scale_fill_manual(values=col_pal, breaks=c("Exposures", "Deaths")) +
  scale_x_continuous(name="Year", breaks=breaks, labels=labels, expand = c(0.005, 0), limits = c(0, 120)) +
  scale_y_continuous(expand = c(0, 0.5), limits=c(0, 45)) +
  labs(y="Exposures & patients") +
  theme_classic() +
  theme(legend.position = "none", legend.title = element_blank(),
        axis.title = element_text(size=14), axis.text = element_text(size=12)) +
  theme( axis.title.x = element_blank(),
           axis.text.x = element_blank()) +
  geom_vline( xintercept = seq(12,144,12), linetype = "dashed", color = "gray", size = 0.3 )

panel_a

# panel_a = ggplot() +
#   geom_col(data=exposures_df, aes(x=month_n, y=n, fill=exp), width=1, alpha=1) +
#   geom_line(data=bite_presentations, aes(x=month_n, y=Bites), size=0.7) +
#   scale_fill_manual(values=col_pal, breaks=c("Exposures", "Deaths")) +
#   scale_x_continuous(name="Year", breaks=breaks, labels=labels) +
#   labs(y="Human exposures") +
#   theme_classic() +
#   theme(legend.position = "none",
#         axis.title = element_text(size=14), axis.text = element_text(size=12))
# panel_a


#----- Panel B -----------------------------------------------------------------

# B) also is formatted like: waningVax_cases_exp_Pemba, but has just animal cases
# and dog vax coverage and has a lower y-axis limit, so the height of panel A is
# higher than of panel B (max of say 15)

# Create dataframes for plotting
cases_df = data.frame("month_n"=1:length(ts_allBites$counts),
                      "Cases"=ts_suspectAnimals$counts)
waning_vax = data.frame("month_n"=1:length(vcWaningTotal_sub),
                        "Waning_cov"=vcWaningTotal_sub)

# Transform waning vaccination coverage to match second axis
waning_vax$Waning_cov = waning_vax$Waning_cov*15

# Set x-axis breaks and labels
breaks = 1:length(ts_allBites$counts)
labels = rep(ys, each=12)
new_labels = every_nth(labels, 12, inverse = T)

# Plot
panel_b = ggplot() +
  geom_col(data=cases_df, aes(x=month_n, y=Cases), fill=col_pal[1], width=1, alpha=1) +
  geom_line(data=waning_vax, aes(x=month_n, y=Waning_cov), size=0.7) +
  scale_x_continuous(name="Year", breaks=breaks, labels=new_labels, expand = c(0,0), limits = c(0, 120)) +
  scale_y_continuous(name="Rabid animal cases", limits=c(0,15), expand = c(0,0),
                     sec.axis = sec_axis(name="Waning vaccination coverage", trans=~./15)) +
  theme_classic() +
  theme(axis.title.y.right = element_text(angle=90),
        axis.title = element_text(size=14), axis.text = element_text(size=12), axis.title.x = element_blank()) +
  geom_vline( xintercept = seq(12,144,12), linetype = "dashed", color = "gray", size = 0.3 )

panel_b

#----- Panel C -----------------------------------------------------------------

# C) is pretty much the Tz map (enlarged from figure_1_raster_map.tiff)

# Extract bounding box for Pemba
Pemba_bb = st_as_sfc(st_bbox(PembaDist))

# Plot in ggplot2
panel_c = ggplot() +
  geom_sf(data=Tz_region, fill="white", color="black") +
  geom_sf(data=Pemba_bb, fill = NA, color = "#D22254", size = 1.2) +
  theme_void()

panel_c

#----- Panel D -----------------------------------------------------------------

# D) is the Pemba dog density map with the clinics overlaid (perhaps as red crosses?)

# Calculate mean dog densities in each year
dogPopMatYear_cell <- matrix(0, nrow=nrow(dogPopWardMat_cell), ncol=ncol(dogPopWardMat_cell)/12)
for(i in 1:ncol(dogPopMatYear_cell)){
  dogPopMatYear_cell[,i] <- rowSums(dogPopWardMat_cell[,((i-1)*12+1):(i*12)])/12
}

# Take only required year (2020)
dogPopMatYear_cell_sub <- dogPopMatYear_cell[,11]

# Create raster grid
dogPopGrid <- PembaGrid
dogPopGrid[which(dogPopGrid[]==0)] <- NA
dogPopGrid[which(!is.na(dogPopGrid[]))] <- dogPopMatYear_cell_sub

# Convert raster and shapefile to dataframe for Pemba map (not Tz inset map)
dogPopGrid_df <- as.data.frame(dogPopGrid, xy = TRUE)

# Plot
panel_d <- ggplot() +
  geom_raster(data=dogPopGrid_df, aes(x=x, y=y, fill=layer)) +
  geom_sf(data=PembaDist, fill=NA, color="black") +
  geom_point(data=clinic_locs, aes(x=eastings, y=northings), size=2, stroke=1, shape=15,
             color="#D22254", fill="#D22254", alpha=1) +
  scale_fill_gradient(name=expression(paste("Dogs/km"^2)), low="white", high="#191919",
                      na.value = "white") +
  theme_void() +
  theme(legend.position = c(.1,.7)) +  # "left",
        #legend.text = element_text(size=12),
        #legend.title = element_text(size=14)) +
  # north(x.min=550000, x.max=600000, y.min=9390000, y.max=9465000,
  #       symbol=9, location="bottomleft", scale=0.1) + # data=PembaDist_df
  # scalebar(x.min=550000, x.max=600000, y.min=9390000, y.max=9465000,
  #          location="bottomleft", dist=5, dist_unit="km", transform=FALSE,
  #          box.fill=c("white", "#323232"), st.size=4, border.size=.5) +
  coord_sf()

panel_d

#----- Panel E -----------------------------------------------------------------

# E) is the maps like in: Pemba/figs/waningVax_cases_exp_Pemba.pdf, keep colours for
# animal cases, human exposures and deaths matching to panels A and B, and grey shading
# for the coverage - much like you explore in the folder: figs/Map_testing

# Create breaks for choropleth map of coverage
chlor_breaks = seq(0, 1, length.out=100)

# Transform waning vaccination coverage to dataframe
vcWaningWard_Yr_df = as.data.frame(vcWaningWard_Yr)
names(vcWaningWard_Yr_df) = paste0("wvc_", ys_vax)

# Add in ward and district info
vcWaningWard_Yr_df$Ward_Name <- sub(".*\\_", "", rownames(vcWaningWard_Yr_df))
vcWaningWard_Yr_df$District_N <- sub("\\_.*", "", rownames(vcWaningWard_Yr_df))

# Merge data into ward shapefile
PembaWard_vax = merge(PembaWard, vcWaningWard_Yr_df, by=c("District_N", "Ward_Name"))
# View(PembaWard_vax)
#' _Mkoani-Shamiani does not have waning vax coverage data_

# Combine cases, exposures and deaths into a single dataframe
animalCases_sub <- dplyr::select(animalCases, Year, Date.bitten, Symptoms.started, "x"=UTM.Easting.jitter, "y"=UTM.Northing.jitter) %>% mutate("Source"="Cases")
humanExposures_sub <- dplyr::select(humanExposures, "Year"=Year.bitten, "x"=UTM.Easting.jitter, "y"=UTM.Northing.jitter) %>% mutate("Source"="Exposures")
humanDeaths_sub <- dplyr::select(humanDeaths, "Year"=Year.bitten, "x"=UTM.Easting.jitter, "y"=UTM.Northing.jitter) %>% mutate("Source"="Deaths")
case_exp_death = bind_rows(animalCases_sub, humanExposures_sub, humanDeaths_sub)

# Identify the earliest case (2016) and last case (2018) identified
arrow_df = data.frame("Year" = c("2014", "2016", "2016", "2018"),
                      "Date"= c(max(case_exp_death$Symptoms.started[which(case_exp_death$Year==2014)], na.rm=T),
                                min(case_exp_death$Symptoms.started[which(case_exp_death$Year==2016)], na.rm=T),
                                min(case_exp_death$Symptoms.started[which(case_exp_death$Year==2016)], na.rm=T),
                                max(case_exp_death$Symptoms.started[which(case_exp_death$Year==2018)], na.rm=T)))

# Extract coordinates for these dates (2x 2016, 1x 2018)
arrow_df$x = c(case_exp_death$x[which(case_exp_death$Symptoms.started == arrow_df$Date[1])],
               case_exp_death$x[which(case_exp_death$Symptoms.started == arrow_df$Date[2])],
               case_exp_death$x[which(case_exp_death$Symptoms.started == arrow_df$Date[4])])
arrow_df$y = c(case_exp_death$y[which(case_exp_death$Symptoms.started == arrow_df$Date[1])],
               case_exp_death$y[which(case_exp_death$Symptoms.started == arrow_df$Date[2])],
               case_exp_death$y[which(case_exp_death$Symptoms.started == arrow_df$Date[4])])

# Set plot standards
pt_alpha = 0.85

# Produce plots for each year
wvc_2010 = ggplot() +
  geom_sf(data=PembaWard_vax, aes(fill=wvc_2010), color=NA) +
  scale_fill_gradient2(name="Vaccination coverage",
                       low="white", mid="#cccccc", high="#191919", midpoint=0.5,
                       breaks=seq(0, 1, 0.25), limits=c(0,1)) +
  geom_sf(data=PembaDist, fill=NA, color="black", size = 0.05) +
  geom_jitter(data=case_exp_death[which(case_exp_death$Year==2010),],
             aes(x=x, y=y, color=Source), size=2, alpha=pt_alpha,
             width=1000, height=1000) +
  scale_color_manual(name="", values=col_pal, labels=c("Rabid animal", "Human exposure", "Human death")) +
  ggtitle("2010") +
  theme_void() +
  theme(legend.title = element_text(size=12), legend.text = element_text(size=10),
        legend.spacing.y = unit(0.1, 'cm')) +
  guides(color = guide_legend(override.aes = list(alpha = 1)))

wvc_2011 = ggplot() +
  geom_sf(data=PembaWard_vax, aes(fill=wvc_2011), color=NA) +
  scale_fill_gradient2(name="Vaccination coverage",
                       low="white", mid="#cccccc", high="#191919", midpoint=0.5,
                       breaks=seq(0, 1, 0.2), limits=c(0,1)) +
  geom_sf(data=PembaDist, fill=NA, color="black", size = 0.05) +
  geom_jitter(data=case_exp_death[which(case_exp_death$Year==2011),],
             aes(x=x, y=y, color=Source), size=2, alpha=pt_alpha,
             width=1000, height=1000) +
  scale_color_manual(name="", values=col_pal, labels=c("Rabid animal", "Human exposure", "Human death")) +
  ggtitle("2011") +
  theme_void() +
  theme(legend.position = "none")

wvc_2012 = ggplot() +
  geom_sf(data=PembaWard_vax, aes(fill=wvc_2012), color=NA) +
  scale_fill_gradient2(name="Vaccination coverage",
                       low="white", mid="#cccccc", high="#191919", midpoint=0.5,
                       breaks=seq(0, 1, 0.2), limits=c(0,1)) +
  geom_sf(data=PembaDist, fill=NA, color="black", size = 0.05) +
  geom_jitter(data=case_exp_death[which(case_exp_death$Year==2012),],
             aes(x=x, y=y, color=Source), size=2, alpha=pt_alpha,
             width=1000, height=1000) +
  scale_color_manual(name="", values=col_pal, labels=c("Rabid animal", "Human exposure", "Human death")) +
  ggtitle("2012") +
  theme_void() +
  theme(legend.position = "none")

wvc_2013 = ggplot() +
  geom_sf(data=PembaWard_vax, aes(fill=wvc_2013), color=NA) +
  scale_fill_gradient2(name="Vaccination coverage",
                       low="white", mid="#cccccc", high="#191919", midpoint=0.5,
                       breaks=seq(0, 1, 0.2), limits=c(0,1)) +
  geom_sf(data=PembaDist, fill=NA, color="black", size = 0.05) +
  geom_jitter(data=case_exp_death[which(case_exp_death$Year==2013),],
             aes(x=x, y=y, color=Source), size=2, alpha=pt_alpha,
             width=1000, height=1000) +
  scale_color_manual(name="", values=col_pal, labels=c("Rabid animal", "Human exposure", "Human death")) +
  ggtitle("2013") +
  theme_void() +
  theme(legend.position = "none")

wvc_2014 = ggplot() +
  geom_sf(data=PembaWard_vax, aes(fill=wvc_2014), color=NA) +
  scale_fill_gradient2(name="Vaccination coverage",
                       low="white", mid="#cccccc", high="#191919", midpoint=0.5,
                       breaks=seq(0, 1, 0.2), limits=c(0,1)) +
  geom_sf(data=PembaDist, fill=NA, color="black", size = 0.05) +
  geom_jitter(data=case_exp_death[which(case_exp_death$Year==2014),],
             aes(x=x, y=y, color=Source), size=2, alpha=pt_alpha,
             width=1000, height=1000) +
  scale_color_manual(name="", values=col_pal, labels=c("Rabid animal", "Human exposure", "Human death")) +
  annotate("segment", x = arrow_df$x[which(arrow_df$Year==2014)]-15000, y = arrow_df$y[which(arrow_df$Year==2014)]+10000,
           xend = arrow_df$x[which(arrow_df$Year==2014)]-1000, yend = arrow_df$y[which(arrow_df$Year==2014)],
           arrow = arrow(length = unit(0.05, "npc"), type = "closed"), size=1) + 
  ggtitle("2014") +
  theme_void() +
  theme(legend.position = "none")

wvc_2015 = ggplot() +
  geom_sf(data=PembaWard_vax, aes(fill=wvc_2015), color=NA) +
  scale_fill_gradient2(name="Vaccination coverage",
                       low="white", mid="#cccccc", high="#191919", midpoint=0.5,
                       breaks=seq(0, 1, 0.2), limits=c(0,1)) +
  geom_sf(data=PembaDist, fill=NA, color="black", size = 0.05) +
  geom_jitter(data=case_exp_death[which(case_exp_death$Year==2015),],
             aes(x=x, y=y, color=Source), size=2, alpha=pt_alpha,
             width=1000, height=1000) +
  scale_color_manual(name="", values=col_pal, labels=c("Rabid animal", "Human exposure", "Human death")) +
  ggtitle("2015") +
  theme_void() +
  theme(legend.position = "none")

wvc_2016 = ggplot() +
  geom_sf(data=PembaWard_vax, aes(fill=wvc_2016), color=NA) +
  scale_fill_gradient2(name="Vaccination coverage",
                       low="white", mid="#cccccc", high="#191919", midpoint=0.5,
                       breaks=seq(0, 1, 0.2), limits=c(0,1)) +
  geom_sf(data=PembaDist, fill=NA, color="black", size = 0.05) +
  geom_jitter(data=case_exp_death[which(case_exp_death$Year==2016),],
             aes(x=x, y=y, color=Source), size=2, alpha=pt_alpha,
             width=1000, height=1000) +
  scale_color_manual(name="", values=col_pal, labels=c("Rabid animal", "Human exposure", "Human death")) +
  annotate("segment", x = arrow_df$x[which(arrow_df$Year==2016)][1]-15000, y = arrow_df$y[which(arrow_df$Year==2016)][1]+10000,
           xend = arrow_df$x[which(arrow_df$Year==2016)][1]-1000, yend = arrow_df$y[which(arrow_df$Year==2016)][1],
           arrow = arrow(length = unit(0.05, "npc"), type = "closed"), size=1) + 
  # annotate("text", label = "detection", 
  #          x = arrow_df$x[which(arrow_df$Year==2016)][1]-20000, 
  #          y = arrow_df$y[which(arrow_df$Year==2016)][2]+10000, 
  #          size = 3) + 
  #  
  ggtitle("2016") +
  theme_void() +
  theme(legend.position = "none")


wvc_2017 = ggplot() +
  geom_sf(data=PembaWard_vax, aes(fill=wvc_2017), color=NA) +
  scale_fill_gradient2(name="Vaccination coverage",
                       low="white", mid="#cccccc", high="#191919", midpoint=0.5,
                       breaks=seq(0, 1, 0.2), limits=c(0,1)) +
  geom_sf(data=PembaDist, fill=NA, color="black", size = 0.05) +
  geom_jitter(data=case_exp_death[which(case_exp_death$Year==2017),],
             aes(x=x, y=y, color=Source), size=2, alpha=pt_alpha,
             width=1000, height=1000) +
  scale_color_manual(name="", values=col_pal, labels=c("Rabid animal", "Human exposure", "Human death")) +
  ggtitle("2017") +
  theme_void() +
  theme(legend.position = "none")

wvc_2018 = ggplot() +
  geom_sf(data=PembaWard_vax, aes(fill=wvc_2018), color=NA) +
  scale_fill_gradient2(name="Vaccination coverage",
                       low="white", mid="#cccccc", high="#191919", midpoint=0.5,
                       breaks=seq(0, 1, 0.2), limits=c(0,1)) +
  geom_sf(data=PembaDist, fill=NA, color="black", size = 0.05) +
  geom_jitter(data=case_exp_death[which(case_exp_death$Year==2018),],
             aes(x=x, y=y, color=Source), size=2, alpha=pt_alpha,
             width=1000, height=1000) +
  scale_color_manual(name="", values=col_pal, labels=c("Rabid animal", "Human exposure", "Human death")) +
  annotate("segment", x = arrow_df$x[which(arrow_df$Year==2018)]-15000, y = arrow_df$y[which(arrow_df$Year==2018)]+10000,
            xend = arrow_df$x[which(arrow_df$Year==2018)]-1000, yend = arrow_df$y[which(arrow_df$Year==2018)],
            arrow = arrow(length = unit(0.05, "npc"), type = "closed"), size=1) + 
  # annotate("text", label = "last case", 
  #          x = arrow_df$x[which(arrow_df$Year==2018)]-20000, 
  #          y = arrow_df$y[which(arrow_df$Year==2018)]+10000, 
  #          size = 3) + 
  
  ggtitle("2018") +
  theme_void() +
  theme(legend.position = "none")

# Save legend from 2010, then remove
wvc_legend <- g_legend(wvc_2010)
wvc_2010 <- wvc_2010 + theme(legend.position = "none")

#----- Combine all panels together ---------------------------------------------

ggarrange(panel_a,
          panel_b,
          ggarrange(panel_c, panel_d, wvc_2010, wvc_2011,
                    wvc_2012, wvc_2013, wvc_2014, wvc_2015,
                    wvc_2016, wvc_2017, wvc_2018, wvc_legend, # grid.draw(
                    nrow=3, ncol=4,
                    labels=c("C", "D", "E", "", "", "", "", "", "", "", "")),
          ncol=1, heights=c(1,1,3), align="v",
          labels=c("A", "B"))
ggsave("figures/Figure_1_Panel.pdf", height=12, width=10)
