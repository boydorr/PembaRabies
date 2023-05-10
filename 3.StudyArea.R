#---------------------------------------------------------
# Author: Kennedy Lushasi
# This script produces the location of the study area Pemba in relation to Tanzania mainland
#---------------------------------------------------------
rm(list=ls())

# Load libraries
library(rgdal)
library(raster)
library(ggplot2)
library(sf)
library(cowplot)
library(ggrepel)
library(RColorBrewer)
library(ggsn)

# Load functions
source("R/shp_to_df.R")

#---------------------------------------------------------
# Load shapefiles
Tz_region <- read_sf("data/GIS/TZ_Region_2012/TZ_Region_2012.shp")
PembaDist <- readOGR("Output/GIS/PembaDist_NBS2012", "PembaDist_NBS2012")

# Load Pemba grid
cell_size <- 1
PembaGrid <- raster(paste("Output/GIS/",cell_size^2,"kmsqGrids/PembaGrid",cell_size^2,"kmsq.grd",sep=""))

# Dog pop by cell/month (1km cell?)
dogPopWardMat_cell <- as.matrix(read.csv("Output/dogPopMatPembaCell.csv",header=F))

# Load clinic locations
clinic_locs <- read.csv("data/Clinic_GPS.csv", stringsAsFactors=FALSE)

#----- Process data ------------------------------------------------------------
# Calculate mean dog densities in each year
dogPopMatYear_cell <- matrix(0, nrow=nrow(dogPopWardMat_cell), ncol=ncol(dogPopWardMat_cell)/12)
for(i in 1:ncol(dogPopMatYear_cell)){
  dogPopMatYear_cell[,i] <- rowSums(dogPopWardMat_cell[,((i-1)*12+1):(i*12)])/12
}
# Take only required year (2020)
dogPopMatYear_cell_sub <- dogPopMatYear_cell[,11]

#----- Create map -------------------------------------------------------------

# Create raster grid
dogPopGrid <- PembaGrid
dogPopGrid[which(dogPopGrid[]==0)] <- NA
dogPopGrid[which(!is.na(dogPopGrid[]))] <- dogPopMatYear_cell_sub

# Extract bounding box for Pemba
Pemba_bb = st_as_sfc(st_bbox(PembaDist))

# Convert raster and shapefile to dataframe for Pemba map (not Tz inset map)
dogPopGrid_df <- as.data.frame(dogPopGrid, xy = TRUE)

# Set log scale breaks
# summary(dogPopMatYear_cell_sub)
# pl_br <- c(0, 1, 5, 20, 50)

# Plot in ggplot2
tz <- ggplot() +
  geom_sf(data=Tz_region, fill="white", color="black") +
  geom_sf(data=Pemba_bb, fill = NA, color = "red", size = 1.2) +
  theme_void()
pemba <- ggplot() +
  geom_raster(data=dogPopGrid_df, aes(x=x, y=y, fill=layer)) +
  geom_polygon(data=PembaDist, aes(x=long, y=lat, group=group), fill=NA, color="black") +
  geom_point(data=clinic_locs, aes(x=eastings, y=northings), shape=21, size=5,
             color="red", fill="red", alpha=0.8) +
  scale_fill_gradient(name=expression(paste("Dogs/km"^2)), low="white", high="#191919",
                      na.value = "white") + # trans="log", breaks=pl_br, labels=pl_br,
  theme_void() +
  theme(legend.position = c(.1,.5), # "left",
        legend.text = element_text(size=12),
        legend.title = element_text(size=14)) +
  north(x.min=550000, x.max=600000, y.min=9390000, y.max=9465000,
        symbol=9, location="bottomleft", scale=0.1) + # data=PembaDist_df
  scalebar(x.min=550000, x.max=600000, y.min=9390000, y.max=9465000,
           location="bottomleft", dist=5, dist_unit="km", transform=FALSE,
           box.fill=c("white", "#323232"), st.size=4, border.size=.5) +
  coord_equal(xlim=c(550000, 600000))

# Combine maps
gg_inset_map1 = ggdraw() +
  draw_plot(pemba) +
  draw_plot(tz, x = 0.01, y = 0.64, width = 0.35, height = 0.35)
gg_inset_map1
ggsave("figures/figure_1_rastermap.tiff", height=8, width=7)

#----- Plot in base R ----------------------------------------------------------

# Create colour palette
# colours <- colorRampPalette(c("white", "black"))(100)
# colours <- colorRampPalette(c("white",brewer.pal(8,"YlOrRd")[2:8]))(100)

# plot(PembaDist)
# plot(dogPopGrid, add=T, col=colours, breaks=seq(0, max(dogPopGrid[], na.rm=T), length.out=100), legend=F)
# plot(PembaDist, add=T)
# plot(dogPopGrid, breaks=seq(0,max(dogPopGrid[],na.rm=T),length.out=100),
#      legend.only=T, add=T,col=colours,
#      legend.args=list(text=expression(paste("Dogs/km"^2)), side=4, font=2, line=3, cex=1.2),
#      axis.args=list(at=seq(0,max(dogPopGrid[],na.rm=T),20),cex.axis=0.8),
#      smallplot=c(0.65,0.66, .25,.75))

