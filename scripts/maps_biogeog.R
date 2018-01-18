#!/usr/bin/env Rscript

# load libs
library("rgbif")
library("rjson")
library("ggmap")
library("celestial")
library("tidyverse")
library("rJava")
library("OpenStreetMap")# needs rJava: first run 'sudo R CMD javareconf', then run 'sudo R' followed by 'install.packages("rJava", dependencies=FALSE)' 
library("maptools")
library("rgdal")
library("sp")
library("raster")
library("rasterVis")
library("rgeos")
library("spatial.tools")
library("RColorBrewer")
library("binr")
library("classInt")

# download Pseudolithoxus records from GBIF
gkey <- name_backbone(name='Pseudolithoxus')$genusKey
max <- occ_count(taxonKey=gkey, georeferenced=TRUE, basisOfRecord="PRESERVED_SPECIMEN")
gdat <- occ_search(taxonKey=gkey, hasCoordinate=TRUE, basisOfRecord="PRESERVED_SPECIMEN", return="data", fields="all", limit=max)
gdat <- data.frame(gdat)

# rename the n. sp.
gdat$specificEpithet[is.na(gdat$specificEpithet)] <- "kinja"
# clean up cat nuns
gdat$catalogNumber <- gsub(".* ", "", gdat$catalogNumber)

# DATA FROM MOL TABLE
ttab <- read_csv(file="../data/mol_samples.csv")
ttab <- ttab %>% dplyr::filter(genus == "Pseudolithoxus")

# Extras from MANUSCRIPT
mtab <- read_csv(file="../data/materials_examined_gps.csv")
# convert refs
mtab$decimalLongitude <- dms2deg(mtab$longitude,sep='dms')
mtab$decimalLatitude <-  dms2deg(mtab$latitude,sep='dms')

# subset and make new DF
# set ccolumns in common
gdat.red <- gdat %>% dplyr::select(institutionCode, catalogNumber, specificEpithet, decimalLatitude, decimalLongitude) %>% dplyr::mutate(catalogNumber=as.character(catalogNumber))
ttab.red <- ttab %>% dplyr::select(institutionCode, catalogNumber, specificEpithet, decimalLatitude, decimalLongitude) %>% dplyr::mutate(catalogNumber=as.character(catalogNumber))
mtab.red <- mtab %>% dplyr::select(institutionCode, catalogNumber, specificEpithet, decimalLatitude, decimalLongitude) %>% dplyr::mutate(catalogNumber=as.character(catalogNumber))

# combine
ff <- dplyr::bind_rows(gdat.red, ttab.red, mtab.red)

# clean
ff <- ff %>% dplyr::mutate(catalogNumber=str_replace_all(catalogNumber, "^0", "")) %>% distinct(catalogNumber, .keep_all=TRUE) %>% filter(!is.na(decimalLatitude))


#### TESTING DOWNLOADED MAPS
# sources:
# http://hydrosheds.cr.usgs.gov/index.php # OLD
# http://www.hydrosheds.org/ # hydrosheds NEW, need to login
# http://www.hydrosheds.org/page/hydrobasins # hydrobasins
# https://www.nceas.ucsb.edu/scicomp/usecases/ReadWriteESRIShapeFiles # opening shape files
# http://www.naturalearthdata.com/downloads/ # good basic maps
# https://en.wikipedia.org/wiki/List_of_GIS_data_sources
# http://fititnt.github.io/gis-dataset-brasil/ # state and municipality boundaries
# http://www.worldwildlife.org/pages/global-lakes-and-wetlands-database # wetlands database # 
# http://www.feow.org/downloads # FW ecoregions


## Prep the Maps
#https://pakillo.github.io/R-GIS-tutorial/
#https://cran.r-project.org/web/packages/raster/vignettes/Raster.pdf
# open raster and shapefiles
dem.ras <- raster("/home/rupert/Dropbox/Projects-temp/ancistrin-barcode-temp/maps/sa_dem_30s.bil")# USGS digital elevation model (void-filled) ?raster
rivs.shp <- readOGR("/home/rupert/Dropbox/Projects-temp/ancistrin-barcode-temp/maps/sa_riv_30s.shp", stringsAsFactors=FALSE)# USGS river network (stream lines)
wat.bod <- readOGR("/home/rupert/Dropbox/Projects-temp/ancistrin-barcode-temp/maps/gis.osm_water_a_free_1.shp")# http://download.geofabrik.de/south-america/brazil.html

## for the main map
# filter the rivers data removing all the small rivers
rivs.shp$UP_CELLS <- as.numeric(rivs.shp$UP_CELLS)
rivs.shp.red <- rivs.shp[rivs.shp$UP_CELLS > 5000, ]

# crop to region of interest
newext <- c(-71, -48.5, -12.5, 10)#c(left, right, bottom, top)
dem.ras.crop <- crop(dem.ras, newext)
rivs.shp.crop <- crop(rivs.shp.red, newext)
wat.bod.crop <- crop(wat.bod, newext)

# make the breaks
ci <- classIntervals(rivs.shp.crop$UP_CELLS, n=6, style="kmeans")
rivs.shp.crop$breaks <- NA
rivs.shp.crop$breaks[which(rivs.shp.crop$UP_CELLS > ci$brks[1] & rivs.shp.crop$UP_CELLS < ci$brks[2])] <- 0.25
rivs.shp.crop$breaks[which(rivs.shp.crop$UP_CELLS > ci$brks[2] & rivs.shp.crop$UP_CELLS < ci$brks[3])] <- 0.75
rivs.shp.crop$breaks[which(rivs.shp.crop$UP_CELLS > ci$brks[3] & rivs.shp.crop$UP_CELLS < ci$brks[4])] <- 1
rivs.shp.crop$breaks[which(rivs.shp.crop$UP_CELLS > ci$brks[4] & rivs.shp.crop$UP_CELLS < ci$brks[5])] <- 1.25
rivs.shp.crop$breaks[which(rivs.shp.crop$UP_CELLS > ci$brks[5] & rivs.shp.crop$UP_CELLS < ci$brks[6])] <- 1.75
rivs.shp.crop$breaks[which(rivs.shp.crop$UP_CELLS > ci$brks[6] & rivs.shp.crop$UP_CELLS < ci$brks[7])] <- 2
br <- rivs.shp.crop$breaks

## plotting data
# make colours  and symbols

# first make a mapping dataframe and join with prev
mapping_df <- data.frame(specificEpithet=unique(ff$specificEpithet), collist=c("gold1", "hotpink1", "chartreuse1", "grey80", "firebrick1", "dodgerblue1"), symb=c(24,21,21,25,22,23), stringsAsFactors=FALSE)
ff <- dplyr::left_join(ff, mapping_df, by="specificEpithet")

# add the types
ff <- ff %>% dplyr::mutate(types=ifelse(catalogNumber=="3220" | catalogNumber=="V-17544" | catalogNumber=="V-17546" | catalogNumber=="28355" | catalogNumber=="51644" | catalogNumber=="355199", "TRUE", "FALSE"))

# subset types
fft <- ff %>% filter(types==TRUE)
ff <- ff %>% filter(types==FALSE)

# col ramp for terrain and rivers
cbr <- rev(colorRampPalette(brewer.pal(n=9, name="YlGn")[1:9])(25))
riv.col <- "#7B9EC8"

# plot the map
pdf(file="../temp2/map_pseudolithoxus.pdf", useDingbats=FALSE, useKerning=FALSE)
plot(dem.ras.crop, add=FALSE, alpha=1, col=cbr, cex.axis=0.75, legend=FALSE)
plot(rivs.shp.crop, add=TRUE, col=riv.col, lwd=br)
plot(wat.bod.crop, add=TRUE, col=riv.col, lty=0)
points(x=ff$decimalLongitude, y=ff$decimalLatitude, col="grey10", bg=ff$collist, pch=ff$symb, lwd=0.5, cex=1.5)
points(x=fft$decimalLongitude, y=fft$decimalLatitude, col="black", bg=fft$collist, pch=fft$symb, lwd=2, cex=1.5)
rect(xleft=-68, ybottom=-1, xright=-66, ytop=1, border="white")
legend(x="topright", legend=unique(mapping_df$specificEpithet), text.font=3, cex=0.9, pt.lwd=0.5, pt.cex=1.25, col="black", pt.bg=mapping_df$collist, pch=mapping_df$symb, bty="n")
dev.off()

### to make a wider map

# crop to region of interest
newext <- c(-85, -35, -57.5, 12.5)#c(left, right, bottom, top)
dem.ras.crop <- crop(dem.ras, newext)
#rivs.shp.crop <- crop(rivs.shp.red, newext)
rivs.shp.red <- rivs.shp[rivs.shp$UP_CELLS > 12000, ]

ci <- classIntervals(rivs.shp.red$UP_CELLS, n=6, style="kmeans")
rivs.shp.red$breaks <- NA
rivs.shp.red$breaks[which(rivs.shp.red$UP_CELLS > ci$brks[1] & rivs.shp.red$UP_CELLS < ci$brks[2])] <- 0.25
rivs.shp.red$breaks[which(rivs.shp.red$UP_CELLS > ci$brks[2] & rivs.shp.red$UP_CELLS < ci$brks[3])] <- 0.75
rivs.shp.red$breaks[which(rivs.shp.red$UP_CELLS > ci$brks[3] & rivs.shp.red$UP_CELLS < ci$brks[4])] <- 1
rivs.shp.red$breaks[which(rivs.shp.red$UP_CELLS > ci$brks[4] & rivs.shp.red$UP_CELLS < ci$brks[5])] <- 1
rivs.shp.red$breaks[which(rivs.shp.red$UP_CELLS > ci$brks[5] & rivs.shp.red$UP_CELLS < ci$brks[6])] <- 2
rivs.shp.red$breaks[which(rivs.shp.red$UP_CELLS > ci$brks[6] & rivs.shp.red$UP_CELLS < ci$brks[7])] <- 2
br <- rivs.shp.red$breaks

# plot
pdf(file="../temp2/map_pseudolithoxus_sa.pdf", useDingbats=FALSE, useKerning=FALSE)
plot(dem.ras.crop, add=FALSE, alpha=0.95, col=cbr, cex.axis=0.75, legend=FALSE, axes=FALSE, frame=FALSE)
plot(rivs.shp.red, add=TRUE, col=rgb(134,194,230,maxColorValue = 255), lwd=br, axes=FALSE, frame=FALSE)#, 
rect(xleft=-70, ybottom=-12.5, xright=-47.5, ytop=10, border="red")
dev.off()


## small maps of upper negro
# crop to region of interest
newext <- c(-68, -66, -1, 1)#c(left, right, bottom, top)
dem.ras.crop <- crop(dem.ras, newext)
rivs.shp.red <- rivs.shp[rivs.shp$UP_CELLS > 400, ]
rivs.shp.crop <- crop(rivs.shp.red, newext)
wat.bod.crop <- crop(wat.bod, newext)

ci <- classIntervals(as.numeric(rivs.shp.crop$UP_CELLS), n=6, style="kmeans")
rivs.shp.crop$breaks <- NA
rivs.shp.crop$breaks[which(rivs.shp.crop$UP_CELLS > ci$brks[1] & rivs.shp.crop$UP_CELLS < ci$brks[2])] <- 0.25
rivs.shp.crop$breaks[which(rivs.shp.crop$UP_CELLS > ci$brks[2] & rivs.shp.crop$UP_CELLS < ci$brks[3])] <- 0.75
rivs.shp.crop$breaks[which(rivs.shp.crop$UP_CELLS > ci$brks[3] & rivs.shp.crop$UP_CELLS < ci$brks[4])] <- 1
rivs.shp.crop$breaks[which(rivs.shp.crop$UP_CELLS > ci$brks[4] & rivs.shp.crop$UP_CELLS < ci$brks[5])] <- 1
rivs.shp.crop$breaks[which(rivs.shp.crop$UP_CELLS > ci$brks[5] & rivs.shp.crop$UP_CELLS < ci$brks[6])] <- 2
rivs.shp.crop$breaks[which(rivs.shp.crop$UP_CELLS > ci$brks[6] & rivs.shp.crop$UP_CELLS < ci$brks[7])] <- 2
br <- rivs.shp.crop$breaks

# subset molecular nicoi
ttab.nic <- ttab[grep("nicoi", ttab$specificEpithet), ]

# plot the map
pdf(file="../temp2/map_pseudolithoxus_negro.pdf", useDingbats=FALSE, useKerning=FALSE)
plot(dem.ras.crop, add=FALSE, alpha=0.95, col=cbr, cex.axis=0.75, legend=FALSE, axes=FALSE, frame=FALSE, bty="n", box=FALSE)
plot(rivs.shp.crop, add=TRUE, col=riv.col, lwd=br, axes=FALSE, frame=FALSE)
plot(wat.bod.crop, add=TRUE, col=riv.col, lty=0, axes=FALSE, frame=FALSE)
points(x=ttab.nic$decimalLongitude, y=ttab.nic$decimalLatitude, col="grey10", bg="gold1", pch=24, lwd=0.5, cex=3)
dev.off()
#https://graphicdesign.stackexchange.com/questions/6419/increment-dynamic-offset-size
