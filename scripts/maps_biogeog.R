#!/usr/bin/env Rscript
require("rgbif")
require("rjson")
require("ggmap")
require("celestial")
require("dplyr")
require("OpenStreetMap")
require("maptools")
require("rgdal")
require("sp")
require("raster")
require("rasterVis")
require("rgeos")
require("spatial.tools")
require("RColorBrewer")
require("binr")
require("classInt")

# DL from GBIF
gkey <- name_backbone(name='Pseudolithoxus')$genusKey
max <- occ_count(taxonKey=gkey, georeferenced=TRUE, basisOfRecord="PRESERVED_SPECIMEN")
gdat <- occ_search(taxonKey=gkey, hasCoordinate=TRUE, basisOfRecord="PRESERVED_SPECIMEN", return="data", fields="all", limit=max)
gdat <- data.frame(gdat)

# rename the n. sp.
gdat$specificEpithet[is.na(gdat$specificEpithet)] <- "viator"
# clean up cat nuns
gdat$catalogNumber <- gsub(".* ", "", gdat$catalogNumber)

# DATA FROM MOL TABLE
ttab <- read.table(file="../data/mol_samples.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
ttab <- ttab[ttab$genus == "Pseudolithoxus", ]
ttab$specificEpithet[which(ttab$identificationQualifier == "n. sp.")] <- "viator" 

# Extras from MANUSCRIPT
mtab <- read.table(file="../data/materials_examined_gps.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
# convert refs
mtab$decimalLongitude <- dms2deg(mtab$longitude,sep='dms')
mtab$decimalLatitude <-  dms2deg(mtab$latitude,sep='dms')

# subset and make new DF
# set ccolumns in common
common_cols <- c("institutionCode", "catalogNumber", "specificEpithet", "decimalLatitude", "decimalLongitude")
# select cols of interest
ff <- rbind(subset(gdat, select=common_cols), subset(ttab, select=common_cols), subset(mtab, select=common_cols))
# sort
ff <- ff[order(ff$specificEpithet), ]
# remove NA
ff <- ff[!is.na(ff$decimalLatitude), ]
# remove duplicates
ff <- ff[!duplicated(ff$catalogNumber), ]



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
#ogrInfo(dsn="usgs_hydrosheds_maps/sa_dem_30s_grid/sa_dem_30s/sa_dem_30s", layer="sa_dem_30s")
#readOGR(dsn="/home/rupert/Downloads/sa_bas_15s_beta", layer="sa_bas_15s_beta")#?readOGR
#counties.mp <- readShapePoly("/home/rupert/Downloads/sa_bas_15s_beta/sa_bas_15s_beta")
#plot(counties.mp, axes=TRUE, border="gray")


## Prep the Maps
#https://pakillo.github.io/R-GIS-tutorial/
#https://cran.r-project.org/web/packages/raster/vignettes/Raster.pdf
# open raster and shapefiles
dem.ras <- raster("../data/maps/sa_dem_30s.bil")# USGS digital elevation model (void-filled)
rivs.shp <- readShapeLines("../data/maps/sa_riv_30s.shp")# USGS river network (stream lines)

# filter the rivers data removing all the 
rivs.shp.red <- rivs.shp[rivs.shp$UP_CELLS > 2500, ]

# crop to region of interest
newext <- c(-70, -47.5, -12.5, 10)#c(left, right, bottom, top)
dem.ras.crop <- crop(dem.ras, newext)
rivs.shp.crop <- crop(rivs.shp.red, newext)


# make the breaks
ci <- classIntervals(rivs.shp.red$UP_CELLS, n=6, style="kmeans")
rivs.shp.red$breaks <- NA
rivs.shp.red$breaks[which(rivs.shp.red$UP_CELLS > ci$brks[1] & rivs.shp.red$UP_CELLS < ci$brks[2])] <- 0.5
rivs.shp.red$breaks[which(rivs.shp.red$UP_CELLS > ci$brks[2] & rivs.shp.red$UP_CELLS < ci$brks[3])] <- 1.5
rivs.shp.red$breaks[which(rivs.shp.red$UP_CELLS > ci$brks[3] & rivs.shp.red$UP_CELLS < ci$brks[4])] <- 2.5
rivs.shp.red$breaks[which(rivs.shp.red$UP_CELLS > ci$brks[4] & rivs.shp.red$UP_CELLS < ci$brks[5])] <- 2.5
rivs.shp.red$breaks[which(rivs.shp.red$UP_CELLS > ci$brks[5] & rivs.shp.red$UP_CELLS < ci$brks[6])] <- 3.5
rivs.shp.red$breaks[which(rivs.shp.red$UP_CELLS > ci$brks[6] & rivs.shp.red$UP_CELLS < ci$brks[7])] <- 4.0
br <- rivs.shp.red$breaks


## Plot the Maps
# make colours (different options)
#cols1 <- brewer.pal(n=6, name="Set1")
#cols1 <- cols1[c(1,2,9,5,8,6)]
cols1 <- c("chartreuse1", "firebrick1", "dodgerblue1", "gold1", "grey80", "hotpink1")

# create cols vector
ff$collist <- NA
for (i in 1:length(cols1)){
    ff$collist[which(ff$specificEpithet %in% unique(ff$specificEpithet)[i])] <- cols1[i]
}

# add the type material column
ff$types <- FALSE
ff$types[which(ff$catalogNumber %in% c("3220", "V-17544", "V-17546", "28355", "51644", "355199.5263750"))] <- TRUE
fft <- ff[ff$types==TRUE, ]
fft <- fft[order(fft$specificEpithet), ]


# MY DATA
pch1 <- c(21,22,23,24,25,21)
ff$symb <- NA
for (i in 1:length(pch1)){
    ff$symb[which(ff$specificEpithet %in% unique(ff$specificEpithet)[i])] <- pch1[i]
}


# col ramp for terrain
cbr <- rev(colorRampPalette(brewer.pal(n=11, name="RdYlGn"))(16))

# plot the map
pdf(file="../temp2/map_pseudolithoxus.pdf", useDingbats=FALSE, useKerning=FALSE)
plot(dem.ras.crop, add=FALSE, alpha=0.9, col=cbr, cex.axis=0.75, legend=FALSE)
plot(rivs.shp.red, add=TRUE, col="grey20", lwd=br)
#points(x=gdat$decimalLongitude, y=gdat$decimalLatitude, col="black", bg=alpha(gdat$collist,1), pch=gdat$symb, lwd=0.5, cex=1.5)
points(x=ff$decimalLongitude, y=ff$decimalLatitude, col="grey10", bg=ff$collist, pch=ff$symb, lwd=0.5, cex=1.5)
points(x=fft$decimalLongitude, y=fft$decimalLatitude, col="black", bg=fft$collist, pch=fft$symb, lwd=2, cex=1.5)
#text(x=ff$decimalLongitude, y=ff$decimalLatitude, labels=ff$catalogNumber, cex=0.5, adj=0.2)# to check
legend(x="topright", legend=unique(ff$specificEpithet), text.font=3, cex=0.9, pt.lwd=0.5, pt.cex=1.25, col="black", pt.bg=cols1, pch=pch1, bty="n")
dev.off()













###### OLD CODE
require("rgbif")
require("rjson")
require("ggmap")
require("celestial")
require("dplyr")
require("OpenStreetMap")
require("maptools")
require("rgdal")
require("sp")
require("raster")
require("rasterVis")
require("rgeos")
require("spatial.tools")
require("RColorBrewer")
require("binr")
require("ape")
require("BioGeoBEARS")#?BioGeoBEARS

# DL from GBIF
gkey <- name_backbone(name='Pseudolithoxus')$genusKey
max <- occ_count(taxonKey=gkey, georeferenced=TRUE, basisOfRecord="PRESERVED_SPECIMEN")
gdat <- occ_search(taxonKey=gkey, hasCoordinate=TRUE, basisOfRecord="PRESERVED_SPECIMEN", return="data", fields="all", limit=max)
# remove dups
gdat <- gdat[grep("Nhamundá", gdat$locality, invert=TRUE), ]
# rename the n. sp.
gdat$specificEpithet[is.na(gdat$specificEpithet)] <- "n. sp."

# load up the data and clean it removing non-Pseudos
ttab <- read.table(file="../data/mol_samples.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
ttab <- ttab[ttab$genus == "Pseudolithoxus", ]
ttab$specificEpithet[which(ttab$identificationQualifier == "n. sp.")] <- "n. sp." 

# set ccolumns in common
common_cols <- c("institutionCode", "catalogNumber", "basisOfRecord", "order", "family", "genus", "specificEpithet", #
"identificationQualifier", "locality", "country", "stateProvince", "decimalLatitude", "decimalLongitude", "otherCatalogNumbers")#"key", "collectionCode", "scientificName", "dateIdentified", "eventDate", "samplingProtocol", "recordedBy", 
# select cols of interest
ff <- rbind(subset(gdat, select=common_cols), subset(ttab, select=common_cols))


# get a base map from google and plot it
map1 <- ggmap(get_map(location=c(lon=-62, lat=3), source="stamen", maptype="terrain", color="color", zoom=6), extent="panel")#?get_map
map1 + geom_point(data=ff, aes(x=decimalLongitude, y=decimalLatitude, colour=specificEpithet), shape=16, size=5, alpha=0.75, na.rm=TRUE) + labs(x=NULL, y=NULL)# + theme(legend.position="none")


ggsave(file="image.svg", width=8, height=8, dpi=150)

# ?openmap
map <- openmap(upperLeft=c(0,-70), lowerRight=c(-6.0,-64.0), type="esri")
m <- autoplot(map)