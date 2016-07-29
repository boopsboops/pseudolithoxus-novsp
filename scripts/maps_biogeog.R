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