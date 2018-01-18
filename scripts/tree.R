#!/usr/bin/env Rscript

# load libs
library("ape")
library("phangorn")
library("phyloch")
library("phytools")
rm(list=ls())


## plot *BEAST species tree

mcc <- read.beast(file="../data/combo_species.tre", digits=2)
lmcc <- ladderize(mcc)
p <- character(length(lmcc$posterior))
co <- c("gray30", "white")
p[lmcc$posterior >= 0.95] <- co[1]
p[lmcc$posterior < 0.95] <- co[2]

# get dims
plot(lmcc, edge.color=0, tip.color=0)
HPDbars(lmcc, col="skyblue", lwd=7)
# plot
pdf(file="../temp2/starbeast_spTree_09-11-16.pdf", width=13, height=10, useDingbats=FALSE)
plot(lmcc, edge.color=0, tip.color=0, x.lim=c(-3.24214, 25.64739))
HPDbars(lmcc, col="skyblue", lwd=7)
plot.phylo.upon(lmcc, cex=1, edge.width=2, font=1, label.offset=0.2, edge.col="gray30", tip.color="grey40")
nodelabels(pch=21, bg=p, cex=1, col="gray30")
data(gradstein04)
data(strat2012)
axisGeo(GTS=strat2012, unit=c("stage", "epoch"), col="yellow", texcol="gray20", ages=TRUE, cex=1, gridty=3, gridcol="gray50")
dev.off()

str(lmcc)
lmcc$"height_95%_HPD_MAX"
lmcc$"height_95%_HPD_MIN"


## plot Gene trees CYTB
# read the tree
cytb.mcc <- read.beast(file="../data/combo_cytb.tre", digits=2)
cytb.lmcc <- ladderize(cytb.mcc)

# make posterior probs
p <- character(length(cytb.lmcc$posterior))
co <- c("gray30", "white")
p[cytb.lmcc$posterior >= 0.95] <- co[1]
p[cytb.lmcc$posterior < 0.95] <- co[2]
# copy tree
ntr <- cytb.lmcc
ttab <- read.table(file="../data/mol_samples.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("NA", ""))

# make the names
noms <- ifelse(test=ttab$taxonRank != "species", yes=paste(ttab$genus, ttab$identificationQualifier), no=paste(ttab$genus, ttab$specificEpithet))
ntr$tip.label <- paste0(ttab$catalogNumber[match(ntr$tip.label, ttab$catalogNumber)], " ", noms[match(ntr$tip.label, ttab$catalogNumber)], " (", ttab$waterBody[match(ntr$tip.label, ttab$catalogNumber)], ")")

# plot
pdf(file="../temp2/CYTB_geneTree_22-07-16.pdf", width=9, height=9, useDingbats=FALSE)
plot.phylo(ntr, cex=0.7, edge.width=2, no.margin=TRUE, font=1, label.offset=0.1, edge.col="gray30", tip.color="grey50")
nodelabels(pch=21, bg=p, cex=0.7, col="gray30")
#add.scale.bar(lwd=2, lcol="gray30", length=1, cex=0.5)
dev.off()


## plot Gene trees RAG1

rag.mcc <- read.beast(file="../data/combo_rag1.tre", digits=2)
rag.lmcc <- ladderize(rag.mcc)

# make posterior probs
p <- character(length(rag.lmcc$posterior))
co <- c("gray30", "white")
p[rag.lmcc$posterior >= 0.95] <- co[1]
p[rag.lmcc$posterior < 0.95] <- co[2]
# copy tree
ntr <- rag.lmcc
ttab <- read.table(file="../data/mol_samples.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("NA", ""))

ntr$tip.label <- gsub("a|b", "", ntr$tip.label)

# make the names
noms <- ifelse(test=ttab$taxonRank != "species", yes=paste(ttab$genus, ttab$identificationQualifier), no=paste(ttab$genus, ttab$specificEpithet))
ntr$tip.label <- paste0(ttab$catalogNumber[match(ntr$tip.label, ttab$catalogNumber)], " ", noms[match(ntr$tip.label, ttab$catalogNumber)], " (", ttab$waterBody[match(ntr$tip.label, ttab$catalogNumber)], ")")

# plot
pdf(file="../temp2/RAG1_geneTree_22-07-16.pdf", width=9, height=12, useDingbats=FALSE)
plot.phylo(ntr, cex=0.5, edge.width=2, no.margin=TRUE, font=1, label.offset=0.1, edge.col="gray30", tip.color="grey50")
nodelabels(pch=21, bg=p, cex=0.7, col="gray30")
#add.scale.bar(lwd=2, lcol="gray30", length=1, cex=0.5)
dev.off()
