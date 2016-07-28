require("ape")
require("phyloch")
require("spider")
require("phangorn")
# rm(list=ls())

## load and process data
# load up nathan extras and extract cytb and rag1 data
nmat <- as.matrix(as.DNAbin(read.nexus.data(file="lujan_extras.nex")))
ncytb <- nmat[,532:1662]# includes trim off the non-cytb parts at the start (originally 516:1663)
nrag <- nmat[,1665:2681]# includes trim off the first and last bases to start the codon correctly (originally 1664:2682)

# write temp files to disk after removing gaps and missing sequences
write.dna(del.gaps(ncytb), file="../temp/cytb_lujan.fasta", format="fasta", colw=9999)
# remove missing seq
nrag <- nrag[-which(labels(nrag) == "T13832"),]
write.dna(del.gaps(nrag), file="../temp/rag1_lujan.fasta", format="fasta", colw=9999)

# load up extracted cytb from this study and nathan's
ncytb <- read.dna(file="../temp/cytb_lujan.fasta", format="fasta", as.matrix=FALSE)
lcytb <- read.dna(file="cytb.fasta", format="fasta", as.matrix=FALSE)
# load up extracted rag1 from this study and nathan's
nrag <- read.dna(file="../temp/rag1_lujan.fasta", format="fasta", as.matrix=FALSE)
lrag <- read.dna(file="rag1.fasta", format="fasta", as.matrix=FALSE)

# cat the two lists
catcytb <- c(ncytb, lcytb)
catrag <- c(nrag, lrag)

# align
catcytbal <- mafft(x=catcytb, path="mafft")
catragal <- mafft(x=catrag, path="mafft")

# write cytb file to disk
catcytbal <- gsub("-", "?", catcytbal)
write.nexus.data(catcytbal, file="../temp/final_alignments/cytb.nex", format="dna", interleaved=FALSE, gap="-", missing="?")
write.dna(catcytbal, file="../temp/final_alignments/cytb.phy", format="sequential", colw=9999)
# write a test dataset with labels to check in geneious
#ttab <- read.table(file="mol_samples.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
#dimnames(catcytbal)[[1]] <- paste(ttab$genus[match(labels(catcytbal), ttab$code)], ttab$species[match(labels(catcytbal), ttab$code)], ttab$code[match(labels(catcytbal), ttab$code)], sep="_")
#write.dna(catcytbal, file="../temp/cytb_all_names.fas", format="fasta", colw=9999)

# write rag1 to disk
catragal <- gsub("-", "?", catragal)
write.nexus.data(catragal, file="../temp/final_alignments/rag1.nex", format="dna", interleaved=FALSE, gap="-", missing="?")
# write a test dataset with labels to check in geneious
#dimnames(catragal)[[1]] <- paste(ttab$genus[match(labels(catragal), ttab$code)], ttab$species[match(labels(catragal), ttab$code)], ttab$code[match(labels(catragal), ttab$code)], sep="_")
#write.dna(catragal, file="../temp/rag_all_names.fas", format="fasta", colw=9999)


## pull out the pseudolithoxus
cytb <- as.DNAbin(read.nexus.data(file="../temp2/final_alignments/cytb.nex"))
rag <- as.DNAbin(read.nexus.data(file="../temp2/final_alignments/rag_phased.nex"))
ttab <- read.table(file="mol_samples.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
cytbpseud <- cytb[na.omit(match(ttab$code[grep("dumus|anthrax|nicoi.gr2|nicoi.gr1|n. sp.", ttab$species)], labels(cytb)))]
rag2 <- rag
names(rag2) <- gsub("a|b", "", names(rag2))
ragpseud <- rag[na.omit(which(labels(rag2) %in% ttab$code[grep("dumus|anthrax|nicoi.gr2|nicoi.gr1|n. sp.", ttab$species)]))]
# write out
write.nexus.data(ragpseud, file="../temp/final_alignments/rag1_pseudolithoxus_phased.nex", format="dna", interleaved=FALSE, gap="-", missing="?")
# split codons
third <- 1:(length(as.matrix(cytbpseud)[1,])/3)*3
thirdchars <- as.matrix(cytbpseud)[, third]
firsec <- as.matrix(cytbpseud)[, -third]
# write out
write.nexus.data(thirdchars, file="../temp/final_alignments/cytb_pseudolithoxus_cp3.nex", format="dna", interleaved=FALSE, gap="-", missing="?")
write.nexus.data(firsec, file="../temp/final_alignments/cytb_pseudolithoxus_cp12.nex", format="dna", interleaved=FALSE, gap="-", missing="?")
write.nexus.data(cytbpseud, file="../temp/final_alignments/cytb_pseudolithoxus.nex", format="dna", interleaved=FALSE, gap="-", missing="?")
write.dna(cytbpseud, file="../temp/final_alignments/cytb_pseudolithoxus.phy", format="sequential", colw=9999)


## export for BP&P
cytb <- as.matrix(as.DNAbin(read.nexus.data(file="../data/cytb_aligned.nex")))
rag <- as.matrix(as.DNAbin(read.nexus.data(file="../data/rag1_phased_aligned.nex")))
ttab <- read.table(file="../data/mol_samples.csv", header=TRUE, sep=",", stringsAsFactors=FALSE, na.strings=c("", "NA"))

cats <- ttab$catalogNumber[which(ttab$specificEpithet == "nicoi" | ttab$identificationQualifier == "n. sp.")]

cattab <- cbind(ttab$catalogNumber[which(ttab$specificEpithet == "nicoi" | ttab$identificationQualifier == "n. sp.")], #
    paste(ttab$waterBody[which(ttab$specificEpithet == "nicoi" | ttab$identificationQualifier == "n. sp.")], #
    ttab$locality[which(ttab$specificEpithet == "nicoi" | ttab$identificationQualifier == "n. sp.")]), seq(from=1, to=length(unique(cats)), by=1))

phasedcats <- c(paste0(cats, "a"), paste0(cats, "b"))

cytb.red <- cytb[match(cats, labels(cytb)), ]
rag.red <- rag[match(phasedcats, labels(rag)), ]

rag.red.tru <- gsub("a|b", "", labels(rag.red))
    
dimnames(rag.red)[[1]] <- paste0(rag.red.tru, "^", cattab[,3][match(rag.red.tru, cattab[,1])])
dimnames(cytb.red)[[1]] <- paste0(labels(cytb.red), "^", cattab[,3][match(labels(cytb.red), cattab[,1])])

write.dna(cytb.red, file="../data/bpp.phy", format="sequential", colw=9999, append=FALSE)
write.dna(rag.red, file="../data/bpp.phy", format="sequential", colw=9999, append=TRUE)

# for the IMAP file
codes <- toupper(rep(letters, length.out=length(unique(cattab[,2]))))

li <- vector(mode="character", length=length(cattab[,2]))
for(i in 1:length(cattab[,2])){
    li[which(cattab[,2] %in% unique(cattab[,2])[i])] <- codes[i]
    }#

# for loci catalogNumber     
(table(li) * 3)

# write out
write.table(cbind(cattab[,3], li), file="../data/bpp_imap.tsv", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)

cbind(cattab[,2], li)

# make a random starting tree
tr <- rtree(n=length(codes), tip.label=codes)
tr$edge.length <- NULL
write.tree(phy=tr, file="../temp2/tr.nwk")

bpptr <- read.tree(text="(H, ((B, C) #0.823470, (G, (D, (E, (A, F) #0.280090) #0.631330) #0.920580) #0.997170) #1.000000) #1.000000;")
bpptr <- read.tree(text="(H, ((B, C) #0.854280, (G, (D, (E, (A, F) #0.299310) #0.660270) #0.939200) #0.998690) #1.000000) #1.000000;")

plot(bpptr)
nodelabels(bpptr$node.label)

## make a concatenated matrix
cytb <- as.matrix(as.DNAbin(read.nexus.data(file="../temp2/final_alignments/cytb.nex")))
rag <- as.matrix(as.DNAbin(read.nexus.data(file="../temp2/final_alignments/rag1.nex")))

# concatenate
li <- list(cytb, rag)
all <- c.genes(single.list=li, match=FALSE)

# save concat matrix
# convert Ns to ?s first
all <- gsub("n", "?", all)
write.nexus.data(all, file="../temp/concat.nex", format="dna", interleaved=FALSE, gap="-", missing="?")


## export some files as phylip for partitionfinder
pcytb <- as.DNAbin(read.nexus.data(file="../temp2/final_alignments/cytb_pseudolithoxus.nex"))
prag <- as.DNAbin(read.nexus.data(file="../temp2/final_alignments/rag1_pseudolithoxus_phased.nex"))
#pall <- as.DNAbin(read.nexus.data(file="../temp/final_alignments/concat.nex"))
#phrag <- as.DNAbin(read.nexus.data(file="../temp/final_alignments/rag_phased.nex"))

write.dna(pcytb, file="../temp2/final_alignments/cytb_pseudolithoxus.phy", format="sequential", colw=9999)
write.dna(prag, file="../temp2/final_alignments/rag1_pseudolithoxus_phased.phy", format="sequential", colw=9999)
#write.dna(pall, file="../temp/final_alignments/concat.phy", format="sequential", colw=9999)
#write.dna(phrag, file="../temp/final_alignments/rag_phased.phy", format="sequential", colw=9999)


## to extract codon positions for checking in jmodeltest
cytb <- as.matrix(as.DNAbin(read.nexus.data(file="../temp2/final_alignments/cytb_pseudolithoxus.nex")))
rag <- as.matrix(as.DNAbin(read.nexus.data(file="../temp2/final_alignments/rag1_phased.nex")))
dat <- rag
dat <- cytb

cp1 <- seq(from=1, to=(dim(dat)[2]-2), by=3)
cp2 <- seq(from=2, to=(dim(dat)[2]-1), by=3)
cp3 <- seq(from=3, to=(dim(dat)[2]-0), by=3)
write.dna(dat[ ,cp1], file="../temp2/final_alignments/cytb_pseudolithoxus_CP1.fasta", format="fasta", colw=9999)
write.dna(dat[ ,cp2], file="../temp2/final_alignments/cytb_pseudolithoxus_CP2.fasta", format="fasta", colw=9999)
write.dna(dat[ ,cp3], file="../temp2/final_alignments/cytb_pseudolithoxus_CP3.fasta", format="fasta", colw=9999)
write.dna(dat[ ,cp1], file="../temp2/final_alignments/rag1_phased_CP1.fasta", format="fasta", colw=9999)
write.dna(dat[ ,cp2], file="../temp2/final_alignments/rag1_phased_CP2.fasta", format="fasta", colw=9999)
write.dna(dat[ ,cp3], file="../temp2/final_alignments/rag1_phased_CP3.fasta", format="fasta", colw=9999)
write.dna(dat[ ,-cp1], file="../temp2/final_alignments/rag1_pseudolithoxus_phased_CP23.fasta", format="fasta", colw=9999)


## splitting into species
# need to convert ? back to - otherwise phase will impute missing values
rag <- as.DNAbin(gsub("\\?", "-", rag))

#detect iupac codes
lk <- lapply(as.list(rag), function(x) grep("r|y|s|w|k|m", x))
tk <- lapply(lk, length)

# get the poly and non-poly spp from the alignment
ttab <- read.table(file="mol_samples.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
poly <- rag[which(as.integer(tk) > 0), ]
polyspp <- unique(ttab$species[match(labels(poly), ttab$code)])
nonpolyspp <- unique(ttab$species)[which(!unique(ttab$species) %in% polyspp)]

# polymorphic spp
for (i in 1:length(polyspp)){#
    write.dna(rag[na.omit(match(ttab$code[which(ttab$species == polyspp[i])], labels(rag))), ], file=paste0("../temp/species_split/", gsub(" |'", "", polyspp[i]), ".poly.fas"), format="fasta", colw=9999)
}#

# non-polymorphic spp (do i need this?)
#for (i in 1:length(nonpolyspp)){#
#    write.dna(rag[na.omit(match(ttab$code[which(ttab$species == nonpolyspp[i])], labels(rag))), ], file=paste0("../temp/species_split/", gsub(" |'", "", nonpolyspp[i]), ".fas"), format="fasta", colw=9999)
#}#
   

# open up the phased file to check it
pop <- read.dna(file="../temp/species_split/poly_phased_all.fas", format="fasta", as.matrix=TRUE)
lapply(as.list(pop), function(x) grep("r|y|s|w|k|m|n", x))
pop2 <- as.DNAbin(gsub("-", "?", pop))
# create copies of indivs without polymorphisms to phase them
nonpolyseqsA <- as.list(rag[na.omit(match(ttab$code[which(ttab$species %in% nonpolyspp)], labels(rag))), ])
nonpolyseqsB <- nonpolyseqsA
names(nonpolyseqsA) <- paste0(names(nonpolyseqsA), "a")
names(nonpolyseqsB) <- paste0(names(nonpolyseqsB), "b")
phased.all <- c(nonpolyseqsA, nonpolyseqsB, as.list(pop2))
# write out
write.nexus.data(phased.all, file="../temp/final_alignments/rag_phased.nex", format="dna", interleaved=FALSE, gap="-", missing="?")


### missing data & stats
cytb <- as.matrix(as.DNAbin(read.nexus.data(file="../data/cytb_aligned.nex")))
rag1 <- as.matrix(as.DNAbin(read.nexus.data(file="../data/rag1_aligned.nex")))

print(cytb)
bstats(data=cytb)
seqStat(cytb)

print(rag1)
bstats(data=rag1)
seqStat(rag1)

# load the function first
bstats <- function(data){#
    # alignment length
    ddat <- dim(data)[2]
    # num parsimony informative sites
    pisnum <- pis(data, abs=TRUE)
    # proportion parsimony informative sites
    pisprop <- round((pisnum / ddat), digits=3)
    # proportion of missing data (Ns)
    propn <- round((length(grep("n|\\?", data)) / length(data)), digits=3)
    # return results
    print(deparse(substitute(data)))
    names(ddat) <- "Alignment length"
    print(ddat)
    names(pisnum) <- "Number parsimony informative sites"
    print(pisnum)
    names(pisprop) <- "Proportion parsimony informative sites"
    print(pisprop)
    names(propn) <- "Proportion of missing data (Ns, ?s)"
    print(propn)
}#
