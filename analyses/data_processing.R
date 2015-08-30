require("ape")
require("phyloch")
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
write.nexus.data(catcytbal, file="cytb.nex", format="dna", interleaved=FALSE, gap="-", missing="?")
# write a test dataset with labels to check in geneious
ttab <- read.table(file="mol_samples.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
dimnames(catcytbal)[[1]] <- paste(ttab$genus[match(labels(catcytbal), ttab$code)], ttab$species[match(labels(catcytbal), ttab$code)], ttab$code[match(labels(catcytbal), ttab$code)], sep="_")
write.dna(catcytbal, file="../temp/cytb_all_names.fas", format="fasta", colw=9999)

# write rag1 to disk
catragal <- gsub("-", "?", catragal)
write.nexus.data(catragal, file="rag1.nex", format="dna", interleaved=FALSE, gap="-", missing="?")
# write a test dataset with labels to check in geneious
dimnames(catragal)[[1]] <- paste(ttab$genus[match(labels(catragal), ttab$code)], ttab$species[match(labels(catragal), ttab$code)], ttab$code[match(labels(catragal), ttab$code)], sep="_")
write.dna(catragal, file="../temp/rag_all_names.fas", format="fasta", colw=9999)


## pull out the pseudolithoxus
cytb <- as.DNAbin(read.nexus.data(file="cytb.nex"))
rag <- as.DNAbin(read.nexus.data(file="rag1.nex"))
# remember to omit those NAs
cytbpseud <- cytb[na.omit(match(ttab$code[ttab$genus == "Pseudolithoxus"], labels(cytb)))]
ragpseud <- rag[na.omit(match(ttab$code[ttab$genus == "Pseudolithoxus"], labels(rag)))]
# write out
write.nexus.data(cytbpseud, file="../temp/cytb_pseudolithoxus.nex", format="dna", interleaved=FALSE, gap="-", missing="?")
write.nexus.data(ragpseud, file="../temp/rag1_pseudolithoxus.nex", format="dna", interleaved=FALSE, gap="-", missing="?")


## make a concatenated matrix
cytb <- as.matrix(as.DNAbin(read.nexus.data(file="cytb.nex")))
rag <- as.matrix(as.DNAbin(read.nexus.data(file="rag1.nex")))

# concatenate
li <- list(cytb, rag)
all <- c.genes(single.list=li, match=FALSE)

# save concat matrix
# convert Ns to ?s first
all <- gsub("n", "?", all)
write.nexus.data(all, file="../temp/concat.nex", format="dna", interleaved=FALSE, gap="-", missing="?")


## export files as phylip for partitionfinder
pcytb <- as.DNAbin(read.nexus.data(file="../temp/partitionfinder_datasets/cytb_pseudolithoxus.nex"))
prag <- as.DNAbin(read.nexus.data(file="../temp/partitionfinder_datasets/rag1_pseudolithoxus.nex"))
pall <- as.DNAbin(read.nexus.data(file="../temp/partitionfinder_datasets/concat.nex"))

write.dna(pcytb, file="../temp/partitionfinder_datasets/cytb_pseudolithoxus.phy", format="sequential", colw=9999)
write.dna(prag, file="../temp/partitionfinder_datasets/rag1_pseudolithoxus.phy", format="sequential", colw=9999)
write.dna(pall, file="../temp/partitionfinder_datasets/concat.phy", format="sequential", colw=9999)


## to extract codon positions
# specify only third position i.e. 1119/3 = 373
cyt.third <- 1:(length(cytb[1,])/3)*3
cyt.third_char <- cytb[, cyt.third]
cyt.12 <- cytb[, -cyt.third]

li <- list(cyt.12, rag)
all <- c.genes(single.list=li, match=FALSE)

# convert "third" to fasta file 
write.dna(all, "../temp/cytb12rag123.fas", format="fasta", colw=10000)
write.dna(cyt.third_char, "../temp/cytb3.fas", format="fasta", colw=10000)
write.dna(cyt.12, "../temp/cytb12.fas", format="fasta", colw=10000)

# write out to nex
write.nexus.data(cyt.12, file="../temp/cytb12.nex", format="dna", interleaved=FALSE, gap="-", missing="?")
write.nexus.data(cyt.third_char, file="../temp/cytb3.nex", format="dna", interleaved=FALSE, gap="-", missing="?")


## splitting into species
#detect iupac codes
lk <- lapply(as.list(rag), function(x) grep("r|y|s|w|k|m", x))
tk <- lapply(lk, length)

# get the poly and non-poly spp from the alignment
poly <- rag[which(as.integer(tk) > 0), ]
polyspp <- unique(ttab$species[match(labels(poly), ttab$code)])
nonpolyspp <- allspp[which(!unique(ttab$species) %in% polyspp)]

# polymorphic spp
for (i in 1:length(polyspp)){#
    write.dna(rag[na.omit(match(ttab$code[which(ttab$species == polyspp[i])], labels(rag))), ], file=paste0("../temp/species_split/", gsub(" |'", "", polyspp[i]), ".poly.fas"), format="fasta", colw=9999)
}#

# non-polymorphic spp (do i need this?)
#for (i in 1:length(nonpolyspp)){#
#    write.dna(rag[na.omit(match(ttab$code[which(ttab$species == nonpolyspp[i])], labels(rag))), ], file=paste0("../temp/species_split/", gsub(" |'", "", nonpolyspp[i]), ".fas"), format="fasta", colw=9999)
#}#
   

# create copies of indivs without polymorphisms
nonpolyseqsA <- as.list(rag[na.omit(match(ttab$code[which(ttab$species %in% nonpolyspp)], labels(rag))), ])
nonpolyseqsB <- nonpolyseqsA
names(nonpolyseqsA) <- paste0(names(nonpolyseqsA), "a")
names(nonpolyseqsB) <- paste0(names(nonpolyseqsB), "b")
write.nexus.data(c(nonpolyseqsA, nonpolyseqsB), file="../temp/rag1_phased.nex", interleaved=FALSE, gap="-", missing="?")
