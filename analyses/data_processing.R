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

