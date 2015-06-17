# start the chunk
## @knitr treecode

require("rentrez")
require("ape")
require("phangorn")
require("phyloch")
require("phytools")

# search entrez
# note: there are other Cramer seqs available. Need to refine search if we want these
rs <- entrez_search(db="nuccore", term="(Pseudolithoxus[Organism] OR Soromonichthys[Organism] OR Lasiancistrus[Organism]) AND (cytb[Gene Name] OR RAG1[Gene Name] OR 16S[All Fields] OR RAG2[Gene Name] OR MyH6[Gene Name])", retmax=100)
# run to check all is okay
es <- entrez_summary(db="nuccore", id=rs$ids)
sapply(es, "[[", "title")
ess <- sapply(es, "[[", "title")

# download the seqs and write to disk
cytb <- entrez_fetch(db="nuccore", rettype='fasta', id=rs$ids[grep("cytb", ess)])
rag1 <- entrez_fetch(db="nuccore", rettype='fasta', id=rs$ids[grep("RAG1", ess)])
r16S <- entrez_fetch(db="nuccore", rettype='fasta', id=rs$ids[grep("16S", ess)])
rag2 <- entrez_fetch(db="nuccore", rettype='fasta', id=rs$ids[grep("RAG2", ess)])
myh6 <- entrez_fetch(db="nuccore", rettype='fasta', id=rs$ids[grep("MyH6", ess)])
write(cytb, "../temp/cytb_raw.fasta")
write(rag1, "../temp/rag1_raw.fasta")
write(r16S, "../temp/16S_raw.fasta")
write(rag2, "../temp/rag2_raw.fasta")
write(myh6, "../temp/myh6_raw.fasta")


# load up from disk
cytbrw <- read.dna("../temp/cytb_raw.fasta", format="fasta", as.matrix=FALSE)
rag1rw <- read.dna("../temp/rag1_raw.fasta", format="fasta", as.matrix=FALSE)
#rag1rw <- rag1rw[-grep("guapore", names(rag1rw))]#remove 'guapore'
#r16Srw <- read.dna("../temp/16S_raw.fasta", format="fasta", as.matrix=FALSE)
#r16Srw <- r16Srw[-grep("12S", names(r16Srw))]#remove MHNG specimens
#rag2rw <- read.dna("../temp/rag2_raw.fasta", format="fasta", as.matrix=FALSE)
#rag2rw <- rag2rw[-grep("guapore", names(rag2rw))]#remove 'guapore'
#myh6rw <- read.dna("../temp/myh6_raw.fasta", format="fasta", as.matrix=FALSE)
#myh6rw <- myh6rw[-grep("guapore", names(myh6rw))]#remove 'guapore'

tab <- read.table(file="mol_samples.csv", header=TRUE, stringsAsFactors=FALSE, sep=",")

nomscytb <- match(labels(cytbrw), tab$code)[!is.na(match(labels(cytbrw), tab$code))]

names(cytbrw) <- paste(tab$genus[nomscytb], tab$species[nomscytb], tab$locality[nomscytb], sep="_")




?na.omit
names(rag1rw) %in% tab$code
# clean names and create a names list

cytbdf <- paste(sapply(strsplit(names(cytbrw), split=" "), function(x) paste(x[c(5,2,3)], collapse=",")), gsub("\\.1", "", sapply(strsplit(names(cytbrw), split="\\|"), function(x) x[4])), sep=",")
write(cytbdf, file="../temp/cytb.csv", sep=",")
rag1df <- paste(sapply(strsplit(names(rag1rw), split=" "), function(x) paste(x[c(5,2,3)], collapse=",")), gsub("\\.1", "", sapply(strsplit(names(rag1rw), split="\\|"), function(x) x[4])), sep=",")
write(rag1df, file="../temp/rag1.csv", sep=",")
#names(r16Srw) <- paste(sapply(strsplit(names(r16Srw), split="\\|"), function(x) x[4]), sapply(strsplit(names(r16Srw), split=" "), function(x) paste(x[2:5], collapse="_")), sep="_")
#names(rag2rw) <- paste(sapply(strsplit(names(rag2rw), split="\\|"), function(x) x[4]), sapply(strsplit(names(rag2rw), split=" "), function(x) paste(x[2:5], collapse="_")), sep="_")
#names(myh6rw) <- paste(sapply(strsplit(names(myh6rw), split="\\|"), function(x) x[4]), sapply(strsplit(names(myh6rw), split=" "), function(x) paste(x[2:5], collapse="_")), sep="_")

names(cytbrw) <- 
names(rag1rw) <- 

# align the seqs
cytbAl <- mafft(x=cytbrw, path="mafft")
rag1Al <- mafft(x=rag1rw, path="mafft")
r16SAl <- mafft(x=r16Srw, path="mafft")
rag2Al <- mafft(x=rag2rw, path="mafft")
myh6Al <- mafft(x=myh6rw, path="mafft")

dat <- cytbAl
dat <- rag1Al
dat <- r16SAl
dat <- rag2Al
dat <- myh6Al


cytb <- read.dna("../temp/cytb_aligned_trimmed.fasta", format="fasta", as.matrix=FALSE)
rag1 <- read.dna("../temp/rag1_aligned_trimmed.fasta", format="fasta", as.matrix=FALSE)

concat <- c.genes(as.matrix(cytb), as.matrix(rag1), match=FALSE)#?c.genes

dat <- cytb
dat <- rag1
dat <- concat

write.dna(concat, file="../temp/concat.phy", format="interleaved", colw=9999)

# make a tree
prat <- pratchet(as.phyDat(dat), rearrangements="SPR")
pars <- acctran(prat, as.phyDat(dat))
mlm <- pml(pars, as.phyDat(dat), k=4, inv=0, model="HKY") 
mlik <- optim.pml(mlm, optNni=TRUE, optGamma=TRUE, optEdge=TRUE, optBf=TRUE, optQ=TRUE, optInv=FALSE, model="HKY")
tr <- mlik$tree# rename tree

noms <- match(tr$tip.label, tab$code)[!is.na(match(tr$tip.label, tab$code))]
tr$tip.label <- paste(tab$code[noms], tab$genus[noms], tab$species[noms], tab$locality[noms], sep="_")

outgr <- grep("Lasiancistrus", tr$tip.label, value=TRUE)
rnod <- fastMRCA(tr, outgr[1], outgr[2])
rtr <- reroot(tr, node.number=rnod, position=0.5*tr$edge.length[which(tr$edge[,2]==rnod)])

ltr <- ladderize(rtr)


code <- sapply(strsplit(ltr$tip.label, split="_"), function(x) x[1])
genus <- sapply(strsplit(ltr$tip.label, split="_"), function(x) x[2])
species <- sapply(strsplit(ltr$tip.label, split="_"), function(x) x[3])
locality <- sapply(strsplit(ltr$tip.label, split="_"), function(x) x[4])
locality <- gsub("NA", "", locality)


ltr$tip.label <- mixedFontLabel(code, genus, species, locality, sep = " ", italic=2:3, bold=2:3, always.upright=NULL)

pdf(file="../temp/pseudolithoxus_tree.pdf", width=11, height=9, useDingbats=FALSE, useKerning=FALSE)
plot.phylo(ltr, edge.color="dodgerblue", edge.width=4, tip.color="gray20", label.offset=0.001, cex=0.9, no.margin=TRUE)
add.scale.bar(length=0.01, lwd=2, lcol="dodgerblue")
dev.off()


cytbtr <- rtr
rag1tr <- rtr
r16Str <- rtr
rag2tr <- rtr
myh6tr <- rtr

# plot and take a look
par(mfrow=c(3,2))
plot.phylo(ladderize(cytbtr), direction="rightwards", main="cytb", cex=1)
plot.phylo(ladderize(rag1tr), direction="rightwards", main="rag1", cex=1)
plot.phylo(ladderize(r16Str), direction="rightwards", main="16S", cex=1)
plot.phylo(ladderize(rag2tr), direction="rightwards", main="rag2", cex=1)
plot.phylo(ladderize(myh6tr), direction="rightwards", main="myh6", cex=1)

# load up Nathan's data and split it

nmat <- as.matrix(as.DNAbin(read.nexus.data(file="../temp/nathan_concat.nex")))
ncytb <- nmat[,516:1663]
n16S <- nmat[,1:515]
nrag1 <- nmat[,1664:2682]
nrag2 <- nmat[,2683:3629]
nmyh <- nmat[,3630:4293]

dat <- nmyh

write.dna(ncytb, file="../temp/ncytb.fas", format="fasta", colw=9999)
write.nexus.data(acytb, file="../temp/ncytb.nex", interleaved=FALSE)
acytb <- read.dna(file="../temp/ncytb.fas", format="fasta")
write(gsub("_", ",", labels(acytb)), file="../temp/cytb.csv", append=TRUE)

write.dna(nrag1, file="../temp/nrag1.fas", format="fasta", colw=9999)
nrag1 <- read.dna("../temp/nrag1.fas", format="fasta", as.matrix=FALSE)

labels(nrag1)[which(labels(nrag1) %in% labels(rag1rw))]

nmat <- as.DNAbin(read.nexus.data(file="../temp/ncytb_sm.nex"))

dat <- nmat
tr <- midpoint(tr)
plot(tr)

### old
# grab data from GenBank
seqs <- read.GB(access.nb=paste("KF", 640094:640157, sep=""))

# make a quick NJ tree and midpoint root it
tr <- midpoint(nj(dist.dna(seqs, model="TN93")))

# convert all negative branches to zero
tr$edge.length[which(tr$edge.length < 0)] <- 0

# convert DNA to phydat format
# to drop a species:
# seqs[grep("atratulus", attr(seqs, "species"), invert=TRUE)]
pdat <- as.phyDat(seqs)

# fit an ML model to the tree
mod <- pml(unroot(tr), pdat, k=4, inv=0, model="HKY")
# optimise ML model params
fmod <- optim.pml(mod, optBf=TRUE, optQ=TRUE, optEdge=TRUE, optGamma=TRUE, model="HKY", optNni=TRUE, optInv=FALSE, control=pml.control(epsilon=1e-08, maxit=10, trace=0))

# reroot the tree again
rtr <- midpoint(fmod$tree)

# drop the outgroup tips 
dtr <- drop.tip(rtr, tip=rtr$tip.label[grep("atratulus", rtr$tip.label, invert=FALSE)])

# plot
plot(ladderize(dtr), edge.width=2, font=1, cex=0.9, adj=0.05, no.margin=TRUE)#?plot.phylo
