# start the chunk
## @knitr treecode

require("rentrez")
require("ape")
require("phangorn")
require("phyloch")
require("phytools")

# search entrez
# note: there are other Cramer seqs available. Need to refine search if we want these
rs <- entrez_search(db="nuccore", term="(Pseudolithoxus[Organism] OR Soromonichthys[Organism] OR Lasiancistrus[Organism]) AND (cytb[Gene Name] OR RAG1[Gene Name])", retmax=100)
# run to check all is okay
es <- entrez_summary(db="nuccore", id=rs$ids)
sapply(es, "[[", "title")
ess <- sapply(es, "[[", "title")

# download the seqs and write to disk
cytb <- entrez_fetch(db="nuccore", rettype='fasta', id=rs$ids[grep("cytb", ess)])
rag1 <- entrez_fetch(db="nuccore", rettype='fasta', id=rs$ids[grep("RAG1", ess)])
write(cytb, "../temp/cytb_raw.fasta")
write(rag1, "../temp/rag1_raw.fasta")

# load up from disk
cytbrw <- read.dna("../temp/cytb_raw.fasta", format="fasta", as.matrix=FALSE)
rag1rw <- read.dna("../temp/rag1_raw.fasta", format="fasta", as.matrix=FALSE)
rag1rw <- rag1rw[-grep("guapore", names(rag1rw))]#remove 'guapore'

# clean names
names(cytbrw) <- paste(sapply(strsplit(names(cytbrw), split="\\|"), function(x) x[4]), sapply(strsplit(names(cytbrw), split=" "), function(x) paste(x[2:3], collapse="_")), sep="_")
names(rag1rw) <- paste(sapply(strsplit(names(rag1rw), split="\\|"), function(x) x[4]), sapply(strsplit(names(rag1rw), split=" "), function(x) paste(x[2:3], collapse="_")), sep="_")

# align the seqs
cytbAl <- mafft(x=cytbrw, path="mafft")
rag1Al <- mafft(x=rag1rw, path="mafft")

dat <- cytbAl
dat <- rag1Al

# make a tree
prat <- pratchet(as.phyDat(dat), rearrangements="SPR")
pars <- acctran(prat, as.phyDat(dat))
mlm <- pml(pars, as.phyDat(dat), k=4, inv=0, model="GTR") 
mlik <- optim.pml(mlm, optNni=TRUE, optGamma=TRUE, optEdge=TRUE, optBf=TRUE, optQ=TRUE, optInv=FALSE)
tr <- mlik$tree# rename tree
outgr <- grep("Lasiancistrus", tr$tip.label, value=TRUE)
rnod <- fastMRCA(tr, outgr[1], outgr[2])
rtr <- reroot(tr, node.number=rnod, position=0.5*tr$edge.length[which(tr$edge[,2]==rnod)])

cytbtr <- rtr
rag1tr <- rtr

# plot and take a look
par(mfrow=c(1,2))
plot.phylo(ladderize(cytbtr), direction="rightwards")
plot.phylo(ladderize(rag1tr), direction="leftwards")



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
