# start the chunk
## @knitr treecode

#require("rentrez")
require("ape")
require("phangorn")
require("phyloch")
require("phytools")
rm(list=ls())


### quick ML tree plot (Feb 16)
cytb <- as.matrix(as.DNAbin(read.nexus.data(file="../temp/final_alignments/cytb.nex")))
rag <- as.matrix(as.DNAbin(read.nexus.data(file="../temp/final_alignments/rag1.nex")))

dat <- cytb
dat <- rag

prat <- pratchet(as.phyDat(dat), rearrangements="SPR")
pars <- acctran(prat, as.phyDat(dat))
mlm <- pml(pars, as.phyDat(dat), k=4, inv=0, model="HKY")
mlik <- optim.pml(mlm, optNni=TRUE, optGamma=TRUE, optInv=FALSE, model="HKY")
tr <- mlik$tree
rtr <- ladderize(reorder(midpoint(tr)))


ttab <- read.table(file="mol_samples.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
rtr$tip.label <- paste(ttab$code[match(rtr$tip.label, ttab$code)], ttab$genus[match(rtr$tip.label, ttab$code)], ttab$species[match(rtr$tip.label, ttab$code)], ttab$locality_drainage[match(rtr$tip.label, ttab$code)], sep="_")
rag.tr <- rtr
cytb.tr <- rtr

pdf(file="../temp/pseudolithoxus_rag1.pdf", width=12, height=15, useDingbats=FALSE)
plot.phylo(rtr, cex=1, edge.width=2, no.margin=TRUE, font=1, label.offset=0.0001, edge.col="gray30", tip.color="grey50")
dev.off()



### 

cytb <- as.matrix(as.DNAbin(read.nexus.data(file="cytb.nex")))
rag <- as.matrix(as.DNAbin(read.nexus.data(file="rag1.nex")))
all <- as.matrix(as.DNAbin(read.nexus.data(file="../temp/concat.nex")))

# convert to phydat
dat <- as.phyDat(cytb)
dat <- as.phyDat(rag)
dat <- as.phyDat(all)

## make parsimony tree
prat <- pratchet(dat, rearrangements="SPR")
pars <- acctran(prat, dat)

## run a modeltest 
mt <- modelTest(dat, tree=pars, G=TRUE, I=FALSE, multicore=TRUE)
mts <- mt[with(mt, order(AICc)), ]# sort by AICc
cbind(mts$Model, mts$AICc- mts$AICc[1])# get the AIC delta values (GTR+G is best)




## make a tree
mlm <- pml(pars, dat, k=4, inv=0, model="GTR") 
mlik <- optim.pml(mlm, optNni=TRUE, optGamma=TRUE, optQ=TRUE, optEdge=TRUE, optBf=TRUE, optInv=FALSE, model="GTR")
tr <- mlik$tree# rename tree
#tr <- pars
#reroot
rnod <- fastMRCA(tr, "B1470", "T13826")
#rnod <- fastMRCA(tr, "Lasiancistrus_schomburgkii_B1470", "Ancistrus_clementinae_T13826")
rtr <- ladderize((reroot(tr, node.number=rnod, position=0.75*tr$edge.length[which(tr$edge[,2]==rnod)])))

#plot.phylo(rtr, no.margin=TRUE, cex=0.8, font=1)
#nodelabels()

# check missing
#ttab$code[which(!ttab$code %in% names(catcytb))]
#names(catcytb)[which(!names(catcytb) %in% ttab$code)]

# bootstaps
btrs <- bootstrap.pml(mlik, bs=80, trees=TRUE, multicore=TRUE, mc.cores=8, optNni=TRUE, optGamma=TRUE, optQ=TRUE, optBf=TRUE, optInv=FALSE, model="GTR")# make bs multiple of 8!
pp <- prop.part(btrs, check.labels=TRUE)
pc <- prop.clades(rtr, part=pp)
bs <- round(pc/80, digits=2)
p <- character(length(bs))
co <- c("gray30", "gray", "white")
p[bs >= 0.95] <- co[1]
p[bs < 0.95 & bs >= 0.70] <- co[2]
p[bs < 0.70] <- co[3]

# copy tree
ntr <- rtr
ttab <- read.table(file="mol_samples.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
ntr$tip.label <- ttab$species[match(ntr$tip.label, ttab$code)]

# interspecies nodes
insp <- c(#
    fastMRCA(ntr, "tentaculatus", "schomburgkii"),
    fastMRCA(ntr, "tentaculatus", "clementinae"),
    fastMRCA(ntr, "ranunculus", "clementinae"),
    fastMRCA(ntr, "ranunculus", "sp.Xingu"),
    fastMRCA(ntr, "sp.Xingu", "macrophthalmus"),
    fastMRCA(ntr, "macrophthalmus", "leucostictus"),
    fastMRCA(ntr, "leucostictus", "megalostomus"),
    fastMRCA(ntr, "megalostomus", "sp.Inambari"),
    fastMRCA(ntr, "sp.Inambari", "bolivianus"),
    fastMRCA(ntr, "kelsorum", "dumus"),
    fastMRCA(ntr, "kelsorum", "tigris"),
    fastMRCA(ntr, "dumus", "stearleyi"),
    fastMRCA(ntr, "dumus", "anthrax"),
    fastMRCA(ntr, "stearleyi", "nicoi.gr2"),
    fastMRCA(ntr, "nicoi.gr2", "nicoi.gr1"),
    fastMRCA(ntr, "nicoi.gr1", "n.sp.")
)#

#plot(ntr)
#nodelabels()

# intraspecies nodes 
want <- c(insp, unlist(lapply(unique(ntr$tip.label), function(x) getMRCA(ntr, x))))
# get numbers
df <- data.frame(cbind(sort(unique(ntr$edge[,1])), p))
names(df) <- c("node", "boot")
wm <- match(want, df$node) 

# add taxon names
rtr$tip.label <- mixedFontLabel(ttab$code[match(rtr$tip.label, ttab$code)], ttab$genus[match(rtr$tip.label, ttab$code)], ttab$species[match(rtr$tip.label, ttab$code)], italic=2:3)

# plot
pdf(file="../temp/pseudolithoxus_tree.pdf", width=8, height=10, useDingbats=FALSE)
plot.phylo(rtr, cex=0.8, edge.width=2, no.margin=TRUE, font=1, label.offset=0.001, edge.col="gray30", tip.color="grey50")
nodelabels(node=want, pch=21, bg=p[wm], cex=1.25, col="gray30")
nodelabels(bs, frame="none", col="red", cex=0.8, adj=c(1.2,1.5))
legend("bottomleft", legend=c("BS > 0.95", "BS < 0.95", "BS < 0.70"), pch=21, cex=0.75, pt.bg=co, pt.cex=1.25, bty="n")
dev.off()

############### *BEAST species tree

mcc <- read.beast(file="../temp/beast_sptree/geocalib_comb.species.tre", digits=2)
lmcc <- ladderize(mcc)
p <- character(length(lmcc$posterior))
co <- c("gray30", "white")
p[lmcc$posterior >= 0.95] <- co[1]
p[lmcc$posterior < 0.95] <- co[2]

pdf(file="../temp/pseudolithoxus_beast_sptree.pdf", width=6, height=6, useDingbats=FALSE)
plot.phylo(lmcc, cex=0.8, edge.width=3, no.margin=TRUE, font=1, label.offset=0.005, edge.col="gray30", tip.color="grey50")
nodelabels(pch=21, bg=p, cex=1.1, col="dodgerblue")
legend("bottomleft", legend=c("Bayesian posterior probability >= 0.95", "Bayesian posterior probability < 0.95"), pch=21, cex=0.75, pt.bg=co, pt.cex=1.25, bty="n", col="dodgerblue")
dev.off()

plot(lmcc, edge.color=0, tip.color=0)
HPDbars(lmcc, col="skyblue", lwd=7)

pdf(file="../temp/pseudolithoxus_beast_sptree.pdf", width=15, height=9, useDingbats=FALSE)
plot(lmcc, edge.color=0, tip.color=0, x.lim=c(-3.40436, 25.61854))
HPDbars(lmcc, col="skyblue", lwd=7)
plot.phylo.upon(lmcc, cex=1, edge.width=2, no.margin=TRUE, font=1, label.offset=0.2, edge.col="gray30", tip.color="grey50")
nodelabels(pch=21, bg=p, cex=0.9, col="gray30")
data(gradstein04)
data(strat2012)
axisGeo(GTS=strat2012, unit=c("stage", "epoch"), col="yellow", texcol="gray20", ages=TRUE, cex=1, gridty=3, gridcol="gray50")
dev.off()

str(lmcc)

lmcc$"height_95%_HPD_MAX"
lmcc$"height_95%_HPD_MIN"

################ read a beast tree
# read the tree
mcc <- read.beast(file="../temp/beast/3parts.tre", digits=2)
lmcc <- ladderize(mcc)
p <- character(length(lmcc$posterior))
co <- c("dodgerblue", "white")
p[lmcc$posterior >= 0.95] <- co[1]
p[lmcc$posterior < 0.95] <- co[2]
# copy tree
ntr <- lmcc
ttab <- read.table(file="mol_samples.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
ntr$tip.label <- ttab$species[match(ntr$tip.label, ttab$code)]
# interspecies nodes
insp <- c(#
    fastMRCA(ntr, "tentaculatus", "schomburgkii"),
    fastMRCA(ntr, "tentaculatus", "clementinae"),
    fastMRCA(ntr, "tentaculatus", "clementinae"),
    fastMRCA(ntr, "ranunculus", "clementinae"),
    fastMRCA(ntr, "ranunculus", "leucostictus"),
    fastMRCA(ntr, "leucostictus", "sp. 'Xingu'"),
    fastMRCA(ntr, "sp. 'Xingu'", "macrophthalmus"),
    fastMRCA(ntr, "sp. 'Xingu'", "megalostomus"),
    fastMRCA(ntr, "megalostomus", "sp. 'Inambari'"),
    fastMRCA(ntr, "sp. 'Inambari'", "bolivianus"),
    fastMRCA(ntr, "kelsorum", "tentaculatus"),
    fastMRCA(ntr, "kelsorum", "tigris"),
    fastMRCA(ntr, "kelsorum", "dumus"),
    fastMRCA(ntr, "dumus", "stearleyi"),
    fastMRCA(ntr, "dumus", "anthrax"),
    fastMRCA(ntr, "stearleyi", "nicoi.gr2"),
    fastMRCA(ntr, "nicoi.gr2", "nicoi.gr1"),
    fastMRCA(ntr, "nicoi.gr1", "n. sp.")
)#
# intraspecies nodes 
want <- c(insp, unlist(lapply(unique(ntr$tip.label), function(x) getMRCA(ntr, x))))
# get numbers from edge matrix
df <- data.frame(cbind(sort(unique(ntr$edge[,1])), p))
names(df) <- c("node", "bpp")
wm <- match(want, df$node) 
# add taxon names
lmcc$tip.label <- mixedFontLabel(ttab$code[match(lmcc$tip.label, ttab$code)], ttab$genus[match(lmcc$tip.label, ttab$code)], #
gsub("\\.gr1|\\.gr2", "", ttab$species[match(lmcc$tip.label, ttab$code)]), italic=2:3, always.upright=c("n. sp.", "sp. 'Inambari'", "sp. 'Xingu'"))
# plot
pdf(file="../temp/pseudolithoxus_beast_tree.pdf", width=8, height=12, useDingbats=FALSE)
plot.phylo(lmcc, cex=0.8, edge.width=3, no.margin=TRUE, font=1, label.offset=0.005, edge.col="gray30", tip.color="grey50")
nodelabels(node=want, pch=21, bg=p[wm], cex=1.1, col="dodgerblue")
legend("bottomleft", legend=c("Bayesian posterior probability >= 0.95", "Bayesian posterior probability < 0.95"), pch=21, cex=0.75, pt.bg=co, pt.cex=1.25, bty="n", col="dodgerblue")
dev.off()









#################### 
# search entrez
# note: there are other Cramer seqs available. Need to refine search if we want these
rs <- entrez_search(db="nuccore", term="(Pseudolithoxus[Organism] OR Ancistrus[Organism] OR Soromonichthys[Organism] OR Lasiancistrus[Organism]) AND (cytb[Gene Name] OR RAG1[Gene Name] OR 16S[All Fields] OR RAG2[Gene Name] OR MyH6[Gene Name])", retmax=100)
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

tab <- read.table(file="mol_samples.csv", header=TRUE, stringsAsFactors=FALSE, sep=",")
noms <- match(tr$tip.label, tab$code)[!is.na(match(tr$tip.label, tab$code))]
tr$tip.label <- paste(tab$code[noms], tab$genus[noms], tab$species[noms], tab$locality[noms], sep="_")

outgr <- grep("ncistrus", tr$tip.label, value=TRUE)
rnod <- fastMRCA(unroot(tr), outgr[3], outgr[9])
rtr <- reroot(tr, node.number=rnod, position=0.5*tr$edge.length[which(tr$edge[,2]==rnod)])

ltr <- ladderize(rtr)


plot(ltr)

# make a starting tree for starbeast
tr2 <- mlik$tree
rnod <- fastMRCA(unroot(tr2), "T09686", "T10383")
rtr <- reroot(tr2, node.number=rnod, position=0.5*tr2$edge.length[which(tr2$edge[,2]==rnod)])
# make a list of the nodes that we want to calibrate
cnod <- fastMRCA(unroot(tr2), "T13829", "T10383")
# set the ages for these nodes and make a dataframe for chronos (penalised likelihood)
ccal <- makeChronosCalib(rtr, node=cnod, age.min=12, soft.bounds=TRUE)
# run chronos using a strict clock (rate cat of 1)
ctr <- chronos(rtr, calibration=ccal, model="discrete", control=chronos.control(nb.rate.cat=1))
plot(ladderize(ctr))
# reclass the tree so it writes correctly
class(ctr) <- "phylo"
write.tree(ladderize(ctr), file="../temp/pseudolithoxus_st_tr.nwk")

keep <- c("CTGA14485", "T09934",  "T09949", "P18073", "V5533", "T09938", "T09898", "T09686", "P6125", "T13829", "B1500", "T08143", "T09397", "B1988", "T10092", "T12872", "T10383")



# reduce to a species tree with correct tips
dtr <- drop.tip(ctr, tip=ctr$tip.label[!ctr$tip.label %in% keep])

noms <- match(dtr$tip.label, tab$code)[!is.na(match(dtr$tip.label, tab$code))]
dtr$tip.label <- tab$species[noms]
plot(dtr)
write.tree(ladderize(dtr), file="../temp/pseudolithoxus_species_tr.nwk")

###

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
