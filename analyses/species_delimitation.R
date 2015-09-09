require("spider")
require("genealogicalSorting")

## load data
cytb <- as.matrix(as.DNAbin(read.nexus.data(file="../temp/final_alignments/cytb.nex")))
rag1 <- as.matrix(as.DNAbin(read.nexus.data(file="../temp/final_alignments/rag1.nex")))
ttab <- read.table(file="mol_samples.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)


## separate the P. nicoi group
cytbpseud <- cytb[na.omit(match(ttab$code[grep("nicoi.gr2|nicoi.gr1|n. sp.", ttab$species)], labels(cytb))), ]
cytbpseudspp <- ttab$species[match(labels(cytbpseud), ttab$code)]
rag1pseud <- rag1[na.omit(match(ttab$code[grep("nicoi.gr2|nicoi.gr1|n. sp.", ttab$species)], labels(rag1))), ]
rag1pseudspp <- ttab$species[match(labels(rag1pseud), ttab$code)]


## seg sites analysis
#To view the nucleotide values 
cytbseg <- nucDiag(cytbpseud, cytbpseudspp)
lapply(cytbseg, length)
rag1seg <- nucDiag(rag1pseud, rag1pseudspp)
lapply(rag1seg, length)


## make some trees
# convert to phydat
dat <- as.phyDat(cytb)
dat <- as.phyDat(rag1)
# make parsimony tree
prat <- pratchet(dat, rearrangements="SPR")
pars <- acctran(prat, dat)
# make an ml tree
mlm <- pml(pars, dat, k=4, inv=0, model="HKY") # cytb
mlik <- optim.pml(mlm, optNni=TRUE, optGamma=TRUE, optEdge=TRUE, optBf=TRUE, optInv=FALSE, model="HKY") # cytb
mlm <- pml(pars, dat, k=4, inv=0, model="K80") # rag1
mlik <- optim.pml(mlm, optNni=TRUE, optGamma=TRUE, optEdge=TRUE, optBf=FALSE, optInv=FALSE, model="K80") # rag1
cytbtr <- midpoint(mlik$tree)
rag1tr <- midpoint(mlik$tree)


## rosenberg analysis
cytbros <- rosenberg(cytbtr)
cytbtr$tip.label <- ttab$species[match(cytbtr$tip.label, ttab$code)]
cytbros[names(cytbros) == getMRCA(cytbtr, tip="n. sp.")]
cytbros[names(cytbros) == getMRCA(cytbtr, tip="nicoi.gr1")]
cytbros[names(cytbros) == getMRCA(cytbtr, tip="nicoi.gr2")]

rag1ros <- rosenberg(rag1tr)
rag1tr$tip.label <- ttab$species[match(rag1tr$tip.label, ttab$code)]
rag1ros[names(rag1ros) == getMRCA(rag1tr, tip="n. sp.")]
# nicoi.gr1 and nicoi.gr2 are not monophyletic


## gsi analysis
# cytb
# drop the non-pseudo tips
dcytbtr <- drop.tip(cytbtr, tip=cytbtr$tip.label[match(ttab$code[grep("ancistrus", ttab$genus, ignore.case=TRUE)], cytbtr$tip.label)])
# prepare the assignments file
write.table(cbind(dcytbtr$tip.label, ttab$species[match(dcytbtr$tip.label, ttab$code)]), file="../temp/cytb_gsi_assignments.txt", sep=" ", row.names=FALSE, col.names=FALSE, quote=TRUE) 
# run the analysis
cytbgsi <- singleTreeAnalysis(tree=dcytbtr, assignmentFile="../temp/cytb_gsi_assignments.txt", nperms=10000, nprocs=1)
names(cytbgsi$gsi) <- cytbgsi$groups
format(round(rbind(cytbgsi$gsi, cytbgsi$pvals), digits=4), scientific=FALSE)

# rag1
# drop the non-pseudo tips
drag1tr <- drop.tip(rag1tr, tip=rag1tr$tip.label[match(ttab$code[grep("ancistrus", ttab$genus, ignore.case=TRUE)], rag1tr$tip.label)])
# prepare the assignments file
write.table(cbind(drag1tr$tip.label, ttab$species[match(drag1tr$tip.label, ttab$code)]), file="../temp/rag1_gsi_assignments.txt", sep=" ", row.names=FALSE, col.names=FALSE, quote=TRUE) 
# run the analysis
rag1gsi <- singleTreeAnalysis(tree=drag1tr, assignmentFile="../temp/rag1_gsi_assignments.txt", nperms=10000, nprocs=1)
names(rag1gsi$gsi) <- rag1gsi$groups
format(round(rbind(rag1gsi$gsi, rag1gsi$pvals), digits=4), scientific=FALSE)