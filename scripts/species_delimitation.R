## load libs
require("spider")
require("genealogicalSorting")
require("phangorn")
# rm(list=ls())

## load data
cytb <- as.matrix(as.DNAbin(read.nexus.data(file="../data/cytb_aligned.nex")))
rag1 <- as.matrix(as.DNAbin(read.nexus.data(file="../data/rag1_aligned.nex")))
ttab <- read.table(file="../data/mol_samples.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)


## set up alternative species delimitation scenarios
nic12 <- ttab # nicoi gr1+gr2
nic12$delimitation <- gsub("nicoi.gr1|nicoi.gr2", "nicoi.gr1+gr2", nic12$delimitation)
nic1nsp <- ttab # nicoi gr1+nsp
nic1nsp$delimitation <- gsub("nicoi.gr1|n.sp.", "nicoi.gr1+n.sp.", nic1nsp$delimitation)
nic12nsp <- ttab # nicoi gr1+gr2+nsp
nic12nsp$delimitation <- gsub("nicoi.gr1|nicoi.gr2|n.sp.", "nicoi.gr1+gr2+n.sp.", nic12nsp$delimitation)


## seg sites analysis
# separate Pseudolithoxus genus and just the nicoi group
#cytbpseud <- cytb[na.omit(match(ttab$catalogNumber[grep("Pseudolithoxus", ttab$genus)], labels(cytb))), ]
#rag1pseud <- rag1[na.omit(match(ttab$catalogNumber[grep("Pseudolithoxus", ttab$genus)], labels(rag1))), ]

cytbnic <- cytb[na.omit(match(ttab$catalogNumber[grep("nicoi.gr2|nicoi.gr1|n.sp.", ttab$delimitation)], labels(cytb))), ]
rag1nic <- rag1[na.omit(match(ttab$catalogNumber[grep("nicoi.gr2|nicoi.gr1|n.sp.", ttab$delimitation)], labels(rag1))), ]

# make the species vectors (note that this writes over prev objects so run one at a time)
#cytbspp <- ttab$delimitation[match(labels(cytbpseud), ttab$catalogNumber)]# all pseudos: 3 spp. delim
#rag1spp <- ttab$delimitation[match(labels(rag1pseud), ttab$catalogNumber)]# all pseudos:3 spp. delim
#
#cytbspp <- nic12$delimitation[match(labels(cytbpseud), nic12$catalogNumber)]# all pseudos: nicoi gr1+gr2
#rag1spp <- nic12$delimitation[match(labels(rag1pseud), nic12$catalogNumber)]# all pseudos: nicoi gr1+gr2
#
#cytbspp <- nic1nsp$delimitation[match(labels(cytbpseud), nic1nsp$catalogNumber)]# all pseudos: nicoi gr1+n.sp.
#rag1spp <- nic1nsp$delimitation[match(labels(rag1pseud), nic1nsp$catalogNumber)]# all pseudos: nicoi gr1+n.sp.
#
#cytbspp <- nic12nsp$delimitation[match(labels(cytbpseud), nic12nsp$catalogNumber)]# all pseudos: nicoi gr1+gr2+n.sp.
#rag1spp <- nic12nsp$delimitation[match(labels(rag1pseud), nic12nsp$catalogNumber)]# all pseudos: nicoi gr1+gr2+n.sp.

#
cytbspp <- ttab$delimitation[match(labels(cytbnic), ttab$catalogNumber)]# just nicoi: 3 spp. delim
rag1spp <- ttab$delimitation[match(labels(rag1nic), ttab$catalogNumber)]# just nicoi: 3 spp. delim
#
cytbspp <- nic12$delimitation[match(labels(cytbnic), nic12$catalogNumber)]# just nicoi: nicoi gr1+gr2
rag1spp <- nic12$delimitation[match(labels(rag1nic), nic12$catalogNumber)]# just nicoi: nicoi gr1+gr2
#
cytbspp <- nic1nsp$delimitation[match(labels(cytbnic), nic1nsp$catalogNumber)]# just nicoi: nicoi gr1+n.sp.
rag1spp <- nic1nsp$delimitation[match(labels(rag1nic), nic1nsp$catalogNumber)]# just nicoi: nicoi gr1+n.sp.
#
#To view the nucleotide values 
#lapply(nucDiag(cytbpseud, cytbspp), length)
#lapply(nucDiag(rag1pseud, rag1spp), length)
#
lapply(nucDiag(cytbnic, cytbspp), length)
lapply(nucDiag(rag1nic, rag1spp), length)


## make some trees
# convert to phydat
dat <- as.phyDat(cytb)
dat <- as.phyDat(rag1)
# make parsimony tree
prat <- pratchet(dat, rearrangements="SPR")
pars <- acctran(prat, dat)
# make an ml tree
mlm <- pml(pars, dat, k=4, inv=0, model="HKY") # cytb
mlik <- optim.pml(mlm, optNni=TRUE, optGamma=TRUE, optInv=FALSE, model="HKY") # cytb
mlm <- pml(pars, dat, k=4, inv=0, model="K80") # rag1
mlik <- optim.pml(mlm, optNni=TRUE, optGamma=TRUE, optInv=FALSE, model="K80") # rag1
# copy trees
# cytb
cytbtr <- midpoint(mlik$tree)
dcytbtr <- drop.tip(cytbtr, tip=cytbtr$tip.label[match(ttab$catalogNumber[grep("ancistrus", ttab$genus, ignore.case=TRUE)], cytbtr$tip.label)])
plot(dcytbtr)
cytbtr.nic12 <- dcytbtr
cytbtr.nic12nsp <- dcytbtr
cytbtr.nic1nsp <- dcytbtr

# rag1
rag1tr <- midpoint(mlik$tree)
drag1tr <- drop.tip(rag1tr, tip=rag1tr$tip.label[match(ttab$catalogNumber[grep("ancistrus", ttab$genus, ignore.case=TRUE)], rag1tr$tip.label)])
plot(drag1tr)
rag1tr.nic12 <- drag1tr
rag1tr.nic12nsp <- drag1tr
rag1tr.nic1nsp <- drag1tr


## rosenberg analysis (cytb)
cytbros <- rosenberg(dcytbtr)
# create the different trees with sp delims
dcytbtr$tip.label <- ttab$delimitation[match(dcytbtr$tip.label, ttab$catalogNumber)]
cytbtr.nic12$tip.label <- nic12$delimitation[match(cytbtr.nic12$tip.label, nic12$catalogNumber)]
cytbtr.nic12nsp$tip.label <- nic12nsp$delimitation[match(cytbtr.nic12nsp$tip.label, nic12nsp$catalogNumber)]
cytbtr.nic1nsp$tip.label <- nic1nsp$delimitation[match(cytbtr.nic1nsp$tip.label, nic1nsp$catalogNumber)]
# extract the rosenberg P values
cytbros[names(cytbros) == getMRCA(dcytbtr, tip="n.sp.")]
cytbros[names(cytbros) == getMRCA(dcytbtr, tip="nicoi.gr1")]
cytbros[names(cytbros) == getMRCA(dcytbtr, tip="nicoi.gr2")]
cytbros[names(cytbros) == getMRCA(cytbtr.nic12, tip="nicoi.gr1+gr2")]
#cytbros[names(cytbros) == getMRCA(cytbtr.nic12nsp, tip="nicoi.gr1+gr2+n.sp.")]
cytbros[names(cytbros) == getMRCA(cytbtr.nic1nsp, tip="nicoi.gr1+n.sp.")]

## rosenberg analysis (rag1)
rag1ros <- rosenberg(drag1tr)
# create the different trees with sp delims
drag1tr$tip.label <- ttab$delimitation[match(drag1tr$tip.label, ttab$catalogNumber)]
rag1tr.nic12$tip.label <- nic12$delimitation[match(rag1tr.nic12$tip.label, nic12$catalogNumber)]
rag1tr.nic12nsp$tip.label <- nic12nsp$delimitation[match(rag1tr.nic12nsp$tip.label, nic12nsp$catalogNumber)]
rag1tr.nic1nsp$tip.label <- nic1nsp$delimitation[match(rag1tr.nic1nsp$tip.label, nic1nsp$catalogNumber)]
# extract the rosenberg P values
rag1ros[names(rag1ros) == getMRCA(drag1tr, tip="n.sp.")]
#rag1ros[names(rag1ros) == getMRCA(rag1tr, tip="nicoi.gr1")] # not monophyletic
#rag1ros[names(rag1ros) == getMRCA(rag1tr, tip="nicoi.gr2")] # not monophyletic
rag1ros[names(rag1ros) == getMRCA(rag1tr.nic12, tip="nicoi.gr1+gr2")]
#rag1ros[names(rag1ros) == getMRCA(rag1tr.nic1nsp, tip="nicoi.gr1+n.sp.")] # not monophyletic

# check for monophyly
plot.phylo(rag1tr.nic12, use.edge.length=FALSE)
is.monophyletic(phy=rag1tr.nic12, tips=which(rag1tr.nic12$tip.label == "nicoi.gr1+gr2"))

### DO BOOTSTRAP or POST PROBS ???

## gsi analysis
# cytb
# remember to run the trees again
# drop the non-pseudo tips
dcytbtr <- drop.tip(cytbtr, tip=cytbtr$tip.label[match(ttab$code[grep("ancistrus", ttab$genus, ignore.case=TRUE)], cytbtr$tip.label)])
# prepare the assignments file
write.table(cbind(dcytbtr$tip.label, ttab$species[match(dcytbtr$tip.label, ttab$code)]), file="../temp/cytb_gsi_assignments.txt", sep=" ", row.names=FALSE, col.names=FALSE, quote=TRUE) 
# run the analysis
cytbgsi <- singleTreeAnalysis(tree=dcytbtr, assignmentFile="../temp/cytb_gsi_assignments.txt", nperms=10000, nprocs=1)
names(cytbgsi$gsi) <- cytbgsi$groups
format(round(rbind(cytbgsi$gsi, cytbgsi$pvals), digits=4), scientific=FALSE)
# for alt delims
write.table(cbind(dcytbtr$tip.label, nic12$species[match(dcytbtr$tip.label, nic12$code)]), file="../temp/cytb_gsi_assignments.txt", sep=" ", row.names=FALSE, col.names=FALSE, quote=TRUE) 
write.table(cbind(dcytbtr$tip.label, nic1nsp$species[match(dcytbtr$tip.label, nic1nsp$code)]), file="../temp/cytb_gsi_assignments.txt", sep=" ", row.names=FALSE, col.names=FALSE, quote=TRUE) 
write.table(cbind(dcytbtr$tip.label, nic12nsp$species[match(dcytbtr$tip.label, nic12nsp$code)]), file="../temp/cytb_gsi_assignments.txt", sep=" ", row.names=FALSE, col.names=FALSE, quote=TRUE) 


# rag1
# drop the non-pseudo tips
# remember to run the trees again
drag1tr <- drop.tip(rag1tr, tip=rag1tr$tip.label[match(ttab$code[grep("ancistrus", ttab$genus, ignore.case=TRUE)], rag1tr$tip.label)])
# prepare the assignments file
write.table(cbind(drag1tr$tip.label, ttab$species[match(drag1tr$tip.label, ttab$code)]), file="../temp/rag1_gsi_assignments.txt", sep=" ", row.names=FALSE, col.names=FALSE, quote=TRUE) 
# run the analysis
rag1gsi <- singleTreeAnalysis(tree=drag1tr, assignmentFile="../temp/rag1_gsi_assignments.txt", nperms=10000, nprocs=1)
names(rag1gsi$gsi) <- rag1gsi$groups
format(round(rbind(rag1gsi$gsi, rag1gsi$pvals), digits=4), scientific=FALSE)
# for alt delims
write.table(cbind(drag1tr$tip.label, nic12$species[match(drag1tr$tip.label, nic12$code)]), file="../temp/rag1_gsi_assignments.txt", sep=" ", row.names=FALSE, col.names=FALSE, quote=TRUE) 
write.table(cbind(drag1tr$tip.label, nic1nsp$species[match(drag1tr$tip.label, nic1nsp$code)]), file="../temp/rag1_gsi_assignments.txt", sep=" ", row.names=FALSE, col.names=FALSE, quote=TRUE) 
write.table(cbind(drag1tr$tip.label, nic12nsp$species[match(drag1tr$tip.label, nic12nsp$code)]), file="../temp/rag1_gsi_assignments.txt", sep=" ", row.names=FALSE, col.names=FALSE, quote=TRUE) 


## other
plot(drag1tr)
?is.binary.tree(drag1tr)
## networks
#nn <- neighborNet(dist.dna(cytbnic, model="raw", pairwise.deletion=TRUE, as.matrix=TRUE))
#nn <- neighborNet(dist.dna(rag1nic, model="raw", pairwise.deletion=TRUE, as.matrix=TRUE))
#plot.networx(nn, "2D", edge.color="black", edge.width=1)