require("spider")
require("genealogicalSorting")
rm(list=ls())

## load data
cytb <- as.matrix(as.DNAbin(read.nexus.data(file="../temp/final_alignments/cytb.nex")))
rag1 <- as.matrix(as.DNAbin(read.nexus.data(file="../temp/final_alignments/rag1.nex")))
ttab <- read.table(file="mol_samples.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)





## set up alternative species delimitation scenarios
# nicoi gr1+gr2
nic12 <- ttab
nic12$species <- gsub("nicoi.gr1|nicoi.gr2", "nicoi.gr1+gr2", nic12$species)
# nicoi gr1+nsp
nic1nsp <- ttab
nic1nsp$species <- gsub("nicoi.gr1|n. sp.", "nicoi.gr1+n.sp.", nic1nsp$species)
# nicoi gr1+gr2+nsp
nic12nsp <- ttab
nic12nsp$species <- gsub("nicoi.gr1|nicoi.gr2|n. sp.", "nicoi.gr1+gr2+n.sp.", nic12nsp$species)


## seg sites analysis
# separate Pseudolithoxus genus and just the nicoi group
cytbpseud <- cytb[na.omit(match(ttab$code[grep("Pseudolithoxus", ttab$genus)], labels(cytb))), ]
rag1pseud <- rag1[na.omit(match(ttab$code[grep("Pseudolithoxus", ttab$genus)], labels(rag1))), ]
cytbnic <- cytb[na.omit(match(ttab$code[grep("nicoi.gr2|nicoi.gr1|n. sp.", ttab$species)], labels(cytb))), ]
rag1nic <- rag1[na.omit(match(ttab$code[grep("nicoi.gr2|nicoi.gr1|n. sp.", ttab$species)], labels(rag1))), ]

# make the species vectors
# note writes over prev objects! run one at a time
cytbspp <- ttab$species[match(labels(cytbpseud), ttab$code)]# all pseudos: 3 spp. delim
rag1spp <- ttab$species[match(labels(rag1pseud), ttab$code)]# all pseudos:3 spp. delim
#
cytbspp <- nic12$species[match(labels(cytbpseud), nic12$code)]# all pseudos: nicoi gr1+gr2
rag1spp <- nic12$species[match(labels(rag1pseud), nic12$code)]# all pseudos: nicoi gr1+gr2
#
cytbspp <- nic1nsp$species[match(labels(cytbpseud), nic1nsp$code)]# all pseudos: nicoi gr1+n.sp.
rag1spp <- nic1nsp$species[match(labels(rag1pseud), nic1nsp$code)]# all pseudos: nicoi gr1+n.sp.
#
cytbspp <- nic12nsp$species[match(labels(cytbpseud), nic12nsp$code)]# all pseudos: nicoi gr1+gr2+n.sp.
rag1spp <- nic12nsp$species[match(labels(rag1pseud), nic12nsp$code)]# all pseudos: nicoi gr1+gr2+n.sp.
#
cytbspp <- ttab$species[match(labels(cytbnic), ttab$code)]# just nicoi: 3 spp. delim
rag1spp <- ttab$species[match(labels(rag1nic), ttab$code)]# just nicoi: 3 spp. delim
#
cytbspp <- nic12$species[match(labels(cytbnic), nic12$code)]# just nicoi: nicoi gr1+gr2
rag1spp <- nic12$species[match(labels(rag1nic), nic12$code)]# just nicoi: nicoi gr1+gr2
#
cytbspp <- nic1nsp$species[match(labels(cytbnic), nic1nsp$code)]# just nicoi: nicoi gr1+n.sp.
rag1spp <- nic1nsp$species[match(labels(rag1nic), nic1nsp$code)]# just nicoi: nicoi gr1+n.sp.
#

#To view the nucleotide values 
lapply(nucDiag(cytbpseud, cytbspp), length)
lapply(nucDiag(rag1pseud, rag1spp), length)

lapply(nucDiag(cytbnic, cytbspp), length)
lapply(nucDiag(rag1nic, rag1spp), length)

## networks
#nn <- neighborNet(dist.dna(cytbnic, model="raw", pairwise.deletion=TRUE, as.matrix=TRUE))
#nn <- neighborNet(dist.dna(rag1nic, model="raw", pairwise.deletion=TRUE, as.matrix=TRUE))
#plot.networx(nn, "2D", edge.color="black", edge.width=1)


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
cytbtr <- midpoint(mlik$tree)
cytbtr.nic12 <- cytbtr
cytbtr.nic12nsp <- cytbtr
cytbtr.nic1nsp <- cytbtr
rag1tr <- midpoint(mlik$tree)
rag1tr.nic12 <- rag1tr
rag1tr.nic12nsp <- rag1tr
rag1tr.nic1nsp <- rag1tr


## rosenberg analysis (cytb)
cytbros <- rosenberg(cytbtr)
# create the different trees with sp delims
cytbtr$tip.label <- ttab$species[match(cytbtr$tip.label, ttab$code)]
cytbtr.nic12$tip.label <- nic12$species[match(cytbtr.nic12$tip.label, nic12$code)]
cytbtr.nic12nsp$tip.label <- nic12nsp$species[match(cytbtr.nic12nsp$tip.label, nic12nsp$code)]
cytbtr.nic1nsp$tip.label <- nic1nsp$species[match(cytbtr.nic1nsp$tip.label, nic1nsp$code)]
# extract the rosenberg P values
cytbros[names(cytbros) == getMRCA(cytbtr, tip="n. sp.")]
cytbros[names(cytbros) == getMRCA(cytbtr, tip="nicoi.gr1")]
cytbros[names(cytbros) == getMRCA(cytbtr, tip="nicoi.gr2")]
cytbros[names(cytbros) == getMRCA(cytbtr.nic12, tip="nicoi.gr1+gr2")]
cytbros[names(cytbros) == getMRCA(cytbtr.nic12nsp, tip="nicoi.gr1+gr2+n.sp.")]
cytbros[names(cytbros) == getMRCA(cytbtr.nic1nsp, tip="nicoi.gr1+n.sp.")]

## rosenberg analysis (rag1)
rag1ros <- rosenberg(rag1tr)
# create the different trees with sp delims
rag1tr$tip.label <- ttab$species[match(rag1tr$tip.label, ttab$code)]
rag1tr.nic12$tip.label <- nic12$species[match(rag1tr.nic12$tip.label, nic12$code)]
rag1tr.nic12nsp$tip.label <- nic12nsp$species[match(rag1tr.nic12nsp$tip.label, nic12nsp$code)]
rag1tr.nic1nsp$tip.label <- nic1nsp$species[match(rag1tr.nic1nsp$tip.label, nic1nsp$code)]
# extract the rosenberg P values
rag1ros[names(rag1ros) == getMRCA(rag1tr, tip="n. sp.")]
#rag1ros[names(rag1ros) == getMRCA(rag1tr, tip="nicoi.gr1")] # not monophyletic
#rag1ros[names(rag1ros) == getMRCA(rag1tr, tip="nicoi.gr2")] # not monophyletic
#rag1ros[names(rag1ros) == getMRCA(rag1tr.nic12, tip="nicoi.gr1+gr2")] # not monophyletic
rag1ros[names(rag1ros) == getMRCA(rag1tr.nic12nsp, tip="nicoi.gr1+gr2+n.sp.")]
rag1ros[names(rag1ros) == getMRCA(rag1tr.nic1nsp, tip="nicoi.gr1+n.sp.")]



plot(rag1tr)

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