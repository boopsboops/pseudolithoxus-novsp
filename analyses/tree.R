# start the chunk
## @knitr treecode

# download seqs from Lujan et al. 2014
# testing with some Rhinichthys

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
