#!/usr/bin/env Rscript
#rm(list=ls())

# read the table
tab <- read.table(file="../analyses/delimitation/18-07-16/bf_results.csv", header=FALSE, sep=",", stringsAsFactors=FALSE)
ord <- read.table(file="../analyses/delimitation/18-07-16/bf_order.csv", header=FALSE, sep=",", stringsAsFactors=FALSE)

method <- rep("", length(ord$V1))
method[grep("pathSamplingAnalysis", ord$V1)] <- "ps"
method[grep("steppingStoneSamplingAnalysis", ord$V1)] <- "ss"

model <- gsub(".mle.log>|.mle.log", "", sapply(strsplit(ord$V1, split="="), function(x) x[2]))

lik <- sapply(strsplit(tab$V1, split=" "), function(x) x[length(x)])

bfdf <- data.frame(as.factor(method), as.factor(model), as.numeric(lik))
names(bfdf) <- c("method", "model", "lik")

# subset the stepping stone results for only the combined runs
ss <- subset(subset(bfdf, method=="ss"), model=="m0run1 m0run2 m0run3" | model=="m1run1 m1run2 m1run3" | model=="m2run1 m2run2 m2run3" | model=="m3run1 m3run2 m3run3")

# add the number of species and the groupings
nsp <- c("1 sp.", "3 spp.", "2 spp.", "2 spp.")
comp <- c("nicoi.gr1+nicoi.gr2+n.sp.", "nicoi.gr1,nicoi.gr2,n.sp.", "nicoi.gr1+n.sp.,nicoi.gr2", "nicoi.gr1+nicoi.gr2,n.sp.")
ss <- cbind(ss, comp, nsp)
ss$model <- as.factor(gsub("run.*", "", as.character(ss$model)))

# rank the models according to likelihood
rank <- as.factor(c(1,2,3,4))
ssr <- cbind(ss[order(ss$lik, decreasing=TRUE), ], rank)

# calc the 2lnBf
Bf <- 2*(ssr$lik[ssr$rank==1] - ssr$lik)
ssr <- cbind(ssr, Bf)

# save the table
write.table(ssr, file="../analyses/delimitation_27-02-16/bf_table.csv", sep=",", row.names=FALSE)
write.table(ssr, file="../analyses/beast_spdelim/test-beast-1.8.4/0.1Prior/bf_table.csv", sep=",", row.names=FALSE)


# calculate Bf for other model we are interested in 
2*(ssr$lik[ssr$model=="m3"] - ssr$lik[ssr$model=="m2"])


# delimiting allopatric spp is hard (fujita 2012)
# BP&P is bad because of polytomies in guide trees (Gummer 2014, Leache 2014)
# Rather than equating gene trees with a species tree or basing species status on some genetic threshold, the relationship between the gene trees and the lineage history is modeled probabilistically using coalescent theory
# the reliability of population demographic results can be questioned when population structure exists (Pybus et al. 2009).