#!/usr/bin/env Rscript
# rm(list = ls())

# load scripts and libs
source("https://raw.githubusercontent.com/legalLab/protocols-scripts/master/scripts/clean_dna.R")
source("genbank_functions.R")
library("tidyverse")
library("traits")
library("ape")
library("stringr")

## script to submit sequences to GenBank

# load samples spreadsheet 
mol.tab <- read_csv(file="../data/mol_samples.csv")
glimpse(mol.tab)

# load nex
cytb.nex <- as.DNAbin(read.nexus.data(file="../data/cytb_aligned.nex"))
length(cytb.nex)
rag1.nex <- as.DNAbin(read.nexus.data(file="../data/rag1_aligned.nex"))
length(rag1.nex)

#cytb.fas <- read.dna(file="../data/cytb.fasta", format="fasta")
#length(cytb.fas)
#rag1.fas <- read.dna(file="../data/rag1.fasta", format="fasta")
#length(rag1.fas)

# check names are the same in both nex files and spreadsheet
setdiff(union(names(cytb.nex),names(rag1.nex)), mol.tab$"otherCatalogNumbers")
setdiff(mol.tab$"otherCatalogNumbers", union(names(cytb.nex),names(rag1.nex)))


# check accessions are available for all seqs
mol.tab.gb <- mol.tab %>% filter(associatedSequences!="TBC")# or filter(remarks=="genbank")
summary(mol.tab.gb$"otherCatalogNumbers" %in% names(cytb.nex))
summary(mol.tab.gb$"otherCatalogNumbers" %in% names(rag1.nex))

# remove the already accessioned data from the DNA files 
cytb.nex.sub <- cytb.nex[which(!names(cytb.nex) %in% mol.tab.gb$"otherCatalogNumbers")]
rag1.nex.sub <- rag1.nex[which(!names(rag1.nex) %in% mol.tab.gb$"otherCatalogNumbers")]


# dealign the sequences leaving Ns not at ends
cytb.nex.de <- rm_missing(dna=cytb.nex.sub, chars=c("?","-","n"), endOnlyNs=TRUE, print=TRUE)
rag1.nex.de <- rm_missing(dna=rag1.nex.sub, chars=c("?","-","n"), endOnlyNs=TRUE, print=TRUE)

# convert back to char 
cytb.nex.ch <- sapply(as.character(as.list(cytb.nex.de)), paste, collapse="")
rag1.nex.ch <- sapply(as.character(as.list(rag1.nex.de)), paste, collapse="")


# filter and add to table
mol.tab.tbc <- mol.tab %>% filter(associatedSequences=="TBC")# or filter(remarks=="genbank")
mol.tab.nuc <- mol.tab.tbc %>% mutate(cytb=cytb.nex.ch[match(otherCatalogNumbers, names(cytb.nex.ch))], rag1=rag1.nex.ch[match(otherCatalogNumbers, names(rag1.nex.ch))])

# filter those with data
mol.tab.nuc.cytb <- mol.tab.nuc %>% filter(!is.na(cytb))
mol.tab.nuc.rag1 <- mol.tab.nuc %>% filter(!is.na(rag1))

# Here we write a GenBank format fasta file containing critical source modifying annotations

# cytb
gb.fas.cytb <- gb_format_fasta_cytb(df=mol.tab.nuc.cytb, gene="CYTB")
write(gb.fas.cytb, file="../temp2/cytb/sequences.fsa", append=FALSE)# write out the fasta file
#
gb.fas.rag1 <- gb_format_fasta_rag1(df=mol.tab.nuc.rag1, gene="RAG1")
write(gb.fas.rag1, file="../temp2/rag1/sequences.fsa", append=FALSE)# write out the fasta file


# https://www.ncbi.nlm.nih.gov/nuccore/NC_005089.1# mus cytb
gb.feat.cytb <- gb_features_cytb(df=mol.tab.nuc.cytb, gene="CYTB", product="cytochrome b")
write(gb.feat.cytb, file="../temp2/cytb/features.tbl", append=FALSE)# write out

# https://www.ncbi.nlm.nih.gov/nuccore/AY215075.1# mus rag1
gb.feat.rag1 <- gb_features_rag1(df=mol.tab.nuc.rag1, gene="RAG1", product="recombination activating 1")
write(gb.feat.rag1, file="../temp2/rag1/features.tbl", append=FALSE)# write out

# NOW edit T07157 by hand - poor data

# to run tbl2asn in terminal
#./tbl2asn -t template.sbt -i sequences.fsa -f features.tbl -a s -V vb -T

## Once done, add the accessions and make a Supplementary table


# add the accessions from GenBank to spreadsheet
cytb.acc <- read_csv(file="../temp2/cytb/CYTB.acc")
rag1.acc <- read_csv(file="../temp2/rag1/RAG1.acc")

# join the cytb and rag1
mol.tab$tempcytb<-rep(NA, length(dim(mol.tab)[1]))
mol.tab$temprag1<-rep(NA, length(dim(mol.tab)[1]))
mol.tab$tempcytb[match(cytb.acc$code, mol.tab$otherCatalogNumbers)] <- cytb.acc$acc
mol.tab$temprag1[match(rag1.acc$code, mol.tab$otherCatalogNumbers)] <- rag1.acc$acc
cyt.rag <- paste(mol.tab$tempcytb, mol.tab$temprag1, sep=" | ")
cyt.rag[which(mol.tab$associatedSequences!="TBC")] <- mol.tab$associatedSequences[which(mol.tab$associatedSequences!="TBC")]
mol.tab$associatedSequences <- cyt.rag

# see cytb accs
str_split_fixed(string=mol.tab$associatedSequences, pattern=" \\| ", n=2)[,1]
# see rag1 accs
str_split_fixed(string=mol.tab$associatedSequences, pattern=" \\| ", n=2)[,2]

# filter the DwC data
suppl.tab <- mol.tab %>% #
    mutate(scientificName=ifelse(test=is.na(identificationQualifier), yes=paste(genus,specificEpithet), no=paste(genus,identificationQualifier))) %>% #
    select(otherCatalogNumbers,associatedSequences,catalogNumber,institutionCode,collectionCode,basisOfRecord,#
    typeStatus,scientificName,kingdom,phylum,class,order,family,genus,specificEpithet,scientificNameAuthorship,identificationQualifier,
    taxonRank,identifiedBy,decimalLatitude,decimalLongitude,country,stateProvince,waterBody,eventDate,locality,recordedBy)
#    
dim(suppl.tab)

# write out 
write_csv(suppl.tab, path="../manuscript/tables/Supplementary_Table_S2.csv")
