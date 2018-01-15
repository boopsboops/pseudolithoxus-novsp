#!/usr/bin/env Rscript
# rm(list = ls())

# load libs
#source("functions_libs.R")
library("tidyverse")
library("traits")
library("ape")


## to make a suppl data table

mol.tab <- read_csv(file="../data/mol_samples.csv")

cytb <- as.matrix(read.nexus.data(file="../data/cytb_aligned.nex"))
rag1 <- as.matrix(read.nexus.data(file="../data/rag1_aligned.nex"))

# check names are the same in both nex file and spreadsheet
setdiff(rownames(cytb), mol.tab$otherCatalogNumbers)
setdiff(rownames(rag1), mol.tab$otherCatalogNumbers)


