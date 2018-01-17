# pseudolithoxus-novsp

---

## This is the GitHub repository for: 

[Collins, R. A., Bifi, A. G., de Oliveira, R. R., Ribeiro, E. D., Lujan, N. K., Rapp Py-Daniel, L. H., Hrbek, T. (2018) **Biogeography and species delimitation of the rheophilic suckermouth-catfish genus _Pseudolithoxus_ (Siluriformes: Loricariidae), with the description of a new species from the Brazilian Amazon**. Systematics and Biodiversity XXX, xxx-xxx, DOI: XXX.](http://dx.doi/XXX)

### Abstract

The rapids-dwelling suckermouth catfish genus _Pseudolithoxus_ was previously only known from the Guiana-Shield-draining Orinoco and Casiquiare river systems of Colombia and Venezuela, but new records expand this range considerably further into the Amazon basin of Brazil, and include occurrences from rivers draining the northern Brazilian Shield. These highly disjunct records are now placed in an evolutionary and phylogeographic context using a dated species tree constructed from mitochondrial (_Cytb_) and nuclear (_RAG1_) gene sequence data. Due to mito-nuclear discordance, we also delimit the putative species using statistical coalescent models and a range of additional metrics. We infer that at least two species of _Pseudolithoxus_ are present in the Amazon basin: _P. nicoi_, previously only recorded from the Casiquiare River but now also reported from the upper rio Negro, and a new species&mdash;which we describe herein&mdash;from south-draining Guiana Shield and north-draining Brazilian Shield drainages. Our data reject a simple model of Miocene vicariance in the group following uplift of the Uaupés Arch separating the Orinoco and Amazon systems, and instead suggest more complex dispersal scenarios through palaeo-connections in the Pliocene and also via the contemporary rio Negro and rio Madeira in the late Pleistocene. 

---

## Organisation

The repository is organised into folders comprising `data`, `manuscript`, and `scripts`. File paths for scripts in the `scripts` folder will be relative to the `data` folder, or a uncommitted `temp` folder (for temporary files).

`data` contains raw data, input files for programs and results tables.

`scripts` contains command files for running the analyses. Mostly `.R` or `.sh`. 

## Files in the `data` folder

* `bf.table.csv`: Results of the Bayes factor species delimitations (Table 1 in paper).
* `combo_cytb.tre`: \*Beast _CYTB_ gene tree for all taxa.
* `combo_rag1.tre`: \*Beast _RAG1_ gene tree for all taxa.
* `combo_species.tre`: \*Beast species tree for all taxa.
* `cytb.fasta`
* `cytb_aligned.nex`: Aligned matrix of _CYTB_ for all taxa (nexus format).
* `cytb_pseudolithoxus_aligned.nex`: Aligned matrix of _CYTB_ for just _Pseudolithoxus_ (nexus format).
* `lujan_extras.nex`
* `lujan_primers.fas`
* `materials_examined_gps.csv`
* `mol_samples.csv`: Metadata for all samples used in the study (Darwin Core archive flatfile format); includes all GenBank accessions.
* `pseudolithoxus_traits_m0.tsv`: Species delimitation models (see paper)
* `pseudolithoxus_traits_m1.tsv`: Species delimitation models (see paper)
* `pseudolithoxus_traits_m2.tsv`: Species delimitation models (see paper)
* `pseudolithoxus_traits_m3.tsv`: Species delimitation models (see paper)
* `rag1.fasta`
* `rag1_aligned.nex`: Aligned matrix of _RAG1_ for all taxa (nexus format).
* `rag1_phased_aligned.nex`: Aligned matrix of phased _RAG1_ haplotypes for all taxa (nexus format).
* `rag1_pseudolithoxus_phased_aligned.nex`: Aligned matrix of phased _RAG1_ haplotypes for just _Pseudolithoxus_ (nexus format).
* `species_delimitation_results.csv`: Table of species delimitation metrics (Table 2 in paper)

