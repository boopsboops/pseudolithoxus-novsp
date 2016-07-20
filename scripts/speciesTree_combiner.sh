#!/usr/bin/env bash
# needs to be bash and not dash!

# $0 = basename of program!
# $1 = name of run1 file
# $2 = name of run2 file
# $3 = name of run3 file
# $4 = name of run4 file
# $5 = name of run5 file
# $6 = name of run6 file
# $7 = name of run7 file
# $8 = name of run8 file
# $9 = name of  combined trees output file 
# $10 = burninTrees value
# $11 = output MCC tree file

#gets the taxon block from one of the files (in this case the first one) and makes a new file.
grep -v -e "tree STATE" ${1} > ${9}
# for species trees, need to remove the last 'End;'
sed -i '$ d' ${9}
#searches for the trees using the term 'tree STATE' in both tree files and dumps them into alternate lines in the combined file 
paste -d"\n" <(grep "tree STATE" ${1}) <(grep "tree STATE" ${2}) <(grep "tree STATE" ${3}) <(grep "tree STATE" ${4}) <(grep "tree STATE" ${5}) <(grep "tree STATE" ${6}) <(grep "tree STATE" ${7}) <(grep "tree STATE" ${8}) >> ${9}
#adds the closing nexus code
echo "End;" >> ${9}
# run tree annotator with 10% burnin
treeannotator -burninTrees ${10} -heights mean ${9} ${11}


# to run
#./speciesTree_combiner.sh ../analyses/species_tree_27-02-16/RUN1.species.trees ../analyses/species_tree_27-02-16/RUN2.species.trees ../analyses/species_tree_27-02-16/RUN3.species.trees ../analyses/species_tree_27-02-16/RUN4.species.trees ../analyses/species_tree_27-02-16/RUN5.species.trees ../analyses/species_tree_27-02-16/RUN6.species.trees ../analyses/species_tree_27-02-16/RUN7.species.trees ../analyses/species_tree_27-02-16/RUN8.species.trees ../analyses/species_tree_27-02-16/COMBO.species.trees 1608 ../analyses/species_tree_27-02-16/COMBO.species.tre

#./speciesTree_combiner.sh ../analyses/species_tree_27-02-16/RUN1.cytb.trees ../analyses/species_tree_27-02-16/RUN2.cytb.trees ../analyses/species_tree_27-02-16/RUN3.cytb.trees ../analyses/species_tree_27-02-16/RUN4.cytb.trees ../analyses/species_tree_27-02-16/RUN5.cytb.trees ../analyses/species_tree_27-02-16/RUN6.cytb.trees ../analyses/species_tree_27-02-16/RUN7.cytb.trees ../analyses/species_tree_27-02-16/RUN8.cytb.trees ../analyses/species_tree_27-02-16/COMBO.cytb.trees 1608 ../analyses/species_tree_27-02-16/COMBO.cytb.tre

#./speciesTree_combiner.sh ../analyses/species_tree_27-02-16/RUN1.rag1_phased.trees ../analyses/species_tree_27-02-16/RUN2.rag1_phased.trees ../analyses/species_tree_27-02-16/RUN3.rag1_phased.trees ../analyses/species_tree_27-02-16/RUN4.rag1_phased.trees ../analyses/species_tree_27-02-16/RUN5.rag1_phased.trees ../analyses/species_tree_27-02-16/RUN6.rag1_phased.trees ../analyses/species_tree_27-02-16/RUN7.rag1_phased.trees ../analyses/species_tree_27-02-16/RUN8.rag1_phased.trees ../analyses/species_tree_27-02-16/COMBO.rag1_phased.trees 1608 ../analyses/species_tree_27-02-16/COMBO.rag1_phased.tre


# run July 2016
#../../scripts/speciesTree_combiner.sh run1_species.trees run2_species.trees run3_species.trees run4_species.trees run5_species.trees run6_species.trees run7_species.trees run8_species.trees combo_species.trees 1608 combo_species.tre

#../../scripts/speciesTree_combiner.sh run1_rag1.trees run2_rag1.trees run3_rag1.trees run4_rag1.trees run5_rag1.trees run6_rag1.trees run7_rag1.trees run8_rag1.trees combo_rag1.trees 1608 combo_rag1.tre

#../../scripts/speciesTree_combiner.sh run1_cytb.trees run2_cytb.trees run3_cytb.trees run4_cytb.trees run5_cytb.trees run6_cytb.trees run7_cytb.trees run8_cytb.trees combo_cytb.trees 1608 combo_cytb.tre
