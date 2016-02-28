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
