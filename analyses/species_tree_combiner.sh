#!/usr/bin/env bash

# $0 = basename of program!
# $1 = name of run1 file
# $2 = name of run2 file
# $3 = name of  combined trees output file 
# $4 = burninTrees value
# $5 = output MCC tree file

#gets the taxon block from one of the files (in this case the first one) and makes a new file.
grep -v -e "tree STATE" "$1" > "$3"
# for species trees, need to remove the last 'End;'
sed -i '$ d' "$3"
#searches for the trees using the term 'tree STATE' in both tree files and dumps them into alternate lines in the combined file 
paste -d"\n" <(grep "tree STATE" "$1") <(grep "tree STATE" "$2") >> "$3"
#adds the closing nexus code
echo "End;" >> "$3"
# run tree annotator with 10% burnin
treeannotator -burninTrees "$4" -heights mean "$3" "$5"
