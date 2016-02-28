#!/bin/sh

# to create the new files
cp RUN1.xml RUN.xml
for i in `seq 2 10`; do cp RUN.xml RUN$i.xml; done
rm RUN.xml

# to change the number gens
# will substitute {} with the filename(s) found.
# the ';' is the end of command, and the '\' is  to escape the interpretation of ; by shell
find . -name "*.xml" -type f -exec sed -i 's/chainLength="500000000"/chainLength="50000000"/g' {} \;


# to change the output file names
# i looked at incrementally changing, but lost patience.
sed -i 's/RUN1/RUN2/g' RUN2.xml
sed -i 's/RUN1/RUN3/g' RUN3.xml
sed -i 's/RUN1/RUN4/g' RUN4.xml
sed -i 's/RUN1/RUN5/g' RUN5.xml
sed -i 's/RUN1/RUN6/g' RUN6.xml
sed -i 's/RUN1/RUN7/g' RUN7.xml
sed -i 's/RUN1/RUN8/g' RUN8.xml
sed -i 's/RUN1/RUN9/g' RUN9.xml
sed -i 's/RUN1/RUN10/g' RUN10.xml
