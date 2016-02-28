#!/usr/bin/env sh

# download and compile/install the program
#git clone https://github.com/ddarriba/jmodeltest2.git
#cd jmodeltest2
#git tag
#git checkout <latest tag>
#ant jar
#wget https://github.com/713/project/raw/master/dependencies/jmodeltest-2.1.3/exe/phyml/PhyML_3.0_linux64
#mv PhyML_3.0_linux64 /home/rupert/Software/jmodeltest2/dist/exe/phyml
#chmod +x /home/rupert/Software/jmodeltest2/dist/exe/phyml
#chmod +x dist/jModelTest.jar

# to run jModelTest
#cd dist
# java -jar jModelTest.jar -help

# can't put a java prog on path without a helper script
java -jar /home/rupert/Software/jmodeltest2/dist/jModelTest.jar -tr 4 -d ../data/rag_phased_aligned.nex -f -g 4 -s 5 -S NNI -t ML -BIC -w -o ../temp2/rag_phased_aligned.console

java -jar /home/rupert/Software/jmodeltest2/dist/jModelTest.jar -tr 4 -d ../data/cytb_aligned.nex -f -g 4 -s 5 -S NNI -t ML -BIC -w -o ../temp2/cytb_aligned.console

