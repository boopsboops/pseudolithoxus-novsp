#!/usr/bin/env bash
# to execute, run './run_phase.sh'
# the 'PHASE' and 'seqphase' programs need to be in your PATH.
# http://seqphase.mpg.de/seqphase/
# http://stephenslab.uchicago.edu/software.html
# your data needs to be split up by populations/species into separate fasta files with the suffix '.poly.fas', but this can be changed below.
# output is a fasta file with all phased individuals of all species

for i in *.poly.fas
    do echo "............... $i ................."; seqphase1.pl -1 $i -p $i
done

for j in *.inp
    do echo "............... $j ................."; PHASE -p0.5 $j $j.out
done

for k in *.inp.out
    do echo "............... $k ................."; seqphase2.pl -c ${k/inp.out/const} -i $k -o $k.phased.fas
done

cat *.phased.fas > poly_phased_all.fas
