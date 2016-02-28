#!/usr/bin/env sh

# run all models thrice

# m0
beast -threads 0 -beagle -beagle_SSE -beagle_CPU -beagle_double -beagle_scaling none -seed 1270216 -overwrite m0run1.xml
beast -threads 0 -beagle -beagle_SSE -beagle_CPU -beagle_double -beagle_scaling none -seed 2270216 -overwrite m0run2.xml
beast -threads 0 -beagle -beagle_SSE -beagle_CPU -beagle_double -beagle_scaling none -seed 3270216 -overwrite m0run3.xml

# m1
beast -threads 0 -beagle -beagle_SSE -beagle_CPU -beagle_double -beagle_scaling none -seed 1270216 -overwrite m1run1.xml
beast -threads 0 -beagle -beagle_SSE -beagle_CPU -beagle_double -beagle_scaling none -seed 2270216 -overwrite m1run2.xml
beast -threads 0 -beagle -beagle_SSE -beagle_CPU -beagle_double -beagle_scaling none -seed 3270216 -overwrite m1run3.xml

# m2
beast -threads 0 -beagle -beagle_SSE -beagle_CPU -beagle_double -beagle_scaling none -seed 1270216 -overwrite m2run1.xml
beast -threads 0 -beagle -beagle_SSE -beagle_CPU -beagle_double -beagle_scaling none -seed 2270216 -overwrite m2run2.xml
beast -threads 0 -beagle -beagle_SSE -beagle_CPU -beagle_double -beagle_scaling none -seed 3270216 -overwrite m2run3.xml

# m3
beast -threads 0 -beagle -beagle_SSE -beagle_CPU -beagle_double -beagle_scaling none -seed 1270216 -overwrite m3run1.xml
beast -threads 0 -beagle -beagle_SSE -beagle_CPU -beagle_double -beagle_scaling none -seed 2270216 -overwrite m3run2.xml
beast -threads 0 -beagle -beagle_SSE -beagle_CPU -beagle_double -beagle_scaling none -seed 3270216 -overwrite m3run3.xml


# run the Bayes factor calculations and extract the results
beast bf_calculations.xml > bf_results.txt
grep "log marginal likelihood" bf_results.txt > bf_results.csv
grep "fileName=" bf_calculations.xml > bf_order.csv