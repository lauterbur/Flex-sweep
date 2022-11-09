#!/bin/bash

echo "simulating"
time singularity exec flex-sweep.sif python3 FlexSweep_simulate.py fullTest ../example/config.txt 1000 20 --numJobs 2 --continue
# singularity exec flex-sweep.sif python3 FlexSweep_simulate.py <outputDir> <config.txt> <numSims> <numChroms> (--numJobs <>) (--continue) (--locusLength <>)

echo "calculating"
time singularity exec flex-sweep.sif python3 FlexSweep_fv.py fullTest 20 --keepSims --keepStats
#singularity exec flex-sweep.sif python3 FlexSweep_fv.py <outputDir> <numChroms> (--numJobs <>) (--continue) (--locusLength <>) (--rMap True) (--keepSims) (--keepStats) (--distCenters <NOT IN USE>) (--minCenter <NOT IN USE>) (--maxCenter <NOT IN USE>)

echo "training"
time singularity exec train.sif python3 FlexSweep_train.py fullTest --fvSplit 80 80
#singularity exec train.sif python3 FlexSweep_train.py <outputDir> (--fvLoc <path>) (--fvSplit 10000 80) (--continue)

#singularity exec flex-sweep.sif python3 FlexSweep_plot.py train ... # plots!

echo "classifying"
time singularity exec flex-sweep.sif python3 FlexSweep_classify.py fullTest fullTestData --hapmap test.hap test.map --keepWindows --keepStats 
time singularity exec train.sif python3 FlexSweep_classify.py fullTest fullTestData --hapmap test.hap test.map --keepWindows --keepStats --continue
#singularity exec train.sif python3 FlexSweep_classify.py <outputdir> <classifyName> (--hapmap <path.hap> <path.map> | --vcf <path.vcf>) (--keepWindows) (--keepStats) (--continue) (--modelLoc <path>) (--normLoc <path>) (--windowStep 10000) (--numJobs <val>) (--locusLength 1200000) (--distCenters 10000) (--minCenter 500000) (--maxCenter 700000) (--threshold 0.5) (--classifyOnly test/classification/fvs) 

#singularity exec flex-sweep.sif python3 FlexSweep_plot.py classify ... # plots!
