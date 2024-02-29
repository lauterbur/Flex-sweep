#!/bin/bash

SIM=$1
SIMPATH=$(dirname "${SIM}")
SIMBASE=$(basename "${SIM}" .out)
LOCUS=$2
HALFLOCUS=$(( (LOCUS + 1) / 2 )) # adding the 1 (really denominator - 1) makes it round
START=$3

# set center point
for CENTER in `seq 500000 10000 700000`
do
	OUTNAME=${SIMBASE}_c${CENTER}

	#echo "Rscript --vanilla calculate_stats/convertHapMap.R -m $SIMPATH/$SIMBASE.map -p $SIMPATH/$SIMBASE.hap -o $SIMPATH/$SIMBASE -l ${LOCUS} -c ${CENTER} -w 1000000 -s $START"
	Rscript --vanilla calculate_stats/convertHapMap.R -m $SIMPATH/$SIMBASE.map -p $SIMPATH/$SIMBASE.hap -o $SIMPATH/$SIMBASE -l ${LOCUS} -c ${CENTER} -w 1000000 -s $START

        # iSAFE - if it works with iSAFE, good, if not, run SAFE
        echo "calculating iSAFE"
        if python2 calculate_stats/iSAFE-1.0/src/isafe.py \
                --format hap -i ${SIMPATH}/${SIMBASE}/${OUTNAME}.hap.iSAFE \
                -o ${SIMPATH}/${SIMBASE}/stats/center_${CENTER}/${OUTNAME} \
                --MaxFreq 1 \
                --MinRegionSize-bp 49000 \
                --MinRegionSize-ps 300 \
                --IgnoreGaps 
        then echo "iSAFE successful" 
        else
                echo "too few SNPs for iSAFE, ran SAFE"
                python2 calculate_stats/iSAFE-1.0/src/isafe.py \
                --format hap -i ${SIMPATH}/${SIMBASE}/${OUTNAME}.hap.iSAFE \
                -o ${SIMPATH}/${SIMBASE}/stats/center_${CENTER}/${OUTNAME} \
                --MaxFreq 1 \
                --MinRegionSize-bp 49000 \
                --MinRegionSize-ps 300 \
                --IgnoreGaps \
                --SAFE
        fi
	mv ${SIMPATH}/${SIMBASE}/stats/center_${CENTER}/${OUTNAME}.iSAFE.out ${SIMPATH}/${SIMBASE}/stats/center_${CENTER}/${OUTNAME}.iSAFE
        mv ${SIMPATH}/${SIMBASE}/stats/center_${CENTER}/${OUTNAME}.SAFE.out ${SIMPATH}/${SIMBASE}/stats/center_${CENTER}/${OUTNAME}.iSAFE

	# nSL
        # using scikit-allel
        python3.6 calculate_stats/nsl_classify.py ${SIMPATH}/${SIMBASE}/${OUTNAME}.hap ${SIMPATH}/${SIMBASE}/${OUTNAME}.map ${SIMPATH}/${SIMBASE}/stats/center_${CENTER}/${OUTNAME}.nsl
	# iHS
	### outputs both standardized and unstandardize scores - default standardization is 2% bins - use undstandardized, will standardize later
	echo "calculating iHS"
	/hapbin/build/ihsbin --hap ${SIMPATH}/${SIMBASE}/${OUTNAME}.hap --map ${SIMPATH}/${SIMBASE}/${OUTNAME}.map --out ${SIMPATH}/${SIMBASE}/stats/center_${CENTER}/${OUTNAME}.ihs.temp --minmaf 0.05 --cutoff 0.05
        awk 'FNR==NR {dict[$2]=$4; next} {$1=($1 in dict) ? dict[$2] : $1}1' ${SIMPATH}/${SIMBASE}/${OUTNAME}.map ${SIMPATH}/${SIMBASE}/stats/center_${CENTER}/${OUTNAME}.ihs.temp > ${SIMPATH}/${SIMBASE}/stats/center_${CENTER}/${OUTNAME}.ihs
done

CENTER=${HALFLOCUS}
OUTNAME=${SIMBASE}_c${CENTER}
        # run hap and map script for full chromosome first
        echo "Rscript --vanilla calculate_stats/convertHapMapLocus.R -m $SIMPATH/$SIMBASE.map -p $SIMPATH/$SIMBASE.hap -o $SIMPATH/$SIMBASE -l ${LOCUS} -c ${HALFLOCUS} -w ${LOCUS} -s $START"
        Rscript --vanilla calculate_stats/convertHapMapLocus.R -m $SIMPATH/$SIMBASE.map -p $SIMPATH/$SIMBASE.hap -o $SIMPATH/$SIMBASE -l ${LOCUS} -c ${HALFLOCUS} -w ${LOCUS} -s $START


	# compress because that's what the perl scripts want
	gzip < ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.map > ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.map.gz
	gzip < ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.hap > ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.hap.gz

        WINDOW=500000 #for H12 and HAF, use default (50000) for all others
        # this is not the same as the Flexsweep window
        perl calculate_stats/calculate_H12.pl ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.hap.gz ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.map.gz ${SIMPATH}/${SIMBASE}/stats/${OUTNAME}.H12 ${LOCUS}
        perl calculate_stats/calculate_HAF.pl ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.hap.gz ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.map.gz ${SIMPATH}/${SIMBASE}/stats/${OUTNAME}.HAF ${LOCUS}
	python3.6 calculate_stats/hapDAF-o.py --pos ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.map.gz --hap ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.hap.gz -o ${SIMPATH}/${SIMBASE}/stats/${OUTNAME}.hDo
	python3.6 calculate_stats/hapDAF-s.py --pos ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.map.gz --hap ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.hap.gz -o ${SIMPATH}/${SIMBASE}/stats/${OUTNAME}.hDs
	python3.6 calculate_stats/DIND.py --pos ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.map.gz --hap ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.hap.gz -o ${SIMPATH}/${SIMBASE}/stats/${OUTNAME}.DIND
	python3.6 calculate_stats/highfreq.py --pos ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.map.gz --hap ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.hap.gz -o ${SIMPATH}/${SIMBASE}/stats/${OUTNAME}.hf
	python3.6 calculate_stats/lowfreq.py --pos ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.map.gz --hap ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.hap.gz -o ${SIMPATH}/${SIMBASE}/stats/${OUTNAME}.lf
	python3.6 calculate_stats/Sratio.py --pos ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.map.gz --hap ${SIMPATH}/${SIMBASE}/${OUTNAME}_full.hap.gz -o ${SIMPATH}/${SIMBASE}/stats/${OUTNAME}.S
