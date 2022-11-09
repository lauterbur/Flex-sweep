#!/bin/bash

SIM=$1
SIMPATH=$(dirname "${SIM}")
SIMBASE=$(basename "${SIM}" .out)
HOMEDIR=$(dirname $(dirname $(dirname ${SIMPATH})))
LOCUS=$2
RECOMB=$3

HALFLOCUS=$(( (LOCUS + 1) / 2 )) # adding the 1 (really denominator - 1) makes it round

OUTNAME=${SIMBASE}
        if [ -s ${SIM} ]
        then
                if [ -z "$RECOMB" ]
                then
                        Rscript --vanilla calculate_stats/discoal2hapmap_fullchromosome.R -f $SIM -o $SIMPATH -l $LOCUS -c $HALFLOCUS -w $LOCUS
                else
                	Rscript --vanilla calculate_stats/discoal2hapmap_fullchromosome.R -f $SIM -o $SIMPATH -l $LOCUS -c $HALFLOCUS -w $LOCUS -r $RECOMB
	        	# center=600000 and window=1200000 makes it use the whole chromosome
                fi
        else
                echo "missing $SIM"
                exit 1
        fi

        #iSAFE
        if python2 calculate_stats/iSAFE-1.0/src/isafe.py \
		--format hap -i ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.hap.iSAFE \
		-o ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS} \
		--MaxFreq 1 \
		--MinRegionSize-bp 49000 \
		--MinRegionSize-ps 300 \
		--IgnoreGaps 
        then echo "iSAFE successful"
        else
                echo "too few SNPs for iSAFE, run SAFE"
                python2 calculate_stats/iSAFE-1.0/src/isafe.py \
		--format hap -i ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.hap.iSAFE \
		-o ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS} \
                --MaxFreq 1 \
                --MinRegionSize-bp 49000 \
                --MinRegionSize-ps 300 \
                --IgnoreGaps \
                --SAFE
        fi

	mv ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.iSAFE.out ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.iSAFE
	mv ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.SAFE.out ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.iSAFE

	# nSL
	### outputs only unstandarized scores
	rm ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.nsl.log
        # using scikit-allel
        python3.6 ${HOMEDIR}/calculate_stats/nsl_scikit.py ${SIM} ${SIMPATH}/${OUTNAME}_c600000_full.map ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.nsl

	# iHS
	/hapbin/build/ihsbin --hap ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.hap --map ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.map --out ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.ihs.temp --minmaf 0.05 --cutoff 0.1
        awk 'FNR==NR {dict[$2]=$4; next} {$1=($1 in dict) ? dict[$2] : $1}1' ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.map ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.ihs.temp > ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.ihs
        rm ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.ihs.temp

        # compress because that's what the perl scripts want

        if [ -s ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.map ]
        then
                gzip < ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.map > ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.map.gz
        fi
        if [ ! -s "${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.map.gz" ]
        then
                echo "missing ${OUTNAME}_c${HALFLOCUS}_full.map.gz"
                exit 1
        fi
        if [ -s ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.hap ]
        then
                gzip < ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.hap > ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.hap.gz
        fi

	WINDOW=500000 # same window for H12 and HAF, different (default = 50000) for all the others
                ### this is different than the pipeline window size
        perl calculate_stats/calculate_H12.pl ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.hap.gz ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.map.gz ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.H12 $WINDOW
        perl calculate_stats/calculate_HAF.pl ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.hap.gz ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.map.gz ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.HAF $WINDOW

        python3.6 calculate_stats/hapDAF-o.py --pos ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.map.gz --hap ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.hap.gz -o ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.hDo
        python3.6 calculate_stats/hapDAF-s.py --pos ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.map.gz --hap ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.hap.gz -o ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.hDs
        python3.6 calculate_stats/DIND.py --pos ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.map.gz --hap ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.hap.gz -o ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.DIND
        python3.6 calculate_stats/highfreq.py --pos ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.map.gz --hap ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.hap.gz -o ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.hf
        python3.6 calculate_stats/lowfreq.py --pos ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.map.gz --hap ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.hap.gz -o ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.lf
        python3.6 calculate_stats/Sratio.py --pos ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.map.gz --hap ${SIMPATH}/${OUTNAME}_c${HALFLOCUS}_full.hap.gz -o ${SIMPATH}/stats/${OUTNAME}_c${HALFLOCUS}.S
