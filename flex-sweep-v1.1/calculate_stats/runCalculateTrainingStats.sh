#!/bin/bash

SIM=$1
SIMPATH=$(dirname "${SIM}")
SIMBASE=$(basename "${SIM}" .out)
HOMEDIR=$(dirname $(dirname $(dirname ${SIMPATH})))

LOCUS=$2
HALFLOCUS=$(( (LOCUS + 1) / 2 )) # adding the 1 (really denominator - 1) makes it round

RECOMB=$3

# set center point
for CENTER in `seq 500000 10000 700000`
do
	OUTNAME=${SIMBASE}_c${CENTER}

        if [ -s ${SIM} ]
        then
                if [ -z "$RECOMB" ]
                then
                        Rscript --vanilla calculate_stats/discoal2hapmap.R -f $SIM -o $SIMPATH -l $LOCUS -c $CENTER -w 1000000           
                else
                        Rscript --vanilla calculate_stats/discoal2hapmap.R -f $SIM -o $SIMPATH -l $LOCUS -c $CENTER -w 1000000 -r $RECOMB
                        # center=600000 and window=1200000 makes it use the whole chromosome
                fi
        else
                echo "missing $SIM"
                exit 1
        fi


	# iSAFE - if it works with iSAFE, good, if not, run SAFE
        if python2 calculate_stats/iSAFE-1.0/src/isafe.py \
		--format hap -i ${SIMPATH}/${OUTNAME}.hap.iSAFE \
		-o ${SIMPATH}/stats/center_${CENTER}/${OUTNAME} \
		--MaxFreq 1 \
		--MinRegionSize-bp 49000 \
		--MinRegionSize-ps 300 \
		--IgnoreGaps 
        then echo "iSAFE successful" 
        else
                echo "too few SNPs for iSAFE, ran SAFE"
        	python2 calculate_stats/iSAFE-1.0/src/isafe.py \
		--format hap -i ${SIMPATH}/${OUTNAME}.hap.iSAFE \
		-o ${SIMPATH}/stats/center_${CENTER}/${OUTNAME} \
		--MaxFreq 1 \
		--MinRegionSize-bp 49000 \
		--MinRegionSize-ps 300 \
		--IgnoreGaps \
                --SAFE
        fi

	mv ${SIMPATH}/stats/center_${CENTER}/${OUTNAME}.iSAFE.out ${SIMPATH}/stats/center_${CENTER}/${OUTNAME}.iSAFE
	mv ${SIMPATH}/stats/center_${CENTER}/${OUTNAME}.SAFE.out ${SIMPATH}/stats/center_${CENTER}/${OUTNAME}.iSAFE
        
	# nSL
        # using scikit-allel
        python3.6 ${HOMEDIR}/calculate_stats/nsl_scikit.py ${SIM} ${SIMPATH}/${OUTNAME}.map ${SIMPATH}/stats/center_${CENTER}/${OUTNAME}.nsl

	# iHS
	### outputs both standardized and unstandardize scores - default standardization is 2% bins - use undstandardized, will standardize later
	/hapbin/build/ihsbin --hap ${SIMPATH}/${OUTNAME}.hap --map ${SIMPATH}/${OUTNAME}.map --out ${SIMPATH}/stats/center_${CENTER}/${OUTNAME}.ihs.temp --minmaf 0.05 --cutoff 0.05
        awk 'FNR==NR {dict[$2]=$4; next} {$1=($1 in dict) ? dict[$2] : $1}1' ${SIMPATH}/${OUTNAME}.map ${SIMPATH}/stats/center_${CENTER}/${OUTNAME}.ihs.temp > ${SIMPATH}/stats/center_${CENTER}/${OUTNAME}.ihs
        rm ${SIMPATH}/stats/center_${CENTER}/${OUTNAME}.ihs.temp
done

if [[ $SIMBASE == "sweep_"* ]] # neutral window stats were created with run_calculate_stats_neutral_norm.sh, no need to redo
then
CENTER=${HALFLOCUS}
OUTNAME=${SIMBASE}_c${CENTER}

        if [ -z "$RECOMB" ]
        then
                Rscript --vanilla calculate_stats/discoal2hapmap_fullchromosome.R -f $SIM -o $SIMPATH -l $LOCUS -c $HALFLOCUS -w $LOCUS
        else
                Rscript --vanilla calculate_stats/discoal2hapmap_fullchromosome.R -f $SIM -o $SIMPATH -l $LOCUS -c $HALFLOCUS -w $LOCUS -r $RECOMB
                # center=600000 and window=1200000 makes it use the whole chromosome
        fi

	# compress because that's what the perl scripts want
        if [ -s ${SIMPATH}/${OUTNAME}_full.map ]
        then
        	gzip < ${SIMPATH}/${OUTNAME}_full.map > ${SIMPATH}/${OUTNAME}_full.map.gz
        fi
        if [ ! -s "${SIMPATH}/${OUTNAME}_full.map.gz" ]
        then
                echo "missing ${OUTNAME}_full.map.gz"
                exit 1
        fi
        if [ -s ${SIMPATH}/${OUTNAME}_full.hap ]
        then
        	gzip < ${SIMPATH}/${OUTNAME}_full.hap > ${SIMPATH}/${OUTNAME}_full.hap.gz
        fi
        if [ ! -s ${SIMPATH}/${OUTNAME}_full.map ]
        then
                echo "missing ${OUTNAME}_full.map.gz"
                exit 1
        fi

        WINDOW=500000 #for H12 and HAF, use default (50000) for all others
        ## this is not the same as the pipeline window
	perl calculate_stats/calculate_H12.pl ${SIMPATH}/${OUTNAME}_full.hap.gz ${SIMPATH}/${OUTNAME}_full.map.gz ${SIMPATH}/stats/${OUTNAME}.H12 ${LOCUS}
	perl calculate_stats/calculate_HAF.pl ${SIMPATH}/${OUTNAME}_full.hap.gz ${SIMPATH}/${OUTNAME}_full.map.gz ${SIMPATH}/stats/${OUTNAME}.HAF ${LOCUS}
	python3.6 calculate_stats/hapDAF-o.py --pos ${SIMPATH}/${OUTNAME}_full.map.gz --hap ${SIMPATH}/${OUTNAME}_full.hap.gz -o ${SIMPATH}/stats/${OUTNAME}.hDo
	python3.6 calculate_stats/hapDAF-s.py --pos ${SIMPATH}/${OUTNAME}_full.map.gz --hap ${SIMPATH}/${OUTNAME}_full.hap.gz -o ${SIMPATH}/stats/${OUTNAME}.hDs
	python3.6 calculate_stats/DIND.py --pos ${SIMPATH}/${OUTNAME}_full.map.gz --hap ${SIMPATH}/${OUTNAME}_full.hap.gz -o ${SIMPATH}/stats/${OUTNAME}.DIND
	python3.6 calculate_stats/highfreq.py --pos ${SIMPATH}/${OUTNAME}_full.map.gz --hap ${SIMPATH}/${OUTNAME}_full.hap.gz -o ${SIMPATH}/stats/${OUTNAME}.hf
	python3.6 calculate_stats/lowfreq.py --pos ${SIMPATH}/${OUTNAME}_full.map.gz --hap ${SIMPATH}/${OUTNAME}_full.hap.gz -o ${SIMPATH}/stats/${OUTNAME}.lf
	python3.6 calculate_stats/Sratio.py --pos ${SIMPATH}/${OUTNAME}_full.map.gz --hap ${SIMPATH}/${OUTNAME}_full.hap.gz -o ${SIMPATH}/stats/${OUTNAME}.S
fi
