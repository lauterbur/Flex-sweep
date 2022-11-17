#!/bin/bash

## this will cut the given hap and map files to regions of 1.2 Mb for use with dl_sweeps

OUT=$1
CLASSNAME=$2
HAP=$3
MAP=$4
WINDOW=1200000
STEP=$5

FIRSTPOS=$( head -1 ${MAP} | awk '{print $2}' )
LASTPOS=$( tail -1 ${MAP} | awk '{print $2}' )
SIZE=$(( $LASTPOS - $FIRSTPOS ))
NUMBER=$(( ${SIZE} / ${STEP}  ))
echo "number of windows: " $NUMBER

STARTPOS=$(( $FIRSTPOS - $WINDOW/2 ))
if [[ $STARTPOS -lt 0 ]]
then
        STARTPOS=0
fi

while [ $(( $STARTPOS + $STEP )) -lt $LASTPOS ]
do
        ENDPOS=$(( $WINDOW + $STARTPOS ))
        awk -v start="$STARTPOS" -v stop="$ENDPOS" '$2 >= start && $2 <= stop {print $0}' $MAP > ${OUT}/${CLASSNAME}_${STARTPOS}-${ENDPOS}.map
        read -r STARTHAP<<<$(awk -v start="$STARTPOS" -v stop="$ENDPOS" '$2 >= start && $2 <= stop {print $2,NR}' $MAP | awk '{print $2}' | head -1)
        read -r STOPHAP<<<$(awk -v start="$STARTPOS" -v stop="$ENDPOS" '$2 >= start && $2 <= stop {print $2,NR}' $MAP | awk '{print $2}' | tail -1)
        sed -n "${STARTHAP},${STOPHAP}p" $HAP > ${OUT}/${CLASSNAME}_${STARTPOS}-${ENDPOS}.hap

        ml=$( wc -l ${OUT}/${CLASSNAME}_${STARTPOS}-${ENDPOS}.map)
        wc -l ${OUT}/${CLASSNAME}_${STARTPOS}-${ENDPOS}.hap

        if [[ $ml -lt 100 ]]
        then
                echo "fewer than 100 variants in region, removing"
                rm ${OUT}/${CLASSNAME}_${STARTPOS}-${ENDPOS}.map
                rm ${OUT}/${CLASSNAME}_${STARTPOS}-${ENDPOS}.hap
        else
                echo "${STARTPOS}-${ENDPOS}" >> ${OUT}/classification_windows.txt
        fi

        STARTPOS=$(( ${ENDPOS} - ${WINDOW} + ${STEP} ))
done
